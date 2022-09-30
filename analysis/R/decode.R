# Copyright 2014 Google Inc. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#
# This library implements the RAPPOR marginal decoding algorithms using LASSO.

suppressPackageStartupMessages(library(glmnet))

# So we don't have to change pwd
source.rappor <- function(rel_path)  {
  abs_path <- paste0(Sys.getenv("RAPPOR_REPO", ""), rel_path)
  source(abs_path)
}

source.rappor('analysis/R/alternative.R')

EstimateBloomCounts <- function(params, obs_counts) {
  # Estimates the number of times each bit in each cohort was set in original
  # Bloom filters.
  #
  # Input:
  #    params: a list of RAPPOR parameters:
  #            k - size of a Bloom filter
  #            h - number of hash functions
  #            m - number of cohorts
  #            p - P(IRR = 1 | PRR = 0)
  #            q - P(IRR = 1 | PRR = 1)
  #            f - Proportion of bits in the Bloom filter that are set randomly
  #                to 0 or 1 regardless of the underlying true bit value
  #    obs_counts: a matrix of size m by (k + 1). Column one contains sample
  #                sizes for each cohort. Other counts indicated how many times
  #                each bit was set in each cohort.
  #
  # Output:
  #    ests: a matrix of size m by k with estimated counts for the probability
  #          of each bit set to 1 in the true Bloom filter.
  #    stds: standard deviation of the estimates.

  p <- params$p
  q <- params$q
  f <- params$f
  m <- params$m
  k <- params$k

  stopifnot(m == nrow(obs_counts), k + 1 == ncol(obs_counts))

  p11 <- q * (1 - f/2) + p * f / 2  # probability of a true 1 reported as 1
  p01 <- p * (1 - f/2) + q * f / 2  # probability of a true 0 reported as 1

  p2 <- p11 - p01  # == (1 - f) * (q - p)

  # When m = 1, obs_counts does not have the right dimensions. Fixing this.
  dim(obs_counts) <- c(m, k + 1)

  ests <- apply(obs_counts, 1, function(cohort_row) {
      N <- cohort_row[1]  # sample size for the cohort -- first column is total
      v <- cohort_row[-1] # counts for individual bits
      (v - p01 * N) / p2  # unbiased estimator for individual bits'
                          # true counts. It can be negative or
                          # exceed the total.
    })

  # NOTE: When k == 1, rows of obs_counts have 2 entries.  Then cohort_row[-1]
  # is a singleton vector, and apply() returns a *vector*.  When rows have 3
  # entries, cohort_row[-1] is a vector of length 2 and apply() returns a
  # *matrix*.
  #
  # Fix this by explicitly setting dimensions.  NOTE: It's k x m, not m x k.
  dim(ests) <- c(k, m)

  total <- sum(obs_counts[,1])

  variances <- apply(obs_counts, 1, function(cohort_row) {
      N <- cohort_row[1]
      v <- cohort_row[-1]
      p_hats <- (v - p01 * N) / (N * p2)  # expectation of a true 1
      p_hats <- pmax(0, pmin(1, p_hats))  # clamp to [0,1]
      r <- p_hats * p11 + (1 - p_hats) * p01  # expectation of a reported 1
      N * r * (1 - r) / p2^2  # variance of the binomial
     })

  dim(variances) <- c(k, m)

  # Transform counts from absolute values to fractional, removing bias due to
  #      variability of reporting between cohorts.
  ests <- apply(ests, 1, function(x) x / obs_counts[,1])
  stds <- apply(variances^.5, 1, function(x) x / obs_counts[,1])

  # Some estimates may be set to infinity, e.g. if f=1. We want to account for
  # this possibility, and set the corresponding counts to 0.
  ests[abs(ests) == Inf] <- 0

  list(estimates = ests, stds = stds)
}

ComputePrivacyGuarantees <- function(params, N) {
  # Compute privacy parameters and guarantees.
  p <- params$p
  q <- params$q
  f <- params$f
  h <- params$h

  q2 <- .5 * f * (p + q) + (1 - f) * q
  p2 <- .5 * f * (p + q) + (1 - f) * p

  exp_e_one <- ((q2 * (1 - p2)) / (p2 * (1 - q2)))^h
  if (exp_e_one < 1) {
    exp_e_one <- 1 / exp_e_one
  }
  e_one <- log(exp_e_one)

  exp_e_inf <- ((1 - .5 * f) / (.5 * f))^(2 * h)
  e_inf <- log(exp_e_inf)

  privacy_names <- c("Effective p", "Effective q", "exp(e_1)",
                     "e_1", "exp(e_inf)", "e_inf")
  privacy_vals <- c(p2, q2, exp_e_one, e_one, exp_e_inf, e_inf)

  privacy <- data.frame(parameters = privacy_names,
                        values = privacy_vals)
  privacy
}

CheckDecodeInputs <- function(counts, map, params) {
  # Returns an error message, or NULL if there is no error.

  if (nrow(map) != (params$m * params$k)) {
    return(sprintf(
        "Map matrix has invalid dimensions: m * k = %d, nrow(map) = %d",
        params$m * params$k, nrow(map)))
  }

  if ((ncol(counts) - 1) != params$k) {
    return(sprintf(paste0(
        "Dimensions of counts file do not match: m = %d, k = %d, ",
        "nrow(counts) = %d, ncol(counts) = %d"), params$m, params$k,
        nrow(counts), ncol(counts)))

  }

  # numerically correct comparison
  if (isTRUE(all.equal((1 - params$f) * (params$p - params$q), 0))) {
    return("Information is lost. Cannot decode.")
  }

  return(NULL)  # no error
}

Decode <- function(counts, map, params, threshold, decision_func) {

  error_msg <- CheckDecodeInputs(counts, map, params)
  if (!is.null(error_msg)) {
    stop(error_msg)
  }

  k <- params$k
  p <- params$p
  q <- params$q
  f <- params$f
  h <- params$h
  m <- params$m

  S <- ncol(map)  # total number of candidates

  N <- sum(counts[, 1])

  filter_cohorts <- which(counts[, 1] != 0)  # exclude cohorts with zero reports

  # stretch cohorts to bits
  filter_bits <- as.vector(matrix(1:nrow(map), ncol = m)[,filter_cohorts, drop = FALSE])

  map_filtered <- map[filter_bits, , drop = FALSE]

  es <- EstimateBloomCounts(params, counts)

  estimates_stds_filtered <-
    list(estimates = es$estimates[filter_cohorts, , drop = FALSE],
         stds = es$stds[filter_cohorts, , drop = FALSE])

  answers <- t(apply(map_filtered, 2, function(all_cohorts) {
    split_cohorts <- split(all_cohorts, ceiling(seq_along(all_cohorts)/k))
    answer <- rep(FALSE, length(split_cohorts))
    temp <- rep(FALSE, length(split_cohorts))
    for (i in 1:length(split_cohorts)) {
      estimate_counts <- estimates_stds_filtered$estimates[i,]
      estimate_std <- estimates_stds_filtered$std[i,]
      cohort <- split_cohorts[i]
      answer[i] <- all(estimate_counts[unlist(cohort)] > threshold) # Todo Dan: show confidence using estimate_std
      temp[i] <- paste(floor((1 - pnorm(threshold, estimate_counts[unlist(cohort)], estimate_std[unlist(cohort)]))*100), collapse = " ")
    }
    answer <- append(paste(answer, temp, sep="\n"), decision_func(answer))
    answer
  }))
  answers <- cbind(colnames(map), data.frame(answers))
  colnames(answers) <- c("Site", paste("Cohort ", filter_cohorts), "Decision")
  answers$Decision <- answers$Decision == "TRUE"
  num_detected <- sum(answers$Decision)

  # Compute summary of the fit.
  parameters <-
      c("Candidate strings", "Detected strings",
        "Sample size (N)")
  values <- c(S, num_detected, N)

  res_summary <- data.frame(parameters = parameters, values = values)

  privacy <- ComputePrivacyGuarantees(params, N)
  params <- data.frame(parameters =
                       c("k", "h", "m", "p", "q", "f", "N"),
                       values = c(k, h, m, p, q, f, N))

  # This is a list of decode stats in a better format than 'summary'.
  metrics <- list(sample_size = N,
                  num_detected = num_detected)

  list(fit = answers, found = answers$Site[answers$Decision], summary = res_summary, privacy = privacy, params = params,
       lasso = NULL, residual = NULL,
       counts = counts[, -1], resid = NULL, metrics = metrics,
       ests = es$estimates  # ests needed by Shiny rappor-sim app      
  )
}

ComputeCounts <- function(reports, cohorts, params) {
  # Counts the number of times each bit in the Bloom filters was set for
  #     each cohort.
  #
  # Args:
  #   reports: A list of N elements, each containing the
  #       report for a given report
  #   cohorts: A list of N elements, each containing the
  #       cohort number for a given report
  #   params: A list of parameters for the problem
  #
  # Returns:
  #   An mx(k+1) array containing the number of times each bit was set
  #       in each cohort.

  # Check that the cohorts are evenly assigned. We assume that if there
  #     are m cohorts, each cohort should have approximately N/m reports.
  #     The constraint we impose here simply says that cohort bins should
  #     each have within N/m reports of one another. Since the most popular
  #     cohort is expected to have about O(logN/loglogN) reports (which we )
  #     approximate as O(logN) bins for practical values of N, a discrepancy of
  #     O(N) bins seems significant enough to alter expected behavior. This
  #     threshold can be changed to be more sensitive if desired.
  N <- length(reports)
  cohort_freqs <- table(factor(cohorts, levels = 1:params$m))
  imbalance_threshold <- N / params$m
  if ((max(cohort_freqs) - min(cohort_freqs)) > imbalance_threshold) {
    cat("\nNote: You are using unbalanced cohort assignments, which can",
        "significantly degrade estimation quality!\n\n")
  }

  # Count the times each bit was set, and add cohort counts to first column
  counts <- lapply(1:params$m, function(i)
                   Reduce("+", reports[which(cohorts == i)]))
  counts[which(cohort_freqs == 0)] <- data.frame(rep(0, params$k))
  cbind(cohort_freqs, do.call("rbind", counts))
}
