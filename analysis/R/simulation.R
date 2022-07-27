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
# RAPPOR simulation library. Contains code for encoding simulated data and
#     creating the map used to encode and decode reports.

library(glmnet)
library(parallel)  # mclapply

SetOfSites <- function(num_sites = 100, proportion_https = 0.7, proportion_hsts_of_https = 0.3) {
  # Generates a set of strings for simulation purposes.
  # urls <- paste0("V_", as.character(1:num_sites), sample(c(".com", ".co.uk", ".fr"), num_sites, replace=TRUE))
  urls <- read.csv('./majestic_million.csv',nrows=num_sites)$Domain
  https <- as.logical(rbinom(n=num_sites, size=1, prob=proportion_https))
  sites <- data.frame(url=urls, https=https)
  sites$hsts <- sites$https & as.logical(rbinom(n=num_sites, size=1, prob=proportion_hsts_of_https))
  sites
}

GetSampleProbs <- function(pop_params) {
  # Generate different underlying distributions for simulations purposes.
  # Args:
  #    - pop_params: a list describing the shape of the true distribution:
  #              c(num_strings, prop_nonzero_strings, decay_type,
  #                rate_exponential).
  nsites <- pop_params[[1]]
  decay <- pop_params[[5]]
  expo <- pop_params[[6]]
  background <- pop_params[[7]]

  probs <- rep(0, nsites)
  ind <- floor(nsites)
  if (decay == "Measured") {
    temp <- read.csv('./majestic_million.csv',nrows=ind)$RefSubNets
    probs[1:ind] <- temp/sum(temp)
  }
  else if (decay == "Zipf") {
    temp <- 1.13/(1:ind) # todo Dan: maybe have a parameter, but 1 for now
    probs[1:ind] <- temp/sum(temp)
  }
  else if (decay == "Linear") {
    probs[1:ind] <- (ind:1) / sum(1:ind)
  } else if (decay == "Constant") {
    probs[1:ind] <- 1 / ind
  } else if (decay == "Exponential") {
    temp <- seq(0, 1, length.out = ind)
    temp <- exp(-temp * expo)
    temp <- temp + background
    temp <- temp / sum(temp)
    probs[1:ind] <- temp
  } else {
    stop('pop_params[[5]] must be in c("Linear", "Exponenential", "Constant", "Measured", "Zipf")')
  }
  probs
}

EncodeAll <- function(x, cohorts, map, params, num_cores = 1) {
  # Encodes the ground truth into RAPPOR reports.
  #
  # Args:
  #   x: Observed strings for each report, Nx1 vector
  #   cohort: Cohort assignment for each report, Nx1 vector
  #   map: list of matrices encoding locations of hashes for each
  #       string, for each cohort
  #   params: System parameters
  #
  # Returns:
  #   RAPPOR reports for each piece of data.

  p <- params$p
  q <- params$q
  f <- params$f
  k <- params$k

  qstar <- (1 - f / 2) * q + (f / 2) * p
  pstar <- (1 - f / 2) * p + (f / 2) * q

  candidates <- colnames(map[[1]])
  if (!all(x %in% candidates)) {
    stop("Some strings are not in the map. set(X) - set(candidates): ",
         paste(setdiff(unique(x), candidates), collapse=" "), "\n")
  }
  bfs <- mapply(function(x, y) y[, x], x, map[cohorts], SIMPLIFY = FALSE,
                USE.NAMES = FALSE)
  reports <- mclapply(bfs, function(x) {
    noise <- sample(0:1, k, replace = TRUE, prob = c(1 - pstar, pstar))
    ind <- which(x)
    noise[ind] <- sample(0:1, length(ind), replace = TRUE,
                         prob = c(1 - qstar, qstar))
    noise
  }, mc.cores = num_cores)

  reports
}

CreateMap <- function(strs, params, generate_pos = TRUE, basic = FALSE) {
  # Creates a list of 0/1 matrices corresponding to mapping between the strs and
  # Bloom filters for each instance of the RAPPOR.
  # Ex. for 3 strings, 2 instances, 1 hash function and Bloom filter of size 4,
  # the result could look this:
  # [[1]]
  #   1 0 0 0
  #   0 1 0 0
  #   0 0 0 1
  # [[2]]
  #   0 1 0 0
  #   0 0 0 1
  #   0 0 1 0
  #
  # Args:
  #    strs: a vector of strings
  #    params: a list of parameters in the following format:
  #         (k, h, m, p, q, f).
  #    generate_pos: Tells whether to generate an object storing the
  #        positions of the nonzeros in the matrix
  #    basic: Tells whether to use basic RAPPOR (only works if h=1).

  M <- length(strs) # number of individual reports
  map_by_cohort <- list()
  k <- params$k # size of the bloom filter instance
  h <- params$h # number of hash functions
  m <- params$m # number of cohorts
  for (i in 1:m) {
    if (basic && (h == 1) && (k == M)) {
      ones <- 1:M
    } else {
      ones <- sample(1:k, M * h, replace = TRUE)
    }
    cols <- rep(1:M, each = h)
    map_by_cohort[[i]] <- sparseMatrix(ones, cols, dims = c(k, M))
    colnames(map_by_cohort[[i]]) <- strs
  }

  all_cohorts_map <- do.call("rbind", map_by_cohort)
  if (generate_pos) {
    map_pos <- t(apply(all_cohorts_map, 2, function(x) {
      ind <- which(x == 1)
      n <- length(ind)
      if (n < h * m) {
        ind <- c(ind, rep(NA, h * m - n))
      }
      ind
    }))
  } else {
    map_pos <- NULL
  }

  list(map_by_cohort = map_by_cohort, all_cohorts_map = all_cohorts_map,
       map_pos = map_pos)
}

GetSample <- function(N, sites, probs) {
  # Sample for the sites population with distribution probs.
  sites[sample(nrow(sites), N, replace = TRUE, prob = probs),]
}

GetTrueBits <- function(samp, map, params) {
  # Convert sample generated by GetSample() to Bloom filters where mapping
  # is defined in map.
  # Output:
  #    - reports: a matrix of size [num_instances x size] where each row
  #               represents the number of times each bit in the Bloom filter
  #               was set for a particular instance.
  # Note: reports[, 1] contains the same size for each instance.

  N <- length(samp)
  k <- params$k
  m <- params$m
  strs <- colnames(map[[1]])
  reports <- matrix(0, m, k + 1)
  inst <- sample(1:m, N, replace = TRUE)
  for (i in 1:m) {
    tab <- table(samp[inst == i])
    tab2 <- rep(0, length(strs))
    tab2[match(names(tab), strs)] <- tab
    counts <- apply(map[[i]], 1, function(x) x * tab2)
    # cat(length(tab2), dim(map[[i]]), dim(counts), "\n")
    reports[i, ] <- c(sum(tab2), apply(counts, 2, sum))
  }
  reports
}

GetNoisyBits <- function(truth, params) {
  # Applies RAPPOR to the Bloom filters.
  # Args:
  #     - truth: a matrix generated by GetTrueBits().

  k <- params$k
  p <- params$p
  q <- params$q
  f <- params$f

  rappors <- apply(truth, 1, function(x) {
    # The following samples considering 4 cases:
    # 1. Signal and we lie on the bit.
    # 2. Signal and we tell the truth.
    # 3. Noise and we lie.
    # 4. Noise and we tell the truth.

    # Lies when signal sampled from the binomial distribution.
    lied_signal <- rbinom(k, x[-1], f)

    # Remaining must be the non-lying bits when signal. Sampled with q.
    truth_signal <- x[-1] - lied_signal

    # Lies when there is no signal which happens x[1] - x[-1] times.
    lied_nosignal <- rbinom(k, x[1] - x[-1], f)

    # Truth when there's no signal. These are sampled with p.
    truth_nosignal <- x[1] - x[-1] - lied_nosignal

    # Total lies and sampling lies with 50/50 for either p or q.
    lied <- lied_signal + lied_nosignal
    lied_p <- rbinom(k, lied, .5)
    lied_q <- lied - lied_p

    # Generating the report where sampling of either p or q occurs.
    rbinom(k, lied_q + truth_signal, q) + rbinom(k, lied_p + truth_nosignal, p)
  })

  cbind(truth[, 1], t(rappors))
}

GenerateMaps <- function(N = 10^5, params, pop_params, prop_missing = 0) {
  # Simulate N reports with pop_params describing the population and
  # params describing the RAPPOR configuration.
  # N - Number of samples
  num_sites <- pop_params[[1]]
  proportion_https <- pop_params[[2]]
  proportion_hsts_of_https <- pop_params[[3]]

  sites <- SetOfSites(num_sites, proportion_https, proportion_hsts_of_https) # Creates a list of sites with HSTS and HTTPS bits
  probs <- GetSampleProbs(pop_params) # creates a probability distribution for sampling sites

  samp <- GetSample(N, sites, probs) # Samples sites according to the distribution

  samp_hsts <- samp[samp$hsts == 1,]$url # Selects sampled sites that have hsts
  samp_nohttps <- samp[samp$https == 0,]$url # Selects sampled sites that have nohttps

  strs_hsts <- sites[sites$hsts == 1,]$url # Selects sites that have hsts
  strs_nohttps <- sites[sites$https == 0,]$url # Selects sites that are nohttps
  strs_https <- sites[sites$https == 1 & sites$hsts == 0,]$url # Selects sites that are https but not hsts


  # hsts
  map_hsts <- CreateMap(strs_hsts, params) # creates a random map of sites
  truth_hsts <- GetTrueBits(samp_hsts, map_hsts$map_by_cohort, params) # creates true bloom filter per cohort; first column counts total per row
  rappors_hsts <- GetNoisyBits(truth_hsts, params) # creates noisy bloom filter per cohort; first column counts total per row BEFORE NOISE

  # nohttps
  map_nohttps <- CreateMap(strs_nohttps, params)
  truth_nohttps <- GetTrueBits(samp_nohttps, map_nohttps$map_by_cohort, params)
  rappors_nohttps <- GetNoisyBits(truth_nohttps, params)

  # nohttps
  map_https <- CreateMap(strs_https, params)

  strs_hsts_apprx <- strs_hsts
  map_hsts_apprx <- map_hsts$all_cohorts_map
  # Remove % of strings to simulate missing variables.
  if (prop_missing > 0) {
    ind <- which(probs > 0)
    removed <- sample(ind, ceiling(prop_missing * length(ind)))
    map_hsts_apprx <- map_hsts$all_cohorts_map[, -removed]
    strs_hsts_apprx <- strs_hsts[-removed]
  }
  strs_nohttps_apprx <- strs_nohttps
  map_nohttps_apprx <- map_nohttps$all_cohorts_map
  # Remove % of strings to simulate missing variables.
  if (prop_missing > 0) {
    ind <- which(probs > 0)
    removed <- sample(ind, ceiling(prop_missing * length(ind)))
    map_nohttps_apprx <- map_nohttps$all_cohorts_map[, -removed]
    strs_nohttps_apprx <- strs_nohttps[-removed]
  }
  strs_https_apprx <- strs_https
  map_https_apprx <- map_https$all_cohorts_map
  # Remove % of strings to simulate missing variables.
  if (prop_missing > 0) {
    ind <- which(probs > 0)
    removed <- sample(ind, ceiling(prop_missing * length(ind)))
    map_https_apprx <- map_https$all_cohorts_map[, -removed]
    strs_https_apprx <- strs_https[-removed]
  }

  # Randomize the columns.
  ind <- sample(1:length(strs_hsts_apprx), length(strs_hsts_apprx))
  map_hsts_apprx <- map_hsts_apprx[, ind]

  ind <- sample(1:length(strs_nohttps_apprx), length(strs_nohttps_apprx))
  map_nohttps_apprx <- map_nohttps_apprx[, ind]

  ind <- sample(1:length(strs_https_apprx), length(strs_https_apprx))
  map_https_apprx <- map_https_apprx[, ind]

  list(sites = sites,
       probs = probs,
       strs_hsts=strs_hsts,
       strs_nohttps=strs_nohttps,
       strs_https=strs_https,
       strs_hsts_apprx = strs_hsts_apprx,
       strs_https_apprx = strs_https_apprx,
       strs_nohttps_apprx = strs_nohttps_apprx,
       truth_hsts=truth_hsts,
       truth_nohttps=truth_nohttps,
       map_hsts = map_hsts,
       map_nohttps = map_nohttps,
       map_https = map_https,
       map_hsts_apprx=map_hsts_apprx,
       map_nohttps_apprx=map_nohttps_apprx,
       map_https_apprx=map_https_apprx,
       rappors_hsts=rappors_hsts,
       rappors_nohttps=rappors_nohttps
  )
}

DecisionFunctions <- function(params, func) {
  decision_func <- function(answer) {TRUE}
  if (func == "any") {
    decision_func <- function(answer) any(answer)
  } else if (func == "all") {
    decision_func <- function(answer) all(answer)
  } else if (func == "all") {
    decision_func <- function(answer) {sum(answer) >  params$m / 2}
  }
  decision_func
}

GenerateSamples <- function(params,
                            decoding_params,
                            sites,
                            probs,
                            strs_hsts, strs_nohttps, strs_https,
                            strs_hsts_apprx, strs_https_apprx, strs_nohttps_apprx,
                            truth_hsts, truth_nohttps,
                            map_hsts, map_nohttps, map_https,
                            map_hsts_apprx, map_nohttps_apprx, map_https_apprx,
                            rappors_hsts, rappors_nohttps) {
  threshold <- decoding_params[[1]]
  primary_decision <- DecisionFunctions(params, decoding_params[[2]])
  secondary_decision <- DecisionFunctions(params, decoding_params[[3]])
  print(primary_decision)
  print(secondary_decision)

  # merge maps map_hsts + map_https_apprx + map_nohttps
  map_merged <- cbind(map_hsts_apprx, map_https_apprx, map_nohttps_apprx)
  # todo Dan: what is wrong with this?
  # ind <- sample(1:length(map_merged), length(map_merged))
  # map_merged <- map_merged[, ind]

  fit_hsts <- Decode(rappors_hsts, map_merged, params, threshold, decision_func = primary_decision)

  hsts_tp <-strs_hsts[na.omit(match(fit_hsts$found, strs_hsts))]
  hsts_fp <-strs_nohttps[na.omit(match(fit_hsts$found, strs_nohttps))]
  hsts_soft_fp <-strs_https[na.omit(match(fit_hsts$found, strs_https))]

  # Add truth column.
  fit_hsts$fit$Truth <- apply(fit_hsts$fit, 1, function(r) {r[1] %in% strs_hsts})

  fit_hsts$summary <- rbind(
    fit_hsts$summary[1:1,],
    c("HSTS strings", length(strs_hsts_apprx)),
    c("HTTPS strings", length(strs_https_apprx)),
    c("No HTTPS Strings", length(strs_nohttps_apprx)),
    fit_hsts$summary[2:2,],
    c("HSTS strings found (True positives) in B1", length(intersect(fit_hsts$found, strs_hsts))),
    c("HTTPS strings found (Soft false positives) in B1", length(intersect(fit_hsts$found, strs_https))),
    c("No HTTPS strings found (Bad false positives) in B1", length(intersect(fit_hsts$found, strs_nohttps))),
    fit_hsts$summary[3:nrow(fit_hsts$summary),]
  )

  fit_hsts$map <- map_hsts$map_by_cohort
  fit_hsts$truth <- truth_hsts
  fit_hsts$strs_full <- sites$url
  fit_hsts$strs <- strs_hsts
  fit_hsts$probs <- probs

  # Can we identify false positives as part of the rappors_nohttps?
  map_tp_apprx <- as(map_hsts_apprx[, hsts_tp, drop=FALSE], "sparseMatrix")
  map_softfp_apprx <- as(map_https_apprx[, hsts_soft_fp, drop=FALSE], "sparseMatrix")
  map_fp_apprx <- as(map_nohttps_apprx[, hsts_fp, drop=FALSE], "sparseMatrix")
  map_merged_positive <- cbind(map_tp_apprx, map_softfp_apprx, map_fp_apprx)
  fit_nohttps <- NULL
  if(ncol(map_merged_positive) != 0) {
    fit_nohttps <- Decode(rappors_nohttps, map_merged_positive, params, threshold, decision_func = secondary_decision)
    fit_nohttps$fit$Truth <- apply(fit_nohttps$fit, 1, function(r) {r[1] %in% strs_nohttps})

    fit_nohttps$map <- map_nohttps$map_by_cohort
    fit_nohttps$truth <- truth_nohttps
    fit_nohttps$strs_nohttps <- strs_nohttps

    fit_hsts$summary <- rbind(
      fit_hsts$summary,
      c("HSTS strings found in B2 (no benefit)", length(intersect(fit_nohttps$found, strs_hsts))),
      c("HTTPS strings found in B2", length(intersect(fit_nohttps$found, strs_https))),
      c("No HTTPS strings not found (disaster) in B2", length(hsts_fp) - length(intersect(fit_nohttps$found, strs_nohttps))),
      c("No HTTPS strings found (disaster averted) in B2", length(intersect(fit_nohttps$found, strs_nohttps))),
      c("Final benefit", length(hsts_tp) - length(intersect(fit_nohttps$found, strs_hsts)))
    )
  }
  fits <- list("hsts" = fit_hsts, "nohttps" = fit_nohttps, "sites" = sites)
  fits
}
# params <- list(k = 128, # size of the bloom filter instance
#          h = 2, # number of hash functions
#          m = 8, # number of cohorts
#          p = 0.5,
#          q = 0.75,
#          f = 0.5)
# pop_params <- list(300, #nsites
#       0.7, #nhttps
#       0.3, #nhsts
#       0.5, #nonzero
#       "Exponential", #decay
#       10, #expo
#       0.05 #background
#       )