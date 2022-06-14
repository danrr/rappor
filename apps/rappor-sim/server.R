library(shiny)
source("../../analysis/R/decode.R")
source("../../analysis/R/simulation.R")
source("../../analysis/R/encode.R")

Plot <- function(x, threshold=NULL, color = "grey") {
  n <- nrow(x)
  if (n < 16) {
    par(mfrow = c(n, 1), mai = c(0, .5, .5, 0))
  } else if (n < 64) {
    par(mfrow = c(n / 2, 2), mai = c(0, .5, .5, 0))
  } else {
    par(mfrow = c(n / 4, 4), mai = c(0, .5, .5, 0))
  }
  for (i in 1:nrow(x)) {
    barplot(x[i, ], main = paste0("Cohort ", i), col = color, border = color)
    if (!is.null(threshold)) {
      abline(h=0)
      abline(h=threshold, lty="dashed")
    }
  }
}

shinyServer(function(input, output) {
  # Example state global variable.
  es <- list()

  # Example buttons states.
  ebs <- rep(0, 3)

  Params <- reactive({
    list(k = as.numeric(input$size),
         h = as.numeric(input$hashes),
         m = as.numeric(input$instances),
         p = as.numeric(input$p),
         q = as.numeric(input$q),
         f = as.numeric(input$f))
  })

  PopParams <- reactive({
    list(as.numeric(input$nsites),
      as.numeric(input$nhttps),
      as.numeric(input$nhsts),
      as.numeric(input$nonzero),
      input$decay,
      as.numeric(input$expo),
      as.numeric(input$background)
      )
  })

  DecodingParams <- reactive({
    as.numeric(input$threshold)
  })

  Sample <- reactive({
    Sample2()$hsts
  })

  Maps <- reactive({
      input$sample
      N <- input$N
      params <- Params()
      pop_params <- PopParams()
      prop_missing <- input$missing
      GenerateMaps(N, params, pop_params, prop_missing = prop_missing)
  })

  Sample2 <- reactive({
    maps <- Maps()
    params <- Params()
    threshold <- DecodingParams()
    fits <- GenerateSamples(
      params,
      sites = maps$sites,
      probs = maps$probs,
      strs_hsts = maps$strs_hsts,
      strs_nohttps = maps$strs_nohttps,
      strs_https = maps$strs_https,
      strs_hsts_apprx = maps$strs_hsts_apprx,
      strs_https_apprx = maps$strs_https_apprx,
      strs_nohttps_apprx = maps$strs_nohttps_apprx,
      truth_hsts = maps$truth_hsts,
      truth_nohttps = maps$truth_nohttps,
      map_hsts = maps$map_hsts,
      map_nohttps = maps$map_nohttps,
      map_https = maps$map_https,
      map_hsts_apprx = maps$map_hsts_apprx,
      map_nohttps_apprx = maps$map_nohttps_apprx,
      map_https_apprx = maps$map_https_apprx,
      rappors_hsts = maps$rappors_hsts,
      rappors_nohttps = maps$rappors_nohttps,
      threshold = threshold)
    fits
  })

  # Results summary.
  output$pr <- renderTable({
    Sample2()$hsts$summary
  },
  include.rownames = FALSE,
  include.colnames = FALSE)

  output$pr2 <- renderTable({
    Sample2()$nohttps$summary
  },
  include.rownames = FALSE,
  include.colnames = FALSE)

  # Results table.
  output$tab <- renderDataTable({
     Sample2()$hsts$fit
   },
  options = list(iDisplayLength = 100))

  output$tab2 <- renderDataTable({
     Sample2()$nohttps$fit
   },
  options = list(iDisplayLength = 100))
  # Epsilon.
  output$epsilon <- renderTable({
    Sample()$privacy
  },
  include.rownames = FALSE,
  include.colnames = FALSE,
  digits = 4)

  # True distribution.
  output$probs <- renderPlot({
    samp <- Sample2()
    PlotPopulation(samp)
  })

  # True bits patterns.
  output$truth <- renderPlot({
    truth <- Sample()$truth
    Plot(truth[, -1, drop = FALSE], threshold = NULL, color = "darkblue")
  })

  # Lasso plot.
  output$lasso <- renderPlot({
    fit <- Sample()$lasso
    if (!is.null(fit)) {
      plot(fit)
    }
  })

  output$resid <- renderPlot({
    resid <- Sample()$residual
    params <- Params()
    plot(resid, xlab = "Bloom filter bits", ylab = "Residuals")
    abline(h = c(-1.96, 1.96), lty = 2, col = 2)
    sq <- qnorm(.025 / length(resid))
    abline(h = c(sq, -sq), lty = 2, col = 3, lwd = 2)
    abline(h = c(-3, 3), lty = 2, col = 4, lwd = 2)
    abline(v = params$k * (0:params$m), lty = 2, col = "blue")
    legend("topright", legend = paste0("SD = ", round(sd(resid), 2)), bty = "n")
  })

  # Estimated bits patterns.
  output$ests <- renderPlot({
    ests <- Sample()$ests
    threshold <- DecodingParams()
    Plot(ests, threshold=threshold, color = "darkred")
  })

  # Estimated vs truth.
  output$ests_truth <- renderPlot({
    plot(unlist(Sample()$ests), unlist(Sample()$truth[, -1]),
         xlab = "Estimates", ylab = "Truth", pch = 19)
    abline(0, 1, lwd = 4, col = "darkred")
  })

  output$example <- renderPlot({
    params <- Params()
    strs <- Sample()$strs
    map <- Sample()$map
    samp <- Sample()

    # First run on app start.
    value <- sample(strs, 1)
    res <- Encode(value, map, strs, params, N = input$N)

    if (input$new_user > ebs[1]) {
      res <- Encode(es$value, map, strs, params, N = input$N)
      ebs[1] <<- input$new_user
    } else if (input$new_value > ebs[2]) {
      res <- Encode(value, map, strs, params, cohort = es$cohort, id = es$id,
                    N = input$N)
      ebs[2] <<- input$new_value
    } else if (input$new_report > ebs[3]) {
      res <- Encode(es$value, map, strs, params, B = es$B,
                    BP = es$BP, cohort = es$cohort, id = es$id, N = input$N)
      ebs[3] <<- input$new_report
    }
    es <<- res
    ExamplePlot(res, params$k, c(ebs, input$new_user, input$new_value, input$new_report))
  })

})
