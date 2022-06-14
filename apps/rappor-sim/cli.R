source("../../analysis/R/decode.R")
source("../../analysis/R/simulation.R")
source("../../analysis/R/encode.R")

main <- function(argv) {
  maps_file <- argv[[1]]
  params_file <- argv[[2]]
  pop_params_file <- argv[[3]]
  decode_file <- argv[[4]]
  params <- read.csv(params_file)

  if (file.exists(maps_file)) {
    load(file=maps_file)
  } else {
    pop_params <- read.csv(pop_params_file)
    N <- pop_params$N
    prop_missing <- pop_params$prop_missing
    maps <- GenerateMaps(N, params, pop_params, prop_missing)
    save(maps, file=maps_file)
  }

  threshold <- read.csv(decode_file)$threshold
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

    print(fits$hsts$summary)
}

main(commandArgs(TRUE))