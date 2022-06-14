CLI usage:
```shell
RAPPOR_REPO=../../ Rscript cli.R <maps>.Rdata params.csv pop_params.csv decode.csv
```

If maps file exists, it will use it. If it doesn't, it will generate a new set of maps based on `pop_params.csv` and 
save it to file handle passed in.

