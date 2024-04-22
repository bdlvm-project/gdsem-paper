This repo holds the code for the simulation study and the case study used in the paper **Gaussian distributional structural equation models: A framework for modeling latent heteroscedasticity**.

Recommended steps:

```r
# 1. Use renv to restore the appropriate package versions to the project environment
# install.packages("renv")
renv::restore()
# Some packages (bdlvm, cmdstanr and SBC) may have to be installed manually
# with devtools::install_github
# 2. Configure the targets package to run e.g. the simulation study:
targets::tar_config_set(
  script = "./simulation-study/_targets.R",
  store = "./simulation-study/_targets"
)
# 3. Run the targets pipeline
targets::tar_make()
# WARNING! The full simulation study takes several days to run and will produce
# many large files. Depending on the computer, it may also exhaust RAM. Check
# the simulation_conditions.R script to modify batch size and total replicates.
```
