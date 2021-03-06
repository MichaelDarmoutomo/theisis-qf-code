library(config)
rm(list=ls())
Sys.setenv(R_CONFIG_ACTIVE = "pension_fund")
config <- config::get()

source("pensionfund/pension_fund.R")
source("pensionfund/adjustment_factor.R")
source("simulation.R")


# Simulate the economy
if (config$load_economy) {
  load(config$save_path)
} else {
  e = simulate_economy(
    paramfile = config$parameters_path, 
    dt = config$dt, 
    T = config$T, 
    nSim = config$nSim, 
    w = config$w, 
    maturities = config$maturities,
    parallel = config$parallel,
    save_economy = config$save_economy,
    save_path = config$save_path)
}

x = c(0.4740, 0.7363)
# x = c(0.4740, 0)
start_time = Sys.time()
print(paste("Starting at",start_time))

cec = optim(x, function(z) pension_fund(z,e,config$nSim), 
            control=list(
              fnscale=-1,
              trace=1
            ))

print(paste("Stopping at",Sys.time()))
print(paste("Total runtime ", Sys.time() - start_time))
