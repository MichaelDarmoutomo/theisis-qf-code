#################################
## -- Load all dependencies -- ##
#################################

rm(list=ls())

if (Sys.info()['sysname'] != 'Linux') {
  if (!('code' %in% getwd())) {
    setwd("C:/Users/Michael.DESKTOP-575L4M7/OneDrive - Erasmus University Rotterdam/Documents/thesis-qf/code")
  }
}


source('optimizer.R')
source('kalman.R')
source('utils/loglikelihood.R')
source('utils/parameters.R')
source('utils/compute_se.R')

# Install dependencies
if("expm" %in% rownames(installed.packages()) == FALSE) {install.packages("expm")}
if("config" %in% rownames(installed.packages()) == FALSE) {install.packages("config")}
if("maxLik" %in% rownames(installed.packages()) == FALSE) {install.packages("maxLik")}
if("trustOptim" %in% rownames(installed.packages()) == FALSE) {install.packages("trustOptim")}
if("numDeriv" %in% rownames(installed.packages()) == FALSE) {install.packages("numDeriv")}

library(expm)
library(readxl)
library(config)
library(maxLik)
library(trustOptim)
library(numDeriv)

Sys.setenv(R_CONFIG_ACTIVE = "default")


config <- config::get()
config_mode = attributes(config)$config

######################
## -- Load data --  ##
######################

load(config$dataset)


if (config$dataset == "../data/data.Rda") {
  # # Get subset of data that is available
  data = DT[410:622, c('X1', 'X5', 'X10', 'X15', 'X20', 'X30', 'HICP', 'MSCI')]
  data$HICP = log(data$HICP)
  data$MSCI = log(data$MSCI)
  data[,1:6] = data[,1:6]/100
  
  data = t(as.matrix(data))
  row.names(data) <- NULL
  
  data[7,] = data[7,] - data[7,1]
  data[8,] = data[8,] - data[8,1]  
} else if (config_mode != "simulate") {
  data = data[325:579, 2:9]
  
  data = t(as.matrix(data))
  row.names(data) <- NULL
  
  # Shift so that HICP and MSCI start with zeros
  data[7,] = data[7,] - data[7,1]
  data[8,] = data[8,] - data[8,1]
  
}

#########################
## -- Kalman filter -- ##
#########################

if (config_mode =="default") {
  res = kalman_optimizer(data)
  param = res$solution
  
  if (config$save_results) {
    save(param, file="results/Parameters.Rdata")
  }
} else if (config_mode == "simulate"){
  for (i in 1:dim(sim_data)[2]) {
    res = kalman_optimizer(sim_data[,i,])
    save(res, file=paste0("results/simulation/", format(Sys.time(), "%Y%m%d_%H%M%S"), ".Rda"))
  }
} else {
  load(config$parameters_path)
  
  # Compute SE
  Hess = maxLik::numericHessian(function(p) loss(p, data), t0=param)
  negInvHess = -solve(Hess)
  diag(negInvHess)
}

