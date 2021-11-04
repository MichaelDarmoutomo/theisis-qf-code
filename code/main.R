#################################
## -- Load all dependencies -- ##
#################################

rm(list=ls())
setwd("~/thesis-qf/code")

source('optimizer.R')
source('kalman.R')
source('utils/loglikelihood.R')
source('utils/parameters.R')

library(ggplot2)
library(expm)

######################
## -- Load data --  ##
######################

load('../data/data.Rda')

# Get subset of data that is available
data = DT[409:622,1:63]
data$MSCI = log(data$MSCI)
data$HICP = log(data$HICP)
data[,4:63] = data[,4:63]/100

data = t(as.matrix(cbind(data[,4:63], data[,2:3])))
# Plot data
# ggplot(data, aes(x=Date)) + 
#   geom_line(aes(x=Date, y=MSCI)) +
#   ggtitle('MSCI Stock Index')
# 
# ggplot(data, aes(x=Date)) + 
#   geom_line(aes(x=Date, y=HICP)) +
#   ggtitle('HICP') 


#########################
## -- Kalman filter -- ##
#########################

res = kalman_optimizer(data)
