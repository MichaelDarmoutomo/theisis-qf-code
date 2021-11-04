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
data = DT[409:622, c('X1', 'X5', 'X10', 'X15', 'X20', 'X30', 'HICP', 'MSCI')]
data$MSCI = log(data$MSCI)
data$HICP = log(data$HICP)
data[,1:6] = data[,1:6]/100

data = t(as.matrix(data))
row.names(data) <- NULL

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
