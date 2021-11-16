#################################
## -- Load all dependencies -- ##
#################################

rm(list=ls())
if (!grep('code', getwd())) {
  setwd("~/thesis-qf/code")
}

source('optimizer.R')
source('kalman.R')
source('utils/loglikelihood.R')
source('utils/parameters.R')

if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}
if("expm" %in% rownames(installed.packages()) == FALSE) {install.packages("expm")}

library(ggplot2)
library(expm)

######################
## -- Load data --  ##
######################

load('../data/data.Rda')

# Get subset of data that is available
data = DT[410:622, c('X1', 'X5', 'X10', 'X15', 'X20', 'X30', 'HICP', 'MSCI')]
data$MSCI = log(data$MSCI)
data$HICP = log(data$HICP)
data[,1:6] = data[,1:6]/100

data = t(as.matrix(data))
row.names(data) <- NULL

ggplot(subset(DT, as.numeric(format.Date(DT$Date, "%Y")) > 2000), aes(x=Date)) +
  geom_line(aes(x=Date, y=MSCI)) +
  ggtitle('MSCI Stock Index')

ggplot(DT, aes(x=Date)) +
  geom_line(aes(x=Date, y=HICP)) +
  ggtitle('HICP')

# Plot data



#########################
## -- Kalman filter -- ##
#########################

res = kalman_optimizer(data)
