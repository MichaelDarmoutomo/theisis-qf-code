#################################
## -- Load all dependencies -- ##
#################################

source('kalman_filter.R')
source('kalman.R')
source('utils/loglikelihood.R')
library(ggplot2)


######################
## -- Load data --  ##
######################

load('../data/data.Rda')

# Get subset of data that is available
data = DT[409:622,1:63]
data$MSCI = log(data$MSCI)
data$HICP = log(data$HICP)

# Plot data
ggplot(data, aes(x=Date)) + 
  geom_line(aes(x=Date, y=MSCI)) +
  ggtitle('MSCI Stock Index')

ggplot(data, aes(x=Date)) + 
  geom_line(aes(x=Date, y=HICP)) +
  ggtitle('HICP') 


#########################
## -- Kalman filter -- ##
#########################

res = kalman_optimizer(data)
