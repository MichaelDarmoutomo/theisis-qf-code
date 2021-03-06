rm(list=ls())
source('utils/data_functions.R')
library(tidyr)
library(dplyr)
library(zoo)
library(data.table)
library(lubridate)

## -- Load and prepare HICP -- ##
file = 'HICP.csv'
HICP = load_data(file)

# Create date
HICP = HICP[c(1,2,3)]
HICP$Year = as.numeric(substr(HICP$V1, 1, 4))
HICP$Month = substr(HICP$V1, 5, 8)
HICP$Month = (match(HICP$Month,month.abb))
HICP$Date = sprintf("%04d-%02d", HICP$Year, HICP$Month)

HICP_ = HICP[,c('Date', 'V2')]
colnames(HICP_) = c('Date', 'HICP')

# Reverse rows
HICP_ = HICP_[(nrow(HICP):1),]

HICP_DF = transform(HICP_, Date = as.Date(as.yearmon(Date)))
HICP_DT = setDT(HICP_DF)

## -- Load and prepare HICP -- ##
file = 'MSCI.xlsx'
MSCI = load_data(file, skip=6)

MSCI = MSCI[1:622,]
MSCI$Year = as.numeric(substr(MSCI$Date, 9, 13))
MSCI$Month = substr(MSCI$Date, 1, 3)
MSCI$Month = (match(MSCI$Month,month.abb))
MSCI$Date = sprintf("%04d-%02d", MSCI$Year, MSCI$Month)

MSCI = transform(MSCI, Date = as.Date(as.yearmon(Date)))
colnames(MSCI) = c('Date', 'MSCI', 'Year', 'Month')
MSCI$MSCI = as.numeric(gsub(",","",MSCI$MSCI))

MSCI_DT = setDT(MSCI[,c(1:2)])

## -- Load and prepare yields -- ##
# Yields 1 (2001 - 2014)
file = 'yields1.xlsx'
yields1 = load_data(file)
yields1$Date = as.character(yields1$Date)
yields1$Date = substr(yields1$Date, 1, 7)
yields1 = transform(yields1, Date = as.Date(as.yearmon(Date)))
Yields1_DT = setDT(yields1)

# Yields 2 (2015- 2021)
file = 'yields2.xlsx'
yields2 = load_data(file)

# Long -> wide
yields2 = reshape(setDT(yields2), idvar = "Periode", timevar = "Looptijd in jaren", direction = "wide")
colnames(yields2) = gsub("waarde.", "X", colnames(yields2))
names(yields2)[names(yields2) == 'Periode'] = "Date"
yields2$Date = substr(yields2$Date, 1, 7)
yields2 = transform(yields2, Date = as.Date(as.yearmon(Date)))
Yields2_DT = setDT(yields2)

Yields_DT = full_join(Yields1_DT, Yields2_DT)
  
## -- MERGE -- ##
DT = full_join(MSCI_DT, HICP_DT, by='Date')
DT = full_join(DT, Yields_DT, by='Date')

# Remove unnessary variables
rm(list=setdiff(ls(), "DT"))
save(DT, file="../data/data.Rda")
