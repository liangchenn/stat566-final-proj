#' README:
#' -------
#' - author: ...
#' - date: 2026-05-25
#' 
#' Desc:
#' -------
#' This file loads the RHC dataset from the ATbounds as the main source of data.
#' 
#' We use Hmisc's RHC data as a source for correction on outcome variable `survival`.
#' We also redefine the death variable to aligned with the mortality rate outcome
#' in the original paper.
#' 
#' Output
#' --------
#' - data/raw/rhc-raw-data.csv

.FILE_NAME <- "01-load-raw-data.R"
cat(sprintf("...Running %s ...", .FILE_NAME))

# Packages ------------------------------------------------------------------------------------

library(data.table)
library(ATbounds) # main data source
library(Hmisc) # secondary data source



# Main ----------------------------------------------------------------------------------------

rhc_data <- ATbounds::RHC |> as.data.table()


#' NOTE: 
#' the survival outcome in ATbounds::RHC is described as 30-day survival indicator,
#' while it is actually the survival indicator within 180-day.
#' We use another data source to recover true 30-day survival/death outcome.

# load secondary data
Hmisc::getHdata("rhc")
setDT(rhc)

# create 30d survival outcome
rhc[, death30d := fifelse(dth30 == "Yes", 1, 0)]
rhc[, survival30d := fifelse(dth30 == "No", 1, 0)]

rhc[, death180d := fifelse(death == "Yes", 1, 0)]
rhc[, survival180d := fifelse(death == "No", 1, 0)]


# append corrected column to the RHC dataset
rhc_data[, death30d := rhc$death30d]
rhc_data[, survival30d := rhc$survival30d]

rhc_data[, death180d := rhc$death180d]
rhc_data[, survival180d := rhc$survival180d]

# save raw data
fwrite(rhc_data, "./data/raw/rhc-raw-data.csv")

cat("...done\n")

# Verifications -------------------------------------------------------------------------------

# # 1. check if RHC data's survival is actually 180d survival
# all(rhc_data$survival == rhc$survival180d)

rm(list = ls()); gc()

# # 2. further refer to original data in "https://hbiostat.org/data/repo/rhc.csv"

