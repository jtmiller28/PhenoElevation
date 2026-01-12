### Title: Calc Phenometrics
### Author: Rob Guralnick
### Editor: JT Miller 
### Date: 12-29-2025

# Purpose: load in prepped data, calculate onset, offset, and duration using Phenesse
setwd("/blue/guralnick/millerjared/PhenoElevation/") # set working-directory 

## Load Libraries
library(sf)
library(dplyr)
library(data.table)
library(phenesse)

## Read in phenometric ready data (see prep-data script)
phenologycv_sf14 <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-ready-data.rds")

## Calc the onset, offset, and median for each species-year-elevbin-hexcell combo
## require at least 7 obs and 3 distinct doys per combination
# Calc OnSet
phenologycv_sf14_onset <-  phenologycv_sf14 %>% 
  group_by(scientific_name,year,bins6,id_cells) %>% 
  filter(n() >= 7) %>%
  filter(n_distinct(doy) > 3) %>%
  group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$doy, percentile = 0.05, bootstraps=250)))

# Calc OffSet
phenologycv_sf14_offset <-  phenologycv_sf14 %>% 
  group_by(scientific_name,year,bins6,id_cells) %>% 
  filter(n() >= 7) %>%
  filter(n_distinct(doy) > 3) %>%
  group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$doy, percentile = 0.95, bootstraps=250)))

# Calc the Median
phenologycv_sf14_median <-  phenologycv_sf14 %>% 
  group_by(scientific_name,year,bins6,id_cells) %>% 
  filter(n() >= 7) %>%
  filter(n_distinct(doy) > 3) %>%
  group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$doy, percentile = 0.50, bootstraps=250)))

## Write out data as an RDS object to be use in sp traits and post pheno cleaning
saveRDS(phenologycv_sf14_onset, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-onset-calc.rds")
saveRDS(phenologycv_sf14_offset, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-offset-calc.rds")
saveRDS(phenologycv_sf14_median, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-median-calc.rds")
