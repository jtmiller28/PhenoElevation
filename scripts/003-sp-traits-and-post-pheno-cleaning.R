### Title: Sp Traits and Post Pheno Cleaning
### Author: Rob Guralnick
### Editor: JT Miller 
### Date: 12-29-2025

### Purpose: load in phenometric calcs, do some post-cleaning and formatting, then add species traits
setwd("/blue/guralnick/millerjared/PhenoElevation/") # set working-directory 

## Load Libraries
library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(data.table)
library(readr)

## Load in phenometric calcs
phenologycv_sf14_onset <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-onset-calc.rds")
phenologycv_sf14_offset <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-offset-calc.rds")
phenologycv_sf14_median <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-median-calc.rds")

## Load in other relevant diagnostics for phenometrics 
phenologycv_sf14_ndistinct <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-distinct-doy-combos.rds")
phenologycv_sf14_obscount <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-obs-combos.rds")

## Post phenometric data assembly steps
# remove nonvital info in OnSet, pivot wider to rotate table. 
phenologycv_sf14_onset_2 <- dplyr::select(phenologycv_sf14_onset, -c(sd, median, min, max, range, skew, trimmed,mad,kurtosis,se))
phenologycv_sf14_onset_3 <- pivot_wider(phenologycv_sf14_onset_2, names_from = column, values_from = mean)

# remove nonvital info in OffSet, pivot wider to rotate table.
phenologycv_sf14_offset_2 <- dplyr::select(phenologycv_sf14_offset, -c(sd, median, min, max, range, skew, trimmed,mad,kurtosis,se))
phenologycv_sf14_offset_3 <- pivot_wider(phenologycv_sf14_offset_2, names_from = column, values_from = mean)

# remove nonvital info the median table, pivot wider to rotate table.
phenologycv_sf14_50_2 <- dplyr::select(phenologycv_sf14_median, -c(sd, median, min, max, range, skew, trimmed,mad,kurtosis,se))
phenologycv_sf14_50_3 <- pivot_wider(phenologycv_sf14_50_2, names_from = column, values_from = mean)

# combine OnSet and OffSet tables for given species, year, spatial cell, and elev bin
phenologycv_sf14_onoff <- merge(phenologycv_sf14_onset_3,phenologycv_sf14_offset_3, by=c("year","bins6","scientific_name","id_cells"))

# combine OnSet-Offset table with the median phenophase table.  
phenologycv_sf14_onoff50 <- merge(phenologycv_sf14_onoff,phenologycv_sf14_50_3, by=c("year","bins6","scientific_name","id_cells"))

# add in number of distinct doys per species x year x elev bin x spatial hex cell
phenologycv_sf14_onoffall <- merge(phenologycv_sf14_onoff50,phenologycv_sf14_ndistinct, by=c("year","bins6","scientific_name","id_cells"))

# add in number of observations per species x year x elev bin x spatial hex cell 
phenologycv_sf14_onoffall2 <- merge(phenologycv_sf14_onoffall,phenologycv_sf14_obscount, by=c("year","bins6","scientific_name","id_cells"))

# rename fields for clarity
phenologycv_sf14_onoffall3 <- phenologycv_sf14_onoffall2 %>%
  dplyr::rename(mean=estimate,mean_low=low_ci, mean_high=high_ci, onset=estimate.x,onset_low=low_ci.x, onset_high=high_ci.x, offset =estimate.y,offset_low=low_ci.y, offset_high=high_ci.y)

## Create centroids for hexcells, join on data
# load in hexed study region
hex_grid_w_ids <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/full-study-region.rds")

# use centroids of hexcells to delimit avg lon-lat values per grouping, join on data
grid_cent <- st_centroid(hex_grid_w_ids)
grid_cent_ll <- st_transform(grid_cent,crs("EPSG:4326"))
grid_cent_ll_point <- grid_cent_ll %>%
  mutate(longitude = st_coordinates(.)[,1],latitude = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  as.data.frame()

# further restrict study region to be less than 60 and more than 25 latitude
grid_cent_ll_point <-  grid_cent_ll_point %>% filter(latitude < 60)
grid_cent_ll_point <- grid_cent_ll_point %>% filter(latitude > 25)

# save centroids for a later delimitation of subregions
saveRDS(grid_cent_ll_point, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/grid_cent_ll_point.rds")
phenologycv_sf14_onoffall4 <- left_join(phenologycv_sf14_onoffall3 , grid_cent_ll_point , by="id_cells")

## Create duration metric based on offset-onset
phenologycv_sf14_onoffall4 <- phenologycv_sf14_onoffall4 %>% 
  mutate(duration = offset-onset)

## scale onset, offset, mean, and duration into z scores. Remove outliers. 
phenologycv_sf14_onoffall5 <-  phenologycv_sf14_onoffall4 %>% 
  group_by(scientific_name) %>% 
  mutate(zonset = scale(onset)) %>% 
  filter(between(zonset,-3.25,+3.25)) # drops 50 

phenologycv_sf14_onoffall6 <-  phenologycv_sf14_onoffall5 %>% 
  group_by(scientific_name) %>% 
  mutate(zoffset = scale(offset)) %>% 
  filter(between(zoffset,-3.25,+3.25)) # drops 58

phenologycv_sf14_onoffall7 <-  phenologycv_sf14_onoffall6 %>% 
  group_by(scientific_name) %>% 
  mutate(zmean = scale(mean)) %>% 
  filter(between(zmean,-3.25,+3.25)) # drops 28

phenologycv_sf14_onoffall8 <-  phenologycv_sf14_onoffall7 %>% 
  group_by(scientific_name) %>% 
  mutate(zdur = scale(duration)) %>% 
  filter(between(zdur,-3.25,+3.25)) # drops 72

## Refilter at this point to remove species that have less than 6 values of phenometrics (any)
phenologycv_sf14_onoffall9 <- phenologycv_sf14_onoffall8 %>% group_by(scientific_name) %>% filter(n()>6) # drops 31

## Assemble Species Traits, Here we're just going to focus on annual/biannual vs perennial. Join traits to pheno data.
# read in trait data, filter for only taxa in current modeled dataset
SpeciesTraits3 <- read_csv("./data/SpeciesTraits3.csv")
species_in_dataset <- phenologycv_sf14_onoffall9$scientific_name
species_in_dataset2 <- gsub(" ", "_", species_in_dataset)
SpeciesTraits3 <- SpeciesTraits3 %>% filter(scientific_name %in% species_in_dataset2)
saveRDS(SpeciesTraits3, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/species-traits-for-filtered-taxa.rds")

# join traits to pheno data
SpeciesTraits3$scientific_name <- gsub("_", " ",SpeciesTraits3$scientific_name)
phenologycv_sf14_onoffall10 <- inner_join(phenologycv_sf14_onoffall9, SpeciesTraits3  , by="scientific_name")

# maintain variable names for clarity
phenologycv_sf14_onoffall11  <- phenologycv_sf14_onoffall10 

## ensure that there is more than 4 distinct sampled doys per species, remove records that have NA for binned elev data. 
phenologycv_sf14_onoffall11 <- phenologycv_sf14_onoffall11 %>% filter(ndistinct>4)
phenologycv_sf14_onoffall11 <- phenologycv_sf14_onoffall11 %>% drop_na(bins6)
phenologycv_sf14_onoffall12 <- as.data.frame(ungroup(phenologycv_sf14_onoffall11))

## Fix elevational bin naming, drop any records that lack elevation data, then find mean timing
phenologycv_sf14_onoffall13 <- phenologycv_sf14_onoffall12 %>% mutate(elev = dplyr::recode(bins6,  "0to.15" = 75, ".15to.3" = 225, ".3to.45" = 375, ".45to.6" = 525, ".6to.75" = 675, ".75to.9" = 825, ".9to1.05" = 975, "1.05to1.2" = 1125, "1.2to1.35" = 1275, "1.35to1.5" = 1425, "1.5to1.65" = 1575, "1.65to1.8" = 1725, "1.8to1.95" = 1875, "1.95to2.1" = 2025, "2.1to2.25" = 2175, "2.25to2.4" = 2325, "2.4to2.55" = 2475, "2.55to2.7" = 2625, "2.7to2.85" = 2775, "2.85to3.0" = 2925))
phenologycv_sf14_onoffall13 <- phenologycv_sf14_onoffall13 %>% drop_na(elev)

# get species mean timing
phenologycv_sf14_sptime <- phenologycv_sf14_onoffall13 %>% group_by(scientific_name) %>% summarize(spmean=mean(mean))
# save sp timebins for use in synchrony analysis 
saveRDS(phenologycv_sf14_sptime, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/species-avg-timings.rds")
# create timebins
phenologycv_sf14_sptime2 <- phenologycv_sf14_sptime %>%  mutate(timebin = ntile(spmean, n=2), timebin2 = ntile(spmean, n=3))
# save sp timebins for use in synchrony analysis 
saveRDS(phenologycv_sf14_sptime2, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/species-avg-timebins.rds")
phenologycv_sf14_onoffall14 <- merge(phenologycv_sf14_onoffall13, phenologycv_sf14_sptime2, by=c("scientific_name"))

## Filter out taxa that are outside our longitude of study region, and those with very short durations
# -95 is min lon 
phenologycv_sf14_onoffall14 <- phenologycv_sf14_onoffall14 %>% filter(longitude < -95)
# require at least 10 doys as duration
phenologycv_sf14_onoffall14 <- phenologycv_sf14_onoffall14 %>% filter(duration > 10)

## Assign traits to species + timebinned data, ensure data is 3000 or less
phenologycv_sf14_onoffall14$scientific_name <- gsub(" ", "_", phenologycv_sf14_onoffall14$scientific_name)
phenologycv_sf14_onoffall15 <- phenologycv_sf14_onoffall14 %>%  mutate(annpertime = paste(geoannper, timebin, sep = '_'))
phenologycv_sf14_onoffall15 <- phenologycv_sf14_onoffall15 %>% filter(elev <= 3000)

## Filter species to ensure there is adequate representation across elevation for these analyses
phenologycv_sf14_onoffall_elevf2 <- phenologycv_sf14_onoffall15 %>% 
  group_by(scientific_name) %>%
  summarize(elev_range = max(elev)-min(elev)) %>%
  filter(elev_range >= 255)

phenologycv_sf14_onoffall_elevf3 <- phenologycv_sf14_onoffall15 %>% 
  group_by(id_cells) %>%
  summarize(elev_range = max(elev)-min(elev)) %>%
  filter(elev_range >= 225)

## Only retain data that meets these thresholds
phenologycv_sf14_onoffall16 <- subset(phenologycv_sf14_onoffall15, scientific_name %in% phenologycv_sf14_onoffall_elevf2$scientific_name)
phenologycv_sf14_onoffall17 <- subset(phenologycv_sf14_onoffall16, id_cells %in% phenologycv_sf14_onoffall_elevf3$id_cells)

## Remove species with low numbers of phenometrics
phenologycv_sf14_onoffall19 <- phenologycv_sf14_onoffall17 %>% 
  select(-notes, -`...11`, -`...12`, -`...13`, -`...14`)
phenologycv_sf14_onoffall20 <- phenologycv_sf14_onoffall19 %>% group_by(scientific_name) %>% filter(n()>=5) %>% ungroup() # make sure to ungroup, otherwise it ruins the rest

## Save as an object to go to run pglmms and synchrony analyses
saveRDS(phenologycv_sf14_onoffall20, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-model-ready-data.rds")