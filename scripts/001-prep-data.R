### Title: Prep data
### Author: Rob Guralnick
### Editor: JT Miller 
### Date: 12-29-2025

# Purpose: load in pheno-phase annotated data for flowering plants, clean to just flowering pheno-phase, delimit to Western North America, clean and ready data for phenometric calcs
setwd("/blue/guralnick/millerjared/PhenoElevation/") # set working-directory 

## Load Libraries
library(arrow)
library(rnaturalearth)
library(sf)
library(dplyr)
library(data.table)
library(ggplot2)
library(readr)
library(terra)

## Set-up phenometric data
# read in elevation data per species
elev_species2 <- readr::read_csv("./data/elev_species2.csv")
elev_species2$species <- gsub("_", " ", elev_species2$species)

# process through large phenovision dataset
phenologycv_csv <- arrow::open_dataset(
  sources = "/blue/guralnick/millerjared/PhenoElevation/data/phenobase-annotations-headers-added-01-2025.csv", 
  col_types = schema(ISBN = string()),
  format = "csv")
# make parquet delimited by year for fast processing
phenologycv_csv |>
  group_by(year) |>
  arrow::write_dataset(path = "./outputs/parquets/phenologycv_presence_data", format = "parquet")

phenologycv_pd <- open_dataset("./outputs/parquets/phenologycv_presence_data")

# retain only data that is 2017 - 2023
phenologycv_pd2 <- phenologycv_pd  |> 
  filter(year >= 2017, year != 2024, trait == "flower", scientific_name %in% elev_species2$species)

# collect into memory
phenologycv_pds3 <- phenologycv_pd2 |> collect()

# remove duplicate records and only retain species with 500 or more unique records
phenologycv_pd4 <- phenologycv_pds3  %>% 
  group_by(scientific_name) %>% 
  distinct(observed_metadata_url, .keep_all = TRUE) %>% 
  filter(n() >= 500)

## Build Study Area
# convert occurrences to spatial 
crs_dd <- 4326
phenologycv_sf5  <- sf::st_as_sf(phenologycv_pd4, coords = c("longitude", "latitude"), crs = crs_dd)

# use naturalearth's country shapes to delimit North America
country_sf <- rnaturalearth::ne_countries(country = c("United States of America","Canada", "Mexico"),returnclass = "sf")
western_country_sf <- st_crop(country_sf, xmin=-171.7911, xmax=-95, ymin=18, ymax=83.23324)
study_region_sf <- st_union(western_country_sf) # combines geoms
study_region_sf <- st_as_sf(st_cast(study_region_sf))  # deals with sinking multiple-shape file issue

## Initial filter by study area and uncertainty
# intersect and filter
phenologycv_sf6 <- st_intersection(phenologycv_sf5, study_region_sf) 
# rename for clarity
phenologycv_sf6 <- rename(phenologycv_sf6, coordinate_uncertainty = coord_uncertainty)
# keep NA uncertainty, if notekeep things with low uncertainty (<6.5km)
phenologycv_sf8 <- phenologycv_sf6 %>% filter( is.na(coordinate_uncertainty) | coordinate_uncertainty < 6500) 

## Build Elevation Bins 
# Load in NA elev data
elevation_1 <-rast("/blue/guralnick/millerjared/PhenoElevation/data/NAelevation4.tif") 
# clamp to only reasonable values for this analysis 
elevation_2 <- clamp(elevation_1, lower = 0, upper = 3000)
# ensure projection is the same 
crs(phenologycv_sf8) == crs(elevation_2) # proceed with extracting...
phenologycv_elev <- as.data.frame(terra::extract(elevation_2, phenologycv_sf8))
phenologycv_sf8$elevation <- phenologycv_elev$NAelevation4

# build many elevation bins
phenologycv_sf8$bins1 <- cut(phenologycv_sf8$elevation, breaks = c(0,500,1000,1500,2000,2500,3000), labels = c("0to.5", ".5to1", "1to1.5", "1.5to2","2to2.5","2.5to3"),include.lowest = TRUE)
phenologycv_sf8$bins2 <- cut(phenologycv_sf8$elevation, breaks = c(0,400,800,1200,1600,2000,2400,2800,3200), labels = c("0to.4", ".4to8", ".8to1.2", "1.2to1.6","1.6to2","2to2.4","2.4to2.8","2.8to3.2"),include.lowest = TRUE)
phenologycv_sf8$bins3 <- cut(phenologycv_sf8$elevation, breaks = c(0,300,600,900,1200,1500,1800,2100,2400,2700,3000), labels = c("0to.3", ".3to6", ".6to.9", ".9to1.2","1.2to1.5","1.5to1.8","1.8to2.1","2.1to2.4","2.4to2.7","2.7to3"),include.lowest = TRUE)
phenologycv_sf8$bins4 <- cut(phenologycv_sf8$elevation, breaks = c(0,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000), labels = c("0to.25", ".25to.5", ".5to.75", ".75to1","1to1.25","1.25to1.5","1.5to1.75","1.75to2","2to2.25","2.25to2.5","2.5to2.75","2.75to3"),include.lowest = TRUE)
phenologycv_sf8$bins5 <- cut(phenologycv_sf8$elevation, breaks = c(0,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000), labels = c("0to.2", ".2to.4", ".4to.6", ".6to.8",".8to1","1to1.2","1.2to1.4","1.4to1.6","1.6to1.8","1.8to2","2to2.2","2.2to2.4","2.4to2.6","2.6to2.8","2.8to3"),include.lowest = TRUE)
phenologycv_sf8$bins6 <- cut(phenologycv_sf8$elevation, breaks = c(0,150,300,450,600,750,900,1050,1200, 1350, 1500, 1650, 1800, 1950, 2100, 2250, 2400, 2550, 2700, 2850, 3000), labels = c("0to.15", ".15to.3", ".3to.45", ".45to.6", ".6to.75", ".75to.9", ".9to1.05", "1.05to1.2", "1.2to1.35", "1.35to1.5", "1.5to1.65", "1.65to1.8", "1.8to1.95", "1.95to2.1", "2.1to2.25", "2.25to2.4", "2.4to2.55", "2.55to2.7", "2.7to2.85", "2.85to3.0"),include.lowest = TRUE)
phenologycv_sf8$bins7 <- cut(phenologycv_sf8$elevation, breaks = c(0,100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000), labels = c("0to.1", ".1to.2", ".2to.3", ".3to.4", ".4to.5", ".5to.6", ".6to.7", ".7to.8", ".8to.9", ".9to1", "1to1.1", "1.1to1.2", "1.2to1.3", "1.3to1.4", "1.4to1.5", "1.5to1.6", "1.6to1.7", "1.7to1.8", "1.8to1.9", "1.9to2", "2to2.1", "2.1to2.2", "2.2to2.3", "2.3to2.4", "2.4to2.5", "2.5to2.6", "2.6to2.7", "2.7to2.8", "2.8to2.9", "2.9to3"),include.lowest = TRUE)

## Create a hex-grid cells over study region, delimit flowering occurrences to these cells
# transform to albers equal area and make grid to study extent
study_region_sf <-st_transform(study_region_sf, sp::CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
hex_grid <- st_make_grid(study_region_sf, cellsize = c(305000, 305000), what = "polygons", square = FALSE)
hex_grid_j <- st_join(st_as_sf(hex_grid), st_as_sf(study_region_sf))
hex_grid_study_region <- hex_grid_j[study_region_sf, col = '#ff000088'] # subset to study region. 
hex_grid_w_ids <- mutate(st_sf(geometry = hex_grid_study_region), id_cells = 1:n())
saveRDS(hex_grid_w_ids, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/full-study-region.rds")

# assign records to grid cells
phenologycv_sf9 <- st_transform(phenologycv_sf8,crs = st_crs(hex_grid))
phenologycv_sf10 <- st_join(phenologycv_sf9, hex_grid_w_ids)
phenologycv_sf11 <- setDT(phenologycv_sf10)

## Down weight abnormally high observation day of years, and remove outliers prior to calculating phenometrics
# down sample doys to be max 3
phenologycv_sf12 <- phenologycv_sf11  %>%
  group_by(scientific_name,year,bins6,id_cells,day_of_year) %>% slice_sample(n = 3)

# remove outliers in day of year that are mid-winter observations (temperate), take distinct obs 
phenologycv_sf13  <- phenologycv_sf12 %>% 
  filter(day_of_year > 20 & day_of_year < 340 ) %>%
  distinct(scientific_name, day_of_year,id_cells,  observed_metadata_url, bins6, .keep_all= TRUE)

# remove day of years that have less than 4 total observations over the course of the data sampling
phenologycv_sf14  <- phenologycv_sf13 %>% group_by(scientific_name,year,id_cells,bins6) %>% 
  filter(n_distinct(day_of_year)>4) 

## Count up flowering occurrences in elevational bands to decide on appropriate interval for analysis
phenologycv_counts_filter400m <- phenologycv_sf14 %>% group_by(scientific_name,year,bins2,id_cells) %>% summarize(num = n()) %>% filter(num >=7)
phenologycv_counts_filter500m <- phenologycv_sf14 %>% group_by(scientific_name,year,bins1,id_cells) %>% summarize(num = n()) %>% filter(num >=7)
phenologycv_counts_filter300m <- phenologycv_sf14 %>% group_by(scientific_name,year,bins3,id_cells) %>% summarize(num = n()) %>% filter(num >=7)
phenologycv_counts_filter250m <- phenologycv_sf14 %>% group_by(scientific_name,year,bins4,id_cells) %>% summarize(num = n()) %>% filter(num >=7)
phenologycv_counts_filter200m <- phenologycv_sf14 %>% group_by(scientific_name,year,bins5,id_cells) %>% summarize(num = n()) %>% filter(num >=7)
phenologycv_counts_filter150m <- phenologycv_sf14 %>% group_by(scientific_name,year,bins6,id_cells) %>% summarize(num = n()) %>% filter(num >=7)
phenologycv_counts_filter100m <- phenologycv_sf14 %>% group_by(scientific_name,year,bins7,id_cells) %>% summarize(num = n()) %>% filter(num >=7)

## We chose 150m bands based on these results. Format data for phenometrics
# get distinct day of years per species, keep only those combos that have at least 7 obs per species-year-elevbin-hexcell grouping
phenologycv_sf14_ndistinct <- phenologycv_sf14 %>%
  group_by(scientific_name,year,bins6,id_cells) %>%
  filter(n() >= 7) %>%
  dplyr::summarise(ndistinct = n_distinct(day_of_year))

# get total number of observations per species, keep only those combos that have at least 7 obs per species-year-elevbin-hexcell grouping
phenologycv_sf14_obscount <- phenologycv_sf14 %>%
  group_by(scientific_name,year,bins6,id_cells) %>%
  filter(n() >= 7) %>%
  dplyr::summarise(obscount = n())

# rename day of year to doy
names(phenologycv_sf14)[4]<-paste("doy")

## Write out data as an RDS object to be use in calc phenometrics script
saveRDS(phenologycv_sf14, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-ready-data.rds")
saveRDS(phenologycv_sf14_ndistinct, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-distinct-doy-combos.rds")
saveRDS(phenologycv_sf14_obscount, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-obs-combos.rds")

