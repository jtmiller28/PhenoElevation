### Title: Extract Climate Data
### Author: JT Miller
### Date: 12-29-2025

### Purpose: Extract Chelsea Climate Data from Study Region
setwd("/blue/guralnick/millerjared/PhenoElevation/") # set working-directory 

## Load Libraries 
library(terra)
library(sf)
library(tidyverse)

## Set up elevation data
elevation_1 <-rast("/blue/guralnick/millerjared/PhenoElevation/data/NAelevation4.tif") 
# clamp to only reasonable values for this analysis 
elevation_2 <- clamp(elevation_1, lower = 0, upper = 3000)

## Set up relevant study region data
# read in full study region
hex_grid_w_ids <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/full-study-region.rds")
# read in pglmm model ready data to filter down hexcells
phenologycv_sf14_onoffall20 <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-model-ready-data.rds")
# filter down hexcells to only include those relevant 
hex_grids_used <- hex_grid_w_ids %>% 
  dplyr::filter(id_cells %in% phenologycv_sf14_onoffall20$id_cells)

## Create a function that goes through by hex-cell and builds summaries...
summarize_hex <- function(hex_geom, hex_id, layer_name, layer_chelsa_name){
  # set up hex cell, crop/mask out our elevation data for this hexcell
  hex <- terra::vect(hex_geom) 
  hex <- terra::project(hex, elevation_2)
  cropped_elev <- terra::crop(elevation_2, hex)
  masked_elev <- terra::mask(cropped_elev, hex)
  # bring in climate data, match elevation data to this coarser resolution & proj, then resample. 
  tif_dir <- "/blue/guralnick/share/Chelsa-climate-data/climatologies/"
  tif_files <- list.files(tif_dir, full.names = TRUE)
  ras <- terra::rast(grep(paste0(layer_chelsa_name, "([^0-9]|$)"), tif_files, value = TRUE))
  masked_elev_rp <- terra::project(masked_elev, ras)
  elev_agg <- terra::resample(masked_elev_rp, ras, method = "bilinear") # interpolation of 3x3 windows
  # bin this newly aggregated elevation data to match our study's elevational bins
  breaks <- c(0,150,300,450,600,750,900,1050,1200,1350,1500,1650,1800,1950,
              2100,2250,2400,2550,2700,2850,3000)
  labels <- c("0to.15", ".15to.3", ".3to.45", ".45to.6", ".6to.75", ".75to.9", 
              ".9to1.05", "1.05to1.2", "1.2to1.35", "1.35to1.5", "1.5to1.65", 
              "1.65to1.8", "1.8to1.95", "1.95to2.1", "2.1to2.25", "2.25to2.4", 
              "2.4to2.55", "2.55to2.7", "2.7to2.85", "2.85to3.0")
  # make a reclassification matrix (rcl) so that each elevation bin gets a new cell value
  rcl <- cbind(
    breaks[-length(breaks)], 
    breaks[-1], 
    seq_along(labels)
  )
  # apply classification
  elev_binned <- classify(elev_agg, rcl = rcl, include.lowest = TRUE, right = FALSE)
  # give a name to clarify these are binned quantities 
  names(elev_binned) <- "elev_binned"
  # bring in climate data, crop/mask to our elev_binned raster, label as the name 
  ras_masked <- terra::mask(terra::crop(ras, elev_binned), elev_binned)
  names(ras_masked) <- layer_name # input layer name
  # combine spatRasters with two layers, climate layer + binned elevation
  stack <- c(ras_masked, elev_binned)
  # extract values as a df, remove NAs
  vals <- as.data.frame(terra::as.data.frame(stack, na.rm = TRUE, xy = FALSE))
  # summarize 
  elev_summaries <- vals %>% 
    dplyr::group_by(elev_binned) %>% 
    dplyr::summarize(mean_value = mean(!!sym(layer_name), na.rm = TRUE), num_cells = dplyr::n(), hex_id = hex_id) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(!!layer_name := mean_value)
  # clean up
  rm(hex, cropped_elev, masked_elev, ras, masked_elev_rp, elev_agg, elev_binned, ras_masked, stack, vals)
  gc()
  return(elev_summaries)
}

# Create a two vectors containing the name of Chelsa variable and an the actual measure to loop through the data and build datasets. 
chelsa_names <- c("bio1", "bio4", "bio5", "bio6", "bio12", "bio15", "scd") 
var_names <- c("mean_temp", "temp_seasonality", "max_temp", "min_temp", "ppt_yr", "ppt_seasonality", "snow_cover_days")

for(j in 1:length(chelsa_names)){
  summary_holding_list <- list()
  for(i in unique(hex_grids_used$id_cells)){
    
    summary_holding_list[[i]] <- summarize_hex(hex_geom = filter(hex_grids_used, id_cells == i),
                                               hex_id = i,
                                               layer_name = var_names[j], 
                                               layer_chelsa_name = chelsa_names[j])
    print(paste("finished cell", i, "of", length(unique(hex_grids_used$id_cells))))
    gc() # clean up for good measure
  } # end of i loop
  summary_table <- data.table::rbindlist(summary_holding_list)
  data.table::fwrite(summary_table, paste0("/blue/guralnick/millerjared/PhenoElevation/data/study-cells-", var_names[j], ".csv" ))
  print(paste("finished intra climate binning for", var_names[j]))
  rm(summary_table)
} # end of j loop