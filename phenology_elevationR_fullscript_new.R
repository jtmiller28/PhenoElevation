library(arrow)
library(rnaturalearth)
library(tidyr)
library(sf)
library(elevatr)
library(dplyr)
library(terra)
library(data.table)
library(readr)
library(rtrees)
library(stringr)
library(lme4)
library(lmerTest)
library(car)
library(performance)
library(brms)
library(sjPlot)
library(ape)
library(TreeTools)
library(brms)
library(DHARMa)
library(car)
library(MuMIn)

setwd("C:/Users/robgu/Downloads")

#filter_list
elev_species2 <- read_csv("elev_species2.csv")
elev_species2$species <- gsub("_", " ", elev_species2$species)

#process_giant_phenologycv_dataset
phenologycv_csv <- open_dataset(
    sources = "annotations_all2.csv", 
    col_types = schema(ISBN = string()),
    format = "csv")

phenologycv_csv |>
  group_by(year) |>
  write_dataset(path = "phenologycv_presence_data", format = "parquet")

phenologycv_pd <- open_dataset("phenologycv_presence_data")

phenologycv_pd2 <- phenologycv_pd  |> 
  filter(year >= 2017, year != 2024, trait == "flower", scientific_name %in% elev_species2$species)

phenologycv_pds3 <- phenologycv_pd2 |> collect()

phenologycv_pd4 <- phenologycv_pds3  %>% 
   group_by(scientific_name) %>% 
   filter(n() >= 500)

#convert to spatial data frame
crs_dd <- 4326
phenologycv_sf5  <- sf::st_as_sf(phenologycv_pd4, coords = c("longitude", "latitude"), crs = crs_dd)

#get area of interest
test <- rnaturalearth::ne_countries(country = c("United States of America","Canada", "Mexico"),returnclass = "sf")
test2 = st_crop(test, xmin=-171.7911, xmax=-95, ymin=18, ymax=83.23324)
test3 <- st_union(test2)

#more filtering of records (spatial)
phenologycv_sf6 <- st_intersection(phenologycv_sf5,test3)
phenologycv_sf8 <- phenologycv_sf6 %>% filter( is.na(coordinate_uncertainty) | coordinate_uncertainty < 6500) 

#get_elevation_data_and_create_elevation_estimate_for_all_records
#elevation_1 <- elevatr::get_elev_raster(locations = test2, z = 9, clip = "locations",override_size_check = TRUE)
elevation_1 <-rast("NAelevation4.tif")
phenologycv_elev <- as.data.frame(terra::extract(elevation_1, phenologycv_sf8))
phenologycv_sf8$elevation <- phenologycv_elev$NAelevation4
#phenologycv_sf8$elevation <- as.numeric(phenologycv_sf8$elevation)

#bin_elevation_I_tried_a_number_of_binnin_options
phenologycv_sf8$bins1 <- cut(phenologycv_sf8$elevation, breaks = c(0,500,1000,1500,2000,2500,3000,3500,4000,4500), labels = c("0to.5", ".5to1", "1to1.5", "1.5to2","2to2.5","2.5to3","3to3.5","3.5to4","4to4.5"),include.lowest = TRUE)
phenologycv_sf8$bins2 <- cut(phenologycv_sf8$elevation, breaks = c(0,400,800,1200,1600,2000,2400,2800,3200,3600,4000,4400), labels = c("0to.4", ".4to8", ".8to1.2", "1.2to1.6","1.6to2","2to2.4","2.4to2.8","2.8to3.2","3.2to3.6","3.6to4","4to4.4"),include.lowest = TRUE)
phenologycv_sf8$bins3 <- cut(phenologycv_sf8$elevation, breaks = c(0,300,600,900,1200,1500,1800,2100,2400,2700,3000,3300,3600,3900,4200,4500), labels = c("0to.3", ".3to6", ".6to.9", ".9to1.2","1.2to1.5","1.5to1.8","1.8to2.1","2.1to2.4","2.4to2.7","2.7to3","3to3.3","3.3to3.6","3.6to3.9","3.9to4.2","4.2to4.5"),include.lowest = TRUE)
phenologycv_sf8$bins4 <- cut(phenologycv_sf8$elevation, breaks = c(0,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,4250,4500), labels = c("0to.25", ".25to.5", ".5to.75", ".75to1","1to1.25","1.25to1.5","1.5to1.75","1.75to2","2to2.25","2.25to2.5","2.5to2.75","2.75to3","3to3.25","3.25to3.5","3.5to3.75","3.75to4","4to4.25","4.25to4.5"),include.lowest = TRUE)
phenologycv_sf8$bins5 <- cut(phenologycv_sf8$elevation, breaks = c(0,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,4200,4400), labels = c("0to.2", ".2to.4", ".4to.6", ".6to.8",".8to1","1to1.2","1.2to1.4","1.4to1.6","1.6to1.8","1.8to2","2to2.2","2.2to2.4","2.4to2.6","2.6to2.8","2.8to3","3to3.2","3.2to3.4","3.4to3.6","3.6to3.8","3.8to4","4to4.2","4.2to4.4"),include.lowest = TRUE)
phenologycv_sf8$bins6 <- cut(phenologycv_sf8$elevation, breaks = c(0,150,300,450,600,750,900,1050,1200, 1350, 1500, 1650, 1800, 1950, 2100, 2250, 2400, 2550, 2700, 2850, 3000, 3150, 3300, 3450, 3600, 3750, 3900, 4050, 4200, 4350, 4500), labels = c("0to.15", ".15to.3", ".3to.45", ".45to.6", ".6to.75", ".75to.9", ".9to1.05", "1.05to1.2", "1.2to1.35", "1.35to1.5", "1.5to1.65", "1.65to1.8", "1.8to1.95", "1.95to2.1", "2.1to2.25", "2.25to2.4", "2.4to2.55", "2.55to2.7", "2.7to2.85", "2.85to3.0", "3.0to3.15", "3.15to3.3", "3.3to3.45", "3.45to3.6", "3.6to3.75", "3.75to3.9", "3.9to4.05", "4.05to4.2", "4.2to4.35", "4.35to4.5"),include.lowest = TRUE)
phenologycv_sf8$bins7 <- cut(phenologycv_sf8$elevation, breaks = c(0,100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500), labels = c("0to.1", ".1to.2", ".2to.3", ".3to.4", ".4to.5", ".5to.6", ".6to.7", ".7to.8", ".8to.9", ".9to1", "1to1.1", "1.1to1.2", "1.2to1.3", "1.3to1.4", "1.4to1.5", "1.5to1.6", "1.6to1.7", "1.7to1.8", "1.8to1.9", "1.9to2", "2to2.1", "2.1to2.2", "2.2to2.3", "2.3to2.4", "2.4to2.5", "2.5to2.6", "2.6to2.7", "2.7to2.8", "2.8to2.9", "2.9to3", "3to3.1", "3.1to3.2", "3.2to3.3", "3.3to3.4", "3.4to3.5", "3.5to3.6", "3.6to3.7", "3.7to3.8", "3.8to3.9", "3.9to4", "4to4.1", "4.1to4.2", "4.2to4.3", "4.3to4.4", "4.4to4.5"),include.lowest = TRUE)

#project_data_and_make_grid_over_region
test3<-st_transform(test2, sp::CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
grids <- st_make_grid(test3, cellsize = c(305000, 305000), what = "polygons", square = FALSE)
grids2 <- st_join(st_as_sf(grids), st_as_sf(test3))
grids3 <- grids2[test3, col = '#ff000088']
grids3 <- mutate(st_sf(geometry = grids3), id_cells = 1:n())

#assign_records_to_grids
phenologycv_sf9 <- st_transform(phenologycv_sf8,crs =st_crs(grids))
phenologycv_sf10 <- st_join(phenologycv_sf9, grids3)
phenologycv_sf11 <- setDT(phenologycv_sf10 )

#account_for_unusual_sampling_days
phenologycv_sf12 <- phenologycv_sf11  %>%
  group_by(scientific_name,year,bins6,id_cells,day_of_year) %>% slice_sample(n = 3)

#clean_outliers_remove_duplicates_and_distinct_days
phenologycv_sf13  <- phenologycv_sf12 %>% 
  filter(day_of_year > 20 & day_of_year < 340 ) %>%
  distinct(scientific_name, day_of_year,id_cells, bins6, observed_metadata_url, .keep_all= TRUE)

#remove_cases_where_too_few_days_of_years_sampled
phenologycv_sf14  <- phenologycv_sf13 %>% group_by(scientific_name,year,id_cells,bins6) %>% 
  filter(n_distinct(day_of_year)>4) 
  
#lets_look_at some_count
phenologycv_counts_filter400m <- phenologycv_sf14 %>% group_by(scientific_name,year,bins2,id_cells) %>% summarize(num = n()) %>% filter(num >=7)
phenologycv_counts_filter500m <- phenologycv_sf14 %>% group_by(scientific_name,year,bins1,id_cells) %>% summarize(num = n()) %>% filter(num >=7)
phenologycv_counts_filter300m <- phenologycv_sf14 %>% group_by(scientific_name,year,bins3,id_cells) %>% summarize(num = n()) %>% filter(num >=7)
phenologycv_counts_filter250m <- phenologycv_sf14 %>% group_by(scientific_name,year,bins4,id_cells) %>% summarize(num = n()) %>% filter(num >=7)
phenologycv_counts_filter200m <- phenologycv_sf14 %>% group_by(scientific_name,year,bins5,id_cells) %>% summarize(num = n()) %>% filter(num >=7)
phenologycv_counts_filter150m <- phenologycv_sf14 %>% group_by(scientific_name,year,bins6,id_cells) %>% summarize(num = n()) %>% filter(num >=7)
phenologycv_counts_filter100m <- phenologycv_sf14 %>% group_by(scientific_name,year,bins7,id_cells) %>% summarize(num = n()) %>% filter(num >=7)

#### ended up using 150m elevation bins for rest of this script
#get_distinct_days
  phenologycv_sf14_ndistinct <- phenologycv_sf14 %>%
  group_by(scientific_name,year,bins6,id_cells) %>%
  filter(n() >= 7) %>%
  dplyr::summarise(ndistinct = n_distinct(day_of_year))


#get_observation_counts
  phenologycv_sf14_obscount <- phenologycv_sf14 %>%
  group_by(scientific_name,year,bins6,id_cells) %>%
  filter(n() >= 7) %>%
  dplyr::summarise(obscount = n())
    
####################################### Run_phenometrics_minimum_record_number_is_7_and_must_have_at_least_4_distinct_days
names(phenologycv_sf14)[4]<-paste("doy")
phenologycv_sf14_onset <-  phenologycv_sf14 %>% 
  group_by(scientific_name,year,bins6,id_cells) %>% 
  filter(n() >= 7) %>%
  filter(n_distinct(doy) > 3) %>%
  group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$doy, percentile = 0.05, bootstraps=250)))

phenologycv_sf14_offset <-  phenologycv_sf14 %>% 
  group_by(scientific_name,year,bins6,id_cells) %>% 
  filter(n() >= 7) %>%
  filter(n_distinct(doy) > 3) %>%
  group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$doy, percentile = 0.95, bootstraps=250)))

phenologycv_sf14_50 <-  phenologycv_sf14 %>% 
  group_by(scientific_name,year,bins6,id_cells) %>% 
  filter(n() >= 7) %>%
  filter(n_distinct(doy) > 3) %>%
  group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$doy, percentile = 0.50, bootstraps=250)))


#######################################post_phenometric_data_assembly_steps
##cleanup_phenometrics_and_assemble_into_formats_compatible_with_other_estimates

phenologycv_sf14_onset_2 <- dplyr::select(phenologycv_sf14_onset, -c(sd, median, min, max, range, skew, trimmed,mad,kurtosis,se))
phenologycv_sf14_onset_3 <- pivot_wider(phenologycv_sf14_onset_2, names_from = column, values_from = mean)

phenologycv_sf14_50_2 <- dplyr::select(phenologycv_sf14_50, -c(sd, median, min, max, range, skew, trimmed,mad,kurtosis,se))
phenologycv_sf14_50_3 <- pivot_wider(phenologycv_sf14_50_2, names_from = column, values_from = mean)

phenologycv_sf14_offset_2 <- dplyr::select(phenologycv_sf14_offset, -c(sd, median, min, max, range, skew, trimmed,mad,kurtosis,se))
phenologycv_sf14_offset_3 <- pivot_wider(phenologycv_sf14_offset_2, names_from = column, values_from = mean)

phenologycv_sf14_onoff <- merge(phenologycv_sf14_onset_3,phenologycv_sf14_offset_3, by=c("year","bins6","scientific_name","id_cells"))
phenologycv_sf14_onoff50 <- merge(phenologycv_sf14_onoff,phenologycv_sf14_50_3, by=c("year","bins6","scientific_name","id_cells"))

phenologycv_sf14_onoffall <- merge(phenologycv_sf14_onoff50,phenologycv_sf14_ndistinct, by=c("year","bins6","scientific_name","id_cells"))
phenologycv_sf14_onoffall2 <- merge(phenologycv_sf14_onoffall,phenologycv_sf14_obscount, by=c("year","bins6","scientific_name","id_cells"))

phenologycv_sf14_onoffall3 <- phenologycv_sf14_onoffall2 %>%
  dplyr::rename(mean=estimate,mean_low=low_ci, mean_high=high_ci, onset=estimate.x,onset_low=low_ci.x, onset_high=high_ci.x, offset =estimate.y,offset_low=low_ci.y, offset_high=high_ci.y)

#####get_lat_lon_for_centroid_of_grid_cells_assign_lat_lon_to_phenometrics
#grids4 <- st_transform(grids3,crs("EPSG:4326"))
grid_cent <- st_centroid(grids3)
grid_cent_ll <- st_transform(grid_cent,crs("EPSG:4326"))
grid_cent_ll_point <- grid_cent_ll %>%
  mutate(longitude = st_coordinates(.)[,1],latitude = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  as.data.frame()
phenologycv_sf14_onoffall4 <- left_join(phenologycv_sf14_onoffall3 , grid_cent_ll_point , by="id_cells")

#add_duration_metric
phenologycv_sf14_onoffall4 <- phenologycv_sf14_onoffall4 %>% 
  mutate(duration = offset-onset)

#remove_outlier and low count species
phenologycv_sf14_onoffall5 <-  phenologycv_sf14_onoffall4 %>% 
  group_by(scientific_name) %>% 
  mutate(zonset = scale(onset)) %>% 
  filter(between(zonset,-3.25,+3.25))

phenologycv_sf14_onoffall6 <-  phenologycv_sf14_onoffall5 %>% 
  group_by(scientific_name) %>% 
  mutate(zoffset = scale(offset)) %>% 
  filter(between(zoffset,-3.25,+3.25))

phenologycv_sf14_onoffall7 <-  phenologycv_sf14_onoffall6 %>% 
  group_by(scientific_name) %>% 
  mutate(zmean = scale(mean)) %>% 
  filter(between(zmean,-3.25,+3.25))

phenologycv_sf14_onoffall8 <-  phenologycv_sf14_onoffall7 %>% 
  group_by(scientific_name) %>% 
  mutate(zdur = scale(duration)) %>% 
  filter(between(zdur,-3.25,+3.25))

phenologycv_sf14_onoffall9 <- phenologycv_sf14_onoffall8 %>% group_by(scientific_name) %>% filter(n()>6)

#assemble_species_traits
SpeciesTraits3 <- read_csv("SpeciesTraits3.csv")
species_in_dataset <- phenologycv_sf14_onoffall9$scientific_name
species_in_dataset2 <- gsub(" ", "_", species_in_dataset)
SpeciesTraits3 <- SpeciesTraits3 %>% filter(scientific_name %in% species_in_dataset2)

#join_traits_to_dataset_do_more_filtering
SpeciesTraits3$scientific_name <- gsub("_", " ",SpeciesTraits3$scientific_name)
phenologycv_sf14_onoffall10 <- left_join(phenologycv_sf14_onoffall9, SpeciesTraits3  , by="scientific_name")
phenologycv_sf14_elevcat <- phenologycv_elev %>%  mutate(mean_elev_bin = ntile(elev_mean, n=4))
phenologycv_sf14_onoffall11 <- left_join(phenologycv_sf14_onoffall10, phenologycv_sf14_elevcat, by="scientific_name")
phenologycv_sf14_onoffall11 <- phenologycv_sf14_onoffall11 %>% filter(ndistinct>4)
phenologycv_sf14_onoffall11 <- phenologycv_sf14_onoffall11 %>% drop_na(bins6)
phenologycv_sf14_onoffall12 <- as.data.frame(ungroup(phenologycv_sf14_onoffall11))

#recode_bins_to_elevation_values_and_drop_any_NAs
phenologycv_sf14_onoffall13 <- phenologycv_sf14_onoffall12 %>% mutate(elev = dplyr::recode(bins6,  "0to.15" = 75, ".15to.3" = 225, ".3to.45" = 375, ".45to.6" = 525, ".6to.75" = 675, ".75to.9" = 825, ".9to1.05" = 975, "1.05to1.2" = 1125, "1.2to1.35" = 1275, "1.35to1.5" = 1425, "1.5to1.65" = 1575, "1.65to1.8" = 1725, "1.8to1.95" = 1875, "1.95to2.1" = 2025, "2.1to2.25" = 2175, "2.25to2.4" = 2325, "2.4to2.55" = 2475, "2.55to2.7" = 2625, "2.7to2.85" = 2775, "2.85to3.0" = 2925, "3.0to3.15" = 3075, "3.15to3.3" = 3225, "3.3to3.45" = 3375, "3.45to3.6" = 3525, "3.6to3.75" = 3675, "3.75to3.9" = 3825, "3.9to4.05" = 3975, "4.05to4.2" = 4125, "4.2to4.35" = 4275, "4.35to4.5" = 4425))
phenologycv_sf14_onoffall13 <- phenologycv_sf14_onoffall13 %>% drop_na(elev)

#get_species_mean_timing
phenologycv_sf14_sptime <- phenologycv_sf14_onoffall13 %>% group_by(scientific_name) %>% summarize(spmean=mean(mean))
phenologycv_sf14_sptime2 <- phenologycv_sf14_sptime %>%  mutate(timebin = ntile(spmean, n=2), timebin2 = ntile(spmean, n=3))
phenologycv_sf14_onoffall14 <- merge(phenologycv_sf14_onoffall13, phenologycv_sf14_sptime2, by=c("scientific_name"))

#more filtering
phenologycv_sf14_onoffall14 <- phenologycv_sf14_onoffall14 %>% filter(longitude < -95)
phenologycv_sf14_onoffall14 <- phenologycv_sf14_onoffall14 %>% filter(duration > 10)
phenologycv_sf14_onoffall14$scientific_name <- gsub(" ", "_", phenologycv_sf14_onoffall14$scientific_name)
phenologycv_sf14_onoffall15 <- phenologycv_sf14_onoffall14 %>%  mutate(annpertime = paste(geoannper, timebin, sep = '_'))
phenologycv_sf14_onoffall15 <- phenologycv_sf14_onoffall15 %>% filter(elev < 3000)

#yet_more_filtering_to_make_sure_theres_variation_across_elevtion_over_cells_and_species
phenologycv_sf14_onoffall_elevf2 <- phenologycv_sf14_onoffall15 %>% 
  group_by(scientific_name) %>%
  summarize(elev_range = max(elev)-min(elev)) %>%
  filter(elev_range >= 255)

phenologycv_sf14_onoffall_elevf3 <- phenologycv_sf14_onoffall15 %>% 
  group_by(id_cells) %>%
  summarize(elev_range = max(elev)-min(elev)) %>%
  filter(elev_range >= 225)

#phenologycv_sf14_onoffall_elevf4 <- phenologycv_sf14_onoffall15 %>% 
#  group_by(id_cells) %>%
#  summarize(elev_cell_mean = mean(elev)) %>%
#  filter(elev_cell_mean >= 150)

phenologycv_sf14_onoffall16 <- subset(phenologycv_sf14_onoffall15, scientific_name %in% phenologycv_sf14_onoffall_elevf2$scientific_name)
phenologycv_sf14_onoffall17 <- subset(phenologycv_sf14_onoffall16, id_cells %in% phenologycv_sf14_onoffall_elevf3$id_cells)
#phenologycv_sf14_onoffall18 <- subset(phenologycv_sf14_onoffall17, id_cells %in% phenologycv_sf14_onoffall_elevf4$id_cells)

#remove_some_extraneous_columns_and_remove_species_with_low#ofphenometrics
phenologycv_sf14_onoffall19 <- phenologycv_sf14_onoffall17[,-c(19:186)]
phenologycv_sf14_onoffall20 <- phenologycv_sf14_onoffall19 %>% group_by(scientific_name) %>% filter(n()>=5)


#scale_and_create_factors
phenologycv_sf14_onoffall20 <- phenologycv_sf14_onoffall20 %>%
  mutate(latitude_sc =scale(latitude), 
         longitude_sc=scale(longitude),
         elev_sc = scale(elev),
         spmean_sc <- scale(spmean)
         geophyte2f = as.factor(geophyte2),
         geoannperf = as.factor(geoannper),
         geoannper2f = as.factor(geopannper2),
         geoannper3f = as.factor(geoannper3),
         mean_elev_binf = as.factor(mean_elev_bin),
         timebinf = as.factor(timebin),
         timebinf2 = as.factor(timebin2),
         annpertimef = as.factor(annpertime),
         id_cellsf = as.factor(id_cells),
         yearf = as.factor(year))


######################### get phylogeny
sp_list <- subset(SpeciesTraits3, select=c('species','genus','family'))
sp_list$scientific_name <- str_to_title(sp_list$species)
sp_list2 <- sp_list %>% filter(sp_list$scientific_name %in% unique(phenologycv_sf14_onoffall20$scientific_name))
plant_tree1 <- get_tree(sp_list = sp_list2, taxon = "plant", scenario = "at_basal_node", show_grafted = FALSE)
A <- ape::vcv.phylo(plant_tree1)


#########################  create_variation_metrics
#  divide region of interest into Northern and Southern parts
grid_cent_ll_point2 <-  grid_cent_ll_point %>% filter(latitude < 60)
grid_cent_ll_point3 <- grid_cent_ll_point2 %>% filter(latitude > 25)
grid_cent_ll_point3$latbin1 <- cut(grid_cent_ll_point3$latitude, breaks = 2, labels = c("NorthLat", "SouthLat"))
grid_cent_ll_point4 <- grid_cent_ll_point3[,169:172] 

# add definition of northern and southern regions to phenometrics
phenologycv_sf14_onoffall21 <- merge(phenologycv_sf14_onoffall20,grid_cent_ll_point4, by=c("id_cells"))

# create coefficient of variation estimates from phenometrics
# community level metrics
phenologycv_sf14_coefvar_sum <- phenologycv_sf14_onoffall20 %>%
  group_by(id_cells,bins6,year) %>%
  summarize(cv_onsetphen=EnvStats::cv(onset), cv_offsetphen=EnvStats::cv(offset), cv_50phen=EnvStats::cv(mean), cv_durphen=EnvStats::cv(duration), count=n(), n_distinct(scientific_name))

#intraspecific metrics
phenologycv_sf14_coefvar_sum_sp <- phenologycv_sf14_onoffall21 %>%
  group_by(scientific_name,bins6,latbin1) %>%
  summarize(cv_onsetphen=EnvStats::cv(onset), cv_offsetphen=EnvStats::cv(offset), cv_50phen=EnvStats::cv(mean), cv_durphen=EnvStats::cv(duration),count=n())

#drop_missing_values_should_work_for_all_metrics
phenologycv_sf14_coefvar_sum <- phenologycv_sf14_coefvar_sum %>% drop_na(cv_50phen)
phenologycv_sf14_coefvar_sum_sp <- phenologycv_sf14_coefvar_sum_sp %>% drop_na(cv_50phen)

#reassign_elevation_values_from_bin_names
phenologycv_sf14_coefvar_sum2 <- phenologycv_sf14_coefvar_sum %>% mutate(elev = dplyr::recode(bins6,  "0to.15" = 75, ".15to.3" = 225, ".3to.45" = 375, ".45to.6" = 525, ".6to.75" = 675, ".75to.9" = 825, ".9to1.05" = 975, "1.05to1.2" = 1125, "1.2to1.35" = 1275, "1.35to1.5" = 1425, "1.5to1.65" = 1575, "1.65to1.8" = 1725, "1.8to1.95" = 1875, "1.95to2.1" = 2025, "2.1to2.25" = 2175, "2.25to2.4" = 2325, "2.4to2.55" = 2475, "2.55to2.7" = 2625, "2.7to2.85" = 2775, "2.85to3.0" = 2925, "3.0to3.15" = 3075, "3.15to3.3" = 3225, "3.3to3.45" = 3375, "3.45to3.6" = 3525, "3.6to3.75" = 3675, "3.75to3.9" = 3825, "3.9to4.05" = 3975, "4.05to4.2" = 4125, "4.2to4.35" = 4275, "4.35to4.5" = 4425))
phenologycv_sf14_coefvar_sum_sp2 <- phenologycv_sf14_coefvar_sum_sp %>% mutate(elev = dplyr::recode(bins6, "0to.15" = 75, ".15to.3" = 225, ".3to.45" = 375, ".45to.6" = 525, ".6to.75" = 675, ".75to.9" = 825, ".9to1.05" = 975, "1.05to1.2" = 1125, "1.2to1.35" = 1275, "1.35to1.5" = 1425, "1.5to1.65" = 1575, "1.65to1.8" = 1725, "1.8to1.95" = 1875, "1.95to2.1" = 2025, "2.1to2.25" = 2175, "2.25to2.4" = 2325, "2.4to2.55" = 2475, "2.55to2.7" = 2625, "2.7to2.85" = 2775, "2.85to3.0" = 2925, "3.0to3.15" = 3075, "3.15to3.3" = 3225, "3.3to3.45" = 3375, "3.45to3.6" = 3525, "3.6to3.75" = 3675, "3.75to3.9" = 3825, "3.9to4.05" = 3975, "4.05to4.2" = 4125, "4.2to4.35" = 4275, "4.35to4.5" = 4425))

#get_rest_of_needed_columns_for_community_metrics_and_filter_to_remove_metrics_where_fewer_than_4_species
phenologycv_sf14_coefvar_sum2$yearf <- as.factor(phenologycv_sf14_coefvar_sum2$year)
phenologycv_sf14_coefvar_sum2 <- left_join(phenologycv_sf14_coefvar_sum2 , grid_cent_ll_point , by=c("id_cells"))
phenologycv_sf14_coefvar_sum3 <- phenologycv_sf14_coefvar_sum2 %>% filter(count >=4)

#get_rest_of_data_for_intraspecific_metrics_and_clean_up
SpeciesTraits3$scientific_name <- gsub(" ", "_",SpeciesTraits3$scientific_name)
phenologycv_sf14_coefvar_sum_sp3 <- merge(phenologycv_sf14_coefvar_sum_sp2 ,SpeciesTraits3 , by=c("scientific_name"))
phenologycv_sf14_sptime$scientific_name <- gsub(" ", "_",phenologycv_sf14_sptime2$scientific_name )
phenologycv_sf14_coefvar_sum_sp4 <- left_join(phenologycv_sf14_coefvar_sum_sp3 ,phenologycv_sf14_sptime2 , by=c("scientific_name"))
phenologycv_sf14_coefvar_sum_sp4 <- left_join(phenologycv_sf14_coefvar_sum_sp3 , phenologycv_sf14_sptime, by=c("scientific_name"))

#set_factors_and_scale_filter_cases_fewer_that_4_species_cell_elevbin_
phenologycv_sf14_coefvar_sum_sp4$elev_sc <- scale(phenologycv_sf14_coefvar_sum_sp4$elev)
phenologycv_sf14_coefvar_sum_sp4$spmean_sc2 <- scale(phenologycv_sf14_coefvar_sum_sp4$spmean)
phenologycv_sf14_coefvar_sum_sp4$geoannperf <- as.factor(phenologycv_sf14_coefvar_sum_sp4$geoannper)
phenologycv_sf14_coefvar_sum_sp4$geoannperf2 <- as.factor(phenologycv_sf14_coefvar_sum_sp4$geopannper2)
phenologycv_sf14_coefvar_sum_sp4$geoannperf3 <- as.factor(phenologycv_sf14_coefvar_sum_sp4$geoannper3)
phenologycv_sf14_coefvar_sum_sp4$latbinf1 <- as.factor(phenologycv_sf14_coefvar_sum_sp4$latbin1)
phenologycv_sf14_coefvar_sum_sp5 <- phenologycv_sf14_coefvar_sum_sp4  %>% filter(count >=4)
phenologycv_sf14_coefvar_sum_sp6 <- phenologycv_sf14_coefvar_sum_sp5  %>% group_by(species) %>% filter(n() >= 4)

#remove_cases_where_cv_is_exactly_zero_and_filter_elevation_where_data_are_missing
phenologycv_sf14_coefvar_sum_sp7 <- phenologycv_sf14_coefvar_sum_sp6 %>% filter(cv_offsetphen>0)
phenologycv_sf14_coefvar_sum_sp8 <- phenologycv_sf14_coefvar_sum_sp7 %>% filter(cv_durphen>0)
phenologycv_sf14_coefvar_sum_sp9 <- phenologycv_sf14_coefvar_sum_sp8 %>% filter(cv_50phen>0)
phenologycv_sf14_coefvar_sum_sp10 <- phenologycv_sf14_coefvar_sum_sp9 %>% filter(elev<2850)

#### models_in_different_r_script
