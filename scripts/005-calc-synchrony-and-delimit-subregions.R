### Title: Calculate Synchrony and Delimit Subregions   
### Author: Rob Guralnick
### Editor: JT Miller 
### Date: 12-29-2025

### Purpose: load in phenometric calcs model ready data, delimit subregions, convert to Covariation metrics, look at interspecific and infraspecific variation
setwd("/blue/guralnick/millerjared/PhenoElevation/") # set working-directory 

## Load Libraries 
library(lme4)
library(lmerTest)
library(performance)
library(dplyr)
library(tidyr)
library(ggplot2)
library(calecopal)
library(sjPlot)
library(patchwork)
library(MuMIn)

## Load previous centroid data & model ready phenology data
grid_cent_ll_point <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/grid_cent_ll_point.rds")
phenologycv_sf14_onoffall20 <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-model-ready-data.rds")
phenologycv_sf14_sptime <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/species-avg-timings.rds")
phenologycv_sf14_sptime2  <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/species-avg-timebins.rds")
SpeciesTraits3 <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/species-traits-for-filtered-taxa.rds")

## divide region of interest into Northern, Middle, Southern parts
grid_cent_ll_point$latbin2 <- cut(grid_cent_ll_point$latitude, breaks = 3, labels = c("SouthLat", "MidLat", "NorthLat"))
grid_cent_ll_point4 <- grid_cent_ll_point %>% select(id_cells, latbin2)

## add these subregions to the phenometrics data
phenologycv_sf14_onoffall21 <- merge(phenologycv_sf14_onoffall20, grid_cent_ll_point4, by=c("id_cells"))

## Create coefficient of variation estimates from phenometrics
# community level metrics
phenologycv_sf14_coefvar_sum <- phenologycv_sf14_onoffall20 %>%
  group_by(id_cells,bins6,year) %>%
  summarize(cv_onsetphen=EnvStats::cv(onset), 
            cv_offsetphen=EnvStats::cv(offset), 
            cv_50phen=EnvStats::cv(mean), # median
            cv_durphen=EnvStats::cv(duration), 
            count=n(), 
            n_distinct(scientific_name)) # redundant, should be the same as count but just as a check.

# IQR v.
phenologyiqr_sf14_coefvar_sum <- phenologycv_sf14_onoffall20 %>%
  group_by(id_cells,bins6,year) %>%
  summarize(
            iqr_onsetphen = EnvStats::iqr(onset, na.rm = TRUE),
            iqr_offsetphen=EnvStats::iqr(offset, na.rm = TRUE),
            iqr_50phen=EnvStats::iqr(mean, na.rm = TRUE), # median
            iqr_durphen=EnvStats::iqr(duration, na.rm = TRUE),
            count=n(), 
            n_distinct(scientific_name))

# infraspecific metrics
phenologycv_sf14_coefvar_sum_sp <- phenologycv_sf14_onoffall21 %>%
  group_by(scientific_name,bins6,latbin2) %>%
  summarize(cv_onsetphen=EnvStats::cv(onset), 
            cv_offsetphen=EnvStats::cv(offset), 
            cv_50phen=EnvStats::cv(mean), # median
            cv_durphen=EnvStats::cv(duration),
            count=n())

# iqr v.
phenologyiqr_sf14_coefvar_sum_sp <- phenologycv_sf14_onoffall21 %>%
  group_by(scientific_name,bins6,latbin2) %>%
  summarize(iqr_onsetphen=EnvStats::iqr(onset, na.rm = TRUE),
            iqr_offsetphen=EnvStats::iqr(offset, na.rm = TRUE), 
            iqr_50phen=EnvStats::iqr(mean, na.rm = TRUE), # median
            iqr_durphen=EnvStats::iqr(duration, na.rm = TRUE),
            count=n())

# drop missing vals 
phenologycv_sf14_coefvar_sum <- phenologycv_sf14_coefvar_sum %>% drop_na(cv_50phen)
phenologycv_sf14_coefvar_sum_sp <- phenologycv_sf14_coefvar_sum_sp %>% drop_na(cv_50phen)

## Reassign elev from bin names
phenologycv_sf14_coefvar_sum2 <- phenologycv_sf14_coefvar_sum %>% mutate(elev = dplyr::recode(bins6,  "0to.15" = 75, ".15to.3" = 225, ".3to.45" = 375, ".45to.6" = 525, ".6to.75" = 675, ".75to.9" = 825, ".9to1.05" = 975, "1.05to1.2" = 1125, "1.2to1.35" = 1275, "1.35to1.5" = 1425, "1.5to1.65" = 1575, "1.65to1.8" = 1725, "1.8to1.95" = 1875, "1.95to2.1" = 2025, "2.1to2.25" = 2175, "2.25to2.4" = 2325, "2.4to2.55" = 2475, "2.55to2.7" = 2625, "2.7to2.85" = 2775, "2.85to3.0" = 2925))
phenologycv_sf14_coefvar_sum_sp2 <- phenologycv_sf14_coefvar_sum_sp %>% mutate(elev = dplyr::recode(bins6, "0to.15" = 75, ".15to.3" = 225, ".3to.45" = 375, ".45to.6" = 525, ".6to.75" = 675, ".75to.9" = 825, ".9to1.05" = 975, "1.05to1.2" = 1125, "1.2to1.35" = 1275, "1.35to1.5" = 1425, "1.5to1.65" = 1575, "1.65to1.8" = 1725, "1.8to1.95" = 1875, "1.95to2.1" = 2025, "2.1to2.25" = 2175, "2.25to2.4" = 2325, "2.4to2.55" = 2475, "2.55to2.7" = 2625, "2.7to2.85" = 2775, "2.85to3.0" = 2925))
phenologyiqr_sf14_coefvar_sum2 <- phenologyiqr_sf14_coefvar_sum %>% mutate(elev = dplyr::recode(bins6,  "0to.15" = 75, ".15to.3" = 225, ".3to.45" = 375, ".45to.6" = 525, ".6to.75" = 675, ".75to.9" = 825, ".9to1.05" = 975, "1.05to1.2" = 1125, "1.2to1.35" = 1275, "1.35to1.5" = 1425, "1.5to1.65" = 1575, "1.65to1.8" = 1725, "1.8to1.95" = 1875, "1.95to2.1" = 2025, "2.1to2.25" = 2175, "2.25to2.4" = 2325, "2.4to2.55" = 2475, "2.55to2.7" = 2625, "2.7to2.85" = 2775, "2.85to3.0" = 2925))
phenologyiqr_sf14_coefvar_sum_sp2 <- phenologyiqr_sf14_coefvar_sum_sp %>% mutate(elev = dplyr::recode(bins6, "0to.15" = 75, ".15to.3" = 225, ".3to.45" = 375, ".45to.6" = 525, ".6to.75" = 675, ".75to.9" = 825, ".9to1.05" = 975, "1.05to1.2" = 1125, "1.2to1.35" = 1275, "1.35to1.5" = 1425, "1.5to1.65" = 1575, "1.65to1.8" = 1725, "1.8to1.95" = 1875, "1.95to2.1" = 2025, "2.1to2.25" = 2175, "2.25to2.4" = 2325, "2.4to2.55" = 2475, "2.55to2.7" = 2625, "2.7to2.85" = 2775, "2.85to3.0" = 2925))

## Grab the rest of the cols for community metrics, join with subregions, filter to remove where there are fewer than 4 sp
# change year to factor
phenologycv_sf14_coefvar_sum2$yearf <- as.factor(phenologycv_sf14_coefvar_sum2$year)
phenologyiqr_sf14_coefvar_sum2$yearf <- as.factor(phenologyiqr_sf14_coefvar_sum2$year)
# join community metrics with subregions
phenologycv_sf14_coefvar_sum2 <- left_join(phenologycv_sf14_coefvar_sum2 , grid_cent_ll_point , by=c("id_cells"))
phenologyiqr_sf14_coefvar_sum2 <- left_join(phenologyiqr_sf14_coefvar_sum2 , grid_cent_ll_point , by=c("id_cells"))
# filter to remove cases where there are fewer than 4 species
phenologycv_sf14_coefvar_sum3 <- phenologycv_sf14_coefvar_sum2 %>% filter(count >=4)
phenologyiqr_sf14_coefvar_sum3 <- phenologyiqr_sf14_coefvar_sum2 %>% filter(count >=4)
# filter to remove cases where there is 0 covariation for median
phenologycv_sf14_coefvar_sum4 <- phenologycv_sf14_coefvar_sum3 %>% filter(cv_50phen>0)
phenologyiqr_sf14_coefvar_sum4 <- phenologyiqr_sf14_coefvar_sum3 %>% filter(iqr_50phen>0)
# drop nas 
phenologycv_sf14_coefvar_sum4c <- phenologycv_sf14_coefvar_sum4 %>% drop_na(cv_onsetphen, elev, latitude,longitude)
phenologyiqr_sf14_coefvar_sum4c <- phenologyiqr_sf14_coefvar_sum4 %>% drop_na(iqr_onsetphen, elev, latitude,longitude)
# save as a rds object for further script 
saveRDS(phenologycv_sf14_coefvar_sum4c, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/cv-community-data.rds")
saveRDS(phenologyiqr_sf14_coefvar_sum4c, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/iqr-community-data.rds")

## Grab the rest of the cols for infraspecific metrics, clean up field names, add in species average metrics calculated prior in sp traits and post pheno cleaning
# format names
SpeciesTraits3$scientific_name <- gsub(" ", "_",SpeciesTraits3$scientific_name)
# merge traits
phenologycv_sf14_coefvar_sum_sp3 <- merge(phenologycv_sf14_coefvar_sum_sp2 ,SpeciesTraits3 , by=c("scientific_name"))
phenologyiqr_sf14_coefvar_sum_sp3 <- merge(phenologyiqr_sf14_coefvar_sum_sp2 ,SpeciesTraits3 , by=c("scientific_name"))
# reformat names
phenologycv_sf14_sptime$scientific_name <- gsub(" ", "_",phenologycv_sf14_sptime2$scientific_name )
# join time bins
phenologycv_sf14_coefvar_sum_sp4 <- left_join(phenologycv_sf14_coefvar_sum_sp3 ,phenologycv_sf14_sptime2 , by=c("scientific_name"))
phenologyiqr_sf14_coefvar_sum_sp4 <- left_join(phenologyiqr_sf14_coefvar_sum_sp3 ,phenologycv_sf14_sptime2 , by=c("scientific_name"))
# join avg species times
phenologycv_sf14_coefvar_sum_sp4 <- left_join(phenologycv_sf14_coefvar_sum_sp3 , phenologycv_sf14_sptime, by=c("scientific_name"))
phenologyiqr_sf14_coefvar_sum_sp4 <- left_join(phenologyiqr_sf14_coefvar_sum_sp3 , phenologycv_sf14_sptime, by=c("scientific_name"))

## For infraspecific, scale, factorize, and filter to remove where there are fewer than 4 species per elevational band
# scale the elev and mean sp timing fields 
phenologycv_sf14_coefvar_sum_sp4$elev_sc <- c(scale(phenologycv_sf14_coefvar_sum_sp4$elev))
phenologyiqr_sf14_coefvar_sum_sp4$elev_sc <- c(scale(phenologyiqr_sf14_coefvar_sum_sp4$elev))
phenologycv_sf14_coefvar_sum_sp4$spmean_sc2 <- c(scale(phenologycv_sf14_coefvar_sum_sp4$spmean))
phenologyiqr_sf14_coefvar_sum_sp4$spmean_sc2 <- c(scale(phenologyiqr_sf14_coefvar_sum_sp4$spmean))
# factorize the plant habit fields & latitude 
phenologycv_sf14_coefvar_sum_sp4$geoannperf <- as.factor(phenologycv_sf14_coefvar_sum_sp4$geoannper)
phenologycv_sf14_coefvar_sum_sp4$geoannperf2 <- as.factor(phenologycv_sf14_coefvar_sum_sp4$geopannper2)
phenologycv_sf14_coefvar_sum_sp4$geoannperf3 <- as.factor(phenologycv_sf14_coefvar_sum_sp4$geoannper3)
phenologyiqr_sf14_coefvar_sum_sp4$geoannperf3 <- as.factor(phenologyiqr_sf14_coefvar_sum_sp4$geoannper3)
phenologycv_sf14_coefvar_sum_sp4$latbinf2 <- as.factor(phenologycv_sf14_coefvar_sum_sp4$latbin2)
phenologyiqr_sf14_coefvar_sum_sp4$latbinf2 <- as.factor(phenologyiqr_sf14_coefvar_sum_sp4$latbin2)
# filter where there are less than 4 species 
phenologycv_sf14_coefvar_sum_sp5 <- phenologycv_sf14_coefvar_sum_sp4  %>% filter(count >=4)
phenologyiqr_sf14_coefvar_sum_sp5 <- phenologyiqr_sf14_coefvar_sum_sp4  %>% filter(count >=4)
phenologycv_sf14_coefvar_sum_sp6 <- phenologycv_sf14_coefvar_sum_sp5  %>% group_by(species) %>% filter(n() >= 4) %>% ungroup()
phenologyiqr_sf14_coefvar_sum_sp6 <- phenologyiqr_sf14_coefvar_sum_sp5  %>% group_by(species) %>% filter(n() >= 4) %>% ungroup()

## Remove cases where cv is exactly zero
phenologycv_sf14_coefvar_sum_sp7 <- phenologycv_sf14_coefvar_sum_sp6 %>% filter(cv_offsetphen>0)
phenologycv_sf14_coefvar_sum_sp8 <- phenologycv_sf14_coefvar_sum_sp7 %>% filter(cv_durphen>0)
phenologycv_sf14_coefvar_sum_sp9 <- phenologycv_sf14_coefvar_sum_sp8 %>% filter(cv_50phen>0)
phenologyiqr_sf14_coefvar_sum_sp7 <- phenologyiqr_sf14_coefvar_sum_sp6 %>% filter(iqr_offsetphen>0)
phenologyiqr_sf14_coefvar_sum_sp8 <- phenologyiqr_sf14_coefvar_sum_sp7 %>% filter(iqr_durphen>0)
phenologyiqr_sf14_coefvar_sum_sp9 <- phenologyiqr_sf14_coefvar_sum_sp8 %>% filter(iqr_50phen>0)
phenologycv_sf14_coefvar_sum_sp10 <- phenologycv_sf14_coefvar_sum_sp9 # removed the filter here...
phenologyiqr_sf14_coefvar_sum_sp10 <- phenologyiqr_sf14_coefvar_sum_sp9
## Check histograms for skew
# without log-transform
ggplot(phenologycv_sf14_coefvar_sum_sp10) + 
  geom_histogram(mapping = aes(x = cv_onsetphen)) + 
  ggtitle("Distribution of CV Onset Infraspecific") # right-skew
ggplot(phenologycv_sf14_coefvar_sum_sp10) + 
  geom_histogram(mapping = aes(x = cv_offsetphen)) + 
  ggtitle("Distribution of CV Offset Infraspecific") # right-skew
ggplot(phenologycv_sf14_coefvar_sum_sp10) + 
  geom_histogram(mapping = aes(x = cv_durphen)) + 
  ggtitle("Distribution of CV Duration Infrapspecific") # right-skew

ggplot(phenologyiqr_sf14_coefvar_sum_sp10) + 
  geom_histogram(mapping = aes(x = iqr_onsetphen)) + 
  ggtitle("Distribution of IQR Onset Infraspecific") # right-skew
ggplot(phenologyiqr_sf14_coefvar_sum_sp10) + 
  geom_histogram(mapping = aes(x = iqr_offsetphen)) + 
  ggtitle("Distribution of IQR Offset Infraspecific") # right-skew
ggplot(phenologyiqr_sf14_coefvar_sum_sp10) + 
  geom_histogram(mapping = aes(x = iqr_durphen)) + 
  ggtitle("Distribution of IQR Duration Infrapspecific") # right-skew

# with log-transform
ggplot(phenologycv_sf14_coefvar_sum_sp10) + 
  geom_histogram(mapping = aes(x = log(cv_onsetphen))) + 
  ggtitle("Distribution of log(CV Onset) Infraspecific") # better
ggplot(phenologycv_sf14_coefvar_sum_sp10) + 
  geom_histogram(mapping = aes(x = log(cv_offsetphen))) + 
  ggtitle("Distribution of log(CV Offset) Infraspecific") # better
ggplot(phenologycv_sf14_coefvar_sum_sp10) + 
  geom_histogram(mapping = aes(x = log(cv_durphen))) + 
  ggtitle("Distribution of log(CV Duration) Infrapspecific") # better

ggplot(phenologyiqr_sf14_coefvar_sum_sp10) + 
  geom_histogram(mapping = aes(x = log(iqr_onsetphen))) + 
  ggtitle("Distribution of log(IQR Onset) Infraspecific") # better
ggplot(phenologyiqr_sf14_coefvar_sum_sp10) + 
  geom_histogram(mapping = aes(x = log(iqr_offsetphen))) + 
  ggtitle("Distribution of log(IQR Offset) Infraspecific") # better
ggplot(phenologyiqr_sf14_coefvar_sum_sp10) + 
  geom_histogram(mapping = aes(x = log(iqr_durphen))) + 
  ggtitle("Distribution of log(IQR Duration) Infraspecific") # better

## Naming stuff for figures
phenologycv_sf14_coefvar_sum_sp10$geoannperf3 <- factor(
  phenologycv_sf14_coefvar_sum_sp10$geoannperf3,
  levels = c("anbien", "herbper"),
  labels = c("Annual & Biannual", "Herbaceous Perennial")
)

phenologyiqr_sf14_coefvar_sum_sp10$geoannperf3 <- factor(
  phenologyiqr_sf14_coefvar_sum_sp10$geoannperf3,
  levels = c("anbien", "herbper"),
  labels = c("Annual & Biannual", "Herbaceous Perennial")
)
phenologycv_sf14_coefvar_sum_sp10$latbinf2 <- factor(
  phenologycv_sf14_coefvar_sum_sp10$latbinf2,
  levels = c("SouthLat", "MidLat", "NorthLat"),
  labels = c("Southern Latitudes", "Middle Latitudes", "Northern Latitudes")
)
phenologyiqr_sf14_coefvar_sum_sp10$latbinf2 <- factor(
  phenologyiqr_sf14_coefvar_sum_sp10$latbinf2,
  levels = c("SouthLat", "MidLat", "NorthLat"),
  labels = c("Southern Latitudes", "Middle Latitudes", "Northern Latitudes")
)

## Get Temp seasonality data for Lat/Elev Bin
temp_seasonality <- data.table::fread("/blue/guralnick/millerjared/PhenoElevation/data/study-cells-temp_seasonality.csv")
data.table::setnames(temp_seasonality, "num_cells", "num_cells_temp_seasonality")
## Match up elevation bin labels with previous labels used
# bring in label names 
elev_labels <- c("0to.15", ".15to.3", ".3to.45", ".45to.6", ".6to.75", ".75to.9", 
                 ".9to1.05", "1.05to1.2", "1.2to1.35", "1.35to1.5", "1.5to1.65", 
                 "1.65to1.8", "1.8to1.95", "1.95to2.1", "2.1to2.25", "2.25to2.4", 
                 "2.4to2.55", "2.55to2.7", "2.7to2.85", "2.85to3.0")
# add median elev labels per interval for modeling
median_elev_labels <- c(75, 225, 375, 525, 675, 825, 975, 1125, 1275, 1425, 1575, 
                        1725, 1875, 2025, 2175, 2325, 2475, 2625, 2775, 2925)
# order we built in extract clim data
elev_bin_labels <- c(1:20) 

# make relational table
elev_relation_table <- data.frame(bins6 = elev_labels, elev_binned = elev_bin_labels, median_elev = median_elev_labels)

# join this with temp seasonality
temp_seasonality <- temp_seasonality %>% 
  left_join(elev_relation_table, by = "elev_binned") %>% 
  rename(id_cells = hex_id)

temp_seasonality_latbinned <- merge(temp_seasonality, grid_cent_ll_point4, by=c("id_cells"))

temp_seasonality_latbinned <- temp_seasonality_latbinned %>% 
  select(latbin2, bins6, temp_seasonality) %>% 
  group_by(latbin2, bins6) %>% 
  summarize(temp_seasonality_avg = mean(temp_seasonality)) %>% 
  ungroup()


# merge temp seasonality with community iqrs with hexcell and elevation bins 
phenologyiqr_sf14_coefvar_sum_sp11  <- phenologyiqr_sf14_coefvar_sum_sp10  %>% 
  left_join(temp_seasonality_latbinned, by = c("bins6", "latbin2"))

phenologyiqr_sf14_coefvar_sum_sp11$temp_seasonality_avg_sc <- c(scale(phenologyiqr_sf14_coefvar_sum_sp11$temp_seasonality_avg))

## Save as RDS object for intraspecific model 
saveRDS(phenologycv_sf14_coefvar_sum_sp10, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/cv-intraspecific-data.rds")
saveRDS(phenologyiqr_sf14_coefvar_sum_sp11, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/iqr-intraspecific-data.rds")

