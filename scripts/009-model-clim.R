### Title: Extract Climate Data
### Author: JT Miller
### Date: 12-29-2025

### Purpose: Model trends of climate variables over this region
setwd("/blue/guralnick/millerjared/PhenoElevation/") # set working-directory 

## Load Libraries 
library(tidyverse)
library(sf)
library(data.table)
library(calecopal)
library(effects)
library(patchwork)

## Set up relevant study region data
# read in full study region
hex_grid_w_ids <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/full-study-region.rds")
# read in pglmm model ready data to filter down hexcells
phenologycv_sf14_onoffall20 <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-model-ready-data.rds")
# filter down hexcells to only include those relevant 
hex_grids_used <- hex_grid_w_ids %>% 
  dplyr::filter(id_cells %in% phenologycv_sf14_onoffall20$id_cells)

## Read in extracted climate data, rename fields for clarity, combine into one singular table
# read in tables made during extract clim data
mean_temp <- data.table::fread("/blue/guralnick/millerjared/PhenoElevation/data/study-cells-mean_temp.csv")
min_temp <- data.table::fread("/blue/guralnick/millerjared/PhenoElevation/data/study-cells-min_temp.csv")
max_temp <- data.table::fread("/blue/guralnick/millerjared/PhenoElevation/data/study-cells-max_temp.csv")
temp_seasonality <- data.table::fread("/blue/guralnick/millerjared/PhenoElevation/data/study-cells-temp_seasonality.csv")
ppt_yr <- data.table::fread("/blue/guralnick/millerjared/PhenoElevation/data/study-cells-ppt_yr.csv")
ppt_seasonality <- data.table::fread("/blue/guralnick/millerjared/PhenoElevation/data/study-cells-ppt_seasonality.csv")
scds <- data.table::fread("/blue/guralnick/millerjared/PhenoElevation/data/study-cells-snow_cover_days.csv")

# rename the number of cells data was extracted from field 
data.table::setnames(mean_temp, "num_cells", "num_cells_mean_temp")
data.table::setnames(min_temp, "num_cells", "num_cells_min_temp")
data.table::setnames(max_temp, "num_cells", "num_cells_max_temp")
data.table::setnames(temp_seasonality, "num_cells", "num_cells_temp_seasonality")
data.table::setnames(ppt_yr, "num_cells", "num_cells_ppt_yr")
data.table::setnames(ppt_seasonality, "num_cells", "num_cells_ppt_seasonality")
data.table::setnames(scds, "num_cells", "num_cells_scds")
# combine 
climate_vars_table <- Reduce(function(x,y) merge(x,y, by = c("elev_binned", "hex_id"), all = TRUE), 
                            list(mean_temp, min_temp, max_temp, temp_seasonality, 
                                 ppt_yr, ppt_seasonality, scds))

# rename hex_id to id_cells for clarity
data.table::setnames(climate_vars_table, "hex_id", "id_cells")

## Load in centroids and divide up with latitude Northern, Middle, and Southern bins as done in calc-sychrony script, combine with climate summary table
grid_cent_ll_point <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/grid_cent_ll_point.rds")

# divide region of interest into Northern, Middle, Southern parts
grid_cent_ll_point$latbin2 <- cut(grid_cent_ll_point$latitude, breaks = 3, labels = c("SouthLat", "MidLat", "NorthLat"))
grid_cent_ll_point4 <- grid_cent_ll_point %>% select(id_cells, latbin2, longitude, latitude)

# merge with climate vars table
climate_vars_table <- climate_vars_table %>% 
  left_join(grid_cent_ll_point4, by = "id_cells")

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

# merge, and select fields for simplicity
grid_w_climate_model <- climate_vars_table %>% 
  left_join(elev_relation_table, by = "elev_binned") %>% 
  select(id_cells, latbin2, mean_temp, min_temp, max_temp, temp_seasonality,
         ppt_yr, ppt_seasonality, snow_cover_days,
         bins6, latitude, longitude, elev_binned)

## scale values for fitting a simple model 
grid_w_climate_model$elev_sc <- scale(grid_w_climate_model$elev)
grid_w_climate_model$lat_sc <- scale(grid_w_climate_model$latitude)
grid_w_climate_model$lon_sc <- scale(grid_w_climate_model$longitude)

## Compute quantiles for elevation (scaled version)
elev_quantiles <- quantile(grid_w_climate_model$elev_sc, probs = c(0.1, 0.5, 0.9))

# create names
low_value  <- elev_quantiles[1] # 10th percentile
mid_value  <- elev_quantiles[2] # 50th percentile (median)
high_value <- elev_quantiles[3] # 90th percentile

## Fit model to predict how climatic variables are interacting with Latitude vs Elevation
mean_temp_mod <- lm(mean_temp ~ lat_sc + elev_sc + lat_sc*elev_sc, 
                     data = grid_w_climate_model)
min_temp_mod <- lm(min_temp ~ lat_sc + elev_sc + lat_sc*elev_sc, 
                      data = grid_w_climate_model)
max_temp_mod <- lm(max_temp ~ lat_sc + elev_sc + lat_sc*elev_sc, 
                     data = grid_w_climate_model)
temp_season_mod <- lm(temp_seasonality ~ lat_sc + elev_sc + lat_sc*elev_sc, 
                      data = grid_w_climate_model)
ppt_yr_mod <- lm(ppt_yr ~ lat_sc + elev_sc + lat_sc*elev_sc, 
                      data = grid_w_climate_model)
ppt_season_mod <- lm(ppt_seasonality ~ lat_sc + elev_sc + lat_sc*elev_sc, 
                      data = grid_w_climate_model)
scds_mod <- lm(snow_cover_days ~ lat_sc + elev_sc + lat_sc*elev_sc, 
                      data = grid_w_climate_model)

## Gather up effects, create a visual 
# get all effects
mean_temp_all_eff <- effects::allEffects(mean_temp_mod, xlevels = list(lat_sc = c(low_value, mid_value, high_value)))
min_temp_all_eff <- effects::allEffects(min_temp_mod,  xlevels = list(lat_sc = c(low_value, mid_value, high_value)))
max_temp_all_eff <- effects::allEffects(max_temp_mod,  xlevels = list(lat_sc = c(low_value, mid_value, high_value)))
temp_season_all_eff <- effects::allEffects(temp_season_mod,  xlevels = list(lat_sc = c(low_value, mid_value, high_value)))
ppt_yr_all_eff <- effects::allEffects(ppt_yr_mod,  xlevels = list(lat_sc = c(low_value, mid_value, high_value)))
ppt_season_all_eff <- effects::allEffects(ppt_season_mod,  xlevels = list(lat_sc = c(low_value, mid_value, high_value)))
scds_all_eff <- effects::allEffects(scds_mod,  xlevels = list(lat_sc = c(low_value, mid_value, high_value)))

# create a df of the interaction 
mean_temp_lat_elev_effect_df <- as.data.frame(mean_temp_all_eff[["lat_sc:elev_sc"]])
min_temp_lat_elev_effect_df <- as.data.frame(min_temp_all_eff[["lat_sc:elev_sc"]])
max_temp_lat_elev_effect_df <- as.data.frame(max_temp_all_eff[["lat_sc:elev_sc"]])
temp_season_lat_elev_effect_df <- as.data.frame(temp_season_all_eff[["lat_sc:elev_sc"]])
ppt_yr_lat_elev_effect_df <- as.data.frame(ppt_yr_all_eff[["lat_sc:elev_sc"]])
ppt_season_lat_elev_effect_df <- as.data.frame(ppt_season_all_eff[["lat_sc:elev_sc"]])
scds_lat_elev_effect_df <- as.data.frame(scds_all_eff[["lat_sc:elev_sc"]])

# label
mean_temp_lat_elev_effect_df <- mean_temp_lat_elev_effect_df %>%
  mutate(lat_label = case_when(
    dplyr::near(lat_sc, low_value)  ~ "Southern",
    dplyr::near(lat_sc, mid_value)  ~ "Middle",
    dplyr::near(lat_sc, high_value) ~ "Northern",
    TRUE ~ NA_character_
  ), 
  lat_label = factor(
    lat_label,
    levels = c(
      "Southern",
      "Middle",
      "Northern"
    )))
min_temp_lat_elev_effect_df <- min_temp_lat_elev_effect_df %>%
  mutate(lat_label = case_when(
    dplyr::near(lat_sc, low_value)  ~ "Southern",
    dplyr::near(lat_sc, mid_value)  ~ "Middle",
    dplyr::near(lat_sc, high_value) ~ "Northern",
    TRUE ~ NA_character_
  ), 
  lat_label = factor(
    lat_label,
    levels = c(
      "Southern",
      "Middle",
      "Northern"
    )))
max_temp_lat_elev_effect_df <- max_temp_lat_elev_effect_df %>%
  mutate(lat_label = case_when(
    dplyr::near(lat_sc, low_value)  ~ "Southern",
    dplyr::near(lat_sc, mid_value)  ~ "Middle",
    dplyr::near(lat_sc, high_value) ~ "Northern",
    TRUE ~ NA_character_
  ), 
  lat_label = factor(
    lat_label,
    levels = c(
      "Southern",
      "Middle",
      "Northern"
    )))
temp_season_lat_elev_effect_df <- temp_season_lat_elev_effect_df %>%
  mutate(lat_label = case_when(
    dplyr::near(lat_sc, low_value)  ~ "Southern",
    dplyr::near(lat_sc, mid_value)  ~ "Middle",
    dplyr::near(lat_sc, high_value) ~ "Northern",
    TRUE ~ NA_character_
  ), 
  lat_label = factor(
    lat_label,
    levels = c(
      "Southern",
      "Middle",
      "Northern"
    )))
ppt_yr_lat_elev_effect_df <- ppt_yr_lat_elev_effect_df %>%
  mutate(lat_label = case_when(
    dplyr::near(lat_sc, low_value)  ~ "Southern",
    dplyr::near(lat_sc, mid_value)  ~ "Middle",
    dplyr::near(lat_sc, high_value) ~ "Northern",
    TRUE ~ NA_character_
  ), 
  lat_label = factor(
    lat_label,
    levels = c(
      "Southern",
      "Middle",
      "Northern"
    )))
ppt_season_lat_elev_effect_df <- ppt_season_lat_elev_effect_df %>%
  mutate(lat_label = case_when(
    dplyr::near(lat_sc, low_value)  ~ "Southern",
    dplyr::near(lat_sc, mid_value)  ~ "Middle",
    dplyr::near(lat_sc, high_value) ~ "Northern",
    TRUE ~ NA_character_
  ), 
  lat_label = factor(
    lat_label,
    levels = c(
      "Southern",
      "Middle",
      "Northern"
    )))
scds_lat_elev_effect_df <- scds_lat_elev_effect_df %>%
  mutate(lat_label = case_when(
    dplyr::near(lat_sc, low_value)  ~ "Southern",
    dplyr::near(lat_sc, mid_value)  ~ "Middle",
    dplyr::near(lat_sc, high_value) ~ "Northern",
    TRUE ~ NA_character_
  ), 
  lat_label = factor(
    lat_label,
    levels = c(
      "Southern",
      "Middle",
      "Northern"
    )))


## Plot 
n_bands <- 20
band_colors <- cal_palette("sierra1", n = n_bands, type = "continuous")
cols <- band_colors[c(1, 8, 20)]

mean_temp_lat_elev_plot <- ggplot(mean_temp_lat_elev_effect_df, mapping = aes(x = elev_sc, y = fit, color = as.factor(lat_label))) +
  geom_line() + 
  geom_ribbon(mapping = aes(ymin = lower, ymax = upper, fill = as.factor(lat_label)), alpha = 0.2, color = NA) + 
  scale_fill_manual(values = cols, guide = guide_legend(title = "Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  labs(x = "Scaled Elevation", 
       y = expression(paste("Predicted Mean Temperature (", degree, "C)")), 
       color = "Relative Latitude\nof Study Region", 
       fill = "Relative Latitude\nof Study Region", 
       title = "Effect of Elevation by Latitude on Mean Temperature") + 
  theme_bw() +
  guides(color = guide_legend(reverse = TRUE),
           fill  = guide_legend(reverse = TRUE))

min_temp_lat_elev_plot <- ggplot(min_temp_lat_elev_effect_df, mapping = aes(x = elev_sc, y = fit, color = as.factor(lat_label))) +
  geom_line() + 
  geom_ribbon(mapping = aes(ymin = lower, ymax = upper, fill = as.factor(lat_label)), alpha = 0.2, color = NA) + 
  scale_fill_manual(values = cols, guide = guide_legend(title = "Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  labs(x = "Scaled Elevation", 
       y = expression(paste("Predicted Minimum Temperature (", degree, "C)")), 
       color = "Relative Latitude\nof Study Region", 
       fill = "Relative Latitude\nof Study Region", 
       title = "Effect of Elevation by Latitude on Minimum Temperature") + 
  theme_bw() + 
  guides(color = guide_legend(reverse = TRUE),
         fill  = guide_legend(reverse = TRUE))

max_temp_lat_elev_plot <- ggplot(max_temp_lat_elev_effect_df, mapping = aes(x = elev_sc, y = fit, color = as.factor(lat_label))) +
  geom_line() + 
  geom_ribbon(mapping = aes(ymin = lower, ymax = upper, fill = as.factor(lat_label)), alpha = 0.2, color = NA) + 
  scale_fill_manual(values = cols, guide = guide_legend(title = "Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  labs(x = "Scaled Elevation", 
       y = expression(paste("Predicted Maximum Temperature (", degree, "C)")), 
       color = "Relative Latitude\nof Study Region", 
       fill = "Relative Latitude\nof Study Region", 
       title = "Effect of Elevation by Latitude on Maximum Temperature") + 
  theme_bw() + 
  guides(color = guide_legend(reverse = TRUE),
         fill  = guide_legend(reverse = TRUE))

temp_season_lat_elev_plot <- ggplot(temp_season_lat_elev_effect_df, mapping = aes(x = elev_sc, y = fit, color = as.factor(lat_label))) +
  geom_line() + 
  geom_ribbon(mapping = aes(ymin = lower, ymax = upper, fill = as.factor(lat_label)), alpha = 0.2, color = NA) + 
  scale_fill_manual(values = cols, guide = guide_legend(title = "Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  labs(x = "Scaled Elevation", 
       y = expression(paste("Predicted Temperature Seasonality (", degree, "C/100)")), 
       color = "Relative Latitude\nof Study Region", 
       fill = "Relative Latitude\nof Study Region", 
       title = "Effect of Elevation by Latitude on Temperature Seasonality") + 
  theme_bw() + 
  guides(color = guide_legend(reverse = TRUE),
         fill  = guide_legend(reverse = TRUE))

ppt_yr_lat_elev_plot <- ggplot(ppt_yr_lat_elev_effect_df, mapping = aes(x = elev_sc, y = fit, color = as.factor(lat_label))) +
  geom_line() + 
  geom_ribbon(mapping = aes(ymin = lower, ymax = upper, fill = as.factor(lat_label)), alpha = 0.2, color = NA) + 
  scale_fill_manual(values = cols, guide = guide_legend(title = "Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  labs(x = "Scaled Elevation", 
       y = expression("Predicted Annual Precipitation (kg m"^{-2}~year^{-1}*")"),
       color = "Relative Latitude\nof Study Region", 
       fill = "Relative Latitude\nof Study Region", 
       title = "Effect of Elevation by Latitude on Annual Precipitation") + 
  theme_bw() + 
  guides(color = guide_legend(reverse = TRUE),
         fill  = guide_legend(reverse = TRUE))

ppt_season_lat_elev_plot <- ggplot(ppt_season_lat_elev_effect_df, mapping = aes(x = elev_sc, y = fit, color = as.factor(lat_label))) +
  geom_line() + 
  geom_ribbon(mapping = aes(ymin = lower, ymax = upper, fill = as.factor(lat_label)), alpha = 0.2, color = NA) + 
  scale_fill_manual(values = cols, guide = guide_legend(title = "Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  labs(x = "Scaled Elevation", 
       y = expression("Predicted Precipitation Seasonality (kg m"^{-2}*")"),
       color = "Relative Latitude\nof Study Region", 
       fill = "Relative Latitude\nof Study Region", 
       title = "Effect of Elevation by Latitude on Precipitation Seasonality") + 
  theme_bw() + 
  guides(color = guide_legend(reverse = TRUE),
         fill  = guide_legend(reverse = TRUE))

scds_lat_elev_plot <- ggplot(scds_lat_elev_effect_df, mapping = aes(x = elev_sc, y = fit, color = as.factor(lat_label))) +
  geom_line() + 
  geom_ribbon(mapping = aes(ymin = lower, ymax = upper, fill = as.factor(lat_label)), alpha = 0.2, color = NA) + 
  scale_fill_manual(values = cols, guide = guide_legend(title = "Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  labs(x = "Scaled Elevation", 
       y = "Predicted Snow Cover Days (number of days)",
       color = "Relative Latitude\nof Study Region", 
       fill = "Relative Latitude\nof Study Region", 
       title = "Effect of Elevation by Latitude on Snow Cover Days") + 
  theme_bw()

## Patchwork these plots into a full figure
# Remove legends from first plots
mean_temp_lat_elev_plot_nolegend <- mean_temp_lat_elev_plot + ggtitle("Mean Temperature") + theme(legend.position = "none")
ppt_yr_lat_elev_plot_nolegend <- ppt_yr_lat_elev_plot + ggtitle("Yearly Precipitation") + theme(legend.position = "none")
temp_season_lat_elev_plot_nolegend <- temp_season_lat_elev_plot + ggtitle("Temperature Seasonality") + theme(legend.position = "none")
ppt_season_lat_elev_plot_legend <- ppt_season_lat_elev_plot + ggtitle("Precipitation Seasonality") # keep legend here

# Combine plots
combined_plot <- wrap_plots(
  mean_temp_lat_elev_plot_nolegend,
  ppt_yr_lat_elev_plot_nolegend,
  temp_season_lat_elev_plot_nolegend,
  ppt_season_lat_elev_plot_legend,
  ncol = 2, nrow = 2
) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Effect of Elevation by Latitude on Climate",
    theme = theme(plot.title = element_text(size = 20))
  ) &
  theme(legend.position = "right", legend.justification = "center")
print(combined_plot)

ggsave("./outputs/figures/elev-by-lat-climate-focus-plots.png", width = 15, height = 10)


