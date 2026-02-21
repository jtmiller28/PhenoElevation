### Title: Hopkins Deviation Post-hoc Analysis
### Author: JT Miller
### Date: 2-1-2026

## Purpose: Calculate species level deviation from Hopkins law, per elevation and latitudnal gain

## Load libraries 
library(data.table)
library(ggplot2)
library(dplyr)
library(terra)
library(purrr)
library(sf)
library(lme4)
library(lmerTest)
library(splines)

## Bring in hexcells, limit hexcells to those use in our phenoestimates
hexcells <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/full-study-region.rds")
pheno_model_rdy_data <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-model-ready-data.rds")
hex_cells_used <- hexcells[hexcells$id_cells %in% unique(pheno_model_rdy_data$id_cells), ]
hex_cells_used <- st_set_geometry(hex_cells_used, "x")

## Load in NA elev data, restrict to 0-3000m, create 150m bins
elevation_1 <-rast("/blue/guralnick/millerjared/PhenoElevation/data/NAelevation4.tif") 
# clamp to only reasonable values for this analysis 
elevation_2 <- clamp(elevation_1, lower = 0, upper = 3000)
# define elev bins
elev_breaks <- seq(0, 3000, by = 150)

## create a function to take each hexcell and calc the latitudnal span within each elevation bin
process_lats_per_hexcell_elev_bin <- function(id_cells, 
                                              hex_geom, 
                                              elevation_rast){
  # crop and mask elevation raster to hexcell
  hex_vect <- vect(hex_geom)
  hex_vect <- terra::project(hex_vect, elevation_2)
  elev_crop <- crop(elevation_rast, hex_vect, mask = TRUE)
  
  # convert to points w/coords
  elev_pts <- as.data.frame(elev_crop, xy=TRUE, na.rm=TRUE)
  colnames(elev_pts) <- c("lon", "lat", "elevation")
  
  # bin elevations
  elev_pts$elev_bin <- cut(elev_pts$elevation, 
                           breaks = elev_breaks, 
                           include.lowest = TRUE, 
                           right = FALSE)
  
  # calc latitude stats per bin
  elev_pts %>% 
    group_by(elev_bin) %>% 
    summarize(
      id_cells = id_cells, 
      lat_min = min(lat), 
      lon_min = min(lon), 
      lat_max = max(lat), 
      lon_max = max(lon), 
      lat_mean = mean(lat),
      lon_mean = mean(lon), 
      lat_median = median(lat), 
      lon_median = median(lon),
      lat_sd = sd(lat),
      lon_sd = sd(lon),
      n_pixels = n(), 
      .groups = "drop"
    )
}

# Apply to all hexcells
bin_summary <- map_df(1:nrow(hex_cells_used), ~{
  hex_row <- hex_cells_used[.x, ]
  process_lats_per_hexcell_elev_bin(
    hex_row$id_cells,
    st_geometry(hex_row),
    elevation_2
  )
})

# save output
saveRDS(bin_summary, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/hexcell_bin_summaries.rds")

bin_summary <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/hexcell_bin_summaries.rds")

# Combine lat summary with phenometric data
breaks <- c("[0,150)","[150,300)","[300,450)","[450,600)", "[600,750)" ,"[750,900)",
            "[900,1.05e+03)","[1.05e+03,1.2e+03)","[1.2e+03,1.35e+03)","[1.35e+03,1.5e+03)",
            "[1.5e+03,1.65e+03)","[1.65e+03,1.8e+03)","[1.8e+03,1.95e+03)", "[1.95e+03,2.1e+03)", 
            "[2.1e+03,2.25e+03)","[2.25e+03,2.4e+03)","[2.4e+03,2.55e+03)","[2.55e+03,2.7e+03)",
            "[2.7e+03,2.85e+03)","[2.85e+03,3e+03]")
labels <- c("0to.15", ".15to.3", ".3to.45", ".45to.6", ".6to.75", ".75to.9", 
            ".9to1.05", "1.05to1.2", "1.2to1.35", "1.35to1.5", "1.5to1.65", 
            "1.65to1.8", "1.8to1.95", "1.95to2.1", "2.1to2.25", "2.25to2.4", 
            "2.4to2.55", "2.55to2.7", "2.7to2.85", "2.85to3.0")
midpoints <- as.integer((elev_breaks[-length(elev_breaks)] + elev_breaks[-1])/2)
elev_mapping <- data.frame(
  elev_bin = breaks, 
  labels = labels, 
  mid_elev = midpoints
)
bin_summary2 <- bin_summary %>% 
 left_join(elev_mapping, by = "elev_bin") %>% 
  rename(bins6 = labels)

## Make spatial delimitations as before, remove species that do not have enough occs across these spatial delimitations
grid_cent_ll_point <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/grid_cent_ll_point.rds")
## divide region of interest into Northern, Middle, Southern parts
grid_cent_ll_point$latbin2 <- cut(grid_cent_ll_point$latitude, breaks = 3, labels = c("SouthLat", "MidLat", "NorthLat"))
grid_cent_ll_point4 <- grid_cent_ll_point %>% select(id_cells, latbin2)

## add these subregions to the phenometrics data
pheno_model_rdy_data <- merge(pheno_model_rdy_data, grid_cent_ll_point4, by=c("id_cells"))

pheno_w_lats <- pheno_model_rdy_data %>% left_join(bin_summary2, by = c("id_cells", "bins6"))

## Remove species with too low of latitude variation
pheno_w_lats$lat_mean_sc <- scale(pheno_w_lats$lat_mean)

overall_lat_range <- range(pheno_w_lats$lat_mean_sc)
overall_lat_width <- diff(overall_lat_range)

# filter for species who surpass this
pheno_w_lats_filtered <- pheno_w_lats %>% 
  group_by(scientific_name) %>% 
  filter((max(lat_mean_sc, na.rm = TRUE) - min(lat_mean_sc, na.rm = TRUE)) >= 0.25 * overall_lat_width) %>% 
  ungroup()


# Use Rob's estimate sp anchors fxn
HOPKINS_LAT_SLOPE <- 4
HOPKINS_ELEV_SLOPE <- 4/122 
estimate_species_anchors <- function(df) {
  df |> group_by(scientific_name) |>
    group_modify(~{
      x <- .x
      slat <- median(x$lat_mean, na.rm = TRUE)
      selev <- median(x$mid_elev, na.rm = TRUE)  
      
      if (nrow(x) >= 10 && length(unique(x$lat_mean[!is.na(x$lat_mean)])) > 2) {
        m <- lm(onset ~ lat_mean, data = x)
        sdoy <- as.numeric(predict(m, newdata = data.frame(lat_mean = slat)))
      } else {
        sdoy <- median(x$onset, na.rm = TRUE)
      }
      tibble(start_lat = slat, start_elev = selev, start_doy = sdoy) 
    }) |>
    ungroup()
}

add_hopkins_predictions <- function(df) {
  df |>
    mutate(
      hopkins_lat = start_doy + (lat_mean - start_lat) * HOPKINS_LAT_SLOPE,
      hopkins_lat_elev = start_doy +
        (lat_mean - start_lat) * HOPKINS_LAT_SLOPE +
        (mid_elev - start_elev) * HOPKINS_ELEV_SLOPE 
    )
}


anchored_sp <- estimate_species_anchors(pheno_w_lats_filtered)

pheno_w_lats_filtered <- pheno_w_lats_filtered %>% 
  left_join(anchored_sp, by = "scientific_name")

pheno_w_pred <- add_hopkins_predictions(pheno_w_lats_filtered)

# pheno_w_dev <- pheno_w_pred %>% 
#   mutate(hopkins_dev = hopkins_lat_elev - onset)

pheno_w_dev <- pheno_w_pred %>%
  mutate(hopkins_dev = onset - hopkins_lat_elev)

### Heatmap of Deviation
hexcell_dev <- pheno_w_dev %>% 
  group_by(id_cells) %>% 
  summarize(mean_abs_dev = mean(abs(hopkins_dev), na.rm = TRUE), 
            n_obs = n()) 

# join back w/hexcells
hex_cells_plot <- hex_cells_used %>% 
  left_join(hexcell_dev, by = "id_cells")

# bring in sf study region object for visualizing
country_sf <- rnaturalearth::ne_countries(country = c("United States of America","Canada", "Mexico"),returnclass = "sf")
western_country_sf <- st_crop(country_sf, xmin=-171.7911, xmax=-95, ymin=18, ymax=83.23324)
study_region_sf <- st_union(western_country_sf) # combines geoms
study_region_sf <- st_as_sf(st_cast(study_region_sf))  # deals with sinking multiple-shape file issue
study_region_sf <- st_transform(study_region_sf, crs = st_crs(hex_cells_used))


# plot 
ggplot(hex_cells_plot) + 
  geom_sf(study_region_sf, mapping = aes()) +
  geom_sf(aes(fill = mean_abs_dev), color = NA) + 
  scale_fill_viridis_c(name = "Mean Absolute\nDeviation (days)", option = "plasma") +
  theme_bw() + 
  labs(title = "Hexcell Deviation from Hopkins Law") + 
  theme(legend.position = "right")

## Species onset by elevation, faceted by lat
lat_breaks <- quantile(pheno_w_dev$lat_mean, c(0, 0.33, 0.67, 1), na.rm = TRUE)


# Function to plot individual species
plot_species_hopkins <- function(species_name, data) {
  sp_data <- data %>% filter(scientific_name == species_name)
  
  ggplot(sp_data, aes(x = mid_elev)) +
    geom_point(aes(y = onset), alpha = 0.5, size = 2) +  # Empirical data
    geom_smooth(aes(y = hopkins_lat_elev), method = "lm", 
                color = "red", se = TRUE, size = 1) +  # Smoothed Hopkins prediction
    facet_wrap(~latbin2, ncol = 3) +
    labs(title = species_name,
         x = "Elevation (m)",
         y = "Onset DOY") +
    theme_bw()

}

# Plot for a single species
plot_species_hopkins("Claytonia_lanceolata", pheno_w_dev)

ggplot() + geom_histogram(pheno_w_dev, mapping = aes(x = hopkins_dev))

# scale predictors, factor categorical vars
pheno_w_dev <- pheno_w_dev %>% 
  mutate(spmean_sc = scale(spmean), 
         elev_mid_sc = scale(mid_elev),
         # lat_mean_sc is already made
         obscount_sc = scale(obscount)) %>% 
  mutate(latbinf2 = factor(latbin2), 
         geoannperf3 = factor(geoannper3))




## Run mixed models on examining deviation in response to latitude, elevation, 
model_onset_dev <- lmer(hopkins_dev ~ elev_mid_sc + lat_mean_sc +
                             geoannperf3 + obscount_sc + 
                             (1 | scientific_name), 
                             data=pheno_w_dev)
model_onset_dev # assumes that there baseline deviation varies by species (intercepts by re), but elevation & latitude effects are consistent across species

model_onset_dev2 <- lmer(hopkins_dev ~ elev_mid_sc*geoannperf3 + 
                          lat_mean_sc*geoannperf3 +
                          obscount_sc*geoannperf3 + 
                          (1 | scientific_name), 
                        data=pheno_w_dev)
model_onset_dev2
lmerTest::step(model_onset_dev2)

# redo with best model 
model_onset_dev3 <- lmer(hopkins_dev ~ elev_mid_sc + 
                           geoannperf3 + 
                           lat_mean_sc + 
                           obscount_sc + 
                           (1 | scientific_name) + 
                           elev_mid_sc:geoannperf3 + 
                           geoannperf3:lat_mean_sc, 
                         data = pheno_w_dev)
model_onset_dev3 # Do elevation and latitude effects differ between herbaceous and perennial plants?

# allows each species to have different interactions with elevation and latitude
model_onset_dev4 <- lmer(hopkins_dev ~ elev_mid_sc + 
                           geoannperf3 + 
                           lat_mean_sc + 
                           obscount_sc + 
                           (1 + elev_mid_sc + lat_mean_sc | scientific_name) + 
                           elev_mid_sc:geoannperf3 + 
                           geoannperf3:lat_mean_sc, 
                         data = pheno_w_dev)


lmerTest::step(model_onset_dev4)
vcov_mat <- vcov(model_onset_dev4)
cov2cor(vcov_mat) # we've got some high correlations, check to see if there is an imbalance in the data
# check for imbalance in perennial vs annual
table(pheno_w_dev$geoannperf3, cut(pheno_w_dev$elev_mid_sc, breaks=3))
table(pheno_w_dev$geoannperf3, cut(pheno_w_dev$lat_mean_sc, breaks=3))

sjPlot::plot_model(model_onset_dev4, type = "int") # interaction plot
# elevation by growth form
sjPlot::plot_model(
  model_onset_dev4,
  type     = "eff",
  terms    = c("elev_mid_sc", "geoannperf3"),
  ci.lvl   = 0.95,
  dot.size = 0,        # no dots, just lines
  line.size = 0.8,
  axis.title = c("Scaled elevation", "Predicted onset deviation"),
  legend.title = "Growth form"
)
sjPlot::plot_model(
  model_onset_dev4,
  type     = "eff",
  terms    = c("lat_mean_sc", "geoannperf3"),
  ci.lvl   = 0.95,
  dot.size = 0,        # no dots, just lines
  line.size = 0.8,
  axis.title = c("Scaled latitude", "Predicted onset deviation"),
  legend.title = "Growth form"
)

sjPlot::plot_model(
  model_onset_dev4,
  type     = "eff",
  terms    = c("elev_mid_sc", "lat_mean_sc", "geoannperf3"),
  ci.lvl   = 0.95,
  dot.size = 0,        # no dots, just lines
  line.size = 0.8#,
  #axis.title = c("Scaled elevation", "Predicted onset deviation"),
  #legend.title = "Growth form"
)

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

# merge temp seasonality with pheno_w_hopkins_dev with hexcell and elevation bins 
pheno_w_dev2 <- pheno_w_dev %>% 
  left_join(temp_seasonality, by = c("bins6", "id_cells"))

## Run scaling 
pheno_w_dev2 <- pheno_w_dev2 %>% 
  mutate(lon_mean_sc = scale(lon_mean), 
         temp_seasonality_sc = scale(temp_seasonality), 
         elev_mid_sc = scale(mid_elev))

## Run a model 
model_temp_seasonality <- lm(temp_seasonality ~ elev_mid_sc + 
                               lat_mean_sc + lon_mean_sc,
                               data = pheno_w_dev2)

model_onset_dev4 <- lmer(hopkins_dev ~ elev_mid_sc + 
                           geoannperf3 + 
                           lat_mean_sc +
                           lon_mean_sc + 
                           obscount_sc + 
                           temp_seasonality_sc +
                           #(1 + elev_mid_sc + lat_mean_sc + lon_mean_sc | scientific_name) + 
                           (1 | scientific_name) + 
                           elev_mid_sc:geoannperf3 +
                           geoannperf3:lon_mean_sc +
                           geoannperf3:lat_mean_sc + 
                           temp_seasonality_sc, 
                         data = pheno_w_dev2)

lmerTest::step(model_onset_dev4)

dat <- pheno_w_dev2

dat <- dat %>% group_by(scientific_name) %>% mutate(inAllBins = all(c("NorthLat", "MidLat", "SouthLat") %in% latbinf2))
dat <- filter(dat, inAllBins == TRUE)
### Robs figure: 
dat$cell_bins6 <- interaction(dat$id_cells, dat$bins6, drop = TRUE)


m_full <- lmer(
  hopkins_dev ~ 
    ns(elev_mid_sc, 2) + ns(lat_mean_sc, 2) + temp_seasonality_sc + 
    ns(elev_mid_sc, 2):ns(lat_mean_sc, 2) + 
    ns(elev_mid_sc, 2):temp_seasonality_sc + 
    ns(elev_mid_sc, 2):geoannperf3 + 
    ns(lat_mean_sc, 2):temp_seasonality_sc + 
    (1 + elev_mid_sc + lat_mean_sc | scientific_name) + 
    (1 | id_cells) + (1| cell_bins6), 
  data = dat,
  weights = 1/sqrt(dat$obscount), 
  REML = TRUE
)

sjPlot::plot_model(
  m_full,
  type     = "eff",
  terms    = c("elev_mid_sc", "lat_mean_sc", "temp_seasonality_sc")
)

### Supplementary figs
test <- fread("/blue/guralnick/millerjared/PhenoElevation/data/pheno_w_hopkins_deviation_w_temp_seasonality.csv")

# Filter for Prunus virginiana and create seasonality bins

dat2 <- dat %>%
  mutate(seasonality_bin = case_when(
    temp_seasonality_sc < -0.5 ~ "Low Seasonality (temp_sc = -1)",
    temp_seasonality_sc >= -0.5 & temp_seasonality_sc < 0.5 ~ "Mid Seasonality (temp_sc = 0)",
    temp_seasonality_sc >= 0.5 ~ "High Seasonality (temp_sc = 1)",
    TRUE ~ NA_character_
  ) %>% factor(levels = c("Low Seasonality (temp_sc = -1)", 
                          "Mid Seasonality (temp_sc = 0)", 
                          "High Seasonality (temp_sc = 1)"))) %>% 
  filter(scientific_name == "Achillea_millefolium")

# Create the plot
ggplot() +
  geom_point(dat2, mapping = aes(x = elev, y = onset, color = latbinf2), alpha = 0.2) +
  geom_smooth(dat2, mapping = aes(x = elev, y = hopkins_lat_elev, color = latbinf2), method = "lm", se = TRUE, linewidth = 1) +
  facet_wrap(~ seasonality_bin, ncol = 3) +
  labs(
    title = "Achillea millefolium's empirical vs expected onset trends",
    x = "Elevation (m)",
    y = "Onset (Day of Year)",
    color = "Latitude Bin"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )
