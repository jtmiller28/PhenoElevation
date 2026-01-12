### Title: Additional Figures
### Author: JT Miller
### Date: 12-29-2025

### Purpose: Wrap up some additional conceptual figures for the paper
setwd("/blue/guralnick/millerjared/PhenoElevation/") # set working-directory 

## Load Libraries 
library(tidyverse)
library(sf)
library(data.table)
library(calecopal)
library(effects)
library(patchwork)
library(ggspatial)

## Study region with latitude subdivisions
# read in full study region
hex_grid_w_ids <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/full-study-region.rds")
# read in pglmm model ready data to filter down hexcells
phenologycv_sf14_onoffall20 <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-model-ready-data.rds")
# filter down hexcells to only include those relevant 
hex_grids_used <- hex_grid_w_ids %>% 
  dplyr::filter(id_cells %in% phenologycv_sf14_onoffall20$id_cells)

# attach subregion delimitations to these hexcells
grid_cent_ll_point <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/grid_cent_ll_point.rds")
grid_cent_ll_point$latbin2 <- cut(grid_cent_ll_point$latitude, breaks = 3, labels = c("SouthLat", "MidLat", "NorthLat"))
grid_cent_ll_point4 <- grid_cent_ll_point %>% select(id_cells, latbin2)
hex_grids_used <- hex_grids_used %>% 
  left_join(grid_cent_ll_point4, by = "id_cells")

# load in a shapefile of North America
country_sf <- rnaturalearth::ne_countries(country = c("United States of America","Canada", "Mexico"),returnclass = "sf")
western_country_sf <- st_crop(country_sf, xmin=-171.7911, xmax=-95, ymin=18, ymax=83.23324)
study_region_sf <- st_union(western_country_sf) # combines geoms
study_region_sf <- st_as_sf(st_cast(study_region_sf))  # deals with sinking multiple-shape file issue

study_region_sf <- st_transform(study_region_sf, crs = st_crs(hex_grids_used))

## Plot
# create cols
n_bands <- 20
band_colors <- cal_palette("sierra1", n = n_bands, type = "continuous")
cols <- band_colors[c(1, 8, 20)]

# change naming
hex_grids_used <- hex_grids_used %>%
  mutate(latbin_label = case_when(
    latbin2 == "SouthLat" ~ "Southern",
    latbin2 == "MidLat"  ~ "Middle",
    latbin2 == "NorthLat" ~ "Northern",
    TRUE ~ NA_character_
  ), 
  latbin_label = factor(
    latbin_label,
    levels = c(
      "Southern",
      "Middle",
      "Northern"
    )))

ggplot() + 
  geom_sf(study_region_sf, mapping = aes()) + 
  geom_sf(hex_grids_used, mapping = aes(fill = as.factor(latbin_label), color = as.factor(latbin_label)), alpha = 0.2) + 
  scale_fill_manual(values = cols, guide = guide_legend(title = "Relative Latitude\nof Study Region")) +
  scale_color_manual(values = cols, guide = "none") + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.10, "in"), pad_y = unit(0.05, "in"), 
                         style = north_arrow_fancy_orienteering()) +
  labs(x = "Longitude", 
       y = "Latitude", 
       color = "Relative Latitude\nof Study Region", 
       fill = "Relative Latitude\nof Study Region", 
       title = "Study Region with Hexagonal Cells") +
  theme_bw() + 
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", linewidth = 0.5), panel.background = element_rect(fill = "transparent")) +
  guides(color = guide_legend(reverse = TRUE),
         fill  = guide_legend(reverse = TRUE)) +
  theme( 
        plot.title = element_text(size = 22, hjust = 0.5),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18), 
        axis.title = element_text(size = 18),
        axis.text  = element_text(size = 16),
        )

ggsave("./outputs/figures/study-region.png", height = 10, width = 10)

## Conceptual Mountain with Elevational Bins

library(dplyr)
library(ggplot2)
library(purrr)
library(calecopal)

# re-parameterize positions and widths
mu_left   <- 1 + 9;  sigma_left   <- 2.2  # narrower, slightly to the left
mu_middle <- 5 + 9; sigma_middle <- 2.8  # tallest, pushed right
mu_right  <- 9 + 9; sigma_right  <- 3.2  # in-between, pushed right

mtn_left   <- data.frame(x = x, y = 1800 * exp(-((x - mu_left  )^2) / (2 * sigma_left^2)))
mtn_middle <- data.frame(x = x, y = 3000 * exp(-((x - mu_middle)^2) / (2 * sigma_middle^2)))
mtn_right  <- data.frame(x = x, y = 2400 * exp(-((x - mu_right )^2) / (2 * sigma_right^2)))

# 150 m elevation bands up to tallest peak
band_breaks <- seq(0, 3000, by = 150)
band_labels <- paste0(band_breaks[-length(band_breaks)], "–", band_breaks[-1], " m")
n_bands <- length(band_labels)
band_colors <- calecopal::cal_palette("sierra1", n = n_bands, type = "continuous")

build_ribbons <- function(mtn) {
  map_df(seq_len(n_bands), function(i) {
    mtn %>%
      mutate(
        ymin = band_breaks[i],
        ymax = pmin(y, band_breaks[i + 1]),
        band = factor(band_labels[i], levels = band_labels)
      ) %>%
      filter(ymax > ymin)
  })
}

ribbons_middle <- build_ribbons(mtn_middle)
ribbons_right  <- build_ribbons(mtn_right)
ribbons_left   <- build_ribbons(mtn_left)

p <- ggplot() +
  geom_ribbon(data = ribbons_middle, aes(x = x, ymin = ymin, ymax = ymax, fill = band), alpha = 1) +
  geom_line(  data = mtn_middle,    aes(x = x, y = y), color = "black", size = 1.1) +
  geom_ribbon(data = ribbons_right, aes(x = x, ymin = ymin, ymax = ymax, fill = band), alpha = 1) +
  geom_line(  data = mtn_right,     aes(x = x, y = y), color = "black", size = 1.1) +
  geom_ribbon(data = ribbons_left,  aes(x = x, ymin = ymin, ymax = ymax, fill = band), alpha = 1) +
  geom_line(  data = mtn_left,      aes(x = x, y = y), color = "black", size = 1.1) +
  scale_fill_manual(values = band_colors, name = "Elevation band") +
  coord_cartesian(xlim = c(0, 30), ylim = c(0, 3050), expand = FALSE, clip = "on") +
  theme_void() +
  ggtitle("Conceptual Diagram of Elevation Bands Across Sample Hexagonal Cell") +
  theme(
    plot.title = element_text(size = 30, hjust = 0.5),
    legend.position = c(0.95, 0.85),         # move higher, inside top-right
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA), 
    legend.title = element_text(size = 27),
    legend.text = element_text(size = 25)
  ) +
  guides(color = guide_legend(reverse = TRUE),
         fill  = guide_legend(reverse = TRUE))
  
#p <- p + coord_cartesian(xlim = c(0, 100), ylim = c(0, 3050), expand = FALSE)
print(p)

ggsave("/blue/guralnick/millerjared/PhenoElevation/outputs/figures/conceptual-mountain-bins.png", height = 10, width = 20)




