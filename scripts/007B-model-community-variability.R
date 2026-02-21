### Title: Model Intraspecific Synchrony  
### Author: Rob Guralnick
### Editor: JT Miller 
### Date: 12-29-2025

### Purpose: load in intraspecific iqr model ready data, fit mixed-effects glms, create visuals 
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

## Build community iqr model
# load community iqr data 
phenologyiqr_sf14_coefvar_sum4c <- readRDS("./outputs/rds-objects/iqr-community-data.rds")

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

# merge temp seasonality with community iqrs with hexcell and elevation bins 
phenologyiqr_sf14_coefvar_sum4c  <- phenologyiqr_sf14_coefvar_sum4c  %>% 
  left_join(temp_seasonality, by = c("bins6", "id_cells"))

# scale predictors 
phenologyiqr_sf14_coefvar_sum4c$elev_sc <- c(scale(phenologyiqr_sf14_coefvar_sum4c$elev))
phenologyiqr_sf14_coefvar_sum4c$lat_sc <- c(scale(phenologyiqr_sf14_coefvar_sum4c$latitude))
phenologyiqr_sf14_coefvar_sum4c$count_sc <- c(scale(phenologyiqr_sf14_coefvar_sum4c$count))
phenologyiqr_sf14_coefvar_sum4c$temp_seasonality_sc <- c(scale(phenologyiqr_sf14_coefvar_sum4c$temp_seasonality))
# log responses
phenologyiqr_sf14_coefvar_sum4c <- phenologyiqr_sf14_coefvar_sum4c %>% 
  mutate(log_iqr_onsetphen = log(iqr_onsetphen), 
         log_iqr_offsetphen = log(iqr_offsetphen), 
         log_iqr_durphen = log(iqr_durphen))


# fit an initial linear model 
model_comm_onset <- lm(log_iqr_onsetphen ~ (elev_sc + lat_sc)^3 + temp_seasonality_sc   + count_sc, data=phenologyiqr_sf14_coefvar_sum4c, na.action=na.fail)
model_comm_offset <- lm(log_iqr_offsetphen ~ (elev_sc + lat_sc)^3  + temp_seasonality_sc + count_sc, data=phenologyiqr_sf14_coefvar_sum4c, na.action=na.fail)
model_comm_duration <- lm(log_iqr_durphen ~ (elev_sc + lat_sc)^3 + temp_seasonality_sc  + count_sc, data=phenologyiqr_sf14_coefvar_sum4c, na.action=na.fail)

# dredge for model selection
model_selection_onset <- MuMIn::dredge(model_comm_onset)
model_selection_offset <- MuMIn::dredge(model_comm_offset)
model_selection_duration <- MuMIn::dredge(model_comm_duration)

# grab best model from dredge
best_model_onset <- MuMIn::get.models(model_selection_onset, subset = 1)[[1]] 
best_model_offset <- MuMIn::get.models(model_selection_offset, subset = 1)[[1]] 
best_model_duration <- MuMIn::get.models(model_selection_duration, subset = 1)[[1]] 

# call best model
best_model_onset$call
best_model_offset$call
best_model_duration$call


# call in best model and fit 
model_comm_onset <- lm(log_iqr_onsetphen ~ count_sc + elev_sc + lat_sc + 
                          temp_seasonality_sc + elev_sc:lat_sc + 1,
                       data = phenologyiqr_sf14_coefvar_sum4c, na.action = na.fail)
model_comm_offset <- lm(log_iqr_offsetphen ~ count_sc + elev_sc + lat_sc + 
                          temp_seasonality_sc + 1,
                       data = phenologyiqr_sf14_coefvar_sum4c, na.action = na.fail)
model_comm_duration <- lm(log_iqr_durphen ~ count_sc + elev_sc + lat_sc + 
                          temp_seasonality_sc + 1, 
                       data = phenologyiqr_sf14_coefvar_sum4c, na.action = na.fail)

# plot community synchrony models
comm_onset_plot <- plot_model(model_comm_onset, type="pred", terms=c("elev_sc","lat_sc","temp_seasonality_sc"))
comm_offset_plot <- plot_model(model_comm_offset, type="pred", terms=c("elev_sc","lat_sc","temp_seasonality_sc"))
comm_duration_plot <- plot_model(model_comm_duration, type="pred", terms=c("elev_sc","lat_sc","temp_seasonality_sc"))

# adjust naming of columns denoting the Western-Eastern Regions
# mapping from sjPlot's default facet strings to your labels
lon_lab <- c(
  `longitude = -121.26` = "Western Study Region",
  `longitude = -113.54` = "Middle Study Region",
  `longitude = -105.83` = "Eastern Study Region"
)

# reapply facet with custom labeller using the correct faceting variable
apply_labeller <- function(p, mapping) {
  fac_var <- names(p$facet$params$facets)[1]  # typically "facet"
  lab <- setNames(list(as_labeller(mapping)), fac_var)
  p + facet_wrap(as.formula(paste("~", fac_var)),
                 labeller = do.call(labeller, lab))
}

comm_onset_plot    <- apply_labeller(comm_onset_plot, lon_lab)
comm_offset_plot   <- apply_labeller(comm_offset_plot, lon_lab)
comm_duration_plot <- apply_labeller(comm_duration_plot, lon_lab)

# grab cols 
n_bands <- 20
band_colors <- cal_palette("sierra1", n = n_bands, type = "continuous")
cols <- band_colors[c(1, 8, 20)]

# adjust labels to make this clearer 
comm_onset_plot <- comm_onset_plot + 
  scale_color_manual(values = cols, labels = c("Southern", "Middle", "Northern")) + 
  scale_fill_manual(values = cols, labels = c("Southern", "Middle", "Northern")) +
  labs(color = "Relative Latitude\nof Study Region", 
       fill = "Relative Latitude\nof Study Region", 
       x = "Elevation (m)", 
       y = "Interquartile Range for Onset") +
  ggtitle("Synchrony of Onset for Community Flowering Phenology") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 10),
    plot.title = element_text(size = 25, hjust = 0.5), 
    strip.text.x = element_text(size = 12)
  )
print(comm_onset_plot)

comm_offset_plot <- comm_offset_plot + 
  scale_color_manual(values = cols, labels = c("Southern", "Middle", "Northern")) + 
  scale_fill_manual(values = cols, labels = c("Southern", "Middle", "Northern")) +
  labs(color = "Relative Latitude\nof Study Region", 
       fill = "Relative Latitude\nof Study Region", 
       x = "Elevation (m)", 
       y = "Interquartile Range for Offset") +
  ggtitle("Synchrony of Offset for Community Flowering Phenology") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 10),
    plot.title = element_text(size = 25, hjust = 0.5), 
    strip.text.x = element_text(size = 12)
  )
print(comm_offset_plot)

comm_duration_plot <- comm_duration_plot + 
  scale_color_manual(values = cols, labels = c("Southern", "Middle", "Northern")) + 
  scale_fill_manual(values = cols, labels = c("Southern", "Middle", "Northern")) +
  labs(color = "Relative Latitude\nof Study Region", 
       fill = "Relative Latitude\nof Study Region", 
       x = "Elevation (m)", 
       y = "Interquartile Range for Duration") +
  ggtitle("Synchrony of Duration for Community Flowering Phenology") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 10),
    plot.title = element_text(size = 25, hjust = 0.5), 
    strip.text.x = element_text(size = 12)
  )
print(comm_duration_plot)

# reverse order of legend 
comm_onset_plot  <- comm_onset_plot  + guides(color = guide_legend(reverse = TRUE),
                                              fill  = guide_legend(reverse = TRUE))
comm_offset_plot <- comm_offset_plot + guides(color = guide_legend(reverse = TRUE),
                                              fill  = guide_legend(reverse = TRUE))
comm_duration_plot <- comm_duration_plot + guides(color = guide_legend(reverse = TRUE),
                                                  fill  = guide_legend(reverse = TRUE))

# Remove legends from the first two plots
comm_onset_plot_nolegend <- comm_onset_plot + ggtitle("Onset") + theme(legend.position = "none")
comm_offset_plot_nolegend <- comm_offset_plot + ggtitle("Offset") + theme(legend.position = "none", ylab = "none")
comm_duration_plot_legend <- comm_duration_plot + ggtitle("Duration") # Keep legend here

# Combine plots
comm_combined_plot <- (comm_onset_plot_nolegend | comm_offset_plot_nolegend | comm_duration_plot_legend) +
  plot_annotation(title = "Synchrony for Community Flowering Phenology", 
                  theme = theme(plot.title = element_text(size = 20))) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(comm_combined_plot)
ggsave("/blue/guralnick/millerjared/PhenoElevation/outputs/figures/iqr-community-combined.png", height = 10, width = 25)

