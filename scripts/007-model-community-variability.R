### Title: Model Intraspecific Synchrony  
### Author: Rob Guralnick
### Editor: JT Miller 
### Date: 12-29-2025

### Purpose: load in intraspecific cv model ready data, fit mixed-effects glms, create visuals 
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

## Build community CV model
# load community cv data 
phenologycv_sf14_coefvar_sum4c <- readRDS("./outputs/rds-objects/cv-community-data.rds")

# fit an initial linear model 
model_comm_onset <- lm(log(cv_onsetphen) ~ (elev + latitude + longitude)^3   + count, data=phenologycv_sf14_coefvar_sum4c, na.action=na.fail)
model_comm_offset <- lm(log(cv_offsetphen) ~ (elev + latitude + longitude)^3   + count, data=phenologycv_sf14_coefvar_sum4c, na.action=na.fail)
model_comm_duration <- lm(log(cv_durphen) ~ (elev + latitude + longitude)^3   + count, data=phenologycv_sf14_coefvar_sum4c, na.action=na.fail)

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
model_comm_onset <- lm(log(cv_onsetphen) ~ count + elev + latitude + longitude + 
                         elev:latitude + elev:longitude + latitude:longitude + elev:latitude:longitude + 
                         1, data = phenologycv_sf14_coefvar_sum4c, na.action = na.fail)
model_comm_offset <- lm(log(cv_offsetphen) ~ count + elev + latitude + longitude + 
                          elev:latitude + elev:longitude + latitude:longitude + elev:latitude:longitude + 
                          1, data = phenologycv_sf14_coefvar_sum4c, na.action = na.fail)
model_comm_duration <- lm(log(cv_durphen) ~ count + elev + latitude + longitude + 
                            elev:latitude + elev:longitude + latitude:longitude + elev:latitude:longitude + 
                            1, data = phenologycv_sf14_coefvar_sum4c, na.action = na.fail)

# plot community synchrony models
comm_onset_plot <- plot_model(model_comm_onset, type="pred", terms=c("elev","latitude","longitude"))
comm_offset_plot <- plot_model(model_comm_offset, type="pred", terms=c("elev","latitude","longitude"))
comm_duration_plot <- plot_model(model_comm_duration, type="pred", terms=c("elev","latitude","longitude"))

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
       y = "Coefficient of Variation for Onset") +
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
       y = "Coefficient of Variation for Offset") +
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
       y = "Coefficient of Variation for Duration") +
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
ggsave("/blue/guralnick/millerjared/PhenoElevation/outputs/figures/cv-community-combined.png", height = 10, width = 25)

