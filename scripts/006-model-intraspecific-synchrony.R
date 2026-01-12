### Title: Model Intraspecific Synchrony  
### Author: Rob Guralnick
### Editor: JT Miller 
### Date: 12-29-2025

### Purpose: load in intraspecific cv model ready data, fit mixed-effects glms, create visuals 
setwd("/blue/guralnick/millerjared/PhenoElevation/") # set working-directory 

## Load Libraries 
library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(ggplot2)
library(calecopal)
library(sjPlot)
library(patchwork)

## Load intraspecific cv data
phenologycv_sf14_coefvar_sum_sp10 <- readRDS("./outputs/rds-objects/cv-intraspecific-data.rds")

## Build Infraspecific CV Model 
# fit an initial mixed-effects model
model_intra_onset <- lmer(log(cv_onsetphen) ~ (elev_sc + spmean_sc2 + latbinf2 + geoannperf3 )^3  # log-transform because of right skew on cv_onsetphen
                          + count  + (1|scientific_name), 
                          data=phenologycv_sf14_coefvar_sum_sp10)
model_intra_offset <- lmer(log(cv_offsetphen) ~ (elev_sc + spmean_sc2 + latbinf2 + geoannperf3 )^3  # log-transform because of right skew on cv_onsetphen
                           + count  + (1|scientific_name), 
                           data=phenologycv_sf14_coefvar_sum_sp10)
model_intra_duration <- lmer(log(cv_durphen) ~ (elev_sc + spmean_sc2 + latbinf2 + geoannperf3 )^3  # log-transform because of right skew on cv_onsetphen
                             + count  + (1|scientific_name), 
                             data=phenologycv_sf14_coefvar_sum_sp10)
# preform stepwise backwards reduction to obtain best model
step(model_intra_onset)
step(model_intra_offset)
step(model_intra_duration)

# take best model found during model selection, fit that 
model_intra_onsetf <- lmer(log(cv_onsetphen) ~ elev_sc + 
                             spmean_sc2 + latbinf2 + 
                             geoannperf3 + count + 
                             (1 | scientific_name) + 
                             elev_sc:latbinf2 + spmean_sc2:latbinf2 + latbinf2:geoannperf3, 
                           data=phenologycv_sf14_coefvar_sum_sp10)
model_intra_offsetf <- lmer(log(cv_offsetphen) ~ elev_sc + 
                              spmean_sc2 + latbinf2 + 
                              geoannperf3 + count + 
                              (1 | scientific_name) + 
                              spmean_sc2:latbinf2 + latbinf2:geoannperf3, 
                            data = phenologycv_sf14_coefvar_sum_sp10)
model_intra_durationf <- lmer(log(cv_durphen) ~ elev_sc + 
                                spmean_sc2 + latbinf2 + 
                                geoannperf3 + count + 
                                (1 | scientific_name) + 
                                elev_sc:spmean_sc2 + elev_sc:latbinf2 + spmean_sc2:latbinf2 + spmean_sc2:geoannperf3 + elev_sc:spmean_sc2:latbinf2, 
                              data = phenologycv_sf14_coefvar_sum_sp10)
# grab cols 
colors <- calecopal::cal_palette("superbloom2", n = 5) 
cols <- colors[c(5, 3, 1)]
# define labels for the facets
geoannperf3.labs <- c(anbien = "Annual & Biannual", herbper = "Herbaceous Perennial")
latbinf2.labs <- c(SouthLat = "Southern Latitudes", MidLat = "Middle Latitudes", NorthLat = "Northern Latitudes")

intra_onset_plot <- sjPlot::plot_model(model_intra_onsetf, type="pred", terms=c("elev_sc","spmean_sc2","geoannperf3","latbinf2"),  bias_correction = TRUE)

intra_onset_plot <- intra_onset_plot + 
  scale_color_manual(values = cols, labels = c("Early", "Middle", "Late")) + 
  scale_fill_manual(values = cols, labels = c("Early", "Middle", "Late")) +
  labs(color = "Seasonal Average\nFlowering Period", 
       fill = "Seasonal Average\nFlowering Period", 
       x = "Scaled Elevation", 
       y = "Coefficient of Variation for Onset") +
  ggtitle("Synchrony of Onset for Intraspecific Flowering Phenology") +
  theme_bw() + 
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 20, hjust = 0.5), 
    strip.text.x = element_text(size =12),
    strip.text.y = element_text(size = 16)
  ) +
  theme(plot.title = element_text(hjust = 0.5)) 

intra_offset_plot <- sjPlot::plot_model(model_intra_offsetf, type="pred", terms=c("elev_sc","spmean_sc2","geoannperf3","latbinf2"),  bias_correction = TRUE)

intra_offset_plot <- intra_offset_plot + 
  scale_color_manual(values = cols, labels = c("Early", "Middle", "Late")) + 
  scale_fill_manual(values = cols, labels = c("Early", "Middle", "Late")) +
  labs(color = "Seasonal Average\nFlowering Period", 
       fill = "Seasonal Average\nFlowering Period", 
       x = "Scaled Elevation", 
       y = "Coefficient of Variation for Offset") +
  ggtitle("Synchrony of Offset for Intraspecific Flowering Phenology") +
  theme_bw() + 
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 20, hjust = 0.5), 
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 16)
  ) +
  theme(plot.title = element_text(hjust = 0.5)) 

intra_duration_plot <- sjPlot::plot_model(model_intra_durationf, type="pred", terms=c("elev_sc","spmean_sc2","geoannperf3","latbinf2"),  bias_correction = TRUE)

intra_duration_plot <- intra_duration_plot + 
  scale_color_manual(values = cols, labels = c("Early", "Middle", "Late")) + 
  scale_fill_manual(values = cols, labels = c("Early", "Middle", "Late")) +
  labs(color = "Seasonal Average\nFlowering Period", 
       fill = "Seasonal Average\nFlowering Period", 
       x = "Scaled Elevation", 
       y = "Coefficient of Variation for Duration") +
  ggtitle("Synchrony of Duration for Intraspecific Flowering Phenology") +
  theme_bw() + 
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 20, hjust = 0.5), 
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 16)
  ) +
  theme(plot.title = element_text(hjust = 0.5)) 


# Remove legends from the first two plots
intra_onset_plot_nolegend <- intra_onset_plot + ggtitle("Onset") + theme(legend.position = "none")
intra_offset_plot_nolegend <- intra_offset_plot + ggtitle("Offset") + theme(legend.position = "none", ylab = "none")
intra_duration_plot_legend <- intra_duration_plot + ggtitle("Duration") # Keep legend here

# Combine plots
combined_plot <- (intra_onset_plot_nolegend | intra_offset_plot_nolegend | intra_duration_plot_legend) +
  plot_annotation(title = "Synchrony for Intraspecific Flowering Phenology", 
                  theme = theme(plot.title = element_text(size = 20)))

print(combined_plot)
ggsave("/blue/guralnick/millerjared/PhenoElevation/outputs/figures/cv-intraspecific-combined.png", height = 10, width = 18)
