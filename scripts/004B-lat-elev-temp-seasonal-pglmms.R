### Title: Latitude Elevation PGLMMs   
### Author: Rob Guralnick
### Editor: JT Miller 
### Date: 12-29-2025

### Purpose: load in phenometric calcs model ready data, scale, run glms and pglmms
setwd("/blue/guralnick/millerjared/PhenoElevation/") # set working-directory 

## Load Libraries
library(lme4)
library(lmerTest)
library(performance)
library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
library(calecopal)
library(rtrees)
library(ape)
library(stringr)
library(posterior)
library(tidybayes)
library(patchwork)

## Load in phenometric calcs model ready data
phenologycv_sf14_onoffall20 <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/pheno-model-ready-data.rds")


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

# merge temp seasonality with pheno estimates with hexcell and elevation bins 
phenologycv_sf14_onoffall20 <- phenologycv_sf14_onoffall20 %>% 
        left_join(temp_seasonality, by = c("bins6", "id_cells"))

## Scale and create factors 
phenologycv_sf14_onoffall21 <- phenologycv_sf14_onoffall20 %>%
  dplyr::mutate(latitude_sc =scale(latitude), 
                longitude_sc=scale(longitude),
                elev_sc = scale(elev),
                spmean_sc = scale(spmean),
                temp_seasonality_sc = scale(temp_seasonality),
                geophyte2f = as.factor(geophyte2),
                geoannperf = as.factor(geoannper),
                geoannper2f = as.factor(geopannper2),
                geoannper3f = as.factor(geoannper3),
                #mean_elev_binf = as.factor(mean_elev_bin), defunct atm
                timebinf = as.factor(timebin),
                timebinf2 = as.factor(timebin2),
                annpertimef = as.factor(annpertime),
                id_cellsf = as.factor(id_cells),
                yearf = as.factor(year))

## Use a generalized linear mixed effect model to set-up model selection
general_trends <- phenologycv_sf14_onoffall21 %>% 
  mutate(elev_sc = scale(elev), 
         latitude_sc = scale(latitude),
         longitude_sc = scale(longitude), 
         id_cells = as.factor(id_cells),
         geoannper3f = as.factor(geoannper3))

mod_onset <- lmer(onset ~ elev_sc + 
                    latitude_sc +
                    temp_seasonality_sc +
                    geoannper3f + 
                    year + 
                    (1 | scientific_name) + 
                    (1 | id_cells) +
                    elev_sc:latitude_sc + 
                    elev_sc:geoannper3f + 
                    latitude_sc:geoannper3f + 
                    latitude_sc:temp_seasonality_sc +
                    elev_sc:temp_seasonality_sc,                          
                  data = general_trends) 
mod_offset <- lmer(offset ~ elev_sc + 
                     latitude_sc +
                     temp_seasonality_sc +
                     geoannper3f + 
                     year + 
                     (1 | scientific_name) + 
                     (1 | id_cells) +
                     elev_sc:latitude_sc + 
                     elev_sc:geoannper3f + 
                     latitude_sc:geoannper3f + 
                     latitude_sc:temp_seasonality_sc +
                     elev_sc:temp_seasonality_sc,                         
                   data = general_trends ) 
mod_duration <- lmer(duration ~ elev_sc + 
                       latitude_sc +
                       temp_seasonality_sc +
                       geoannper3f + 
                       year + 
                       (1 | scientific_name) + 
                       (1 | id_cells) +
                       elev_sc:latitude_sc + 
                       elev_sc:geoannper3f + 
                       latitude_sc:geoannper3f + 
                       latitude_sc:temp_seasonality_sc +
                       elev_sc:temp_seasonality_sc,
                     data = general_trends ) 

## use model selection on these mixed-effects models
print("Check onset glmms model selection")
step(mod_onset)
print("Check offset glmms model selection")
step(mod_offset)
print("Check duration glmms model selection")
step(mod_duration)

## check preformance and inflation
print("Check onset glmms inflation")
performance::check_collinearity(mod_onset) # low 
print("Check offset glmms inflation")
performance::check_collinearity(mod_offset) # low
print("Check duration glmms inflation")
performance::check_collinearity(mod_duration) # low

## visually evaluate whether transforming the data would be helpful
ggplot() + geom_histogram(general_trends, mapping = aes(x = onset))
ggplot() + geom_histogram(general_trends, mapping = aes(x = duration))
ggplot() + geom_histogram(general_trends, mapping = aes(x = offset))
ggplot() + geom_histogram(general_trends, mapping = aes(x = onset^2)) # no
ggplot() + geom_histogram(general_trends, mapping = aes(x = log(duration))) # yes
ggplot() + geom_histogram(general_trends, mapping = aes(x = log(offset))) # no


## Apply transformations where applicable 
general_trends <- general_trends %>% 
  mutate(log_duration = log(duration))

## Rerun model for transformations + check model selection and inflation
mod_duration <- lmer(log_duration ~ elev_sc + 
                       latitude_sc +
                       temp_seasonality_sc +
                       geoannper3f + 
                       year + 
                       (1 | scientific_name) + 
                       (1 | id_cells) +
                       elev_sc:latitude_sc + 
                       elev_sc:geoannper3f + 
                       latitude_sc:geoannper3f + 
                       latitude_sc:temp_seasonality_sc +
                       elev_sc:temp_seasonality_sc,
                     data = general_trends ) 
print("Check log-transformed duration glmm model selection")
step(mod_duration)
print("Check log-transformed duration glmm inflation")
performance::check_collinearity(mod_duration)

## Read in traits, derive species list, gather up a phylogeny for these taxa
SpeciesTraits3 <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/species-traits-for-filtered-taxa.rds")
sp_list <- subset(SpeciesTraits3, select=c('species','genus','family'))
sp_list$scientific_name <- str_to_title(sp_list$species)
sp_list2 <- sp_list %>% filter(sp_list$scientific_name %in% unique(phenologycv_sf14_onoffall20$scientific_name))
plant_tree1 <- rtrees::get_tree(sp_list = sp_list2, taxon = "plant", scenario = "at_basal_node", show_grafted = FALSE)
A <- ape::vcv.phylo(plant_tree1)
# remove data for names in general trends that fail to show up in the phylogeny 
levels_names_in_data <- sort(unique(general_trends$scientific_name))
levels_in_A <- sort(rownames(A))
names_missing <- setdiff(levels_names_in_data, levels_in_A)
general_trends2 <- general_trends %>% 
  filter(!scientific_name %in% names_missing)
## Use best models founds during mixed-effects modeling to fit a phylogenetic glmm using brms. 
## Note that the best model was the same for onset, offset, and log duration. 
## We will drop one mixed effect of latitude to plant habit as its not directly relevant to our question and simplifies the model convergence.
brm_onset <- brm(
  onset ~ elev_sc + 
    latitude_sc + 
    temp_seasonality_sc +
    geoannper3f +
    elev_sc:latitude_sc + 
    elev_sc:geoannper3f + 
    elev_sc:temp_seasonality_sc +
    # latitude_sc:geoannper3f + # dropping to increase bulk sampling density...
    (1|gr(scientific_name, cov = A)) + 
    (1 | id_cells) , 
  data = general_trends2, 
  data2 = list(A = A), 
  family = gaussian(), 
  chains = 4, cores = 4, iter = 15000, control = list(adapt_delta = 0.99, max_treedepth = 15)
)
saveRDS(brm_onset, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/onset-brms-phylo-model-w-temp-seasonality.rds") # 12 hours est 
brm_offset <- brm(
  offset ~ elev_sc + 
    latitude_sc + 
    temp_seasonality_sc +
    geoannper3f +
    elev_sc:latitude_sc + 
    elev_sc:geoannper3f + 
    elev_sc:temp_seasonality_sc +
    # latitude_sc:geoannper3f + # dropping to increase bulk sampling density...
    (1|gr(scientific_name, cov = A)) + 
    (1 | id_cells) ,
  data = general_trends2, 
  data2 = list(A = A), 
  family = gaussian(), 
  chains = 4, cores = 4, iter = 15000, control = list(adapt_delta = 0.99, max_treedepth = 15)
)
saveRDS(brm_offset, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/offset-brms-phylo-model-w-temp-seasonality.rds")
brm_duration <- brm(
  log_duration ~ elev_sc + 
    latitude_sc + 
    temp_seasonality_sc +
    geoannper3f +
    elev_sc:latitude_sc + 
    elev_sc:geoannper3f + 
    elev_sc:temp_seasonality_sc +
    # latitude_sc:geoannper3f + # dropping to increase bulk sampling density...
    (1|gr(scientific_name, cov = A)) + 
    (1 | id_cells) ,
  data = general_trends2, 
  data2 = list(A = A), 
  family = gaussian(), 
  chains = 4, cores = 4, iter = 15000, control = list(adapt_delta = 0.99, max_treedepth = 15)
)
saveRDS(brm_duration, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/duration-brms-phylo-model-w-temp-seasonality.rds")

## Assess model fit 
# check overall model summary
print("Check onset pglmms diagnostics")
summary(brm_onset) # diagnositics look good, Rhats are all 1.00, chains converged well, ESS are relatively high. 
# check posterior predictive check
print("Check offset pglmms diagnostics")
pp_check(brm_onset) # generally looks pretty good for goodness of fit
summary(brm_offset) # diagnositics look good, Rhats are all 1.00, chains converged well, ESS are relatively high. 
# check posterior predictive check
pp_check(brm_offset) # generally looks pretty good for goodness of fit
print("Check duration pglmms diagnostics")
summary(brm_duration) # diagnositics look good, Rhats are all 1.00, chains converged well, ESS are relatively high. 
# check posterior predictive check
pp_check(brm_duration) # generally looks pretty good for goodness of fit

## Read in models (redundant unless your saving above modeling time)
brm_onset <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/onset-brms-phylo-model-w-temp-seasonality.rds")
brm_offset <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/offset-brms-phylo-model-w-temp-seasonality.rds")
brm_duration <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/duration-brms-phylo-model-w-temp-seasonality.rds")

## Visualize the posterior estimates for latitude and elevation interaction + elevation
# construct new data for predictions
lat_vals <- quantile(general_trends2$latitude_sc, c(0.1, 0.5, 0.9)) # low, median, and high latitude
temp_seas_vals <- quantile(general_trends2$temp_seasonality_sc, c(0.1, 0.5, 0.9)) # low, median, high temp seasonality
elev_seq <- seq(from = min(general_trends2$elev_sc),
                to = max(general_trends2$elev_sc),
                length.out = 80)

newdf <- expand.grid(
  latitude_sc = lat_vals,
  elev_sc = elev_seq,
  temp_seasonality_sc = temp_seas_vals,
  geoannper3f = c("anbien", "herbper"),
  scientific_name = NA,
  id_cells = NA
)

epreds_onset <- add_epred_draws(brm_onset, newdata = newdf, value = "onset_pred", re_formula = ~0)
epreds_offset <- add_epred_draws(brm_offset, newdata = newdf, value = "offset_pred", re_formula = ~0)
epreds_duration <- add_epred_draws(brm_duration, newdata = newdf, value = "duration_pred", re_formula = ~0)

plot_df_onset <- epreds_onset %>%
  group_by(elev_sc, latitude_sc, temp_seasonality_sc) %>%
  mean_qi(onset_pred, .width = 0.89) %>%
  mutate(
    lat_label = case_when(
      dplyr::near(latitude_sc, lat_vals[1]) ~ "Southern (10th percentile)",
      dplyr::near(latitude_sc, lat_vals[2]) ~ "Middle (50th percentile)",
      dplyr::near(latitude_sc, lat_vals[3]) ~ "Northern (90th percentile)",
      TRUE ~ NA_character_
    ),
    lat_label = factor(
      lat_label,
      levels = c(
        "Northern (90th percentile)",
        "Middle (50th percentile)",
        "Southern (10th percentile)"
      )
    ), 
    temp_seas_label = case_when(
      dplyr::near(temp_seasonality_sc, temp_seas_vals[1]) ~ "Low Temp Seasonality (10th percentile)",
      dplyr::near(temp_seasonality_sc, temp_seas_vals[2]) ~ "Medium Temp Seasonality (50th percentile)",
      dplyr::near(temp_seasonality_sc, temp_seas_vals[3]) ~ "High Temp Seasonality (90th percentile)",
      TRUE ~ NA_character_
    ),
    temp_seas_label = factor(
      temp_seas_label,
      levels = c(
        "Low Temp Seasonality (10th percentile)",
        "Medium Temp Seasonality (50th percentile)",
        "High Temp Seasonality (90th percentile)"
      )
    )
  )
plot_df_offset <- epreds_offset %>%
  group_by(elev_sc, latitude_sc, temp_seasonality_sc) %>%
  mean_qi(offset_pred, .width = 0.89) %>%
  mutate(
    lat_label = case_when(
      dplyr::near(latitude_sc, lat_vals[1]) ~ "Southern (10th percentile)",
      dplyr::near(latitude_sc, lat_vals[2]) ~ "Middle (50th percentile)",
      dplyr::near(latitude_sc, lat_vals[3]) ~ "Northern (90th percentile)",
      TRUE ~ NA_character_
    ),
    lat_label = factor(
      lat_label,
      levels = c(
        "Northern (90th percentile)",
        "Middle (50th percentile)",
        "Southern (10th percentile)"
      )
    ), 
    temp_seas_label = case_when(
      dplyr::near(temp_seasonality_sc, temp_seas_vals[1]) ~ "Low Temp Seasonality (10th percentile)",
      dplyr::near(temp_seasonality_sc, temp_seas_vals[2]) ~ "Medium Temp Seasonality (50th percentile)",
      dplyr::near(temp_seasonality_sc, temp_seas_vals[3]) ~ "High Temp Seasonality (90th percentile)",
      TRUE ~ NA_character_
    ),
    temp_seas_label = factor(
      temp_seas_label,
      levels = c(
        "Low Temp Seasonality (10th percentile)",
        "Medium Temp Seasonality (50th percentile)",
        "High Temp Seasonality (90th percentile)"
      )
    )
  )
plot_df_duration <- epreds_duration %>%
  group_by(elev_sc, latitude_sc, temp_seasonality_sc) %>%
  mean_qi(duration_pred, .width = 0.89) %>%
  mutate(
    lat_label = case_when(
      dplyr::near(latitude_sc, lat_vals[1]) ~ "Southern (10th percentile)",
      dplyr::near(latitude_sc, lat_vals[2]) ~ "Middle (50th percentile)",
      dplyr::near(latitude_sc, lat_vals[3]) ~ "Northern (90th percentile)",
      TRUE ~ NA_character_
    ),
    lat_label = factor(
      lat_label,
      levels = c(
        "Northern (90th percentile)",
        "Middle (50th percentile)",
        "Southern (10th percentile)"
      )
    ), 
    temp_seas_label = case_when(
      dplyr::near(temp_seasonality_sc, temp_seas_vals[1]) ~ "Low Temp Seasonality (10th percentile)",
      dplyr::near(temp_seasonality_sc, temp_seas_vals[2]) ~ "Medium Temp Seasonality (50th percentile)",
      dplyr::near(temp_seasonality_sc, temp_seas_vals[3]) ~ "High Temp Seasonality (90th percentile)",
      TRUE ~ NA_character_
    ),
    temp_seas_label = factor(
      temp_seas_label,
      levels = c(
        "Low Temp Seasonality (10th percentile)",
        "Medium Temp Seasonality (50th percentile)",
        "High Temp Seasonality (90th percentile)"
      )
    ),
    # transform out of log
    duration_pred_unlogged = exp(duration_pred),
    .lower_unlogged        = exp(.lower),
    .upper_unlogged        = exp(.upper)
  )

  

## Plot
n_bands <- 20
band_colors <- cal_palette("sierra1", n = n_bands, type = "continuous")
cols <- band_colors[c(20, 8, 1)]

onset_elev_lat_plot <- ggplot(plot_df_onset, aes(x = elev_sc, y = onset_pred, color = as.factor(lat_label), fill = as.factor(lat_label))) +
  scale_fill_manual(values = cols, guide = guide_legend(title = "Relative Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, linetype = 0) +
  facet_wrap(~temp_seas_label) +
  labs(
    x = "Scaled Elevation",
    y = "Predicted Onset (posterior mean)",
    color = "Elevation\n(quantiles)",
    fill = "Elevation\n(quantiles)",
    title = "Predicted Onset across Latitude and Elevation"
  ) +
  theme_bw() + 
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 24, hjust = 0.5)
  )

print(onset_elev_lat_plot)

offset_elev_lat_plot <- ggplot(plot_df_offset, aes(x = elev_sc, y = offset_pred, color = as.factor(lat_label), fill = as.factor(lat_label))) +
  scale_fill_manual(values = cols, guide = guide_legend(title = "Relative Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, linetype = 0) +
  facet_wrap(~temp_seas_label) +
  labs(
    x = "Scaled Elevation",
    y = "Predicted offset (posterior mean)",
    color = "Elevation\n(quantiles)",
    fill = "Elevation\n(quantiles)",
    title = "Predicted Offset across Latitude and Elevation"
  ) +
  theme_bw() + 
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 24, hjust = 0.5)
  )

print(offset_elev_lat_plot)


duration_elev_lat_plot <- ggplot(plot_df_duration, aes(x = elev_sc, y = duration_pred_unlogged, color = as.factor(lat_label), fill = as.factor(lat_label))) +
  scale_fill_manual(values = cols, guide = guide_legend(title = "Relative Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  geom_line() +
  geom_ribbon(aes(ymin = .lower_unlogged, ymax = .upper_unlogged), alpha = 0.2, linetype = 0) +
  facet_wrap(~temp_seas_label) +
  labs(
    x = "Scaled Elevation",
    y = "Predicted duration (posterior mean)",
    color = "Elevation\n(quantiles)",
    fill = "Elevation\n(quantiles)",
    title = "Predicted Duration across Latitude and Elevation"
  ) +
  theme_bw() + 
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 24, hjust = 0.5)
  )

print(duration_elev_lat_plot)


## Patchwork these plots into a full figure
# Remove legends from the first two plots
onset_elev_lat_plot_nolegend <- onset_elev_lat_plot + ggtitle("Onset") +labs(x = NULL) +theme(legend.position = "none")
offset_elev_lat_plot_nolegend <- offset_elev_lat_plot + ggtitle("Offset") +labs(x = NULL) + theme(legend.position = "none")
duration_elev_lat_plot_legend <- duration_elev_lat_plot + ggtitle("Duration")+ theme(legend.position = "none")

# Combine plots
combined_plot <- (onset_elev_lat_plot_nolegend / offset_elev_lat_plot_nolegend / duration_elev_lat_plot_legend) +
  plot_annotation(title = "Predicted Phenometrics across Latitude and Elevation", 
                  theme = theme(plot.title = element_text(size = 20))) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", 
        strip.text = element_text(size = 12)) # adjusts the size of the facet wrap titles

print(combined_plot)

ggsave(plot = combined_plot, file = "/blue/guralnick/millerjared/PhenoElevation/outputs/figures/conditional-effects-elev-lat-temp-seas.png", height = 17, width = 12)



### Visualize plant habit effects
## Visualize the posterior estimates for latitude and elevation interaction + elevation
# construct new data for predictions
lat_vals <- quantile(general_trends2$latitude_sc, c(0.1, 0.5, 0.9)) # low, median, and high latitude
temp_seas_vals <- quantile(general_trends2$temp_seasonality_sc, c(0.1, 0.5, 0.9)) # low, median, high temp seasonality
elev_seq <- seq(from = min(general_trends2$elev_sc),
                to = max(general_trends2$elev_sc),
                length.out = 80)

newdf_combined <- expand.grid(
  latitude_sc = lat_vals,
  elev_sc = elev_seq,
  temp_seasonality_sc = temp_seas_vals,
  geoannper3f = c("anbien", "herbper"),
  scientific_name = NA,
  id_cells = NA
)

## Get predictions for both life forms
epreds_onset <- add_epred_draws(brm_onset, newdata = newdf_combined, value = "onset_pred", re_formula = ~0, ndraws = 4000)
epreds_offset <- add_epred_draws(brm_offset, newdata = newdf_combined, value = "offset_pred", re_formula = ~0, ndraws = 4000)
epreds_duration <- add_epred_draws(brm_duration, newdata = newdf_combined, value = "duration_pred", re_formula = ~0, ndraws = 4000)

## Create plotting dataframes with labels
create_plot_df <- function(epreds, pred_var) {
  epreds %>%
    group_by(elev_sc, latitude_sc, temp_seasonality_sc, geoannper3f) %>%
    mean_qi(!!sym(pred_var), .width = 0.89) %>%
    mutate(
      lat_label = case_when(
        dplyr::near(latitude_sc, lat_vals[1]) ~ "Southern (10th percentile)",
        dplyr::near(latitude_sc, lat_vals[2]) ~ "Middle (50th percentile)",
        dplyr::near(latitude_sc, lat_vals[3]) ~ "Northern (90th percentile)",
        TRUE ~ NA_character_
      ),
      lat_label = factor(
        lat_label,
        levels = c(
          "Northern (90th percentile)",
          "Middle (50th percentile)",
          "Southern (10th percentile)"
        )
      ),
      temp_seas_label = case_when(
        dplyr::near(temp_seasonality_sc, temp_seas_vals[1]) ~ "Low Temp Seasonality (10th percentile)",
        dplyr::near(temp_seasonality_sc, temp_seas_vals[2]) ~ "Medium Temp Seasonality (50th percentile)",
        dplyr::near(temp_seasonality_sc, temp_seas_vals[3]) ~ "High Temp Seasonality (90th percentile)",
        TRUE ~ NA_character_
      ),
      temp_seas_label = factor(
        temp_seas_label,
        levels = c(
          "Low Temp Seasonality (10th percentile)",
          "Medium Temp Seasonality (50th percentile)",
          "High Temp Seasonality (90th percentile)"
        )
      ),
      life_form = factor(geoannper3f, 
                         levels = c("anbien", "herbper"),
                         labels = c("Annuals", "Perennials"))
    )
}


plot_df_onset <- create_plot_df(epreds_onset, "onset_pred")
plot_df_offset <- create_plot_df(epreds_offset, "offset_pred")
plot_df_duration <- create_plot_df(epreds_duration, "duration_pred") %>%
  mutate(
    duration_pred_unlogged = exp(duration_pred),
    .lower_unlogged = exp(.lower),
    .upper_unlogged = exp(.upper)
  )

## Plot with life form as facet variable
n_bands <- 20
band_colors <- cal_palette("sierra1", n = n_bands, type = "continuous")
cols <- band_colors[c(20, 8, 1)]


onset_elev_lat_plot <- ggplot(plot_df_onset, 
                              aes(x = elev_sc, y = onset_pred, 
                                  color = as.factor(lat_label), 
                                  fill = as.factor(lat_label))) +
  scale_fill_manual(values = cols, guide = guide_legend(title = "Relative Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, linetype = 0) +
  facet_grid(life_form ~ temp_seas_label) +  # Grid with life form in rows
  labs(
    x = "Scaled Elevation",
    y = "Predicted Onset (posterior mean)",
    title = "Onset"
  ) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 24, hjust = 0.5),
    strip.text = element_text(size = 12)
  )


offset_elev_lat_plot <- ggplot(plot_df_offset, 
                               aes(x = elev_sc, y = offset_pred, 
                                   color = as.factor(lat_label), 
                                   fill = as.factor(lat_label))) +
  scale_fill_manual(values = cols, guide = guide_legend(title = "Relative Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, linetype = 0) +
  facet_grid(life_form ~ temp_seas_label) +
  labs(
    x = "Scaled Elevation",
    y = "Predicted Offset (posterior mean)",
    title = "Offset"
  ) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 24, hjust = 0.5),
    strip.text = element_text(size = 12)
  )

duration_elev_lat_plot <- ggplot(plot_df_duration, 
                                 aes(x = elev_sc, y = duration_pred_unlogged, 
                                     color = as.factor(lat_label), 
                                     fill = as.factor(lat_label))) +
  scale_fill_manual(values = cols, guide = guide_legend(title = "Relative Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  geom_line() +
  geom_ribbon(aes(ymin = .lower_unlogged, ymax = .upper_unlogged), alpha = 0.2, linetype = 0) +
  facet_grid(life_form ~ temp_seas_label) +
  labs(
    x = "Scaled Elevation",
    y = "Predicted Duration (posterior mean)",
    title = "Duration"
  ) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 24, hjust = 0.5),
    strip.text = element_text(size = 12)
  )

## Combine plots
onset_elev_lat_plot_nolegend <- onset_elev_lat_plot + labs(x = NULL) + theme(legend.position = "none")
offset_elev_lat_plot_nolegend <- offset_elev_lat_plot + labs(x = NULL) + theme(legend.position = "none")
duration_elev_lat_plot_legend <- duration_elev_lat_plot + theme(legend.position = "none")

combined_plot <- (onset_elev_lat_plot_nolegend / offset_elev_lat_plot_nolegend / duration_elev_lat_plot_legend) +
  plot_annotation(title = "Predicted Phenometrics across Latitude and Elevation",
                  theme = theme(plot.title = element_text(size = 20))) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12))

print(combined_plot)
ggsave(plot = combined_plot, file = "/blue/guralnick/millerjared/PhenoElevation/outputs/figures/conditional-effects-elev-lat-temp-seas-plant-habit.png", height = 17, width = 12)








## Reorient for Rob: 
# lat_vals <- quantile(general_trends2$latitude_sc, c(0.1, 0.5, 0.9))
# elev_vals <- quantile(general_trends2$elev_sc, c(0.1, 0.5, 0.9))
# temp_seas_seq <- seq(
#   from = min(general_trends2$temp_seasonality_sc),
#   to   = max(general_trends2$temp_seasonality_sc),
#   length.out = 80
# )
# 
# newdf <- expand.grid(
#   latitude_sc = lat_vals,
#   elev_sc = elev_vals,
#   temp_seasonality_sc = temp_seas_seq,
#   geoannper3f = "anbien"
# )
# 
# epreds_onset <- add_epred_draws(
#   brm_onset, newdata = newdf,
#   value = "onset_pred",
#   re_formula = ~0
# )
# 
# epreds_offset <- add_epred_draws(
#   brm_offset, newdata = newdf,
#   value = "offset_pred",
#   re_formula = ~0
# )
# 
# epreds_duration <- add_epred_draws(
#   brm_duration, newdata = newdf,
#   value = "duration_pred",
#   re_formula = ~0
# )
# 
# plot_df_onset <- epreds_onset %>%
#   group_by(temp_seasonality_sc, latitude_sc, elev_sc) %>%
#   mean_qi(onset_pred, .width = 0.89) %>%
#   mutate(
#     lat_label = case_when(
#       dplyr::near(latitude_sc, lat_vals[1]) ~ "Southern (10th percentile)",
#       dplyr::near(latitude_sc, lat_vals[2]) ~ "Middle (50th percentile)",
#       dplyr::near(latitude_sc, lat_vals[3]) ~ "Northern (90th percentile)",
#       TRUE ~ NA_character_
#     ),
#     lat_label = factor(
#       lat_label,
#       levels = c("Northern (90th percentile)", "Middle (50th percentile)", "Southern (10th percentile)")
#     ),
#     elev_label = case_when(
#       dplyr::near(elev_sc, elev_vals[1]) ~ "Low Elevation (10th percentile)",
#       dplyr::near(elev_sc, elev_vals[2]) ~ "Medium Elevation (50th percentile)",
#       dplyr::near(elev_sc, elev_vals[3]) ~ "High Elevation (90th percentile)",
#       TRUE ~ NA_character_
#     ),
#     elev_label = factor(
#       elev_label,
#       levels = c("Low Elevation (10th percentile)", "Medium Elevation (50th percentile)", "High Elevation (90th percentile)")
#     )
#   )
# 
# plot_df_offset <- epreds_offset %>%
#   group_by(temp_seasonality_sc, latitude_sc, elev_sc) %>%
#   mean_qi(offset_pred, .width = 0.89) %>%
#   mutate(
#     lat_label = case_when(
#       dplyr::near(latitude_sc, lat_vals[1]) ~ "Southern (10th percentile)",
#       dplyr::near(latitude_sc, lat_vals[2]) ~ "Middle (50th percentile)",
#       dplyr::near(latitude_sc, lat_vals[3]) ~ "Northern (90th percentile)",
#       TRUE ~ NA_character_
#     ),
#     lat_label = factor(
#       lat_label,
#       levels = c("Northern (90th percentile)", "Middle (50th percentile)", "Southern (10th percentile)")
#     ),
#     elev_label = case_when(
#       dplyr::near(elev_sc, elev_vals[1]) ~ "Low Elevation (10th percentile)",
#       dplyr::near(elev_sc, elev_vals[2]) ~ "Medium Elevation (50th percentile)",
#       dplyr::near(elev_sc, elev_vals[3]) ~ "High Elevation (90th percentile)",
#       TRUE ~ NA_character_
#     ),
#     elev_label = factor(
#       elev_label,
#       levels = c("Low Elevation (10th percentile)", "Medium Elevation (50th percentile)", "High Elevation (90th percentile)")
#     )
#   )
# 
# plot_df_duration <- epreds_duration %>%
#   group_by(temp_seasonality_sc, latitude_sc, elev_sc) %>%
#   mean_qi(duration_pred, .width = 0.89) %>%
#   mutate(
#     lat_label = case_when(
#       dplyr::near(latitude_sc, lat_vals[1]) ~ "Southern (10th percentile)",
#       dplyr::near(latitude_sc, lat_vals[2]) ~ "Middle (50th percentile)",
#       dplyr::near(latitude_sc, lat_vals[3]) ~ "Northern (90th percentile)",
#       TRUE ~ NA_character_
#     ),
#     lat_label = factor(
#       lat_label,
#       levels = c("Northern (90th percentile)", "Middle (50th percentile)", "Southern (10th percentile)")
#     ),
#     elev_label = case_when(
#       dplyr::near(elev_sc, elev_vals[1]) ~ "Low Elevation (10th percentile)",
#       dplyr::near(elev_sc, elev_vals[2]) ~ "Medium Elevation (50th percentile)",
#       dplyr::near(elev_sc, elev_vals[3]) ~ "High Elevation (90th percentile)",
#       TRUE ~ NA_character_
#     ),
#     elev_label = factor(
#       elev_label,
#       levels = c("Low Elevation (10th percentile)", "Medium Elevation (50th percentile)", "High Elevation (90th percentile)")
#     )
#   )
# 
# n_bands <- 20
# band_colors <- cal_palette("sierra1", n = n_bands, type = "continuous")
# cols <- band_colors[c(20, 8, 1)]
# 
# onset_temp_lat_facetElev_plot <- ggplot(
#   plot_df_onset,
#   aes(x = temp_seasonality_sc, y = onset_pred, color = lat_label, fill = lat_label)
# ) +
#   geom_line(linewidth = 1) +
#   geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, linetype = 0) +
#   facet_wrap(~ elev_label) +
#   scale_fill_manual(values = cols, guide = guide_legend(title = "Relative Latitude")) +
#   scale_color_manual(values = cols, guide = "none") +
#   labs(
#     x = "Scaled Temperature Seasonality",
#     y = "Predicted mean onset",
#     title = "Predicted onset across temperature seasonality (faceted by elevation)"
#   ) +
#   theme_bw() +
#   theme(
#     legend.title = element_text(size = 16),
#     legend.text  = element_text(size = 14),
#     axis.title   = element_text(size = 18),
#     axis.text    = element_text(size = 14),
#     plot.title   = element_text(size = 20, hjust = 0.5)
#   )
# 
# print(onset_temp_lat_facetElev_plot)
# ggsave(plot = onset_temp_lat_facetElev_plot, "/blue/guralnick/millerjared/PhenoElevation/outputs/figures/")