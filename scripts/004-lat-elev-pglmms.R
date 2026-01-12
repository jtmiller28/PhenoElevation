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

## Scale and create factors 
phenologycv_sf14_onoffall21 <- phenologycv_sf14_onoffall20 %>%
  dplyr::mutate(latitude_sc =scale(latitude), 
                longitude_sc=scale(longitude),
                elev_sc = scale(elev),
                spmean_sc = scale(spmean),
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
                    geoannper3f + 
                    year + 
                    (1 | scientific_name) + 
                    (1 | id_cells) +
                    + elev_sc:latitude_sc + elev_sc:geoannper3f + latitude_sc:geoannper3f,                          
                  data = general_trends) 
mod_offset <- lmer(offset ~ elev_sc + 
                     latitude_sc +
                     geoannper3f + 
                     year + 
                     (1 | scientific_name) + 
                     (1 | id_cells) +
                     + elev_sc:latitude_sc + elev_sc:geoannper3f + latitude_sc:geoannper3f,                          
                  data = general_trends ) 
mod_duration <- lmer(duration ~ elev_sc + 
                     latitude_sc +
                     geoannper3f + 
                     year + 
                     (1 | scientific_name) + 
                     (1 | id_cells) +
                     + elev_sc:latitude_sc + elev_sc:geoannper3f + latitude_sc:geoannper3f,                          
                   data = general_trends ) 

## use model selection on these mixed-effects models
print("Check onset glmms model selection")
step(mod_onset)
print("Check offset glmms model selection")
step(mod_offset)
print("Check duration glmms model selection")
step(mod_duration)

## check inflation using the package performance for mixed-effects
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
                       geoannper3f + 
                       year + 
                       (1 | scientific_name) + 
                       (1 | id_cells) +
                       + elev_sc:latitude_sc + elev_sc:geoannper3f + latitude_sc:geoannper3f,                          
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
    geoannper3f +
    elev_sc:latitude_sc + 
    elev_sc:geoannper3f + 
    # latitude_sc:geoannper3f + # dropping to increase bulk sampling density...
    (1|gr(scientific_name, cov = A)) + 
    (1 | id_cells) , 
  data = general_trends2, 
  data2 = list(A = A), 
  family = gaussian(), 
  chains = 4, cores = 4, iter = 15000, control = list(adapt_delta = 0.99, max_treedepth = 15)
)
saveRDS(brm_onset, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/onset-brms-phylo-model.rds") # 12 hours est 
brm_offset <- brm(
  offset ~ elev_sc + 
    latitude_sc + 
    geoannper3f +
    elev_sc:latitude_sc + 
    elev_sc:geoannper3f + 
    # latitude_sc:geoannper3f + # dropping to increase bulk sampling density...
    (1|gr(scientific_name, cov = A)) + 
    (1 | id_cells) , 
  data = general_trends2, 
  data2 = list(A = A), 
  family = gaussian(), 
  chains = 4, cores = 4, iter = 15000, control = list(adapt_delta = 0.99, max_treedepth = 15)
)
saveRDS(brm_offset, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/offset-brms-phylo-model.rds")
brm_duration <- brm(
  log_duration ~ elev_sc + 
    latitude_sc + 
    geoannper3f +
    elev_sc:latitude_sc + 
    elev_sc:geoannper3f + 
    # latitude_sc:geoannper3f + # dropping to increase bulk sampling density...
    (1|gr(scientific_name, cov = A)) + 
    (1 | id_cells) , 
  data = general_trends2, 
  data2 = list(A = A), 
  family = gaussian(), 
  chains = 4, cores = 4, iter = 15000, control = list(adapt_delta = 0.99, max_treedepth = 15)
)
saveRDS(brm_duration, "/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/duration-brms-phylo-model.rds")

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
brm_onset <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/onset-brms-phylo-model.rds")
brm_offset <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/offset-brms-phylo-model.rds")
brm_duration <- readRDS("/blue/guralnick/millerjared/PhenoElevation/outputs/rds-objects/duration-brms-phylo-model.rds")

## Visualize the posterior estimates for latitude and elevation interaction
# construct new data for predictions
lat_vals <- quantile(general_trends2$latitude_sc, c(0.1, 0.5, 0.9)) # low, median, and high elev
elev_seq <- seq(from = min(general_trends2$elev_sc), 
                to = max(general_trends2$elev_sc), 
                length.out = 80)
newdf <- expand.grid(
  latitude_sc = lat_vals, 
  elev_sc = elev_seq, 
  geoannper3f = "anbien", # use annual/biannual for this
  scientific_name = NA, 
  id_cells = NA
)

# get posterior predictions for all combos, disclude re 
epreds_onset <- add_epred_draws(brm_onset, newdata = newdf, value = "onset_pred", re_formula = NA)
epreds_offset <- add_epred_draws(brm_offset, newdata = newdf, value = "offset_pred", re_formula = NA)
epreds_duration <- add_epred_draws(brm_duration, newdata = newdf, value = "duration_pred", re_formula = NA)

# summarize for plotting 
plot_df_onset <- epreds_onset %>%
  group_by(elev_sc, latitude_sc) %>%
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
    )
  )
plot_df_offset <- epreds_offset %>%
  group_by(elev_sc, latitude_sc) %>%
  mean_qi(offset_pred, .width = 0.89) %>%
  mutate(
    lat_label = case_when(
      dplyr::near(latitude_sc, lat_vals[1]) ~ "Northern (90th percentile)",
      dplyr::near(latitude_sc, lat_vals[2]) ~ "Middle (50th percentile)",
      dplyr::near(latitude_sc, lat_vals[3]) ~ "Southern (10th percentile)",
      TRUE ~ NA_character_
    ),
    lat_label = factor(
      lat_label,
      levels = c(
        "Northern (90th percentile)",
        "Middle (50th percentile)",
        "Southern (10th percentile)"
      )
    )
  )

plot_df_duration <- epreds_duration %>%
  group_by(elev_sc, latitude_sc) %>%
  mean_qi(duration_pred, .width = 0.89) %>%
  mutate(
    lat_label = case_when(
      dplyr::near(latitude_sc, lat_vals[1]) ~ "Northern (90th percentile)",
      dplyr::near(latitude_sc, lat_vals[2]) ~ "Middle (50th percentile)",
      dplyr::near(latitude_sc, lat_vals[3]) ~ "Southern (10th percentile)",
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
    plot.title = element_text(size = 20, hjust = 0.5)
  )

print(onset_elev_lat_plot)

offset_elev_lat_plot <- ggplot(plot_df_offset, aes(x = elev_sc, y = offset_pred, color = as.factor(lat_label), fill = as.factor(lat_label))) +
  scale_fill_manual(values = cols, guide = guide_legend(title = "Relative Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, linetype = 0) +
  labs(
    x = "Scaled Elevation",
    y = "Predicted Offset (posterior mean)",
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
    plot.title = element_text(size = 20, hjust = 0.5)
  )

print(offset_elev_lat_plot)

duration_elev_lat_plot <- ggplot(plot_df_duration, aes(x = elev_sc, y = duration_pred_unlogged, color = as.factor(lat_label), fill = as.factor(lat_label))) +
  scale_fill_manual(values = cols, guide = guide_legend(title = "Relative Latitude")) +
  scale_color_manual(values = cols, guide = "none") +
  geom_line() +
  geom_ribbon(aes(ymin = .lower_unlogged, ymax = .upper_unlogged), alpha = 0.2, linetype = 0) +
  labs(
    x = "Scaled Elevation",
    y = "Predicted Duration (posterior mean)",
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
    plot.title = element_text(size = 20, hjust = 0.5)
  )

print(duration_elev_lat_plot)

## Patchwork these plots into a full figure
# Remove legends from the first two plots
onset_elev_lat_plot_nolegend <- onset_elev_lat_plot + ggtitle("Onset") + theme(legend.position = "none")
offset_elev_lat_plot_nolegend <- offset_elev_lat_plot + ggtitle("Offset") + theme(legend.position = "none")
duration_elev_lat_plot_legend <- duration_elev_lat_plot + ggtitle("Duration") # keep legend here

# Combine plots
combined_plot <- (onset_elev_lat_plot_nolegend | offset_elev_lat_plot_nolegend | duration_elev_lat_plot_legend) +
  plot_annotation(title = "Predicted Phenometrics across Latitude and Elevation", 
                  theme = theme(plot.title = element_text(size = 20)))

print(combined_plot)

ggsave("./outputs/figures/lat-elev-phenometric-pred-plot.png", width = 15, height = 10)


## Evaluate Phylogenetic Signal on Onset, Offset, and Duration
onset_post <- as_draws_df(brm_onset)
onset_draws <- as_draws_df(brm_onset)
onset_lambda <- with(onset_draws, sd_scientific_name__Intercept^2 / (sd_scientific_name__Intercept^2 + sigma^2))
print("Onset lambda estimate")
c(mean = mean(onset_lambda),
  quantile(onset_lambda, c(0.025, 0.5, 0.975)),
  p_gt0 = mean(onset_lambda > 0))

offset_post <- as_draws_df(brm_offset)
offset_draws <- as_draws_df(brm_offset)
offset_lambda <- with(offset_draws, sd_scientific_name__Intercept^2 / (sd_scientific_name__Intercept^2 + sigma^2))
print("Offset lambda estimate")
c(mean = mean(offset_lambda),
  quantile(offset_lambda, c(0.025, 0.5, 0.975)),
  p_gt0 = mean(offset_lambda > 0))

duration_post <- as_draws_df(brm_duration)
duration_draws <- as_draws_df(brm_duration)
duration_lambda <- with(duration_draws, sd_scientific_name__Intercept^2 / (sd_scientific_name__Intercept^2 + sigma^2))
print("Duration lambda estimate")
c(mean = mean(duration_lambda),
  quantile(duration_lambda, c(0.025, 0.5, 0.975)),
  p_gt0 = mean(duration_lambda > 0))

