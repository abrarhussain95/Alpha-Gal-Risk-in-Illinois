

library(dplyr)
library(readr)
library(stringr)

## Read tick data
file_path <- "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Data/alpha_gal_data.csv"

data_aa <- read_csv(file_path)

#View(data_aa)


##Shape File

library(sf)

# Import shapefile
IL_shape <- st_read("C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/TBD Risk Score Illinois/Illinois_shape_file/IL_BNDY_County_Py.shp")

# Print column names
print(colnames(IL_shape))


# Keep only county name + geometry
IL_shape <- IL_shape %>%
  dplyr::select(COUNTY_NAM, geometry)

# Fix county name ONLY for merge
IL_shape <- IL_shape %>%
  mutate(County_std = COUNTY_NAM |> str_to_lower() |> str_replace_all(" ", ""))

data_aa <- data_aa %>%
  mutate(County_std = County    |> str_to_lower() |> str_replace_all(" ", ""))

# Merge shapefile 
data_aa <- IL_shape %>%
  left_join(data_aa, by = "County_std")

# check
nrow(data_aa) 

# Drop helper columns
data_aa <- data_aa %>%
  dplyr::select(-COUNTY_NAM, -County_std)

# Confirm final structure
print(colnames(data_aa))

# View
#View(data_aa)


median(data_aa$tick_aa); range(data_aa$tick_aa)
median(data_aa$tbd_cases_aa); range(data_aa$tbd_cases_aa)






## Build adjacency matrix (queen contiguity)
library(spdep)

# poly2nb needs an sp-style object, not sf
data_aa_sp <- as(data_aa, "Spatial")

nb <- poly2nb(data_aa_sp, queen = TRUE)

# does every county have at least one neighbor?
no_neighbors <- which(card(nb) == 0)
if (length(no_neighbors) > 0) {
  cat("WARNING: these counties have NO neighbors detected:\n")
  print(data_aa$County[no_neighbors])
} else {
  cat("All counties have at least one neighbor. Good.\n")
}

# Convert neighbor list to a binary adjacency matrix
W <- nb2mat(nb, style = "B", zero.policy = TRUE)

cat("Adjacency matrix dimensions:", dim(W)[1], "x", dim(W)[2], "\n")
cat("Average neighbors per county:", round(mean(rowSums(W)), 2), "\n")
cat("Min neighbors:", min(rowSums(W)), " Max neighbors:", max(rowSums(W)), "\n")




# code tick establishment status into numeric weight
data_aa <- data_aa %>%
  dplyr::mutate(
    status_weight = dplyr::case_when(
      Status_aa == "Established"          ~ 1.0,
      Status_aa == "Reported"              ~ 0.5,
      Status_aa == "No Known Population"   ~ 0.0,
      TRUE ~ NA_real_
    )
  )

# Sanity check — every county should have a value, no NAs, no unexpected categories
table(data_aa$Status_aa, useNA = "always")
sum(is.na(data_aa$status_weight))





library(CARBayes)

# Log-population offset for the case model
data_aa$log_pop <- log(data_aa$population)

model_data <- data_aa %>% st_drop_geometry()

## TBD

set.seed(123)

model_cases <- S.CARleroux(
  formula = tbd_cases_aa ~ 1 + offset(log_pop),
  family = "poisson",
  data = model_data,
  W = W,
  burnin = 5000,
  n.sample = 20000,
  thin = 10,
  verbose = TRUE
)

cat("\n=== CASE MODEL SUMMARY ===\n")
print(model_cases$summary.results)

## Tick

set.seed(123)

model_ticks <- S.CARleroux(
  formula = tick_aa ~ 1,
  family = "poisson",
  data = model_data,
  W = W,
  burnin = 5000,
  n.sample = 20000,
  thin = 10,
  verbose = TRUE
)

cat("\n=== TICK MODEL SUMMARY ===\n")
print(model_ticks$summary.results)


# Extract posterior mean fitted values (the smoothed risk surface) from each model
data_aa$case_risk_raw <- model_cases$fitted.values
data_aa$tick_risk_raw <- model_ticks$fitted.values

data_aa %>% 
  sf::st_drop_geometry() %>% 
  dplyr::select(County, tbd_cases_aa, case_risk_raw, tick_aa, tick_risk_raw) %>% 
  head(10)


# Convert fitted case counts into a rate per 100,000 — comparable across counties
data_aa$case_risk_rate <- (data_aa$case_risk_raw / data_aa$population) * 100000

# Check it
data_aa %>% 
  sf::st_drop_geometry() %>% 
  dplyr::select(County, population, tbd_cases_aa, case_risk_raw, case_risk_rate) %>% 
  head(10)







min_max_norm <- function(x) (x - min(x)) / (max(x) - min(x))
data_aa$case_risk_norm <- min_max_norm(data_aa$case_risk_rate)
data_aa$tick_risk_norm <- min_max_norm(data_aa$tick_risk_raw)

w_case   <- 0.5
w_tick   <- 0.45
w_status <- 0.05

data_aa$combined_risk_score <- w_case * data_aa$case_risk_norm +
  w_tick * data_aa$tick_risk_norm +
  w_status * data_aa$status_weight

data_aa %>%
  sf::st_drop_geometry() %>%
  dplyr::select(County, Status_aa, status_weight, case_risk_norm, tick_risk_norm, combined_risk_score) %>%
  dplyr::arrange(dplyr::desc(combined_risk_score)) %>%
  head(15)






library(sf)
library(dplyr)
library(tmap)

# ONE-TIME STEP: only run this block if shawnee_boundary.shp doesn't exist yet  
#shawnee_url <- paste0(
# "https://apps.fs.usda.gov/arcx/rest/services/EDW/EDW_ForestSystemBoundaries_01/",
# "MapServer/0/query?where=FORESTNAME+LIKE+'%25Shawnee%25'",
# "&outFields=*&f=geojson"
#)
#shawnee <- st_read(shawnee_url)

#dir.create("C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Data/forest shape file",
# showWarnings = FALSE, recursive = TRUE)

#st_write(shawnee,
# "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Data/forest shape file/shawnee_boundary.shp",
#     delete_layer = TRUE)


shawnee <- st_read("C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Data/forest shape file/shawnee_boundary.shp")
shawnee <- st_transform(shawnee, st_crs(data_aa))

map_combined <- tm_shape(data_aa) +
  tm_fill("combined_risk_score",
          fill.scale = tm_scale_continuous(
            values = c("#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026"),
            trans = "sqrt"),
          fill.legend = tm_legend(title = "Risk score", show.na = FALSE)) +
  tm_borders(lwd = 0.3, col = "grey40") +
  tm_shape(shawnee) +
  tm_borders(lwd = 1.5, col = "darkgreen", lty = "solid") +
  tm_add_legend(type = "lines",
                labels = "Shawnee National Forest",
                col = "darkgreen",
                lwd = 1.5) +
  tm_layout(legend.outside = TRUE,
            frame = FALSE,
            legend.frame = FALSE)

print(map_combined)

tmap_save(map_combined,
          filename = "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Output/illinois_combined_risk_map.tiff",
          width = 8, height = 6, units = "in", dpi = 300,
          compression = "lzw")








case_samples <- model_cases$samples$fitted
tick_samples <- model_ticks$samples$fitted

pop_vec <- data_aa$population
case_rate_samples <- sweep(case_samples, 2, pop_vec, FUN = "/") * 100000

normalize_row <- function(row) {
  (row - min(row)) / (max(row) - min(row))
}

case_norm_samples <- t(apply(case_rate_samples, 1, normalize_row))
tick_norm_samples  <- t(apply(tick_samples,      1, normalize_row))

status_matrix <- matrix(rep(data_aa$status_weight, each = 1500),
                        nrow = 1500, ncol = 102)

combined_samples <- w_case * case_norm_samples +
  w_tick * tick_norm_samples +
  w_status * status_matrix

data_aa$combined_mean  <- apply(combined_samples, 2, mean)
data_aa$combined_lower <- apply(combined_samples, 2, quantile, probs = 0.025)
data_aa$combined_upper <- apply(combined_samples, 2, quantile, probs = 0.975)
data_aa$combined_width <- data_aa$combined_upper - data_aa$combined_lower

data_aa %>%
  sf::st_drop_geometry() %>%
  dplyr::select(County, Status_aa, status_weight, combined_mean, combined_lower, combined_upper) %>%
  dplyr::arrange(dplyr::desc(combined_mean)) %>%
  head(15)








plot_data <- data_aa %>%
  sf::st_drop_geometry() %>%
  dplyr::arrange(dplyr::desc(combined_mean)) %>%
  head(20) %>%
  dplyr::mutate(County = stringr::str_to_title(County))

plot_data$County <- factor(plot_data$County, levels = rev(plot_data$County))

library(ggplot2)

forest_plot <- ggplot(plot_data, aes(x = combined_mean, y = County)) +
  geom_col(fill = "orangered4", width = 0.65, alpha = 0.85) +
  geom_errorbar(aes(xmin = combined_lower, xmax = combined_upper),
                width = 0.25, linewidth = 0.5, color = "#4A1B0C") +
  scale_x_continuous(limits = c(0, 1.05), expand = c(0.01, 0.01)) +
  labs(x = "Risk score", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 11, color = "grey20"),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(size = 11, color = "black", margin = margin(t = 10)),
    plot.margin = margin(15, 20, 15, 15)
  )

forest_plot

ggsave(
  plot = forest_plot,
  filename = "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Output/illinois_risk_forest_plot.tiff",
  width = 8, height = 6, units = "in", dpi = 300,
  compression = "lzw"
)






library(dplyr)
library(knitr)

table_full <- data_aa %>%
  sf::st_drop_geometry() %>%
  dplyr::arrange(dplyr::desc(combined_mean)) %>%
  dplyr::select(County, combined_mean, combined_lower, combined_upper) %>%
  dplyr::mutate(
    County = stringr::str_to_title(County),
    combined_mean  = round(combined_mean, 3),
    combined_lower = round(combined_lower, 3),
    combined_upper = round(combined_upper, 3)
  )

# View it in the console
print(dplyr::as_tibble(table_full), n = 102)

# Render as a formatted table in your R Markdown output
table_full %>%
  kable(
    col.names = c("County", "Risk score (mean)", "95% CrI lower", "95% CrI upper"),
    caption = "Table S2. Combined alpha-gal syndrome proxy risk score with 95% credible intervals for all 102 Illinois counties, ranked from highest to lowest."
  )

write.csv(table_full,
          "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Output/Supplementary file.csv",
          row.names = FALSE)




library(ggplot2)
library(tidyr)
library(dplyr)

#Extract MCMC samples for both models
case_trace <- as.data.frame(model_cases$samples$beta)
colnames(case_trace) <- "Intercept"
case_trace$tau2 <- as.vector(model_cases$samples$tau2)
case_trace$rho  <- as.vector(model_cases$samples$rho)
case_trace$iteration <- 1:nrow(case_trace)
case_trace$Model <- "Ehrlichiosis case model"

tick_trace <- as.data.frame(model_ticks$samples$beta)
colnames(tick_trace) <- "Intercept"
tick_trace$tau2 <- as.vector(model_ticks$samples$tau2)
tick_trace$rho  <- as.vector(model_ticks$samples$rho)
tick_trace$iteration <- 1:nrow(tick_trace)
tick_trace$Model <- "Tick abundance model"

# Combine and reshape to long format for faceting
trace_all <- bind_rows(case_trace, tick_trace) %>%
  pivot_longer(cols = c(Intercept, tau2, rho),
               names_to = "Parameter", values_to = "Value")

trace_all$Parameter <- factor(trace_all$Parameter,
                              levels = c("Intercept", "tau2", "rho"),
                              labels = c("Intercept", "tau\u00b2", "rho"))

# Plot
trace_plot <- ggplot(trace_all, aes(x = iteration, y = Value, color = Model)) +
  geom_line(linewidth = 0.3, alpha = 0.8) +
  facet_grid(Parameter ~ Model, scales = "free_y") +
  scale_color_manual(values = c("Ehrlichiosis case model" = "#993C1D",
                                "Tick abundance model" = "#0F6E56")) +
  labs(x = "MCMC sample index", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(15, 15, 15, 15)
  )

trace_plot

ggsave(
  plot = trace_plot,
  filename = "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Output/figure_S1_trace_plots.tiff",
  width = 9, height = 8, units = "in", dpi = 300,
  compression = "lzw"
)










## Model validation: posterior predictive checks

library(dplyr)
library(tidyr)

# Fitted expected counts from CARBayes posterior samples
case_lambda <- model_cases$samples$fitted
tick_lambda <- model_ticks$samples$fitted

# Observed counts
obs_cases <- model_data$tbd_cases_aa
obs_ticks <- model_data$tick_aa

# Dimensions
n_samp <- nrow(case_lambda)
n_county <- ncol(case_lambda)

set.seed(123)

# Simulate replicated datasets from posterior predictive distribution
case_rep <- matrix(
  rpois(n_samp * n_county, lambda = as.vector(case_lambda)),
  nrow = n_samp,
  ncol = n_county
)

tick_rep <- matrix(
  rpois(n_samp * n_county, lambda = as.vector(tick_lambda)),
  nrow = n_samp,
  ncol = n_county
)

# Function to calculate summary statistics for each replicated dataset
ppc_stats <- function(rep_matrix) {
  data.frame(
    zero_count  = apply(rep_matrix, 1, function(x) sum(x == 0)),
    mean_count  = apply(rep_matrix, 1, mean),
    max_count   = apply(rep_matrix, 1, max),
    total_count = apply(rep_matrix, 1, sum)
  )
}

case_ppc <- ppc_stats(case_rep)
tick_ppc <- ppc_stats(tick_rep)

# Observed summary statistics
case_obs_stats <- data.frame(
  outcome = "Ehrlichiosis cases",
  statistic = c("zero_count", "mean_count", "max_count", "total_count"),
  observed = c(
    sum(obs_cases == 0),
    mean(obs_cases),
    max(obs_cases),
    sum(obs_cases)
  )
)

tick_obs_stats <- data.frame(
  outcome = "Tick abundance",
  statistic = c("zero_count", "mean_count", "max_count", "total_count"),
  observed = c(
    sum(obs_ticks == 0),
    mean(obs_ticks),
    max(obs_ticks),
    sum(obs_ticks)
  )
)

# Summarize posterior predictive distributions
summarize_ppc <- function(ppc_df, outcome_name) {
  ppc_df %>%
    pivot_longer(
      cols = everything(),
      names_to = "statistic",
      values_to = "replicated"
    ) %>%
    group_by(statistic) %>%
    summarise(
      ppc_median = median(replicated),
      ppc_lower  = quantile(replicated, 0.025),
      ppc_upper  = quantile(replicated, 0.975),
      .groups = "drop"
    ) %>%
    mutate(outcome = outcome_name)
}

case_ppc_summary <- summarize_ppc(case_ppc, "Ehrlichiosis cases")
tick_ppc_summary <- summarize_ppc(tick_ppc, "Tick abundance")

ppc_summary <- bind_rows(case_ppc_summary, tick_ppc_summary) %>%
  left_join(
    bind_rows(case_obs_stats, tick_obs_stats),
    by = c("outcome", "statistic")
  ) %>%
  dplyr::select(outcome, statistic, observed, ppc_median, ppc_lower, ppc_upper)

print(ppc_summary)

write.csv(
  ppc_summary,
  "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Output/posterior_predictive_check_summary.csv",
  row.names = FALSE
)


### FIGURES  ###

## Posterior predictive check plot

library(ggplot2)
library(dplyr)
library(tidyr)

case_ppc_long <- case_ppc %>%
  dplyr::mutate(draw = dplyr::row_number(), outcome = "Ehrlichiosis cases") %>%
  tidyr::pivot_longer(
    cols = c(zero_count, mean_count, max_count, total_count),
    names_to = "statistic",
    values_to = "replicated"
  )

tick_ppc_long <- tick_ppc %>%
  dplyr::mutate(draw = dplyr::row_number(), outcome = "Tick abundance") %>%
  tidyr::pivot_longer(
    cols = c(zero_count, mean_count, max_count, total_count),
    names_to = "statistic",
    values_to = "replicated"
  )

ppc_long <- dplyr::bind_rows(case_ppc_long, tick_ppc_long)
obs_long <- dplyr::bind_rows(case_obs_stats, tick_obs_stats)

# Clean facet labels
stat_labels <- c(
  zero_count  = "Zero count",
  mean_count  = "Mean count",
  max_count   = "Max. count",
  total_count = "Total count"
)

ppc_long <- ppc_long %>%
  dplyr::mutate(
    statistic = factor(
      statistic,
      levels = c("max_count", "mean_count", "total_count", "zero_count"),
      labels = c("Max. count", "Mean count", "Total count", "Zero count")
    ),
    outcome = factor(
      outcome,
      levels = c("Ehrlichiosis cases", "Tick abundance")
    )
  )

obs_long <- obs_long %>%
  dplyr::mutate(
    statistic = factor(
      statistic,
      levels = c("max_count", "mean_count", "total_count", "zero_count"),
      labels = c("Max. count", "Mean count", "Total count", "Zero count")
    ),
    outcome = factor(
      outcome,
      levels = c("Ehrlichiosis cases", "Tick abundance")
    )
  )

# Posterior predictive median for black dotted line
ppc_medians <- ppc_long %>%
  dplyr::group_by(outcome, statistic) %>%
  dplyr::summarise(
    ppc_median = median(replicated),
    .groups = "drop"
  )

ppc_plot <- ggplot(ppc_long, aes(x = replicated)) +
  geom_histogram(
    bins = 30,
    fill = "grey75",
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    data = ppc_medians,
    aes(xintercept = ppc_median),
    color = "black",
    linetype = "dotted",
    linewidth = 0.9
  ) +
  geom_vline(
    data = obs_long,
    aes(xintercept = observed),
    color = "#08306B",
    linewidth = 1.2
  ) +
  facet_grid(outcome ~ statistic, scales = "free") +
  labs(
    x = "Posterior predictive statistic",
    y = "Frequency",
    caption = "Dark blue line = observed value; black dotted line = posterior predictive median."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold", color = "black", size = 11),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", face = "bold", size = 11),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.line = element_line(color = "black", linewidth = 0.3),
    plot.caption = element_text(color = "black", size = 9, hjust = 0),
    panel.spacing = unit(1.0, "lines")
  )

ppc_plot

ggsave(
  plot = ppc_plot,
  filename = "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Output/posterior_predictive_checks.tiff",
  width = 9,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)







library(dplyr)
library(ggplot2)
library(tmap)
library(knitr)

# Scenario 1 weights: no establishment status
# Equal weighting of ehrlichiosis and tick abundance only
w_case_s1 <- 0.5
w_tick_s1 <- 0.5

combined_samples_s1 <- w_case_s1 * case_norm_samples +
  w_tick_s1 * tick_norm_samples

# Save posterior mean and 95% credible intervals
data_aa$combined_mean_s1  <- apply(combined_samples_s1, 2, mean)
data_aa$combined_lower_s1 <- apply(combined_samples_s1, 2, quantile, probs = 0.025)
data_aa$combined_upper_s1 <- apply(combined_samples_s1, 2, quantile, probs = 0.975)
data_aa$combined_width_s1 <- data_aa$combined_upper_s1 - data_aa$combined_lower_s1

# Check top counties
data_aa %>%
  sf::st_drop_geometry() %>%
  dplyr::select(County, Status_aa, combined_mean_s1, combined_lower_s1, combined_upper_s1) %>%
  dplyr::arrange(dplyr::desc(combined_mean_s1)) %>%
  head(15)

# CSV output

table_s1 <- data_aa %>%
  sf::st_drop_geometry() %>%
  dplyr::arrange(dplyr::desc(combined_mean_s1)) %>%
  dplyr::select(County, combined_mean_s1, combined_lower_s1, combined_upper_s1) %>%
  dplyr::mutate(
    County = stringr::str_to_title(County),
    combined_mean_s1  = round(combined_mean_s1, 3),
    combined_lower_s1 = round(combined_lower_s1, 3),
    combined_upper_s1 = round(combined_upper_s1, 3)
  )

print(dplyr::as_tibble(table_s1), n = 102)


table_s1 %>%
  kable(
    col.names = c("County", "Risk score mean", "95% CrI lower", "95% CrI upper"),
    caption = "Supplementary file. Sensitivity analysis 1: equal weighting of ehrlichiosis, tick abundance, and establishment status."
  )

write.csv(
  table_s1,
  "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Output/supplementary file (sensitivity analysis 1).csv",
  row.names = FALSE
)


# Map output

# Load Shawnee National Forest boundary
if (!exists("shawnee")) {
  shawnee <- sf::st_read(
    "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Data/forest shape file/shawnee_boundary.shp"
  )
  shawnee <- sf::st_transform(shawnee, sf::st_crs(data_aa))
}

map_s1 <- tm_shape(data_aa) +
  tm_fill(
    "combined_mean_s1",
    fill.scale = tm_scale_continuous(
      values = c("#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026"),
      trans = "sqrt"
    ),
    fill.legend = tm_legend(title = "Risk score", show.na = FALSE)
  ) +
  tm_borders(lwd = 0.3, col = "grey40") +
  tm_shape(shawnee) +
  tm_borders(lwd = 1.5, col = "darkgreen", lty = "solid") +
  tm_add_legend(
    type = "lines",
    labels = "Shawnee National Forest",
    col = "darkgreen",
    lwd = 1.5
  ) +
  tm_layout(
    legend.outside = TRUE,
    frame = FALSE,
    legend.frame = FALSE
  )

print(map_s1)

tmap_save(
  map_s1,
  filename = "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Output/illinois_combined_risk_map (sensitivity analysis 1).tiff",
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)


# Forest plot output


plot_data_s1 <- data_aa %>%
  sf::st_drop_geometry() %>%
  dplyr::arrange(dplyr::desc(combined_mean_s1)) %>%
  head(20) %>%
  dplyr::mutate(County = stringr::str_to_title(County))

plot_data_s1$County <- factor(plot_data_s1$County, levels = rev(plot_data_s1$County))

forest_plot_s1 <- ggplot(plot_data_s1, aes(x = combined_mean_s1, y = County)) +
  geom_col(fill = "orangered4", width = 0.65, alpha = 0.85) +
  geom_errorbar(
    aes(xmin = combined_lower_s1, xmax = combined_upper_s1),
    width = 0.25,
    linewidth = 0.5,
    color = "#4A1B0C"
  ) +
  scale_x_continuous(limits = c(0, 1.05), expand = c(0.01, 0.01)) +
  labs(x = "Risk score", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 11, color = "grey20"),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(size = 11, color = "black", margin = margin(t = 10)),
    plot.margin = margin(15, 20, 15, 15)
  )

forest_plot_s1

ggsave(
  plot = forest_plot_s1,
  filename = "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Output/illinois_risk_forest_plot (sensitivity analysis 1).tiff",
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)










library(dplyr)
library(ggplot2)
library(tmap)
library(knitr)

# Scenario 2 weights: Disease-heavy
w_case_s2   <- 0.7
w_tick_s2   <- 0.25
w_status_s2 <- 0.05

# Status is already coded as:
# Established = 1, Reported = 0.5, No Known Population = 0
status_matrix_s2 <- matrix(
  rep(data_aa$status_weight, each = nrow(case_norm_samples)),
  nrow = nrow(case_norm_samples),
  ncol = ncol(case_norm_samples)
)

combined_samples_s2 <- w_case_s2 * case_norm_samples +
  w_tick_s2 * tick_norm_samples +
  w_status_s2 * status_matrix_s2

# Save posterior mean and 95% credible intervals
data_aa$combined_mean_s2  <- apply(combined_samples_s2, 2, mean)
data_aa$combined_lower_s2 <- apply(combined_samples_s2, 2, quantile, probs = 0.025)
data_aa$combined_upper_s2 <- apply(combined_samples_s2, 2, quantile, probs = 0.975)
data_aa$combined_width_s2 <- data_aa$combined_upper_s2 - data_aa$combined_lower_s2

# Check top counties
data_aa %>%
  sf::st_drop_geometry() %>%
  dplyr::select(County, Status_aa, combined_mean_s2, combined_lower_s2, combined_upper_s2) %>%
  dplyr::arrange(dplyr::desc(combined_mean_s2)) %>%
  head(15)


# CSV output

table_s2 <- data_aa %>%
  sf::st_drop_geometry() %>%
  dplyr::arrange(dplyr::desc(combined_mean_s2)) %>%
  dplyr::select(County, combined_mean_s2, combined_lower_s2, combined_upper_s2) %>%
  dplyr::mutate(
    County = stringr::str_to_title(County),
    combined_mean_s2  = round(combined_mean_s2, 3),
    combined_lower_s2 = round(combined_lower_s2, 3),
    combined_upper_s2 = round(combined_upper_s2, 3)
  )

print(dplyr::as_tibble(table_s2), n = 102)

table_s2 %>%
  kable(
    col.names = c("County", "Risk score mean", "95% CrI lower", "95% CrI upper"),
    caption = "Supplementary file. Sensitivity analysis 2: disease-heavy weighting."
  )

write.csv(
  table_s2,
  "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Output/supplementary file (sensitivity analysis 2).csv",
  row.names = FALSE
)

# Map output

# Load Shawnee National Forest boundary
if (!exists("shawnee")) {
  shawnee <- sf::st_read(
    "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Data/forest shape file/shawnee_boundary.shp"
  )
  shawnee <- sf::st_transform(shawnee, sf::st_crs(data_aa))
}

map_s2 <- tm_shape(data_aa) +
  tm_fill(
    "combined_mean_s2",
    fill.scale = tm_scale_continuous(
      values = c("#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026"),
      trans = "sqrt"
    ),
    fill.legend = tm_legend(title = "Risk score", show.na = FALSE)
  ) +
  tm_borders(lwd = 0.3, col = "grey40") +
  tm_shape(shawnee) +
  tm_borders(lwd = 1.5, col = "darkgreen", lty = "solid") +
  tm_add_legend(
    type = "lines",
    labels = "Shawnee National Forest",
    col = "darkgreen",
    lwd = 1.5
  ) +
  tm_layout(
    legend.outside = TRUE,
    frame = FALSE,
    legend.frame = FALSE
  )

print(map_s2)

tmap_save(
  map_s2,
  filename = "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Output/illinois_combined_risk_map (sensitivity analysis 2).tiff",
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

# Forest plot output

plot_data_s2 <- data_aa %>%
  sf::st_drop_geometry() %>%
  dplyr::arrange(dplyr::desc(combined_mean_s2)) %>%
  head(20) %>%
  dplyr::mutate(County = stringr::str_to_title(County))

plot_data_s2$County <- factor(plot_data_s2$County, levels = rev(plot_data_s2$County))

forest_plot_s2 <- ggplot(plot_data_s2, aes(x = combined_mean_s2, y = County)) +
  geom_col(fill = "orangered4", width = 0.65, alpha = 0.85) +
  geom_errorbar(
    aes(xmin = combined_lower_s2, xmax = combined_upper_s2),
    width = 0.25,
    linewidth = 0.5,
    color = "#4A1B0C"
  ) +
  scale_x_continuous(limits = c(0, 1.05), expand = c(0.01, 0.01)) +
  labs(x = "Risk score", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 11, color = "grey20"),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(size = 11, color = "black", margin = margin(t = 10)),
    plot.margin = margin(15, 20, 15, 15)
  )

forest_plot_s2

ggsave(
  plot = forest_plot_s2,
  filename = "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Output/illinois_risk_forest_plot (sensitivity analysis 2).tiff",
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)








library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)

# Output folder
output_dir <- "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Output"

# File paths
primary_file <- file.path(output_dir, "Supplementary file.csv")
s1_file      <- file.path(output_dir, "supplementary file (sensitivity analysis 1).csv")
s2_file      <- file.path(output_dir, "supplementary file (sensitivity analysis 2).csv")

# Function to read risk-score CSV and standardize score column name
read_risk_file <- function(file, scenario_name) {
  
  df <- readr::read_csv(file, show_col_types = FALSE)
  
  # Detect the risk-score mean column automatically
  score_col <- grep("combined_mean|Risk score", names(df), value = TRUE)[1]
  
  df %>%
    dplyr::select(County, score = all_of(score_col)) %>%
    dplyr::mutate(
      County = stringr::str_to_title(County),
      scenario = scenario_name
    )
}

# Read primary and sensitivity files
primary_df <- read_risk_file(primary_file, "Primary")
s1_df      <- read_risk_file(s1_file, "Sensitivity 1")
s2_df      <- read_risk_file(s2_file, "Sensitivity 2")

# Combine into long and wide format
risk_long <- dplyr::bind_rows(primary_df, s1_df, s2_df)

risk_wide <- risk_long %>%
  tidyr::pivot_wider(
    names_from = scenario,
    values_from = score
  )

# Spearman correlation using all 102 counties

cor_all <- risk_wide %>%
  dplyr::select(Primary, `Sensitivity 1`, `Sensitivity 2`) %>%
  cor(method = "spearman", use = "complete.obs")

print(cor_all)


# Add confidence intervals for the two key comparisons reported in text
library(RVAideMemoire)

# All-counties comparisons (top row of your figure)
ci_all_p_s1 <- spearman.ci(risk_wide$Primary, risk_wide$`Sensitivity 1`, nrep = 1000)
ci_all_p_s2 <- spearman.ci(risk_wide$Primary, risk_wide$`Sensitivity 2`, nrep = 1000)
ci_all_s1_s2 <- spearman.ci(risk_wide$`Sensitivity 1`, risk_wide$`Sensitivity 2`, nrep = 1000)

# Top-10 comparisons (bottom row of your figure) — uses the same top10 county filter you already built
ci_top10_p_s1 <- spearman.ci(risk_top10$Primary, risk_top10$`Sensitivity 1`, nrep = 1000)
ci_top10_p_s2 <- spearman.ci(risk_top10$Primary, risk_top10$`Sensitivity 2`, nrep = 1000)
ci_top10_s1_s2 <- spearman.ci(risk_top10$`Sensitivity 1`, risk_top10$`Sensitivity 2`, nrep = 1000)

cat("=== All 102 counties ===\n")
print(ci_all_p_s1); print(ci_all_p_s2); print(ci_all_s1_s2)

cat("\n=== Top 10 counties ===\n")
print(ci_top10_p_s1); print(ci_top10_p_s2); print(ci_top10_s1_s2)

# Convert correlation matrix to long format for plotting
cor_all_long <- as.data.frame(as.table(cor_all)) %>%
  dplyr::rename(
    Scenario_1 = Var1,
    Scenario_2 = Var2,
    Spearman_rho = Freq
  )


# Top-10 overlap with primary model


top10_primary <- risk_wide %>%
  dplyr::arrange(dplyr::desc(Primary)) %>%
  dplyr::slice(1:10) %>%
  dplyr::pull(County)

top10_s1 <- risk_wide %>%
  dplyr::arrange(dplyr::desc(`Sensitivity 1`)) %>%
  dplyr::slice(1:10) %>%
  dplyr::pull(County)

top10_s2 <- risk_wide %>%
  dplyr::arrange(dplyr::desc(`Sensitivity 2`)) %>%
  dplyr::slice(1:10) %>%
  dplyr::pull(County)

top10_overlap <- data.frame(
  Comparison = c("Primary vs Sensitivity 1", "Primary vs Sensitivity 2"),
  Top10_overlap_count = c(
    length(intersect(top10_primary, top10_s1)),
    length(intersect(top10_primary, top10_s2))
  ),
  Top10_overlap_percent = c(
    length(intersect(top10_primary, top10_s1)) / 10 * 100,
    length(intersect(top10_primary, top10_s2)) / 10 * 100
  ),
  Overlapping_counties = c(
    paste(intersect(top10_primary, top10_s1), collapse = ", "),
    paste(intersect(top10_primary, top10_s2), collapse = ", ")
  )
)

print(top10_overlap)

# Spearman correlation restricted to union of top-10 counties

top10_union <- unique(c(top10_primary, top10_s1, top10_s2))

risk_top10_union <- risk_wide %>%
  dplyr::filter(County %in% top10_union)

cor_top10_union <- risk_top10_union %>%
  dplyr::select(Primary, `Sensitivity 1`, `Sensitivity 2`) %>%
  cor(method = "spearman", use = "complete.obs")

print(cor_top10_union)

cor_top10_long <- as.data.frame(as.table(cor_top10_union)) %>%
  dplyr::rename(
    Scenario_1 = Var1,
    Scenario_2 = Var2,
    Spearman_rho = Freq
  ) %>%
  dplyr::mutate(Analysis = "Top-10 counties")

cor_all_long <- cor_all_long %>%
  dplyr::mutate(Analysis = "All counties")

cor_plot_data <- dplyr::bind_rows(cor_all_long, cor_top10_long) %>%
  dplyr::mutate(
    Scenario_1 = factor(Scenario_1, levels = c("Primary", "Sensitivity 1", "Sensitivity 2")),
    Scenario_2 = factor(Scenario_2, levels = c("Primary", "Sensitivity 1", "Sensitivity 2"))
  )


# Correlation heatmap plot

cor_plot <- ggplot(cor_plot_data, aes(x = Scenario_1, y = Scenario_2, fill = Spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = round(Spearman_rho, 2)), color = "black", fontface = "bold", size = 4) +
  facet_wrap(~ Analysis) +
  scale_fill_gradient(
    low = "#FEE8C8",
    high = "#7F0000",
    limits = c(0, 1),
    name = "Spearman\nrho"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", angle = 35, hjust = 1),
    axis.text.y = element_text(color = "black"),
    strip.text = element_text(face = "bold", color = "black", size = 11),
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "black"),
    plot.margin = margin(10, 15, 10, 10)
  )

cor_plot

ggsave(
  plot = cor_plot,
  filename = file.path(output_dir, "sensitivity_spearman_correlation_plot.tiff"),
  width = 8,
  height = 4.8,
  units = "in",
  dpi = 300,
  compression = "lzw"
)



####_____________

## Pairwise sensitivity comparison scatter plots


library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)

# Output folder
output_dir <- "C:/Users/abrar/OneDrive - University of Illinois - Urbana/Documents/Abrar Hussain/Ph.D/Ph.D. Papers/As Authors/Alpha Gal Risk in Illinois/Output"

# File paths
primary_file <- file.path(output_dir, "Supplementary file.csv")
s1_file      <- file.path(output_dir, "supplementary file (sensitivity analysis 1).csv")
s2_file      <- file.path(output_dir, "supplementary file (sensitivity analysis 2).csv")

# Function to read risk-score CSV and standardize score column
read_risk_file <- function(file, scenario_name) {
  
  df <- readr::read_csv(file, show_col_types = FALSE)
  
  score_col <- grep("combined_mean|Risk score", names(df), value = TRUE)[1]
  
  df %>%
    dplyr::select(County, score = all_of(score_col)) %>%
    dplyr::mutate(
      County = stringr::str_to_title(County),
      scenario = scenario_name
    )
}

# Read files
primary_df <- read_risk_file(primary_file, "Primary")
s1_df      <- read_risk_file(s1_file, "Sensitivity 1")
s2_df      <- read_risk_file(s2_file, "Sensitivity 2")

# Combine and reshape
risk_long <- bind_rows(primary_df, s1_df, s2_df)

risk_wide <- risk_long %>%
  pivot_wider(
    names_from = scenario,
    values_from = score
  )

# Top 10 counties from the primary model
top10_primary <- risk_wide %>%
  arrange(desc(Primary)) %>%
  slice(1:10) %>%
  pull(County)

# Helper function to create pairwise data
make_pair_data <- function(df, xvar, yvar, comparison_name, analysis_name) {
  out <- df %>%
    dplyr::select(County, all_of(xvar), all_of(yvar)) %>%
    dplyr::rename(x = all_of(xvar), y = all_of(yvar)) %>%
    dplyr::mutate(
      Comparison = comparison_name,
      Analysis = analysis_name
    )
  
  rho_val <- suppressWarnings(cor(out$x, out$y, method = "spearman", use = "complete.obs"))
  
  out <- out %>%
    dplyr::mutate(
      rho = rho_val
    )
  
  return(out)
}

# All 102 counties
pair_all_1 <- make_pair_data(risk_wide, "Primary", "Sensitivity 1",
                             "Primary vs Sensitivity 1", "All 102 counties")

pair_all_2 <- make_pair_data(risk_wide, "Primary", "Sensitivity 2",
                             "Primary vs Sensitivity 2", "All 102 counties")

pair_all_3 <- make_pair_data(risk_wide, "Sensitivity 1", "Sensitivity 2",
                             "Sensitivity 1 vs Sensitivity 2", "All 102 counties")

# Top 10 counties from primary model
risk_top10 <- risk_wide %>%
  filter(County %in% top10_primary)

pair_top10_1 <- make_pair_data(risk_top10, "Primary", "Sensitivity 1",
                               "Primary vs Sensitivity 1", "Top 10 counties")

pair_top10_2 <- make_pair_data(risk_top10, "Primary", "Sensitivity 2",
                               "Primary vs Sensitivity 2", "Top 10 counties")

pair_top10_3 <- make_pair_data(risk_top10, "Sensitivity 1", "Sensitivity 2",
                               "Sensitivity 1 vs Sensitivity 2", "Top 10 counties")

# Combine all panels
pair_plot_data <- bind_rows(
  pair_all_1, pair_all_2, pair_all_3,
  pair_top10_1, pair_top10_2, pair_top10_3
) %>%
  mutate(
    Comparison = factor(
      Comparison,
      levels = c("Primary vs Sensitivity 1",
                 "Primary vs Sensitivity 2",
                 "Sensitivity 1 vs Sensitivity 2")
    ),
    Analysis = factor(
      Analysis,
      levels = c("All 102 counties", "Top 10 counties")
    )
  )

# Annotation data for Spearman rho
rho_labels <- pair_plot_data %>%
  group_by(Comparison, Analysis) %>%
  summarise(
    rho = unique(rho)[1],
    x_pos = min(x, na.rm = TRUE) + 0.00 * (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
    y_pos = max(y, na.rm = TRUE) - 0.02 * (max(y, na.rm = TRUE) - min(y, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("Spearman rho = ", round(rho, 2)))


# Plot
pair_plot <- ggplot(pair_plot_data, aes(x = x, y = y)) +
  geom_point(
    size = 2.1,
    alpha = 0.85,
    color = "#1F4E79"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey40"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8,
    color = "#B22222"
  ) +
  geom_text(
    data = rho_labels,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.6,
    color = "black",
    fontface = "bold"
  ) +
  facet_grid(Analysis ~ Comparison, scales = "free") +
  labs(
    x = "Risk score",
    y = "Risk score"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", face = "bold", size = 11),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold", color = "black", size = 11),
    plot.margin = margin(12, 12, 12, 12)
  )

pair_plot

ggsave(
  plot = pair_plot,
  filename = file.path(output_dir, "sensitivity_pairwise_comparison_plot.tiff"),
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 300,
  compression = "lzw"
)
