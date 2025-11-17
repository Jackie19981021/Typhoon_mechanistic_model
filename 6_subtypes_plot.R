# ============================================================================
# TYPHOON IMPACT ANALYSIS - WITH ATTACK RATE PERIODS
# COMPLETE ANALYSIS INCLUDING NEW ATTACK RATE VISUALIZATIONS
# ============================================================================
Sys.setenv('R_MAX_VSIZE' = 64 * 1024^3)

library(rstan)
library(ggplot2)
library(dplyr)
library(lubridate)
library(bayesplot)
library(tidyr)
library(reshape2)
library(patchwork)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# ============================================================================
# DATA LOADING AND PREPROCESSING
# ============================================================================
cat("=== Loading and preprocessing data ===\n")

data_path <- "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/HK_ILI_COVID_Sep.csv"
raw_data <- read.csv(data_path)
raw_data$date <- as.Date(raw_data$date, format = "%Y/%m/%d")

start_date <- as.Date("2023-01-07")
typhoon_start_date <- as.Date("2025-09-20")
end_date <- as.Date("2025-11-01")

fitting_data <- raw_data[raw_data$date >= start_date & raw_data$date < typhoon_start_date, ]
validation_data <- raw_data[raw_data$date >= typhoon_start_date & raw_data$date <= end_date, ]

fitting_data <- fitting_data[complete.cases(fitting_data), ]
validation_data <- validation_data[complete.cases(validation_data), ]

strain_names <- c("B", "H3", "H1", "COVID", "RSV", "HFMD")
fitting_data <- fitting_data[, c("date", strain_names)]
validation_data <- validation_data[, c("date", strain_names)]

cat("Fitting period:", as.character(range(fitting_data$date)), "\n")
cat("Validation period:", as.character(range(validation_data$date)), "\n")
cat("Total fitting weeks:", nrow(fitting_data), "\n")
cat("Total validation weeks:", nrow(validation_data), "\n")
cat("Strains:", paste(strain_names, collapse=", "), "\n\n")

# ============================================================================
# DEFINE ATTACK RATE PERIODS
# ============================================================================
cat("=== Defining Attack Rate Periods ===\n")

# Period 1: 2025-07-12 to 2025-09-27 (before typhoon)
period1_start_date <- as.Date("2025-07-12")
period1_end_date <- as.Date("2025-09-27")

# Period 2: 2025-07-12 to 2025-10-18 (including some forecast)
period2_start_date <- period1_start_date
period2_end_date <- as.Date("2025-10-18")

# Calculate week indices
period1_start_week <- which(fitting_data$date == period1_start_date)
period1_end_week <- nrow(fitting_data)  # This is 2025-09-27
period2_weeks_forecast <- as.numeric(difftime(period2_end_date, period1_end_date, units = "weeks"))

cat("Period 1 (Before Typhoon):\n")
cat("  Start:", as.character(period1_start_date), "(week", period1_start_week, ")\n")
cat("  End:", as.character(period1_end_date), "(week", period1_end_week, ")\n")
cat("  Duration:", period1_end_week - period1_start_week + 1, "weeks\n\n")

cat("Period 2 (Including Typhoon + Recovery):\n")
cat("  Start:", as.character(period2_start_date), "(week", period1_start_week, ")\n")
cat("  End:", as.character(period2_end_date), "\n")
cat("  Forecast weeks needed:", period2_weeks_forecast, "\n\n")

# ============================================================================
# STAN DATA PREPARATION
# ============================================================================
T_weeks <- nrow(fitting_data)
T_weeks_validation <- nrow(validation_data)
T_weeks_forecast <- 8
typhoon_weeks <- 1
N_strains <- 6
T_weeks_total <- T_weeks + T_weeks_forecast

cases_matrix <- as.matrix(fitting_data[, strain_names])
cases_matrix[cases_matrix < 0] <- 0
cases_matrix <- round(cases_matrix)

validation_matrix <- as.matrix(validation_data[, strain_names])
validation_matrix[validation_matrix < 0] <- NA
validation_matrix <- round(validation_matrix)

stan_data <- list(
  T_weeks = T_weeks,
  T_weeks_forecast = T_weeks_forecast,
  N_strains = N_strains,
  cases = cases_matrix,
  num_knots = 20,
  spline_degree = 3,
  population = 7524000,
  typhoon_weeks = typhoon_weeks,
  
  # NEW: Attack rate period parameters
  period1_start_week = period1_start_week,
  period1_end_week = period1_end_week,
  period2_weeks_forecast = period2_weeks_forecast
)

cat("Stan data prepared:\n")
cat("  period1_start_week:", stan_data$period1_start_week, "\n")
cat("  period1_end_week:", stan_data$period1_end_week, "\n")
cat("  period2_weeks_forecast:", stan_data$period2_weeks_forecast, "\n\n")

# ============================================================================
# MODEL COMPILATION AND FITTING
# ============================================================================
cat("=== Compiling Stan model ===\n")
model_file <- "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/6_subtypes_attack_rate.stan"
model <- stan_model(model_file)

cat("\n=== Starting MCMC sampling ===\n")
cat("Note: Computing 47 main scenarios + 281 extended scenarios\n")
cat("Plus attack rates for two specific time periods\n\n")

# Uncomment to run new fit:
fit <- sampling(
  model,
  data = stan_data,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  thin = 1,
  cores = 4,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)
saveRDS(fit, file = "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/typhoon_model_fit_attack_rate.rds")

# Or load previous fit:
# cat("=== LOADING PREVIOUS FIT ===\n")
# fit <- readRDS("/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/typhoon_model_fit_attack_rate.rds")

# ============================================================================
# EXTRACT RESULTS
# ============================================================================
cat("\n=== Extracting results ===\n")

pred_cases <- rstan::extract(fit, pars = "pred_cases")$pred_cases
forecast_cases <- rstan::extract(fit, pars = "forecast_cases")$forecast_cases
forecast_cases_typhoon <- rstan::extract(fit, pars = "forecast_cases_typhoon")$forecast_cases_typhoon
R_eff <- rstan::extract(fit, pars = "R_eff")$R_eff
R_eff_scenarios <- rstan::extract(fit, pars = "R_eff_scenarios")$R_eff_scenarios
R0_t <- rstan::extract(fit, pars = "R0_t")$R0_t
reduction_typhoon <- rstan::extract(fit, pars = "reduction_typhoon_period")$reduction_typhoon_period
avg_weekly_typhoon <- rstan::extract(fit, pars = "avg_weekly_typhoon")$avg_weekly_typhoon
avg_weekly_recovery <- rstan::extract(fit, pars = "avg_weekly_recovery")$avg_weekly_recovery

# Extended scenarios
cases_extended <- rstan::extract(fit, pars = "cases_extended_typhoon")$cases_extended_typhoon
incidence_per10k_extended <- rstan::extract(fit, pars = "incidence_per10k_extended")$incidence_per10k_extended
reduction_extended <- rstan::extract(fit, pars = "reduction_extended_typhoon")$reduction_extended_typhoon

# Original attack rates
attack_rate_baseline <- rstan::extract(fit, pars = "attack_rate_baseline")$attack_rate_baseline
attack_rate_scenarios <- rstan::extract(fit, pars = "attack_rate_scenarios")$attack_rate_scenarios
attack_rate_diff <- rstan::extract(fit, pars = "attack_rate_diff")$attack_rate_diff

# NEW: Attack rates for specific periods
attack_rate_period1_historical <- rstan::extract(fit, pars = "attack_rate_period1_historical")$attack_rate_period1_historical
attack_rate_period1_baseline <- rstan::extract(fit, pars = "attack_rate_period1_baseline")$attack_rate_period1_baseline
attack_rate_period1_scenarios <- rstan::extract(fit, pars = "attack_rate_period1_scenarios")$attack_rate_period1_scenarios
attack_rate_period1_extended <- rstan::extract(fit, pars = "attack_rate_period1_extended")$attack_rate_period1_extended

attack_rate_period2_historical <- rstan::extract(fit, pars = "attack_rate_period2_historical")$attack_rate_period2_historical
attack_rate_period2_baseline <- rstan::extract(fit, pars = "attack_rate_period2_baseline")$attack_rate_period2_baseline
attack_rate_period2_scenarios <- rstan::extract(fit, pars = "attack_rate_period2_scenarios")$attack_rate_period2_scenarios
attack_rate_period2_extended <- rstan::extract(fit, pars = "attack_rate_period2_extended")$attack_rate_period2_extended

cat("\n=== Checking extracted dimensions ===\n")
cat("pred_cases dimensions:", dim(pred_cases), "\n")
cat("forecast_cases dimensions:", dim(forecast_cases), "\n")
cat("forecast_cases_typhoon dimensions:", dim(forecast_cases_typhoon), "\n")
cat("R_eff dimensions:", dim(R_eff), "\n")
cat("R_eff_scenarios dimensions:", dim(R_eff_scenarios), "\n")
cat("R0_t dimensions:", dim(R0_t), "\n")
cat("reduction_typhoon dimensions:", dim(reduction_typhoon), "\n")
cat("cases_extended dimensions:", dim(cases_extended), "\n")
cat("reduction_extended dimensions:", dim(reduction_extended), "\n")
cat("attack_rate_baseline dimensions:", dim(attack_rate_baseline), "\n")
cat("attack_rate_scenarios dimensions:", dim(attack_rate_scenarios), "\n")
cat("attack_rate_diff dimensions:", dim(attack_rate_diff), "\n")
cat("attack_rate_period1_historical dimensions:", dim(attack_rate_period1_historical), "\n")
cat("attack_rate_period1_scenarios dimensions:", dim(attack_rate_period1_scenarios), "\n")
cat("attack_rate_period1_extended dimensions:", dim(attack_rate_period1_extended), "\n")
cat("attack_rate_period2_historical dimensions:", dim(attack_rate_period2_historical), "\n")
cat("attack_rate_period2_scenarios dimensions:", dim(attack_rate_period2_scenarios), "\n")
cat("attack_rate_period2_extended dimensions:", dim(attack_rate_period2_extended), "\n")

# ============================================================================
# DEFINE SCENARIO LABELS
# ============================================================================
typhoon_scenarios <- c(
  "0% (Baseline)",
  # 3 days
  "3d -50%", "3d -30%", "3d -20%", "3d -10%", "3d -5%",
  "3d +5%", "3d +10%", "3d +20%", "3d +30%", "3d +50%",
  # 7 days
  "7d -50%", "7d -30%", "7d -20%", "7d -10%", "7d -5%",
  "7d +5%", "7d +10%", "7d +20%", "7d +30%", "7d +50%",
  # 10 days
  "10d -50%", "10d -30%", "10d -20%", "10d -10%", "10d -5%",
  "10d +5%", "10d +10%", "10d +20%", "10d +30%", "10d +50%",
  # 14 days
  "14d -50%", "14d -30%", "14d -20%", "14d -10%", "14d -5%",
  "14d +5%", "14d +10%", "14d +20%", "14d +30%", "14d +50%",
  # Shifted timing (7 days, -50%)
  "7d -50% E7", "7d -50% E5", "7d -50% E3",
  "7d -50% D3", "7d -50% D5", "7d -50% D7"
)

# ============================================================================
# CALCULATE STATISTICS
# ============================================================================
cat("\n=== Calculating statistics ===\n")

pred_mean <- apply(pred_cases, c(2,3), mean)
pred_median <- apply(pred_cases, c(2,3), median)
pred_lower <- apply(pred_cases, c(2,3), quantile, probs = 0.025)
pred_upper <- apply(pred_cases, c(2,3), quantile, probs = 0.975)

forecast_mean <- apply(forecast_cases, c(2,3), mean)
forecast_median <- apply(forecast_cases, c(2,3), median)
forecast_lower <- apply(forecast_cases, c(2,3), quantile, probs = 0.025)
forecast_upper <- apply(forecast_cases, c(2,3), quantile, probs = 0.975)

# 47 typhoon scenarios
forecast_typhoon_mean <- array(NA, dim = c(47, T_weeks_forecast, N_strains))
forecast_typhoon_lower <- array(NA, dim = c(47, T_weeks_forecast, N_strains))
forecast_typhoon_upper <- array(NA, dim = c(47, T_weeks_forecast, N_strains))

for (typhoon_idx in 1:47) {
  forecast_typhoon_mean[typhoon_idx,,] <- apply(forecast_cases_typhoon[,typhoon_idx,,], c(2,3), mean)
  forecast_typhoon_lower[typhoon_idx,,] <- apply(forecast_cases_typhoon[,typhoon_idx,,], c(2,3), quantile, probs = 0.025)
  forecast_typhoon_upper[typhoon_idx,,] <- apply(forecast_cases_typhoon[,typhoon_idx,,], c(2,3), quantile, probs = 0.975)
}

R_eff_scenarios_mean <- array(NA, dim = c(47, T_weeks + T_weeks_forecast, N_strains))
R_eff_scenarios_lower <- array(NA, dim = c(47, T_weeks + T_weeks_forecast, N_strains))
R_eff_scenarios_upper <- array(NA, dim = c(47, T_weeks + T_weeks_forecast, N_strains))

for (typhoon_idx in 1:47) {
  R_eff_scenarios_mean[typhoon_idx,,] <- apply(R_eff_scenarios[,typhoon_idx,,], c(2,3), mean)
  R_eff_scenarios_lower[typhoon_idx,,] <- apply(R_eff_scenarios[,typhoon_idx,,], c(2,3), quantile, probs = 0.025)
  R_eff_scenarios_upper[typhoon_idx,,] <- apply(R_eff_scenarios[,typhoon_idx,,], c(2,3), quantile, probs = 0.975)
}

# Original attack rate statistics
attack_rate_baseline_mean <- apply(attack_rate_baseline, 2, mean)
attack_rate_baseline_lower <- apply(attack_rate_baseline, 2, quantile, probs = 0.025)
attack_rate_baseline_upper <- apply(attack_rate_baseline, 2, quantile, probs = 0.975)

attack_rate_scenarios_mean <- apply(attack_rate_scenarios, c(2,3), mean)
attack_rate_scenarios_lower <- apply(attack_rate_scenarios, c(2,3), quantile, probs = 0.025)
attack_rate_scenarios_upper <- apply(attack_rate_scenarios, c(2,3), quantile, probs = 0.975)

attack_rate_diff_mean <- apply(attack_rate_diff, c(2,3), mean)
attack_rate_diff_lower <- apply(attack_rate_diff, c(2,3), quantile, probs = 0.025)
attack_rate_diff_upper <- apply(attack_rate_diff, c(2,3), quantile, probs = 0.975)

# NEW: Attack rate statistics for periods
# Period 1
attack_rate_period1_historical_mean <- apply(attack_rate_period1_historical, 2, mean)
attack_rate_period1_baseline_mean <- apply(attack_rate_period1_baseline, 2, mean)
attack_rate_period1_scenarios_mean <- apply(attack_rate_period1_scenarios, c(2,3), mean)
attack_rate_period1_extended_mean <- apply(attack_rate_period1_extended, c(2,3), mean)

# Period 2
attack_rate_period2_historical_mean <- apply(attack_rate_period2_historical, 2, mean)
attack_rate_period2_baseline_mean <- apply(attack_rate_period2_baseline, 2, mean)
attack_rate_period2_scenarios_mean <- apply(attack_rate_period2_scenarios, c(2,3), mean)
attack_rate_period2_extended_mean <- apply(attack_rate_period2_extended, c(2,3), mean)

# Extended scenarios mean
incidence_per10k_mean <- apply(incidence_per10k_extended, c(2,3), mean)
reduction_extended_mean <- apply(reduction_extended, c(2,3), mean)

cat("Statistics calculated successfully\n\n")

# ============================================================================
# PREPARE VISUALIZATION DATA
# ============================================================================
theme_set(theme_minimal(base_size = 11))

forecast_dates <- seq(typhoon_start_date, by = "week", length.out = T_weeks_forecast)
all_dates <- c(fitting_data$date, forecast_dates)

# Date markers for visualization
typhoon_start_date_line <- max(fitting_data$date)
typhoon_end_date_line <- max(fitting_data$date) + 7
typhoon_period_end <- max(fitting_data$date) + 7
recovery_start_date <- max(fitting_data$date) + 7

cat("\n=== Date markers ===\n")
cat("Last fitting date:", as.character(max(fitting_data$date)), "\n")
cat("First forecast date:", as.character(min(forecast_dates)), "\n")
cat("Typhoon start date (data split):", as.character(typhoon_start_date), "\n")
cat("Typhoon start line (visualization):", as.character(typhoon_start_date_line), "\n")
cat("Typhoon end (red line):", as.character(typhoon_end_date_line), "\n")
cat("Recovery start:", as.character(recovery_start_date), "\n\n")

# Prepare data frames
results_df <- data.frame(
  week = rep(1:T_weeks, N_strains),
  date = rep(fitting_data$date, N_strains),
  strain = factor(rep(strain_names, each = T_weeks), levels = strain_names),
  observed = as.vector(cases_matrix),
  predicted = as.vector(pred_median),
  pred_mean = as.vector(pred_mean),
  pred_lower = as.vector(pred_lower),
  pred_upper = as.vector(pred_upper)
)

last_fit_point <- results_df %>%
  group_by(strain) %>%
  slice_tail(n = 1) %>%
  ungroup()

forecast_baseline_raw <- data.frame(
  week = rep((T_weeks + 1):(T_weeks + T_weeks_forecast), N_strains),
  date = rep(forecast_dates, N_strains),
  strain = factor(rep(strain_names, each = T_weeks_forecast), levels = strain_names),
  observed = NA,
  predicted = as.vector(forecast_median),
  pred_mean = as.vector(forecast_mean),
  pred_lower = as.vector(forecast_lower),
  pred_upper = as.vector(forecast_upper)
)

forecast_baseline <- rbind(
  last_fit_point %>% select(week, date, strain, observed, predicted, pred_mean, pred_lower, pred_upper),
  forecast_baseline_raw
)

typhoon_forecast_list <- list()
for (typhoon_idx in 1:47) {
  for (strain_idx in 1:N_strains) {
    last_point <- last_fit_point %>% filter(strain == strain_names[strain_idx])
    
    scenario_forecast <- data.frame(
      week = (T_weeks + 1):(T_weeks + T_weeks_forecast),
      date = forecast_dates,
      strain = strain_names[strain_idx],
      observed = NA,
      predicted = forecast_typhoon_mean[typhoon_idx, , strain_idx],
      pred_mean = forecast_typhoon_mean[typhoon_idx, , strain_idx],
      pred_lower = forecast_typhoon_lower[typhoon_idx, , strain_idx],
      pred_upper = forecast_typhoon_upper[typhoon_idx, , strain_idx],
      scenario = typhoon_scenarios[typhoon_idx]
    )
    
    last_point$scenario <- typhoon_scenarios[typhoon_idx]
    scenario_with_overlap <- rbind(
      last_point %>% select(week, date, strain, observed, predicted, pred_mean, pred_lower, pred_upper, scenario),
      scenario_forecast
    )
    
    typhoon_forecast_list[[length(typhoon_forecast_list) + 1]] <- scenario_with_overlap
  }
}

typhoon_forecast_df <- do.call(rbind, typhoon_forecast_list)
typhoon_forecast_df$strain <- factor(typhoon_forecast_df$strain, levels = strain_names)

validation_plot_df <- data.frame(
  date = rep(validation_data$date, N_strains),
  strain = factor(rep(strain_names, each = T_weeks_validation), levels = strain_names),
  observed = as.vector(validation_matrix)
)
validation_plot_df <- validation_plot_df[!is.na(validation_plot_df$observed), ]

cat("Visualization data prepared\n\n")

# ============================================================================
# FIGURE 1: Historical Fit + Baseline Forecast + Validation Data
# ============================================================================
cat("\n=== Creating Figure 1: Fit, Baseline Forecast, and Validation ===\n")

fit_and_forecast_df <- rbind(
  cbind(results_df, period = "Fitted"),
  cbind(forecast_baseline, period = "Forecast")
)

p1_fit_forecast <- ggplot(fit_and_forecast_df, aes(x = date)) +
  annotate("rect",
           xmin = typhoon_start_date_line,
           xmax = typhoon_end_date_line,
           ymin = -Inf, ymax = Inf,
           fill = "pink", alpha = 0.15) +
  annotate("rect",
           xmin = recovery_start_date,
           xmax = max(forecast_dates),
           ymin = -Inf, ymax = Inf,
           fill = "lightgreen", alpha = 0.1) +
  geom_ribbon(data = filter(fit_and_forecast_df, period == "Fitted"),
              aes(ymin = pred_lower, ymax = pred_upper),
              alpha = 0.25, fill = "#377EB8") +
  geom_line(data = filter(fit_and_forecast_df, period == "Fitted"),
            aes(y = predicted), color = "#377EB8", size = 1.3) +
  geom_ribbon(data = filter(fit_and_forecast_df, period == "Forecast"),
              aes(ymin = pred_lower, ymax = pred_upper),
              alpha = 0.25, fill = "#E41A1C") +
  geom_line(data = filter(fit_and_forecast_df, period == "Forecast"),
            aes(y = predicted), color = "#E41A1C", size = 1.3) +
  geom_point(data = filter(fit_and_forecast_df, !is.na(observed), period == "Fitted"),
             aes(y = observed), color = "black", size = 1.2, alpha = 0.6) +
  geom_point(data = validation_plot_df,
             aes(x = date, y = observed),
             color = "gold", fill = "yellow", shape = 23, size = 1.2, stroke = 0.8) +
  geom_vline(xintercept = typhoon_start_date_line,
             linetype = "dashed", color = "gray40", size = 0.8) +
  geom_vline(xintercept = typhoon_end_date_line,
             linetype = "dashed", color = "red", size = 0.8) +
  facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
  scale_y_continuous(trans = "sqrt") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "4 months") +
  labs(
    title = "Typhoon Impact Simulation: Historical Fit and Baseline Forecast",
    subtitle = "Black dots: Observed (fitted) | Yellow diamonds: Validation data | Blue: Model fit | Red: Baseline forecast\nGray dashed: Typhoon start | Red dashed: Typhoon end",
    x = NULL,
    y = "Weekly Cases/Hospitalizations (√ scale)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "gray30")
  )

print(p1_fit_forecast)
ggsave("figure1_typhoon_baseline_forecast_validation.png", p1_fit_forecast,
       width = 14, height = 9, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 2: Recent Period + Selected Typhoon Scenarios (Multiple Duration Plots)
# ============================================================================
cat("\n=== Creating Figure 2: Recent Period + Scenarios ===\n")

cutoff_date <- as.Date("2025-09-01")
recent_fit_df <- results_df %>%
  filter(date >= cutoff_date) %>%
  mutate(scenario = "Historical Fit")

recent_and_scenarios <- rbind(recent_fit_df, typhoon_forecast_df)
recent_and_scenarios$scenario <- factor(
  recent_and_scenarios$scenario,
  levels = c("Historical Fit", typhoon_scenarios)
)

last_observed <- recent_fit_df %>%
  group_by(strain) %>%
  slice_tail(n = 1) %>%
  ungroup()

# Function to create scenario plot
create_scenario_plot <- function(selected_scenarios, plot_title, colors, linetypes) {
  df_filtered <- recent_and_scenarios %>% filter(scenario %in% selected_scenarios)
  df_filtered$scenario <- factor(df_filtered$scenario, levels = selected_scenarios)
  
  p <- ggplot(df_filtered, aes(x = date)) +
    annotate("rect",
             xmin = typhoon_start_date_line,
             xmax = typhoon_end_date_line,
             ymin = -Inf, ymax = Inf,
             fill = "pink", alpha = 0.15) +
    annotate("rect",
             xmin = recovery_start_date,
             xmax = max(forecast_dates),
             ymin = -Inf, ymax = Inf,
             fill = "lightgreen", alpha = 0.1) +
    geom_ribbon(data = filter(df_filtered, scenario == "Historical Fit"),
                aes(ymin = pred_lower, ymax = pred_upper),
                alpha = 0.2, fill = "#377EB8") +
    geom_line(data = filter(df_filtered, scenario == "Historical Fit"),
              aes(y = predicted), color = "#377EB8", size = 1.4) +
    geom_point(data = filter(df_filtered, scenario == "Historical Fit", !is.na(observed)),
               aes(y = observed), color = "black", size = 1.3, alpha = 0.6) +
    geom_point(data = last_observed, aes(y = observed),
               color = "black", size = 3, shape = 21, fill = "white", stroke = 1.3) +
    geom_point(data = validation_plot_df,
               aes(x = date, y = observed),
               color = "black", fill = "yellow", shape = 23, size = 2.0, stroke = 1.2) +
    geom_ribbon(data = filter(df_filtered, scenario != "Historical Fit"),
                aes(ymin = pred_lower, ymax = pred_upper, fill = scenario),
                alpha = 0.18) +
    geom_line(data = filter(df_filtered, scenario != "Historical Fit"),
              aes(y = predicted, color = scenario, linetype = scenario),
              size = 1.2) +
    geom_vline(xintercept = typhoon_start_date_line,
               linetype = "dashed", color = "gray40", size = 0.8) +
    geom_vline(xintercept = typhoon_end_date_line,
               linetype = "dashed", color = "red", size = 0.8) +
    facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
    scale_color_manual(name = "Scenarios:", values = colors) +
    scale_fill_manual(name = "Scenarios:", values = colors) +
    scale_linetype_manual(name = "Scenarios:", values = linetypes) +
    scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
    labs(
      title = plot_title,
      subtitle = "Blue: Model fit | Yellow diamonds: Validation data | Colored lines: Forecast scenarios\nGray dashed: Typhoon start | Red dashed: Typhoon end | Pink: Typhoon (1 week) | Green: Recovery (7 weeks)",
      x = NULL,
      y = "Weekly Cases/Hospitalizations"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 8.5, color = "gray30")
    ) +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE),
      fill = guide_legend(nrow = 2, byrow = TRUE),
      linetype = guide_legend(nrow = 2, byrow = TRUE)
    )
  
  return(p)
}

# Define colors and linetypes for each duration
colors_3d <- c("0% (Baseline)" = "#999999",
               "3d -50%" = "#1B4F72", "3d -30%" = "#21618C", "3d -20%" = "#3498DB",
               "3d -10%" = "#5DADE2", "3d -5%" = "#AED6F1",
               "3d +5%" = "#FADBD8", "3d +10%" = "#F5B7B1",
               "3d +20%" = "#E74C3C", "3d +30%" = "#C0392B", "3d +50%" = "#922B21")

linetypes_3d <- c("0% (Baseline)" = "solid",
                  "3d -50%" = "solid", "3d -30%" = "dashed", "3d -20%" = "dotted",
                  "3d -10%" = "dotdash", "3d -5%" = "longdash",
                  "3d +5%" = "longdash", "3d +10%" = "dotdash",
                  "3d +20%" = "dotted", "3d +30%" = "dashed", "3d +50%" = "solid")

selected_3d <- c("Historical Fit", "0% (Baseline)", 
                 "3d -50%", "3d -30%", "3d -20%", "3d -10%", "3d -5%",
                 "3d +5%", "3d +10%", "3d +20%", "3d +30%", "3d +50%")

p2_3days <- create_scenario_plot(selected_3d, 
                                 "Typhoon Impact: Recent Period and 3 Days Duration Scenarios", 
                                 colors_3d, linetypes_3d)
print(p2_3days)
ggsave("figure2_typhoon_scenarios_3days.png", p2_3days, width = 14, height = 10, dpi = 300, bg = "white")

# 7 days
colors_7d <- c("0% (Baseline)" = "#999999",
               "7d -50%" = "#1B4F72", "7d -30%" = "#21618C", "7d -20%" = "#3498DB",
               "7d -10%" = "#5DADE2", "7d -5%" = "#AED6F1",
               "7d +5%" = "#FADBD8", "7d +10%" = "#F5B7B1",
               "7d +20%" = "#E74C3C", "7d +30%" = "#C0392B", "7d +50%" = "#922B21")

linetypes_7d <- c("0% (Baseline)" = "solid",
                  "7d -50%" = "solid", "7d -30%" = "dashed", "7d -20%" = "dotted",
                  "7d -10%" = "dotdash", "7d -5%" = "longdash",
                  "7d +5%" = "longdash", "7d +10%" = "dotdash",
                  "7d +20%" = "dotted", "7d +30%" = "dashed", "7d +50%" = "solid")

selected_7d <- c("Historical Fit", "0% (Baseline)", 
                 "7d -50%", "7d -30%", "7d -20%", "7d -10%", "7d -5%",
                 "7d +5%", "7d +10%", "7d +20%", "7d +30%", "7d +50%")

p2_7days <- create_scenario_plot(selected_7d, 
                                 "Typhoon Impact: Recent Period and 7 Days Duration Scenarios", 
                                 colors_7d, linetypes_7d)
print(p2_7days)
ggsave("figure2_typhoon_scenarios_7days.png", p2_7days, width = 14, height = 10, dpi = 300, bg = "white")

# 10 days
colors_10d <- c("0% (Baseline)" = "#999999",
                "10d -50%" = "#1B4F72", "10d -30%" = "#21618C", "10d -20%" = "#3498DB",
                "10d -10%" = "#5DADE2", "10d -5%" = "#AED6F1",
                "10d +5%" = "#FADBD8", "10d +10%" = "#F5B7B1",
                "10d +20%" = "#E74C3C", "10d +30%" = "#C0392B", "10d +50%" = "#922B21")

linetypes_10d <- c("0% (Baseline)" = "solid",
                   "10d -50%" = "solid", "10d -30%" = "dashed", "10d -20%" = "dotted",
                   "10d -10%" = "dotdash", "10d -5%" = "longdash",
                   "10d +5%" = "longdash", "10d +10%" = "dotdash",
                   "10d +20%" = "dotted", "10d +30%" = "dashed", "10d +50%" = "solid")

selected_10d <- c("Historical Fit", "0% (Baseline)", 
                  "10d -50%", "10d -30%", "10d -20%", "10d -10%", "10d -5%",
                  "10d +5%", "10d +10%", "10d +20%", "10d +30%", "10d +50%")

p2_10days <- create_scenario_plot(selected_10d, 
                                  "Typhoon Impact: Recent Period and 10 Days Duration Scenarios", 
                                  colors_10d, linetypes_10d)
print(p2_10days)
ggsave("figure2_typhoon_scenarios_10days.png", p2_10days, width = 14, height = 10, dpi = 300, bg = "white")

# 14 days
colors_14d <- c("0% (Baseline)" = "#999999",
                "14d -50%" = "#1B4F72", "14d -30%" = "#21618C", "14d -20%" = "#3498DB",
                "14d -10%" = "#5DADE2", "14d -5%" = "#AED6F1",
                "14d +5%" = "#FADBD8", "14d +10%" = "#F5B7B1",
                "14d +20%" = "#E74C3C", "14d +30%" = "#C0392B", "14d +50%" = "#922B21")

linetypes_14d <- c("0% (Baseline)" = "solid",
                   "14d -50%" = "solid", "14d -30%" = "dashed", "14d -20%" = "dotted",
                   "14d -10%" = "dotdash", "14d -5%" = "longdash",
                   "14d +5%" = "longdash", "14d +10%" = "dotdash",
                   "14d +20%" = "dotted", "14d +30%" = "dashed", "14d +50%" = "solid")

selected_14d <- c("Historical Fit", "0% (Baseline)", 
                  "14d -50%", "14d -30%", "14d -20%", "14d -10%", "14d -5%",
                  "14d +5%", "14d +10%", "14d +20%", "14d +30%", "14d +50%")

p2_14days <- create_scenario_plot(selected_14d, 
                                  "Typhoon Impact: Recent Period and 14 Days Duration Scenarios", 
                                  colors_14d, linetypes_14d)
print(p2_14days)
ggsave("figure2_typhoon_scenarios_14days.png", p2_14days, width = 14, height = 10, dpi = 300, bg = "white")

# Shifted scenarios
selected_shifted <- c("Historical Fit", "0% (Baseline)", 
                      "7d -50% E7", "7d -50% E5", "7d -50% E3",
                      "7d -50% D3", "7d -50% D5", "7d -50% D7")

colors_shifted <- c("0% (Baseline)" = "#999999",
                    "7d -50% E7" = "#1B4F72", "7d -50% E5" = "#21618C", "7d -50% E3" = "#3498DB",
                    "7d -50% D3" = "#E74C3C", "7d -50% D5" = "#C0392B", "7d -50% D7" = "#922B21")

linetypes_shifted <- c("0% (Baseline)" = "solid",
                       "7d -50% E7" = "solid", "7d -50% E5" = "dashed", "7d -50% E3" = "dotted",
                       "7d -50% D3" = "dotted", "7d -50% D5" = "dashed", "7d -50% D7" = "solid")

p2_shifted <- create_scenario_plot(selected_shifted, 
                                   "Typhoon Impact: Recent Period and Shifted Timing Scenarios (7d -50%)", 
                                   colors_shifted, linetypes_shifted)
print(p2_shifted)
ggsave("figure2_typhoon_scenarios_shifted.png", p2_shifted, width = 14, height = 10, dpi = 300, bg = "white")

cat("Figures 1-2 complete\n\n")






# FIGURES 3-6: Baseline vs Typhoon Forecasts (WITH CONTINUITY)
# ============================================================================

# FIGURE 3: 3 days duration
cat("\n=== Creating Figure 3: Baseline vs Typhoon Forecasts - 3 days (WITH CONTINUITY) ===\n")

forecast_3days <- typhoon_forecast_df %>%
  filter(scenario %in% c("0% (Baseline)", "3d -50%", "3d -30%", "3d -20%", "3d -10%", "3d -5%",
                         "3d +5%", "3d +10%", "3d +20%", "3d +30%", "3d +50%"))

forecast_3days$scenario <- factor(forecast_3days$scenario, 
                                  levels = c("0% (Baseline)", "3d -50%", "3d -30%", "3d -20%", "3d -10%", "3d -5%",
                                             "3d +5%", "3d +10%", "3d +20%", "3d +30%", "3d +50%"))

p3_3days <- ggplot(forecast_3days, aes(x = date, y = predicted, color = scenario, fill = scenario, linetype = scenario)) +
  annotate("rect",
           xmin = typhoon_start_date_line,
           xmax = typhoon_end_date_line,
           ymin = -Inf, ymax = Inf,
           fill = "pink", alpha = 0.15) +
  annotate("rect",
           xmin = recovery_start_date,
           xmax = max(forecast_dates),
           ymin = -Inf, ymax = Inf,
           fill = "lightgreen", alpha = 0.1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper), alpha = 0.18) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = typhoon_start_date_line,
             linetype = "dashed", color = "gray40", size = 0.8) +
  geom_vline(xintercept = typhoon_end_date_line,
             linetype = "dashed", color = "red", size = 0.8) +
  geom_point(data = validation_plot_df,
             aes(x = date, y = observed),
             color = "black", fill = "yellow", shape = 23, size = 2.0, stroke = 1.2, inherit.aes = FALSE) +
  facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
  scale_color_manual(values = colors_3d) +
  scale_fill_manual(values = colors_3d) +
  scale_linetype_manual(values = linetypes_3d) +
  scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
  labs(
    title = "Typhoon Impact Forecasts: 3 Days Duration Scenarios",
    subtitle = "Yellow diamonds: Validation data | Continuous lines from last fitted point | Gray dashed: Typhoon start | Red dashed: Typhoon end\nPink: Typhoon period | Green: Recovery period",
    x = NULL,
    y = "Weekly Cases/Hospitalizations"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "gray30")
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )

print(p3_3days)
ggsave("figure3_typhoon_forecasts_3days.png", p3_3days, width = 14, height = 10, dpi = 300, bg = "white")

# FIGURE 4: 7 days duration
cat("\n=== Creating Figure 4: Baseline vs Typhoon Forecasts - 7 days (WITH CONTINUITY) ===\n")

forecast_7days <- typhoon_forecast_df %>%
  filter(scenario %in% c("0% (Baseline)", "7d -50%", "7d -30%", "7d -20%", "7d -10%", "7d -5%",
                         "7d +5%", "7d +10%", "7d +20%", "7d +30%", "7d +50%"))

forecast_7days$scenario <- factor(forecast_7days$scenario, 
                                  levels = c("0% (Baseline)", "7d -50%", "7d -30%", "7d -20%", "7d -10%", "7d -5%",
                                             "7d +5%", "7d +10%", "7d +20%", "7d +30%", "7d +50%"))

p4_7days <- ggplot(forecast_7days, aes(x = date, y = predicted, color = scenario, fill = scenario, linetype = scenario)) +
  annotate("rect",
           xmin = typhoon_start_date_line,
           xmax = typhoon_end_date_line,
           ymin = -Inf, ymax = Inf,
           fill = "pink", alpha = 0.15) +
  annotate("rect",
           xmin = recovery_start_date,
           xmax = max(forecast_dates),
           ymin = -Inf, ymax = Inf,
           fill = "lightgreen", alpha = 0.1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper), alpha = 0.18) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = typhoon_start_date_line,
             linetype = "dashed", color = "gray40", size = 0.8) +
  geom_vline(xintercept = typhoon_end_date_line,
             linetype = "dashed", color = "red", size = 0.8) +
  geom_point(data = validation_plot_df,
             aes(x = date, y = observed),
             color = "black", fill = "yellow", shape = 23, size = 2.0, stroke = 1.2, inherit.aes = FALSE) +
  facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
  scale_color_manual(values = colors_7d) +
  scale_fill_manual(values = colors_7d) +
  scale_linetype_manual(values = linetypes_7d) +
  scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
  labs(
    title = "Typhoon Impact Forecasts: 7 Days Duration Scenarios",
    subtitle = "Yellow diamonds: Validation data | Continuous lines from last fitted point | Gray dashed: Typhoon start | Red dashed: Typhoon end\nPink: Typhoon period | Green: Recovery period",
    x = NULL,
    y = "Weekly Cases/Hospitalizations"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "gray30")
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )

print(p4_7days)
ggsave("figure4_typhoon_forecasts_7days.png", p4_7days, width = 14, height = 10, dpi = 300, bg = "white")

# FIGURE 5: 10 days duration
cat("\n=== Creating Figure 5: Baseline vs Typhoon Forecasts - 10 days (WITH CONTINUITY) ===\n")

forecast_10days <- typhoon_forecast_df %>%
  filter(scenario %in% c("0% (Baseline)", "10d -50%", "10d -30%", "10d -20%", "10d -10%", "10d -5%",
                         "10d +5%", "10d +10%", "10d +20%", "10d +30%", "10d +50%"))

forecast_10days$scenario <- factor(forecast_10days$scenario, 
                                   levels = c("0% (Baseline)", "10d -50%", "10d -30%", "10d -20%", "10d -10%", "10d -5%",
                                              "10d +5%", "10d +10%", "10d +20%", "10d +30%", "10d +50%"))

p5_10days <- ggplot(forecast_10days, aes(x = date, y = predicted, color = scenario, fill = scenario, linetype = scenario)) +
  annotate("rect",
           xmin = typhoon_start_date_line,
           xmax = typhoon_end_date_line,
           ymin = -Inf, ymax = Inf,
           fill = "pink", alpha = 0.15) +
  annotate("rect",
           xmin = recovery_start_date,
           xmax = max(forecast_dates),
           ymin = -Inf, ymax = Inf,
           fill = "lightgreen", alpha = 0.1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper), alpha = 0.18) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = typhoon_start_date_line,
             linetype = "dashed", color = "gray40", size = 0.8) +
  geom_vline(xintercept = typhoon_end_date_line,
             linetype = "dashed", color = "red", size = 0.8) +
  geom_point(data = validation_plot_df,
             aes(x = date, y = observed),
             color = "black", fill = "yellow", shape = 23, size = 2.0, stroke = 1.2, inherit.aes = FALSE) +
  facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
  scale_color_manual(values = colors_10d) +
  scale_fill_manual(values = colors_10d) +
  scale_linetype_manual(values = linetypes_10d) +
  scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
  labs(
    title = "Typhoon Impact Forecasts: 10 Days Duration Scenarios",
    subtitle = "Yellow diamonds: Validation data | Continuous lines from last fitted point | Gray dashed: Typhoon start | Red dashed: Typhoon end\nPink: Typhoon period | Green: Recovery period",
    x = NULL,
    y = "Weekly Cases/Hospitalizations"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "gray30")
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )

print(p5_10days)
ggsave("figure5_typhoon_forecasts_10days.png", p5_10days, width = 14, height = 10, dpi = 300, bg = "white")

# FIGURE 6: 14 days duration
cat("\n=== Creating Figure 6: Baseline vs Typhoon Forecasts - 14 days (WITH CONTINUITY) ===\n")

forecast_14days <- typhoon_forecast_df %>%
  filter(scenario %in% c("0% (Baseline)", "14d -50%", "14d -30%", "14d -20%", "14d -10%", "14d -5%",
                         "14d +5%", "14d +10%", "14d +20%", "14d +30%", "14d +50%"))

forecast_14days$scenario <- factor(forecast_14days$scenario, 
                                   levels = c("0% (Baseline)", "14d -50%", "14d -30%", "14d -20%", "14d -10%", "14d -5%",
                                              "14d +5%", "14d +10%", "14d +20%", "14d +30%", "14d +50%"))

p6_14days <- ggplot(forecast_14days, aes(x = date, y = predicted, color = scenario, fill = scenario, linetype = scenario)) +
  annotate("rect",
           xmin = typhoon_start_date_line,
           xmax = typhoon_end_date_line,
           ymin = -Inf, ymax = Inf,
           fill = "pink", alpha = 0.15) +
  annotate("rect",
           xmin = recovery_start_date,
           xmax = max(forecast_dates),
           ymin = -Inf, ymax = Inf,
           fill = "lightgreen", alpha = 0.1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper), alpha = 0.18) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = typhoon_start_date_line,
             linetype = "dashed", color = "gray40", size = 0.8) +
  geom_vline(xintercept = typhoon_end_date_line,
             linetype = "dashed", color = "red", size = 0.8) +
  geom_point(data = validation_plot_df,
             aes(x = date, y = observed),
             color = "black", fill = "yellow", shape = 23, size = 2.0, stroke = 1.2, inherit.aes = FALSE) +
  facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
  scale_color_manual(values = colors_14d) +
  scale_fill_manual(values = colors_14d) +
  scale_linetype_manual(values = linetypes_14d) +
  scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
  labs(
    title = "Typhoon Impact Forecasts: 14 Days Duration Scenarios",
    subtitle = "Yellow diamonds: Validation data | Continuous lines from last fitted point | Gray dashed: Typhoon start | Red dashed: Typhoon end\nPink: Typhoon period | Green: Recovery period",
    x = NULL,
    y = "Weekly Cases/Hospitalizations"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "gray30")
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )

print(p6_14days)
ggsave("figure6_typhoon_forecasts_14days.png", p6_14days, width = 14, height = 10, dpi = 300, bg = "white")

# FIGURE 6b: Shifted scenarios
cat("\n=== Creating Figure 6b: Baseline vs Typhoon Forecasts - Shifted (WITH CONTINUITY) ===\n")

forecast_shifted <- typhoon_forecast_df %>%
  filter(scenario %in% c("0% (Baseline)", "7d -50% E7", "7d -50% E5", "7d -50% E3",
                         "7d -50% D3", "7d -50% D5", "7d -50% D7"))

forecast_shifted$scenario <- factor(forecast_shifted$scenario, 
                                    levels = c("0% (Baseline)", "7d -50% E7", "7d -50% E5", "7d -50% E3",
                                               "7d -50% D3", "7d -50% D5", "7d -50% D7"))

p6b_shifted <- ggplot(forecast_shifted, aes(x = date, y = predicted, color = scenario, fill = scenario, linetype = scenario)) +
  annotate("rect",
           xmin = typhoon_start_date_line,
           xmax = typhoon_end_date_line,
           ymin = -Inf, ymax = Inf,
           fill = "pink", alpha = 0.15) +
  annotate("rect",
           xmin = recovery_start_date,
           xmax = max(forecast_dates),
           ymin = -Inf, ymax = Inf,
           fill = "lightgreen", alpha = 0.1) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper), alpha = 0.18) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = typhoon_start_date_line,
             linetype = "dashed", color = "gray40", size = 0.8) +
  geom_vline(xintercept = typhoon_end_date_line,
             linetype = "dashed", color = "red", size = 0.8) +
  geom_point(data = validation_plot_df,
             aes(x = date, y = observed),
             color = "black", fill = "yellow", shape = 23, size = 2.0, stroke = 1.2, inherit.aes = FALSE) +
  facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
  scale_color_manual(values = colors_shifted) +
  scale_fill_manual(values = colors_shifted) +
  scale_linetype_manual(values = linetypes_shifted) +
  scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
  labs(
    title = "Typhoon Impact Forecasts: Shifted Timing Scenarios (7d -50%)",
    subtitle = "Yellow diamonds: Validation data | Continuous lines from last fitted point | Gray dashed: Typhoon start | Red dashed: Typhoon end\nPink: Typhoon period | Green: Recovery period",
    x = NULL,
    y = "Weekly Cases/Hospitalizations"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "gray30")
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )

print(p6b_shifted)
ggsave("figure6b_typhoon_forecasts_shifted.png", p6b_shifted, width = 14, height = 10, dpi = 300, bg = "white")

cat("\n=== CONTINUITY FIX SUMMARY ===\n")
cat("✓ Figure 1: Baseline forecast now continuous with historical fit\n")
cat("✓ Figure 2: All scenario forecasts now continuous with historical fit\n")
cat("✓ Figures 3-6: All typhoon scenario forecasts now continuous from last fitted point\n")
cat("✓ Lines and confidence intervals connect smoothly at typhoon transition point\n\n")


# Due to length, I'll create a streamlined version that generates all
# For complete code, these would be similar to Figure 2 but forecast-only

# [Figures 3-6 would be similar structure - omitting for length but would be included]

cat("Figures 3-6 placeholder (similar to Figure 2 structure)\n\n")

# ============================================================================
# FIGURE 7: Effective Reproduction Number (Baseline)
# ============================================================================
cat("\n=== Creating Figure 7: R_eff Baseline ===\n")

Reff_mean <- R_eff_scenarios_mean[1,,]
Reff_lower <- R_eff_scenarios_lower[1,,]
Reff_upper <- R_eff_scenarios_upper[1,,]

Reff_df <- data.frame(
  week = rep(1:T_weeks_total, N_strains),
  date = rep(all_dates, N_strains),
  strain = factor(rep(strain_names, each = T_weeks_total), levels = strain_names),
  R_eff = as.vector(Reff_mean),
  R_lower = as.vector(Reff_lower),
  R_upper = as.vector(Reff_upper),
  period = rep(c(rep("Historical", T_weeks), rep("Forecast", T_weeks_forecast)), N_strains)
)

strain_colors <- c("B" = "#E41A1C", "H3" = "#377EB8", "H1" = "#4DAF4A",
                   "COVID" = "#984EA3", "RSV" = "#FF7F00", "HFMD" = "#8B4513")

p7_reff_baseline <- ggplot(Reff_df, aes(x = date)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.5, size = 0.8) +
  annotate("rect",
           xmin = typhoon_start_date_line,
           xmax = typhoon_end_date_line,
           ymin = -Inf, ymax = Inf,
           fill = "pink", alpha = 0.15) +
  annotate("rect",
           xmin = recovery_start_date,
           xmax = max(forecast_dates),
           ymin = -Inf, ymax = Inf,
           fill = "lightgreen", alpha = 0.1) +
  geom_ribbon(data = filter(Reff_df, period == "Historical"),
              aes(ymin = R_lower, ymax = R_upper, fill = strain), alpha = 0.2) +
  geom_line(data = filter(Reff_df, period == "Historical"),
            aes(y = R_eff, color = strain), size = 1.2) +
  geom_ribbon(data = filter(Reff_df, period == "Forecast"),
              aes(ymin = R_lower, ymax = R_upper, fill = strain), alpha = 0.15) +
  geom_line(data = filter(Reff_df, period == "Forecast"),
            aes(y = R_eff, color = strain), size = 1, linetype = "dashed") +
  geom_vline(xintercept = typhoon_start_date_line,
             linetype = "dashed", color = "gray40", size = 0.8) +
  geom_vline(xintercept = typhoon_end_date_line,
             linetype = "dashed", color = "red", size = 0.8) +
  facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
  scale_color_manual(values = strain_colors) +
  scale_fill_manual(values = strain_colors) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "4 months") +
  scale_y_continuous(breaks = seq(0, 10, 1)) +
  coord_cartesian(ylim = c(0, 5)) +
  labs(
    title = expression(paste("Weekly Effective Reproduction Number (", R[eff], ") - Baseline Scenario")),
    subtitle = expression(paste("Solid: Historical | Dashed: Forecast | Red line: ", R[eff], "=1 threshold")),
    x = NULL,
    y = expression(R[eff])
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", size = 0.2),
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "gray30")
  )

print(p7_reff_baseline)
ggsave("figure7_reproduction_number_baseline.png", p7_reff_baseline,
       width = 14, height = 9, dpi = 300, bg = "white")

cat("Figure 7 complete\n\n")

# ============================================================================
# FIGURES 8-9: R_eff by Scenario and Changes (all durations)
# ============================================================================
# ============================================================================
# FIGURE 8: Effective Reproduction Number (Rt) by Scenario (Recent Period)
# ============================================================================
cat("\n=== Creating Figure 8: R_eff by Scenario (Recent Period) - Separate for each duration ===\n")

reff_start_week <- T_weeks - 3
reff_scenarios_list <- list()

for (typhoon_idx in 1:47) {
  for (strain_idx in 1:N_strains) {
    weeks_subset <- reff_start_week:T_weeks_total
    dates_subset <- all_dates[weeks_subset]
    
    reff_scenarios_list[[length(reff_scenarios_list) + 1]] <- data.frame(
      week = weeks_subset,
      date = dates_subset,
      strain = strain_names[strain_idx],
      R_eff = R_eff_scenarios_mean[typhoon_idx, weeks_subset, strain_idx],
      R_lower = R_eff_scenarios_lower[typhoon_idx, weeks_subset, strain_idx],
      R_upper = R_eff_scenarios_upper[typhoon_idx, weeks_subset, strain_idx],
      scenario = typhoon_scenarios[typhoon_idx]
    )
  }
}

reff_scenarios_df <- do.call(rbind, reff_scenarios_list)
reff_scenarios_df$strain <- factor(reff_scenarios_df$strain, levels = strain_names)
reff_scenarios_df$scenario <- factor(reff_scenarios_df$scenario, levels = typhoon_scenarios)

# Function to create reff plot for selected scenarios
create_reff_plot <- function(selected_scenarios, plot_title, colors, linetypes) {
  df_filtered <- reff_scenarios_df %>% filter(scenario %in% selected_scenarios)
  df_filtered$scenario <- factor(df_filtered$scenario, levels = selected_scenarios)
  
  p <- ggplot(df_filtered, aes(x = date)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.5, size = 0.8) +
    geom_vline(xintercept = max(fitting_data$date),
               linetype = "dashed", color = "gray40", size = 0.6, alpha = 0.6) +
    annotate("rect",
             xmin = typhoon_start_date_line,
             xmax = typhoon_end_date_line,
             ymin = -Inf, ymax = Inf,
             fill = "pink", alpha = 0.15) +
    annotate("rect",
             xmin = recovery_start_date,
             xmax = max(forecast_dates),
             ymin = -Inf, ymax = Inf,
             fill = "lightgreen", alpha = 0.1) +
    geom_vline(xintercept = typhoon_end_date_line,
               linetype = "dashed", color = "red", size = 0.7, alpha = 0.7) +
    geom_ribbon(aes(ymin = R_lower, ymax = R_upper, fill = scenario), alpha = 0.15) +
    geom_line(aes(y = R_eff, color = scenario, linetype = scenario), size = 1.1) +
    facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
    scale_color_manual(name = "Scenarios:", values = colors) +
    scale_fill_manual(name = "Scenarios:", values = colors) +
    scale_linetype_manual(name = "Scenarios:", values = linetypes) +
    scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
    scale_y_continuous(breaks = seq(0, 10, 1)) +
    coord_cartesian(ylim = c(0, 5)) +
    labs(
      title = plot_title,
      subtitle = "Comparing R_eff trajectories under different transmission modification scenarios\nPink: Typhoon | Green: Recovery",
      x = NULL,
      y = "R_eff"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", size = 0.2),
      panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      panel.spacing = unit(0.5, "lines"),
      strip.text = element_text(size = 11, face = "bold", margin = margin(b = 5)),
      strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
      axis.text = element_text(size = 9, color = "gray20"),
      axis.title.y = element_text(size = 10, face = "bold", margin = margin(r = 8)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.background = element_rect(fill = "white", color = "gray40", size = 0.3),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8.5),
      legend.margin = margin(t = 5, b = 5),
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 8.5, color = "gray30")
    ) +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE),
      fill = guide_legend(nrow = 2, byrow = TRUE),
      linetype = guide_legend(nrow = 2, byrow = TRUE)
    )
  
  return(p)
}

# Create R_eff plots for each duration
selected_3d_reff <- c("0% (Baseline)", "3d -50%", "3d -30%", "3d -20%", "3d -10%", "3d -5%",
                      "3d +5%", "3d +10%", "3d +20%", "3d +30%", "3d +50%")
p8_3days <- create_reff_plot(selected_3d_reff, 
                             "Effective Reproduction Number (R_eff) by Typhoon Scenario - 3 Days Duration", 
                             colors_3d, linetypes_3d)
print(p8_3days)
ggsave("figure8_reff_scenarios_3days.png", p8_3days, width = 14, height = 10, dpi = 300, bg = "white")

selected_7d_reff <- c("0% (Baseline)", "7d -50%", "7d -30%", "7d -20%", "7d -10%", "7d -5%",
                      "7d +5%", "7d +10%", "7d +20%", "7d +30%", "7d +50%")
p8_7days <- create_reff_plot(selected_7d_reff, 
                             "Effective Reproduction Number (R_eff) by Typhoon Scenario - 7 Days Duration", 
                             colors_7d, linetypes_7d)
print(p8_7days)
ggsave("figure8_reff_scenarios_7days.png", p8_7days, width = 14, height = 10, dpi = 300, bg = "white")

selected_10d_reff <- c("0% (Baseline)", "10d -50%", "10d -30%", "10d -20%", "10d -10%", "10d -5%",
                       "10d +5%", "10d +10%", "10d +20%", "10d +30%", "10d +50%")
p8_10days <- create_reff_plot(selected_10d_reff, 
                              "Effective Reproduction Number (R_eff) by Typhoon Scenario - 10 Days Duration", 
                              colors_10d, linetypes_10d)
print(p8_10days)
ggsave("figure8_reff_scenarios_10days.png", p8_10days, width = 14, height = 10, dpi = 300, bg = "white")

selected_14d_reff <- c("0% (Baseline)", "14d -50%", "14d -30%", "14d -20%", "14d -10%", "14d -5%",
                       "14d +5%", "14d +10%", "14d +20%", "14d +30%", "14d +50%")
p8_14days <- create_reff_plot(selected_14d_reff, 
                              "Effective Reproduction Number (R_eff) by Typhoon Scenario - 14 Days Duration", 
                              colors_14d, linetypes_14d)
print(p8_14days)
ggsave("figure8_reff_scenarios_14days.png", p8_14days, width = 14, height = 10, dpi = 300, bg = "white")

selected_shifted_reff <- c("0% (Baseline)", "7d -50% E7", "7d -50% E5", "7d -50% E3",
                           "7d -50% D3", "7d -50% D5", "7d -50% D7")
p8_shifted <- create_reff_plot(selected_shifted_reff, 
                               "Effective Reproduction Number (R_eff) by Typhoon Scenario - Shifted Timing Scenarios", 
                               colors_shifted, linetypes_shifted)
print(p8_shifted)
ggsave("figure8_reff_scenarios_shifted.png", p8_shifted, width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 9: Change in R_eff from Baseline - Separate for each duration
# ============================================================================
cat("\n=== Creating Figure 9: Change in R_eff from Baseline ===\n")

reff_change_start_week <- T_weeks - 4
reff_change_list <- list()

for (typhoon_idx in 2:47) {
  for (strain_idx in 1:N_strains) {
    weeks_subset <- reff_change_start_week:T_weeks_total
    dates_subset <- all_dates[weeks_subset]
    
    # Get baseline (scenario 1) posterior samples
    base_reff_samples <- R_eff_scenarios[, 1, weeks_subset, strain_idx]
    
    # Get current scenario posterior samples
    scen_reff_samples <- R_eff_scenarios[, typhoon_idx, weeks_subset, strain_idx]
    
    # Calculate difference distribution
    diff_dist <- scen_reff_samples - base_reff_samples
    
    # Calculate statistics across samples
    R_eff_change_mean <- apply(diff_dist, 2, mean)
    R_change_lower <- apply(diff_dist, 2, quantile, probs = 0.025)
    R_change_upper <- apply(diff_dist, 2, quantile, probs = 0.975)
    
    reff_change_list[[length(reff_change_list) + 1]] <- data.frame(
      week = weeks_subset,
      date = dates_subset,
      strain = strain_names[strain_idx],
      R_eff_change = R_eff_change_mean,
      R_change_lower = R_change_lower,
      R_change_upper = R_change_upper,
      scenario = typhoon_scenarios[typhoon_idx]
    )
  }
}

reff_change_df <- do.call(rbind, reff_change_list)
reff_change_df$strain <- factor(reff_change_df$strain, levels = strain_names)
reff_change_df$scenario <- factor(reff_change_df$scenario, levels = typhoon_scenarios[2:47])

# Function to create reff change plot for selected scenarios
create_reff_change_plot <- function(selected_scenarios, plot_title, colors, linetypes) {
  df_filtered <- reff_change_df %>% filter(scenario %in% selected_scenarios)
  df_filtered$scenario <- factor(df_filtered$scenario, levels = selected_scenarios)
  
  p <- ggplot(df_filtered, aes(x = date)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5, size = 0.8) +
    geom_vline(xintercept = max(fitting_data$date),
               linetype = "dashed", color = "gray40", size = 0.6, alpha = 0.6) +
    annotate("rect",
             xmin = typhoon_start_date_line,
             xmax = typhoon_end_date_line,
             ymin = -Inf, ymax = Inf,
             fill = "pink", alpha = 0.15) +
    annotate("rect",
             xmin = recovery_start_date,
             xmax = max(forecast_dates),
             ymin = -Inf, ymax = Inf,
             fill = "lightgreen", alpha = 0.1) +
    geom_vline(xintercept = typhoon_end_date_line,
               linetype = "dashed", color = "red", size = 0.7, alpha = 0.7) +
    geom_ribbon(aes(ymin = R_change_lower, ymax = R_change_upper, fill = scenario), alpha = 0.15) +
    geom_line(aes(y = R_eff_change, color = scenario, linetype = scenario), size = 1.1) +
    facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
    scale_color_manual(name = "Scenarios:", values = colors) +
    scale_fill_manual(name = "Scenarios:", values = colors) +
    scale_linetype_manual(name = "Scenarios:", values = linetypes) +
    scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
    labs(
      title = plot_title,
      subtitle = "Negative values indicate reduction in transmission | Positive values indicate increase | ΔR_eff = R_eff(scenario) - R_eff(baseline)",
      x = NULL,
      y = "ΔR_eff"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", size = 0.2),
      panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      panel.spacing = unit(0.5, "lines"),
      strip.text = element_text(size = 11, face = "bold", margin = margin(b = 5)),
      strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
      axis.text = element_text(size = 9, color = "gray20"),
      axis.title.y = element_text(size = 10, face = "bold", margin = margin(r = 8)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.background = element_rect(fill = "white", color = "gray40", size = 0.3),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8.5),
      legend.margin = margin(t = 5, b = 5),
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 8.5, color = "gray30")
    ) +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE),
      fill = guide_legend(nrow = 2, byrow = TRUE),
      linetype = guide_legend(nrow = 2, byrow = TRUE)
    )
  
  return(p)
}

# Define colors for change plots (excluding baseline)
colors_3d_change <- colors_3d[names(colors_3d) != "0% (Baseline)"]
linetypes_3d_change <- linetypes_3d[names(linetypes_3d) != "0% (Baseline)"]
selected_3d_change <- selected_3d[selected_3d != "0% (Baseline)" & selected_3d != "Historical Fit"]

p9_3days <- create_reff_change_plot(selected_3d_change, 
                                    "Change in R_eff from Baseline - 3 Days Duration", 
                                    colors_3d_change, linetypes_3d_change)
print(p9_3days)
ggsave("figure9_reff_change_3days.png", p9_3days, width = 14, height = 10, dpi = 300, bg = "white")

colors_7d_change <- colors_7d[names(colors_7d) != "0% (Baseline)"]
linetypes_7d_change <- linetypes_7d[names(linetypes_7d) != "0% (Baseline)"]
selected_7d_change <- selected_7d[selected_7d != "0% (Baseline)" & selected_7d != "Historical Fit"]

p9_7days <- create_reff_change_plot(selected_7d_change, 
                                    "Change in R_eff from Baseline - 7 Days Duration", 
                                    colors_7d_change, linetypes_7d_change)
print(p9_7days)
ggsave("figure9_reff_change_7days.png", p9_7days, width = 14, height = 10, dpi = 300, bg = "white")

colors_10d_change <- colors_10d[names(colors_10d) != "0% (Baseline)"]
linetypes_10d_change <- linetypes_10d[names(linetypes_10d) != "0% (Baseline)"]
selected_10d_change <- selected_10d[selected_10d != "0% (Baseline)" & selected_10d != "Historical Fit"]

p9_10days <- create_reff_change_plot(selected_10d_change, 
                                     "Change in R_eff from Baseline - 10 Days Duration", 
                                     colors_10d_change, linetypes_10d_change)
print(p9_10days)
ggsave("figure9_reff_change_10days.png", p9_10days, width = 14, height = 10, dpi = 300, bg = "white")

colors_14d_change <- colors_14d[names(colors_14d) != "0% (Baseline)"]
linetypes_14d_change <- linetypes_14d[names(linetypes_14d) != "0% (Baseline)"]
selected_14d_change <- selected_14d[selected_14d != "0% (Baseline)" & selected_14d != "Historical Fit"]

p9_14days <- create_reff_change_plot(selected_14d_change, 
                                     "Change in R_eff from Baseline - 14 Days Duration", 
                                     colors_14d_change, linetypes_14d_change)
print(p9_14days)
ggsave("figure9_reff_change_14days.png", p9_14days, width = 14, height = 10, dpi = 300, bg = "white")

colors_shifted_change <- colors_shifted[names(colors_shifted) != "0% (Baseline)"]
linetypes_shifted_change <- linetypes_shifted[names(linetypes_shifted) != "0% (Baseline)"]
selected_shifted_change <- selected_shifted[selected_shifted != "0% (Baseline)" & selected_shifted != "Historical Fit"]

p9_shifted <- create_reff_change_plot(selected_shifted_change, 
                                      "Change in R_eff from Baseline - Shifted Timing Scenarios", 
                                      colors_shifted_change, linetypes_shifted_change)
print(p9_shifted)

# ============================================================================
# FIGURE 10: Typhoon Impact Effectiveness (Reduction During Typhoon Period)
# ============================================================================
cat("\n=== Creating Figure 10: Typhoon Impact Effectiveness ===\n")

reduction_typhoon_mean <- apply(reduction_typhoon, c(2,3), mean)
reduction_typhoon_lower <- apply(reduction_typhoon, c(2,3), quantile, probs = 0.025)
reduction_typhoon_upper <- apply(reduction_typhoon, c(2,3), quantile, probs = 0.975)

create_reduction_plot <- function(scenario_indices, plot_title, colors) {
  selected_scenarios <- typhoon_scenarios[scenario_indices]
  
  reduction_df <- data.frame(
    scenario = rep(selected_scenarios, N_strains),
    strain = rep(strain_names, each = length(selected_scenarios)),
    mean = as.vector(reduction_typhoon_mean[scenario_indices,]),
    lower = as.vector(reduction_typhoon_lower[scenario_indices,]),
    upper = as.vector(reduction_typhoon_upper[scenario_indices,])
  )
  
  reduction_df$scenario <- factor(reduction_df$scenario, levels = selected_scenarios)
  reduction_df$strain <- factor(reduction_df$strain, levels = strain_names)
  
  p <- ggplot(reduction_df, aes(x = scenario, y = mean, fill = scenario)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = position_dodge(0.8), width = 0.25, size = 0.6) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
    facet_wrap(~ strain, scales = "free_y") +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title = plot_title,
      subtitle = "Error bars: 95% credible intervals | Typhoon Period ONLY\nNegative = increase | Positive = reduction",
      x = "Typhoon Scenario",
      y = "Change During Typhoon Period (%)",
      fill = "Scenario"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom",
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 8.5, color = "gray30")
    )
  
  return(p)
}

selected_3d_red_indices <- c(1, 2:11)
p10_3days <- create_reduction_plot(selected_3d_red_indices,
                                   "Percentage Change in Cases (Typhoon Period) - 3 Days Duration",
                                   colors_3d)
print(p10_3days)
ggsave("figure10_typhoon_reduction_3days.png", p10_3days,
       width = 12, height = 12, dpi = 300, bg = "white")

selected_7d_red_indices <- c(1, 12:21)
p10_7days <- create_reduction_plot(selected_7d_red_indices,
                                   "Percentage Change in Cases (Typhoon Period) - 7 Days Duration",
                                   colors_7d)
print(p10_7days)
ggsave("figure10_typhoon_reduction_7days.png", p10_7days,
       width = 12, height = 12, dpi = 300, bg = "white")

selected_10d_red_indices <- c(1, 22:31)
p10_10days <- create_reduction_plot(selected_10d_red_indices,
                                    "Percentage Change in Cases (Typhoon Period) - 10 Days Duration",
                                    colors_10d)
print(p10_10days)
ggsave("figure10_typhoon_reduction_10days.png", p10_10days,
       width = 12, height = 12, dpi = 300, bg = "white")

selected_14d_red_indices <- c(1, 32:41)
p10_14days <- create_reduction_plot(selected_14d_red_indices,
                                    "Percentage Change in Cases (Typhoon Period) - 14 Days Duration",
                                    colors_14d)
print(p10_14days)
ggsave("figure10_typhoon_reduction_14days.png", p10_14days,
       width = 12, height = 12, dpi = 300, bg = "white")

selected_shifted_red_indices <- c(1, 42:47)
p10_shifted <- create_reduction_plot(selected_shifted_red_indices,
                                     "Percentage Change in Cases (Typhoon Period) - Shifted Timing",
                                     colors_shifted)
print(p10_shifted)
ggsave("figure10_typhoon_reduction_shifted.png", p10_shifted,
       width = 12, height = 12, dpi = 300, bg = "white")

cat("Figure 10 complete\n\n")

# ============================================================================
# FIGURE 11: Average Weekly Cases - Typhoon vs Recovery
# ============================================================================
cat("\n=== Creating Figure 11: Average Weekly Cases by Period ===\n")

avg_typhoon_mean <- apply(avg_weekly_typhoon, c(2,3), mean)
avg_typhoon_lower <- apply(avg_weekly_typhoon, c(2,3), quantile, probs = 0.025)
avg_typhoon_upper <- apply(avg_weekly_typhoon, c(2,3), quantile, probs = 0.975)

avg_recovery_mean <- apply(avg_weekly_recovery, c(2,3), mean)
avg_recovery_lower <- apply(avg_weekly_recovery, c(2,3), quantile, probs = 0.025)
avg_recovery_upper <- apply(avg_weekly_recovery, c(2,3), quantile, probs = 0.975)

create_avg_cases_plot <- function(selected_scenario_indices, plot_title) {
  selected_scenario_names <- typhoon_scenarios[selected_scenario_indices]
  
  avg_df_typhoon <- data.frame(
    scenario = rep(selected_scenario_names, N_strains),
    strain = rep(strain_names, each = length(selected_scenario_names)),
    mean = as.vector(avg_typhoon_mean[selected_scenario_indices,]),
    lower = as.vector(avg_typhoon_lower[selected_scenario_indices,]),
    upper = as.vector(avg_typhoon_upper[selected_scenario_indices,]),
    period = "Typhoon Period (Avg Weekly)"
  )
  
  avg_df_recovery <- data.frame(
    scenario = rep(selected_scenario_names, N_strains),
    strain = rep(strain_names, each = length(selected_scenario_names)),
    mean = as.vector(avg_recovery_mean[selected_scenario_indices,]),
    lower = as.vector(avg_recovery_lower[selected_scenario_indices,]),
    upper = as.vector(avg_recovery_upper[selected_scenario_indices,]),
    period = "Recovery Period (Avg Weekly)"
  )
  
  avg_df <- rbind(avg_df_typhoon, avg_df_recovery)
  avg_df$scenario <- factor(avg_df$scenario, levels = selected_scenario_names)
  avg_df$strain <- factor(avg_df$strain, levels = strain_names)
  avg_df$period <- factor(avg_df$period, levels = c("Typhoon Period (Avg Weekly)", "Recovery Period (Avg Weekly)"))
  
  p <- ggplot(avg_df, aes(x = scenario, y = mean, fill = period)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.8), width = 0.25, size = 0.6) +
    facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
    scale_fill_manual(values = c("Typhoon Period (Avg Weekly)" = "#E41A1C", "Recovery Period (Avg Weekly)" = "#4DAF4A")) +
    labs(
      title = plot_title,
      subtitle = "Error bars: 95% credible intervals",
      x = "Typhoon Scenario",
      y = "Average Weekly Cases",
      fill = "Period"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom",
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 8.5, color = "gray30")
    )
  
  return(p)
}

p11_3days <- create_avg_cases_plot(selected_3d_red_indices, 
                                   "Average Weekly Cases: Typhoon vs Recovery - 3 Days")
print(p11_3days)
ggsave("figure11_avg_weekly_cases_3days.png", p11_3days, width = 14, height = 10, dpi = 300, bg = "white")

p11_7days <- create_avg_cases_plot(selected_7d_red_indices, 
                                   "Average Weekly Cases: Typhoon vs Recovery - 7 Days")
print(p11_7days)
ggsave("figure11_avg_weekly_cases_7days.png", p11_7days, width = 14, height = 10, dpi = 300, bg = "white")

cat("Figures 10-11 complete\n\n")

# ============================================================================
# FIGURE 12: Validation Comparison
# ============================================================================
cat("\n=== Creating Figure 12: Validation Comparison ===\n")

# 解决 'validation_comparison' not found 错误
# 将观测数据 (validation_plot_df) 与基线预测数据 (forecast_baseline_raw) 合并
validation_comparison <- validation_plot_df %>%
  left_join(
    forecast_baseline_raw %>% 
      # 选择并重命名预测所需的列
      select(date, strain, baseline_pred = pred_mean, baseline_lower = pred_lower, baseline_upper = pred_upper),
    by = c("date", "strain")
  ) %>%
  # 仅保留观测值和预测值均存在的行
  filter(!is.na(observed) & !is.na(baseline_pred))

if (nrow(validation_comparison) > 0) {
  strain_colors <- c("B" = "#E41A1C", "H3" = "#377EB8", "H1" = "#4DAF4A",
                     "COVID" = "#984EA3", "RSV" = "#FF7F00", "HFMD" = "#FFFF33")
  
  p12_validation <- ggplot(validation_comparison, aes(x = baseline_pred, y = observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1) +
    geom_errorbarh(aes(xmin = baseline_lower, xmax = baseline_upper),
                   height = 0, alpha = 0.3, color = "gray50") +
    geom_point(aes(color = strain), size = 3, alpha = 0.7) +
    facet_wrap(~strain, scales = "free", ncol = 2, nrow = 3) +
    scale_color_manual(values = strain_colors) +
    labs(
      title = "Out-of-Sample Validation: Observed vs Baseline Predicted Cases",
      subtitle = "Red line = perfect predictions | Error bars = 95% credible intervals | Data NOT used in fitting",
      x = "Baseline Predicted Cases",
      y = "Observed Cases (Validation Data)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_line(color = "gray90", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = "gray30"),
      legend.position = "none",
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray30")
    )
  
  print(p12_validation)
  ggsave("figure12_validation_comparison.png", p12_validation,
         width = 12, height = 10, dpi = 300, bg = "white")
}

# ============================================================================
# NEW: FIGURE 13 - ATTACK RATE ANALYSIS
# ============================================================================
cat("\n=== Creating Figure 13: Attack Rate Analysis (Up to T_weeks - 2) ===\n")

# Prepare attack rate data
attack_rate_df <- data.frame()

for (strain_idx in 1:N_strains) {
  # Baseline
  attack_rate_df <- rbind(attack_rate_df, data.frame(
    scenario = "0% (Baseline)",
    strain = strain_names[strain_idx],
    attack_rate = attack_rate_baseline_mean[strain_idx],
    lower = attack_rate_baseline_lower[strain_idx],
    upper = attack_rate_baseline_upper[strain_idx],
    stringsAsFactors = FALSE
  ))
  
  # All 47 scenarios
  for (typhoon_idx in 1:47) {
    attack_rate_df <- rbind(attack_rate_df, data.frame(
      scenario = typhoon_scenarios[typhoon_idx],
      strain = strain_names[strain_idx],
      attack_rate = attack_rate_scenarios_mean[typhoon_idx, strain_idx],
      lower = attack_rate_scenarios_lower[typhoon_idx, strain_idx],
      upper = attack_rate_scenarios_upper[typhoon_idx, strain_idx],
      stringsAsFactors = FALSE
    ))
  }
}

attack_rate_df$scenario <- factor(attack_rate_df$scenario, levels = typhoon_scenarios)
attack_rate_df$strain <- factor(attack_rate_df$strain, levels = strain_names)

# Difference from baseline
attack_rate_diff_df <- data.frame()

for (strain_idx in 1:N_strains) {
  for (typhoon_idx in 1:47) {
    attack_rate_diff_df <- rbind(attack_rate_diff_df, data.frame(
      scenario = typhoon_scenarios[typhoon_idx],
      strain = strain_names[strain_idx],
      diff = attack_rate_diff_mean[typhoon_idx, strain_idx],
      diff_lower = attack_rate_diff_lower[typhoon_idx, strain_idx],
      diff_upper = attack_rate_diff_upper[typhoon_idx, strain_idx],
      stringsAsFactors = FALSE
    ))
  }
}

attack_rate_diff_df$scenario <- factor(attack_rate_diff_df$scenario, levels = typhoon_scenarios)
attack_rate_diff_df$strain <- factor(attack_rate_diff_df$strain, levels = strain_names)

# Function to create attack rate comparison plot
create_attack_rate_plot <- function(selected_scenarios, plot_title, colors) {
  df_filtered <- attack_rate_df %>% filter(scenario %in% selected_scenarios)
  df_filtered$scenario <- factor(df_filtered$scenario, levels = selected_scenarios)
  
  p <- ggplot(df_filtered, aes(x = scenario, y = attack_rate, fill = scenario)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = position_dodge(0.8), width = 0.25, size = 0.6) +
    facet_wrap(~ strain, scales = "free_y") +
    scale_fill_manual(values = colors) +
    labs(
      title = plot_title,
      subtitle = "Error bars: 95% credible intervals | Cumulative incidence up to T_weeks - 2 (before typhoon)",
      x = "Scenario",
      y = "Cumulative Attack Rate (Proportion)",
      fill = "Scenario"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom",
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 8.5, color = "gray30")
    )
  
  return(p)
}

# Function to create attack rate difference plot
create_attack_rate_diff_plot <- function(selected_scenarios, plot_title, colors) {
  df_filtered <- attack_rate_diff_df %>% filter(scenario %in% selected_scenarios)
  df_filtered$scenario <- factor(df_filtered$scenario, levels = selected_scenarios)
  
  p <- ggplot(df_filtered, aes(x = scenario, y = diff, fill = scenario)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = diff_lower, ymax = diff_upper),
                  position = position_dodge(0.8), width = 0.25, size = 0.6) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
    facet_wrap(~ strain, scales = "free_y") +
    scale_fill_manual(values = colors) +
    labs(
      title = plot_title,
      subtitle = "Error bars: 95% credible intervals | Difference = Scenario - Baseline\nPositive = higher attack rate | Negative = lower attack rate",
      x = "Scenario",
      y = "Difference in Attack Rate from Baseline",
      fill = "Scenario"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom",
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 8.5, color = "gray30")
    )
  
  return(p)
}

# Create attack rate plots for each duration
p13a_3days <- create_attack_rate_plot(selected_3d, "Attack Rate Comparison - 3 Days Duration", colors_3d)
print(p13a_3days)
ggsave("figure13a_attack_rate_3days.png", p13a_3days, width = 12, height = 12, dpi = 300, bg = "white")

p13b_3days <- create_attack_rate_diff_plot(selected_3d[selected_3d != "0% (Baseline)"], 
                                           "Attack Rate Difference from Baseline - 3 Days Duration", 
                                           colors_3d[names(colors_3d) != "0% (Baseline)"])
print(p13b_3days)
ggsave("figure13b_attack_rate_diff_3days.png", p13b_3days, width = 12, height = 12, dpi = 300, bg = "white")

p13a_7days <- create_attack_rate_plot(selected_7d, "Attack Rate Comparison - 7 Days Duration", colors_7d)
print(p13a_7days)
ggsave("figure13a_attack_rate_7days.png", p13a_7days, width = 12, height = 12, dpi = 300, bg = "white")

p13b_7days <- create_attack_rate_diff_plot(selected_7d[selected_7d != "0% (Baseline)"], 
                                           "Attack Rate Difference from Baseline - 7 Days Duration", 
                                           colors_7d[names(colors_7d) != "0% (Baseline)"])
print(p13b_7days)
ggsave("figure13b_attack_rate_diff_7days.png", p13b_7days, width = 12, height = 12, dpi = 300, bg = "white")

p13a_10days <- create_attack_rate_plot(selected_10d, "Attack Rate Comparison - 10 Days Duration", colors_10d)
print(p13a_10days)
ggsave("figure13a_attack_rate_10days.png", p13a_10days, width = 12, height = 12, dpi = 300, bg = "white")

p13b_10days <- create_attack_rate_diff_plot(selected_10d[selected_10d != "0% (Baseline)"], 
                                            "Attack Rate Difference from Baseline - 10 Days Duration", 
                                            colors_10d[names(colors_10d) != "0% (Baseline)"])
print(p13b_10days)
ggsave("figure13b_attack_rate_diff_10days.png", p13b_10days, width = 12, height = 12, dpi = 300, bg = "white")

p13a_14days <- create_attack_rate_plot(selected_14d, "Attack Rate Comparison - 14 Days Duration", colors_14d)
print(p13a_14days)
ggsave("figure13a_attack_rate_14days.png", p13a_14days, width = 12, height = 12, dpi = 300, bg = "white")

p13b_14days <- create_attack_rate_diff_plot(selected_14d[selected_14d != "0% (Baseline)"], 
                                            "Attack Rate Difference from Baseline - 14 Days Duration", 
                                            colors_14d[names(colors_14d) != "0% (Baseline)"])
print(p13b_14days)
ggsave("figure13b_attack_rate_diff_14days.png", p13b_14days, width = 12, height = 12, dpi = 300, bg = "white")

p13a_shifted <- create_attack_rate_plot(selected_shifted, "Attack Rate Comparison - Shifted Timing", colors_shifted)
print(p13a_shifted)
ggsave("figure13a_attack_rate_shifted.png", p13a_shifted, width = 12, height = 12, dpi = 300, bg = "white")

p13b_shifted <- create_attack_rate_diff_plot(selected_shifted[selected_shifted != "0% (Baseline)"], 
                                             "Attack Rate Difference from Baseline - Shifted Timing", 
                                             colors_shifted[names(colors_shifted) != "0% (Baseline)"])
print(p13b_shifted)
ggsave("figure13b_attack_rate_diff_shifted.png", p13b_shifted, width = 12, height = 12, dpi = 300, bg = "white")

# ============================================================================
# Create attack rate summary table
# ============================================================================
cat("\n=== Creating Attack Rate Summary Table ===\n")

attack_rate_table <- data.frame()

for (typhoon_idx in 1:47) {
  for (strain_idx in 1:N_strains) {
    attack_rate_table <- rbind(attack_rate_table, data.frame(
      Scenario = typhoon_scenarios[typhoon_idx],
      Strain = strain_names[strain_idx],
      Baseline_Attack_Rate = attack_rate_baseline_mean[strain_idx],
      Scenario_Attack_Rate = attack_rate_scenarios_mean[typhoon_idx, strain_idx],
      Difference = attack_rate_diff_mean[typhoon_idx, strain_idx],
      Diff_Lower_CI = attack_rate_diff_lower[typhoon_idx, strain_idx],
      Diff_Upper_CI = attack_rate_diff_upper[typhoon_idx, strain_idx],
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(attack_rate_table,
          file = "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/attack_rate_summary.csv",
          row.names = FALSE)

cat("Attack rate summary table saved to: attack_rate_summary.csv\n")

# Print summary statistics
cat("\n=== Attack Rate Summary Statistics ===\n")
cat("Baseline attack rates (up to T_weeks - 2):\n")
for (i in 1:N_strains) {
  cat(sprintf("%s: %.6f [%.6f, %.6f]\n", 
              strain_names[i],
              attack_rate_baseline_mean[i],
              attack_rate_baseline_lower[i],
              attack_rate_baseline_upper[i]))
}

# ============================================================================
# Extract and Print Key Parameters
# ============================================================================
cat("\n=== KEY PARAMETERS SUMMARY ===\n")

params <- list(
  R0 = apply(rstan::extract(fit, "log_R0_spline_coeff")$log_R0_spline_coeff,
             c(1,2), function(x) exp(mean(x))),
  sigma = rstan::extract(fit, "sigma")$sigma,
  gamma = rstan::extract(fit, "gamma")$gamma,
  mu = rstan::extract(fit, "mu")$mu,
  detection = rstan::extract(fit, "detection_rate")$detection_rate,
  hosp = rstan::extract(fit, "hospitalization_rate")$hospitalization_rate,
  phi = rstan::extract(fit, "phi")$phi,
  cross_flu = rstan::extract(fit, "cross_immunity_flu")$cross_immunity_flu,
  cross_flu_covid = rstan::extract(fit, "cross_immunity_flu_covid")$cross_immunity_flu_covid,
  cross_flu_rsv = rstan::extract(fit, "cross_immunity_flu_rsv")$cross_immunity_flu_rsv,
  cross_covid_rsv = rstan::extract(fit, "cross_immunity_covid_rsv")$cross_immunity_covid_rsv,
  init_state = rstan::extract(fit, "init_state")$init_state
)

print_param <- function(x, name, strains = NULL) {
  if(is.null(dim(x)) || length(dim(x)) == 1) {
    stats <- c(Mean = mean(x), SD = sd(x),
               quantile(x, c(0.025, 0.5, 0.975)))
    cat("\n", name, ":\n", sep="")
    print(round(stats, 4))
  } else {
    stats <- t(apply(x, 2, function(v) c(Mean = mean(v), SD = sd(v),
                                         quantile(v, c(0.025, 0.5, 0.975)))))
    if(!is.null(strains)) rownames(stats) = strains
    cat("\n", name, ":\n", sep="")
    print(round(stats, 4))
  }
}

print_param(params$R0, "R0 (Basic Reproduction Number)", strain_names)
print_param(params$sigma, "Sigma (incubation rate, 1/days)", strain_names)
print_param(params$gamma, "Gamma (recovery rate, 1/days)", strain_names)
print_param(params$mu, "Mu (immunity waning, 1/days)", strain_names)

cat("\nDerived periods (days):\n")
print(data.frame(
  Strain = strain_names,
  Incubation = round(1/colMeans(params$sigma), 1),
  Infectious = round(1/colMeans(params$gamma), 1),
  Immunity = round(1/colMeans(params$mu), 0)
))

print_param(params$cross_flu, "Cross-immunity: Flu-Flu")
print_param(params$cross_flu_covid, "Cross-immunity: Flu-COVID")
print_param(params$cross_flu_rsv, "Cross-immunity: Flu-RSV")
print_param(params$cross_covid_rsv, "Cross-immunity: COVID-RSV")
print_param(params$detection, "Detection rates", strain_names[1:5])
print_param(params$hosp, "Hospitalization rate (HFMD)")
print_param(params$phi, "Phi (overdispersion)", strain_names)

# ============================================================================
# FIGURE 14: Visualize Posterior Distributions
# ============================================================================
cat("\n=== Creating Figure 14: Posterior Distribution Figure ===\n")

make_df <- function(mat, name, labs) {
  if(is.null(dim(mat))) mat <- matrix(mat, ncol=1)
  data.frame(
    value = as.vector(mat),
    param = name,
    group = rep(labs, each=nrow(mat))
  )
}

post_data <- rbind(
  make_df(params$R0, "R0", strain_names),
  make_df(params$sigma, "Sigma", strain_names),
  make_df(params$gamma, "Gamma", strain_names),
  make_df(params$mu, "Mu", strain_names),
  make_df(params$detection, "Detection", strain_names[1:5]),
  make_df(params$hosp, "Hosp.Rate", "HFMD"),
  make_df(params$cross_flu, "Cross:Flu", "Flu"),
  make_df(params$cross_flu_covid, "Cross:F-C", "F-C"),
  make_df(params$cross_flu_rsv, "Cross:F-R", "F-R"),
  make_df(params$cross_covid_rsv, "Cross:C-R", "C-R"),
  make_df(params$phi, "Phi", strain_names)
)

post_data$param <- factor(post_data$param,
                          levels = c("R0", "Sigma", "Gamma", "Mu", "Detection", "Hosp.Rate",
                                     "Cross:Flu", "Cross:F-C", "Cross:F-R", "Cross:C-R", "Phi"))

p14_post <- ggplot(post_data, aes(x = value, fill = group)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~param, scales = "free", ncol = 3) +
  scale_fill_manual(values = c(B="#E41A1C", H3="#377EB8", H1="#4DAF4A",
                               COVID="#984EA3", RSV="#FF7F00", HFMD="#8B4513",
                               Flu="#E41A1C", `F-C`="#984EA3", `F-R`="#FF7F00",
                               `C-R`="#66C2A5")) +
  labs(title = "Posterior Distributions of Model Parameters",
       x = "Value", y = "Density", fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "gray95", color = "gray30"),
    panel.border = element_rect(fill = NA, color = "gray30"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12)
  )

print(p14_post)
ggsave("figure14_posterior_distributions.png", p14_post,
       width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================================
# Print Initial States
# ============================================================================
cat("\n", rep("=", 80), "\n", sep="")
cat("INITIAL COMPARTMENT STATES (SEIR)\n")
cat(rep("=", 80), "\n", sep="")

init_summary <- t(apply(params$init_state, 2, function(x) {
  c(Mean = mean(x), SD = sd(x), quantile(x, c(0.025, 0.5, 0.975)))
}))

E_indices <- seq(2, 19, by=3)
I_indices <- seq(3, 19, by=3)
R_indices <- seq(4, 19, by=3)

cat("\nS (Susceptible):\n")
print(round(init_summary[1,], 6))

cat("\nE (Exposed) by strain:\n")
E_summary <- init_summary[E_indices,]
rownames(E_summary) <- strain_names
print(round(E_summary, 6))

cat("\nI (Infected) by strain:\n")
I_summary <- init_summary[I_indices,]
rownames(I_summary) <- strain_names
print(round(I_summary, 6))

cat("\nR (Recovered) by strain:\n")
R_summary <- init_summary[R_indices,]
rownames(R_summary) <- strain_names
print(round(R_summary, 6))

cat("\nVerification - Sum of all compartments = ", round(sum(init_summary[,1]), 6), "\n")

# ============================================================================
# Print summary statistics to verify correctness
# ============================================================================
cat("\n=== VERIFICATION: Reduction Statistics (TYPHOON PERIOD ONLY) ===\n")
cat("\n=== Typhoon Period Changes (1 week only) ===\n")

for (typhoon_idx in 2:47) {
  cat(sprintf("\n%s:\n", typhoon_scenarios[typhoon_idx]))
  for (strain_idx in 1:N_strains) {
    red_dist <- reduction_typhoon[, typhoon_idx, strain_idx]
    cat(sprintf("  %s: %.1f%% [%.1f%%, %.1f%%]\n",
                strain_names[strain_idx],
                mean(red_dist),
                quantile(red_dist, 0.025),
                quantile(red_dist, 0.975)))
  }
}

# ============================================================================
# Calculate mean values for extended scenarios
# ============================================================================
cat("\n=== Calculating extended scenario statistics ===\n")

# Calculate mean across posterior samples
incidence_per10k_mean <- apply(incidence_per10k_extended, c(2,3), mean)
reduction_extended_mean <- apply(reduction_extended, c(2,3), mean)

cat("incidence_per10k_mean dimensions:", dim(incidence_per10k_mean), "\n")
cat("reduction_extended_mean dimensions:", dim(reduction_extended_mean), "\n")

# ============================================================================
# FIGURE 15: Sunburst and Heatmap (Original Incidence)
# ============================================================================
cat("\n=== Creating Figure 15: Original Sunburst and Heatmap ===\n")

# Build scenario labels for extended scenarios
durations <- c(3, 7, 10, 14)
intensities <- c(-0.50, -0.30, -0.20, -0.10, -0.05, 0.05, 0.10, 0.20, 0.30, 0.50)
shifts <- c(-7, -5, -3, 0, 3, 5, 7)

# Calculate mean values for extended scenarios
incidence_per10k_mean <- apply(incidence_per10k_extended, c(2,3), mean)
reduction_extended_mean <- apply(reduction_extended, c(2,3), mean)

# Function to create sunburst for original incidence
create_incidence_sunburst <- function(strain_idx, incidence_data, plot_title) {
  sunburst_data <- data.frame()
  
  for (dur_idx in 1:4) {
    dur_days <- durations[dur_idx]
    
    for (int_idx in 1:10) {
      int_val <- intensities[int_idx]
      
      temp_incidences <- numeric(7)
      for (shift_idx in 1:7) {
        scenario_idx <- 1 + (dur_idx-1)*70 + (int_idx-1)*7 + shift_idx
        temp_incidences[shift_idx] <- incidence_data[scenario_idx, strain_idx]
      }
      
      max_inc_group <- max(temp_incidences)
      min_inc_group <- min(temp_incidences)
      
      for (shift_idx in 1:7) {
        shift_val <- shifts[shift_idx]
        scenario_idx <- 1 + (dur_idx-1)*70 + (int_idx-1)*7 + shift_idx
        incidence_val <- incidence_data[scenario_idx, strain_idx]
        
        if (shift_val < 0) {
          shift_label <- paste0("E", abs(shift_val))
        } else if (shift_val > 0) {
          shift_label <- paste0("D", shift_val)
        } else {
          shift_label <- "On"
        }
        
        segment_idx <- (dur_idx-1)*70 + (int_idx-1)*7 + shift_idx
        
        if (max_inc_group > min_inc_group) {
          extension_length <- 0.3 + 0.9 * (incidence_val - min_inc_group) /
            (max_inc_group - min_inc_group)
        } else {
          extension_length <- 0.75
        }
        
        theta_start <- (segment_idx - 1) * (360 / 280)
        theta_end <- segment_idx * (360 / 280)
        theta_mid <- (theta_start + theta_end) / 2
        
        sunburst_data <- rbind(sunburst_data, data.frame(
          duration = dur_days,
          intensity = int_val * 100,
          shift = shift_val,
          shift_label = shift_label,
          incidence = incidence_val,
          segment_idx = segment_idx,
          dur_idx = dur_idx,
          int_idx = int_idx,
          shift_idx = shift_idx,
          extension_length = extension_length,
          theta_start = theta_start,
          theta_end = theta_end,
          theta_mid = theta_mid,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  duration_data <- sunburst_data %>%
    group_by(dur_idx) %>%
    summarise(
      theta_start = min(theta_start),
      theta_end = max(theta_end),
      theta_mid = (min(theta_start) + max(theta_end)) / 2,
      duration = first(duration),
      .groups = "drop"
    )
  
  intensity_data <- sunburst_data %>%
    group_by(dur_idx, int_idx) %>%
    summarise(
      theta_start = min(theta_start),
      theta_end = max(theta_end),
      theta_mid = (min(theta_start) + max(theta_end)) / 2,
      intensity = first(intensity),
      .groups = "drop"
    )
  
  p <- ggplot() +
    geom_rect(
      data = sunburst_data,
      aes(xmin = 2.5,
          xmax = 2.5 + extension_length * 1.5,
          ymin = theta_start,
          ymax = theta_end,
          fill = incidence),
      color = "white", size = 0.15
    ) +
    geom_rect(
      data = intensity_data,
      aes(xmin = 1.5, xmax = 2.5,
          ymin = theta_start, ymax = theta_end),
      fill = "gray85", color = "white", size = 0.25
    ) +
    geom_rect(
      data = duration_data,
      aes(xmin = 0.5, xmax = 1.5,
          ymin = theta_start, ymax = theta_end),
      fill = "gray70", color = "white", size = 0.4
    ) +
    geom_text(
      data = duration_data,
      aes(x = 1, y = theta_mid,
          label = paste0(duration, "d"),
          angle = ifelse(theta_mid > 180, theta_mid - 90, theta_mid + 90)),
      size = 2.5, fontface = "bold", color = "black"
    ) +
    geom_text(
      data = intensity_data,
      aes(x = 2, y = theta_mid,
          label = paste0(intensity, "%"),
          angle = ifelse(theta_mid > 180, theta_mid - 90, theta_mid + 90)),
      size = 1.2, color = "black"
    ) +
    geom_text(
      data = sunburst_data,
      aes(x = 2.5 + (extension_length * 1.5) * 0.5,
          y = theta_mid,
          label = shift_label,
          angle = ifelse(theta_mid > 180, theta_mid - 90, theta_mid + 90)),
      size = 0.6, color = "white", fontface = "bold"
    ) +
    scale_fill_gradientn(
      colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
      values = scales::rescale(c(
        min(sunburst_data$incidence),
        quantile(sunburst_data$incidence, 0.25),
        quantile(sunburst_data$incidence, 0.40),
        quantile(sunburst_data$incidence, 0.50),
        quantile(sunburst_data$incidence, 0.60),
        quantile(sunburst_data$incidence, 0.75),
        quantile(sunburst_data$incidence, 0.90),
        max(sunburst_data$incidence)
      )),
      limits = c(min(sunburst_data$incidence), max(sunburst_data$incidence)),
      name = "Incidence\n(per 10K)",
      guide = guide_colorbar(
        barwidth = 0.6,
        barheight = 6,
        title.position = "top",
        title.hjust = 0.5,
        title.theme = element_text(size = 7, face = "bold"),
        label.theme = element_text(size = 6)
      )
    ) +
    coord_polar(theta = "y") +
    xlim(0, 5.5) +
    labs(title = plot_title) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 3)),
      legend.position = "right",
      legend.justification = "center",
      plot.margin = margin(5, 5, 5, 5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

p15a_B <- create_incidence_sunburst(1, incidence_per10k_mean, "B")
p15a_H3 <- create_incidence_sunburst(2, incidence_per10k_mean, "H3")
p15a_H1 <- create_incidence_sunburst(3, incidence_per10k_mean, "H1")
p15a_COVID <- create_incidence_sunburst(4, incidence_per10k_mean, "COVID")
p15a_RSV <- create_incidence_sunburst(5, incidence_per10k_mean, "RSV")
p15a_HFMD <- create_incidence_sunburst(6, incidence_per10k_mean, "HFMD")

p15a_combined <- (p15a_B | p15a_H3 | p15a_H1) / (p15a_COVID | p15a_RSV | p15a_HFMD) +
  plot_annotation(
    title = "Typhoon Impact Scenarios: Case Incidence per 10,000 Population",
    subtitle = "Inner: Duration | Middle: Intensity | Outer: Timing | Each strain has independent color scale",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30", lineheight = 1.2)
    )
  )

print(p15a_combined)
ggsave("figure15a_sunburst_combined_6strains.png", p15a_combined,
       width = 20, height = 13, dpi = 300, bg = "white")

cat("Figure 15 complete\n\n")

cat("All original figures (1-15) now complete\n\n")

# ============================================================================
# NOW ADD NEW ATTACK RATE VISUALIZATIONS
# ============================================================================

cat("\n", rep("=", 100), "\n", sep="")
cat("CREATING NEW ATTACK RATE VISUALIZATIONS (FIGURES 16A-D)\n")
cat(rep("=", 100), "\n", sep="")

# ============================================================================
# FIGURE 16A: SUNBURST for Period 1 Attack Rates
# (2025-07-12 to 2025-09-27 - Before Typhoon)
# ============================================================================
cat("\n=== Creating Figure 16A: Attack Rate Sunburst - Period 1 (Before Typhoon) ===\n")

# Convert attack rates to per 10k population
attack_rate_period1_per10k <- attack_rate_period1_extended_mean * 10000

# Define scenario parameters
durations <- c(3, 7, 10, 14)
intensities <- c(-0.50, -0.30, -0.20, -0.10, -0.05, 0.05, 0.10, 0.20, 0.30, 0.50)
shifts <- c(-7, -5, -3, 0, 3, 5, 7)

create_attack_rate_sunburst <- function(strain_idx, attack_rate_data, plot_title) {
  sunburst_data <- data.frame()
  
  # Build hierarchical structure
  for (dur_idx in 1:4) {
    dur_days <- durations[dur_idx]
    
    for (int_idx in 1:10) {
      int_val <- intensities[int_idx]
      
      # Extract all 7 timing attack rates for this group
      temp_attack_rates <- numeric(7)
      for (shift_idx in 1:7) {
        scenario_idx <- 1 + (dur_idx-1)*70 + (int_idx-1)*7 + shift_idx
        temp_attack_rates[shift_idx] <- attack_rate_data[scenario_idx, strain_idx]
      }
      
      # Calculate group-wise min/max for normalization
      max_ar_group <- max(temp_attack_rates)
      min_ar_group <- min(temp_attack_rates)
      
      # Now process each of the 7 timing scenarios
      for (shift_idx in 1:7) {
        shift_val <- shifts[shift_idx]
        scenario_idx <- 1 + (dur_idx-1)*70 + (int_idx-1)*7 + shift_idx
        attack_rate_val <- attack_rate_data[scenario_idx, strain_idx]
        
        # Timing labels
        if (shift_val < 0) {
          shift_label <- paste0("E", abs(shift_val))
        } else if (shift_val > 0) {
          shift_label <- paste0("D", shift_val)
        } else {
          shift_label <- "On"
        }
        
        segment_idx <- (dur_idx-1)*70 + (int_idx-1)*7 + shift_idx
        
        # Group-wise normalization for EXTENSION LENGTH
        if (max_ar_group > min_ar_group) {
          extension_length <- 0.3 + 0.9 * (attack_rate_val - min_ar_group) /
            (max_ar_group - min_ar_group)
        } else {
          extension_length <- 0.75
        }
        
        # Calculate angular positions (280 equal segments)
        theta_start <- (segment_idx - 1) * (360 / 280)
        theta_end <- segment_idx * (360 / 280)
        theta_mid <- (theta_start + theta_end) / 2
        
        sunburst_data <- rbind(sunburst_data, data.frame(
          duration = dur_days,
          intensity = int_val * 100,
          shift = shift_val,
          shift_label = shift_label,
          attack_rate = attack_rate_val,
          segment_idx = segment_idx,
          dur_idx = dur_idx,
          int_idx = int_idx,
          shift_idx = shift_idx,
          extension_length = extension_length,
          theta_start = theta_start,
          theta_end = theta_end,
          theta_mid = theta_mid,
          group_id = paste(dur_idx, int_idx, sep = "_"),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Aggregated data for inner and middle rings
  duration_data <- sunburst_data %>%
    group_by(dur_idx) %>%
    summarise(
      theta_start = min(theta_start),
      theta_end = max(theta_end),
      theta_mid = (min(theta_start) + max(theta_end)) / 2,
      duration = first(duration),
      .groups = "drop"
    )
  
  intensity_data <- sunburst_data %>%
    group_by(dur_idx, int_idx) %>%
    summarise(
      theta_start = min(theta_start),
      theta_end = max(theta_end),
      theta_mid = (min(theta_start) + max(theta_end)) / 2,
      intensity = first(intensity),
      .groups = "drop"
    )
  
  # Create plot
  p <- ggplot() +
    # Outer ring: variable xmax for extension
    geom_rect(
      data = sunburst_data,
      aes(xmin = 2.5,
          xmax = 2.5 + extension_length * 1.5,
          ymin = theta_start,
          ymax = theta_end,
          fill = attack_rate),
      color = "white", size = 0.15
    ) +
    
    # Middle ring: Intensity levels
    geom_rect(
      data = intensity_data,
      aes(xmin = 1.5, xmax = 2.5,
          ymin = theta_start, ymax = theta_end),
      fill = "gray85", color = "white", size = 0.25
    ) +
    
    # Inner ring: Duration levels
    geom_rect(
      data = duration_data,
      aes(xmin = 0.5, xmax = 1.5,
          ymin = theta_start, ymax = theta_end),
      fill = "gray70", color = "white", size = 0.4
    ) +
    
    # Labels
    geom_text(
      data = duration_data,
      aes(x = 1, y = theta_mid,
          label = paste0(duration, "d"),
          angle = ifelse(theta_mid > 180, theta_mid - 90, theta_mid + 90)),
      size = 2.5, fontface = "bold", color = "black"
    ) +
    
    geom_text(
      data = intensity_data,
      aes(x = 2, y = theta_mid,
          label = paste0(intensity, "%"),
          angle = ifelse(theta_mid > 180, theta_mid - 90, theta_mid + 90)),
      size = 1.2, color = "black"
    ) +
    
    geom_text(
      data = sunburst_data,
      aes(x = 2.5 + (extension_length * 1.5) * 0.5,
          y = theta_mid,
          label = shift_label,
          angle = ifelse(theta_mid > 180, theta_mid - 90, theta_mid + 90)),
      size = 0.6, color = "white", fontface = "bold"
    )
  
  # === BEGIN FIX FOR CONSTANT DATA ERROR ===
  min_ar <- min(sunburst_data$attack_rate)
  max_ar <- max(sunburst_data$attack_rate)
  
  if (abs(max_ar - min_ar) < 1e-6) {
    # 数据是恒定的（Period 1）。创建极小的非零范围和固定颜色集，以避免插值错误。
    fill_limits <- c(min_ar, min_ar + 1e-5)
    fill_values <- c(0, 1)
    # 使用中间颜色，确保图表有一致的颜色。
    fill_colors <- c("#D1E5F0", "#D1E5F0") 
    
  } else {
    # 数据是可变的。使用原始的、基于分位数的颜色刻度。
    fill_limits <- c(min_ar, max_ar)
    fill_values <- scales::rescale(c(
      min_ar,
      quantile(sunburst_data$attack_rate, 0.25),
      quantile(sunburst_data$attack_rate, 0.40),
      quantile(sunburst_data$attack_rate, 0.50),
      quantile(sunburst_data$attack_rate, 0.60),
      quantile(sunburst_data$attack_rate, 0.75),
      quantile(sunburst_data$attack_rate, 0.90),
      max_ar
    ))
    fill_colors <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
  }
  
  p <- p + scale_fill_gradientn(
    colors = fill_colors,
    values = fill_values,
    limits = fill_limits,
    name = "Attack Rate\n(per 10K)",
    guide = guide_colorbar(
      barwidth = 0.6,
      barheight = 6,
      title.position = "top",
      title.hjust = 0.5,
      title.theme = element_text(size = 7, face = "bold"),
      label.theme = element_text(size = 6)
    )
  )
  # === END FIX FOR CONSTANT DATA ERROR ===
  
  p <- p + coord_polar(theta = "y") +
    xlim(0, 5.5) +
    
    labs(title = plot_title) +
    
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 3)),
      legend.position = "right",
      legend.justification = "center",
      plot.margin = margin(5, 5, 5, 5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}
# Create sunburst plots for Period 1
p16a_B <- create_attack_rate_sunburst(1, attack_rate_period1_per10k, "B")
p16a_H3 <- create_attack_rate_sunburst(2, attack_rate_period1_per10k, "H3")
p16a_H1 <- create_attack_rate_sunburst(3, attack_rate_period1_per10k, "H1")
p16a_COVID <- create_attack_rate_sunburst(4, attack_rate_period1_per10k, "COVID")
p16a_RSV <- create_attack_rate_sunburst(5, attack_rate_period1_per10k, "RSV")
p16a_HFMD <- create_attack_rate_sunburst(6, attack_rate_period1_per10k, "HFMD")

# Combine using patchwork
p16a_combined <- (p16a_B | p16a_H3 | p16a_H1) / (p16a_COVID | p16a_RSV | p16a_HFMD) +
  plot_annotation(
    title = "Attack Rate by Scenario - Period 1 (2025-07-12 to 2025-09-27, Before Typhoon)",
    subtitle = "Inner: Duration | Middle: Intensity | Outer: Timing (extension = relative attack rate within group)\nHistorical data only - Same across scenarios (no typhoon effect yet) | Each strain has independent color scale",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30", lineheight = 1.2)
    )
  )

print(p16a_combined)
ggsave(
  filename = "figure16a_attack_rate_period1_sunburst.png",
  plot = p16a_combined,
  width = 20,
  height = 13,
  dpi = 300,
  bg = "white"
)

cat("Figure 16A saved: figure16a_attack_rate_period1_sunburst.png\n\n")

# ============================================================================
# FIGURE 16B: SUNBURST for Period 2 Attack Rates
# (2025-07-12 to 2025-10-18 - Including Typhoon + Recovery)
# ============================================================================
cat("\n=== Creating Figure 16B: Attack Rate Sunburst - Period 2 (Including Typhoon) ===\n")

# Convert attack rates to per 10k population
attack_rate_period2_per10k <- attack_rate_period2_extended_mean * 10000

# Create sunburst plots for Period 2
p16b_B <- create_attack_rate_sunburst(1, attack_rate_period2_per10k, "B")
p16b_H3 <- create_attack_rate_sunburst(2, attack_rate_period2_per10k, "H3")
p16b_H1 <- create_attack_rate_sunburst(3, attack_rate_period2_per10k, "H1")
p16b_COVID <- create_attack_rate_sunburst(4, attack_rate_period2_per10k, "COVID")
p16b_RSV <- create_attack_rate_sunburst(5, attack_rate_period2_per10k, "RSV")
p16b_HFMD <- create_attack_rate_sunburst(6, attack_rate_period2_per10k, "HFMD")

# Combine using patchwork
p16b_combined <- (p16b_B | p16b_H3 | p16b_H1) / (p16b_COVID | p16b_RSV | p16b_HFMD) +
  plot_annotation(
    title = "Attack Rate by Scenario - Period 2 (2025-07-12 to 2025-10-18, Including Typhoon + Recovery)",
    subtitle = "Inner: Duration | Middle: Intensity | Outer: Timing (extension = relative attack rate within group)\nIncludes historical + forecast with typhoon scenarios | Each strain has independent color scale",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30", lineheight = 1.2)
    )
  )

print(p16b_combined)
ggsave(
  filename = "figure16b_attack_rate_period2_sunburst.png",
  plot = p16b_combined,
  width = 20,
  height = 13,
  dpi = 300,
  bg = "white"
)

cat("Figure 16B saved: figure16b_attack_rate_period2_sunburst.png\n\n")

# ============================================================================
# FIGURE 16C: HEATMAP for Period 1 Attack Rates
# (2025-07-12 to 2025-09-27 - Before Typhoon)
# ============================================================================
cat("\n=== Creating Figure 16C: Attack Rate Heatmap - Period 1 (Before Typhoon) ===\n")

# Prepare heatmap data for Period 1
heatmap_period1_list <- list()

for (dur_idx in 1:4) {
  for (int_idx in 1:10) {
    for (shift_idx in 1:7) {
      scenario_idx <- 1 + (dur_idx-1)*70 + (int_idx-1)*7 + shift_idx
      
      for (strain_idx in 1:N_strains) {
        heatmap_period1_list[[length(heatmap_period1_list) + 1]] <- data.frame(
          duration = factor(durations[dur_idx], levels = durations),
          intensity = factor(paste0(as.integer(intensities[int_idx]*100), "%")),
          shift = factor(shifts[shift_idx]),
          attack_rate = attack_rate_period1_per10k[scenario_idx, strain_idx],
          strain = strain_names[strain_idx],
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

heatmap_period1_df <- do.call(rbind, heatmap_period1_list)
heatmap_period1_df$strain <- factor(heatmap_period1_df$strain, levels = strain_names)
heatmap_period1_df$intensity <- factor(heatmap_period1_df$intensity, 
                                       levels = c("-50%", "-30%", "-20%", "-10%", "-5%", 
                                                  "5%", "10%", "20%", "30%", "50%"))

p16c_heatmap <- ggplot(heatmap_period1_df, aes(x = shift, y = intensity, fill = attack_rate)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%.1f", attack_rate)), size = 1.5, color = "black") +
  facet_grid(duration ~ strain,
             labeller = labeller(duration = function(x) paste(x, "days"))) +
  scale_fill_gradientn(
    colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
    name = "Attack Rate\n(per 10K)",
    limits = c(min(heatmap_period1_df$attack_rate), max(heatmap_period1_df$attack_rate))
  ) +
  scale_x_discrete(labels = c("-7" = "E7", "-5" = "E5", "-3" = "E3", "0" = "On",
                              "3" = "D3", "5" = "D5", "7" = "D7")) +
  labs(
    title = "Attack Rate Heatmap - Period 1 (2025-07-12 to 2025-09-27, Before Typhoon)",
    subtitle = "Values show attack rate per 10,000 population | Historical data only (same across scenarios)\nColumns: Strains | Rows: Duration | X: Timing | Y: Intensity",
    x = "Timing Shift (days)",
    y = "Transmission Change Intensity"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
    strip.text.x = element_text(size = 9, face = "bold"),
    strip.text.y = element_text(size = 9, face = "bold", angle = 0),
    strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
    axis.text.x = element_text(angle = 0, size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 8, color = "gray30")
  )

print(p16c_heatmap)
ggsave("figure16c_attack_rate_period1_heatmap.png", p16c_heatmap,
       width = 16, height = 10, dpi = 300, bg = "white")

cat("Figure 16C saved: figure16c_attack_rate_period1_heatmap.png\n\n")

# ============================================================================
# FIGURE 16D: HEATMAP for Period 2 Attack Rates
# (2025-07-12 to 2025-10-18 - Including Typhoon + Recovery)
# ============================================================================
cat("\n=== Creating Figure 16D: Attack Rate Heatmap - Period 2 (Including Typhoon) ===\n")

# Prepare heatmap data for Period 2
heatmap_period2_list <- list()

for (dur_idx in 1:4) {
  for (int_idx in 1:10) {
    for (shift_idx in 1:7) {
      scenario_idx <- 1 + (dur_idx-1)*70 + (int_idx-1)*7 + shift_idx
      
      for (strain_idx in 1:N_strains) {
        heatmap_period2_list[[length(heatmap_period2_list) + 1]] <- data.frame(
          duration = factor(durations[dur_idx], levels = durations),
          intensity = factor(paste0(as.integer(intensities[int_idx]*100), "%")),
          shift = factor(shifts[shift_idx]),
          attack_rate = attack_rate_period2_per10k[scenario_idx, strain_idx],
          strain = strain_names[strain_idx],
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

heatmap_period2_df <- do.call(rbind, heatmap_period2_list)
heatmap_period2_df$strain <- factor(heatmap_period2_df$strain, levels = strain_names)
heatmap_period2_df$intensity <- factor(heatmap_period2_df$intensity, 
                                       levels = c("-50%", "-30%", "-20%", "-10%", "-5%", 
                                                  "5%", "10%", "20%", "30%", "50%"))

p16d_heatmap <- ggplot(heatmap_period2_df, aes(x = shift, y = intensity, fill = attack_rate)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%.1f", attack_rate)), size = 1.5, color = "black") +
  facet_grid(duration ~ strain,
             labeller = labeller(duration = function(x) paste(x, "days"))) +
  scale_fill_gradientn(
    colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
    name = "Attack Rate\n(per 10K)",
    limits = c(min(heatmap_period2_df$attack_rate), max(heatmap_period2_df$attack_rate))
  ) +
  scale_x_discrete(labels = c("-7" = "E7", "-5" = "E5", "-3" = "E3", "0" = "On",
                              "3" = "D3", "5" = "D5", "7" = "D7")) +
  labs(
    title = "Attack Rate Heatmap - Period 2 (2025-07-12 to 2025-10-18, Including Typhoon + Recovery)",
    subtitle = "Values show attack rate per 10,000 population | Includes historical + forecast with scenarios\nColumns: Strains | Rows: Duration | X: Timing | Y: Intensity",
    x = "Timing Shift (days)",
    y = "Transmission Change Intensity"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
    strip.text.x = element_text(size = 9, face = "bold"),
    strip.text.y = element_text(size = 9, face = "bold", angle = 0),
    strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
    axis.text.x = element_text(angle = 0, size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 8, color = "gray30")
  )

print(p16d_heatmap)
ggsave("figure16d_attack_rate_period2_heatmap.png", p16d_heatmap,
       width = 16, height = 10, dpi = 300, bg = "white")

cat("Figure 16D saved: figure16d_attack_rate_period2_heatmap.png\n\n")

# ============================================================================
# CREATE SUMMARY TABLES FOR ATTACK RATES
# ============================================================================
cat("\n=== Creating Attack Rate Summary Tables ===\n")

# Period 1 Summary
period1_summary <- data.frame(
  Strain = rep(strain_names, each = 281),
  Scenario_Index = rep(1:281, times = 6),
  Attack_Rate_Historical = as.vector(t(replicate(281, attack_rate_period1_historical_mean))),
  Attack_Rate_Per_10K = as.vector(t(attack_rate_period1_per10k)),
  stringsAsFactors = FALSE
)

# Add scenario labels
get_scenario_label <- function(idx) {
  if (idx == 1) return("Baseline")
  
  dur_idx <- ((idx - 2) %/% 70) + 1
  int_idx <- (((idx - 2) %% 70) %/% 7) + 1
  shift_idx <- ((idx - 2) %% 7) + 1
  
  dur <- durations[dur_idx]
  int <- intensities[int_idx]
  shift <- shifts[shift_idx]
  
  int_label <- paste0(as.integer(int * 100), "%")
  if (shift < 0) {
    shift_label <- paste0("E", abs(shift))
  } else if (shift > 0) {
    shift_label <- paste0("D", shift)
  } else {
    shift_label <- "On"
  }
  
  return(paste0(dur, "d_", int_label, "_", shift_label))
}

period1_summary$Scenario <- sapply(period1_summary$Scenario_Index, get_scenario_label)

write.csv(period1_summary, 
          file = "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/attack_rate_period1_summary.csv",
          row.names = FALSE)

# Period 2 Summary
period2_summary <- data.frame(
  Strain = rep(strain_names, each = 281),
  Scenario_Index = rep(1:281, times = 6),
  Attack_Rate_Historical = as.vector(t(replicate(281, attack_rate_period2_historical_mean))),
  Attack_Rate_Baseline = as.vector(t(replicate(281, attack_rate_period2_baseline_mean))),
  Attack_Rate_Per_10K = as.vector(t(attack_rate_period2_per10k)),
  stringsAsFactors = FALSE
)

period2_summary$Scenario <- sapply(period2_summary$Scenario_Index, get_scenario_label)

# Calculate difference from baseline
period2_summary$Diff_From_Baseline <- period2_summary$Attack_Rate_Per_10K - 
  period2_summary$Attack_Rate_Baseline * 10000

write.csv(period2_summary, 
          file = "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/attack_rate_period2_summary.csv",
          row.names = FALSE)

cat("Attack rate summary tables saved:\n")
cat("  - attack_rate_period1_summary.csv\n")
cat("  - attack_rate_period2_summary.csv\n\n")

# ============================================================================
# PRINT SUMMARY STATISTICS
# ============================================================================
cat("\n", rep("=", 100), "\n", sep="")
cat("ATTACK RATE SUMMARY STATISTICS\n")
cat(rep("=", 100), "\n", sep="")

cat("\nPeriod 1 (2025-07-12 to 2025-09-27, Before Typhoon):\n")
cat("All scenarios are identical (historical data only)\n\n")
for (i in 1:N_strains) {
  cat(sprintf("%s: %.6f (%.3f per 10K)\n", 
              strain_names[i],
              attack_rate_period1_historical_mean[i],
              attack_rate_period1_historical_mean[i] * 10000))
}

cat("\nPeriod 2 (2025-07-12 to 2025-10-18, Including Typhoon + Recovery):\n")
cat("Baseline Attack Rates:\n")
for (i in 1:N_strains) {
  cat(sprintf("%s: %.6f (%.3f per 10K)\n", 
              strain_names[i],
              attack_rate_period2_baseline_mean[i],
              attack_rate_period2_baseline_mean[i] * 10000))
}

cat("\nPeriod 2 - Selected Scenario Comparisons:\n")
# Show a few key scenarios (e.g., 7d -50%, 7d 0%, 7d +50%)
scenario_indices_to_show <- c(1, 12, 17, 21)  # Baseline, 7d -50%, 7d 0%, 7d +50%
scenario_names_to_show <- c("Baseline", "7d -50% On-time", "7d 0% On-time", "7d +50% On-time")

for (sc_idx in 1:length(scenario_indices_to_show)) {
  scenario_extended_idx <- 1 + (1)*70 + (scenario_indices_to_show[sc_idx] - 12)*7 + 4  # Shift=0
  if (scenario_indices_to_show[sc_idx] == 1) scenario_extended_idx <- 1
  
  cat(sprintf("\n%s:\n", scenario_names_to_show[sc_idx]))
  for (i in 1:N_strains) {
    ar_val <- attack_rate_period2_per10k[scenario_extended_idx, i]
    diff_from_baseline <- ar_val - attack_rate_period2_baseline_mean[i] * 10000
    cat(sprintf("  %s: %.3f per 10K (%.3f from baseline)\n", 
                strain_names[i],
                ar_val,
                diff_from_baseline))
  }
}

cat("\n", rep("=", 100), "\n", sep="")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 100), "\n", sep="")

cat("\nGenerated New Figures:\n")
cat("  Figure 16A: Attack Rate Sunburst - Period 1 (Before Typhoon)\n")
cat("  Figure 16B: Attack Rate Sunburst - Period 2 (Including Typhoon)\n")
cat("  Figure 16C: Attack Rate Heatmap - Period 1 (Before Typhoon)\n")
cat("  Figure 16D: Attack Rate Heatmap - Period 2 (Including Typhoon)\n\n")

cat("Generated Summary Tables:\n")
cat("  - attack_rate_period1_summary.csv\n")
cat("  - attack_rate_period2_summary.csv\n\n")

cat("All visualizations and analyses complete!\n")
cat(rep("=", 100), "\n\n")

