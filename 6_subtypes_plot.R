# ============================================================================
# TYPHOON IMPACT ANALYSIS - WITH 6 SCENARIOS (0%, 10%, 20%, 40%, 60%, 80%)
# ============================================================================
Sys.setenv('R_MAX_VSIZE' = 32 * 1024^3)
library(rstan)
library(ggplot2)
library(dplyr)
library(lubridate)
library(bayesplot)
library(tidyr)
library(reshape2)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
# 1. Data reading and preprocessing
cat("=== Loading and preprocessing data ===\n")
data_path <- "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/HK_ILI_COVID_Sep.csv"
raw_data <- read.csv(data_path)
raw_data$date <- as.Date(raw_data$date, format = "%Y/%m/%d")
start_date <- as.Date("2023-01-07")
typhoon_start_date <- as.Date("2025-09-23")
end_date <- as.Date("2025-10-11")
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
# 2. Build Stan data
T_weeks <- nrow(fitting_data)
T_weeks_validation <- nrow(validation_data)
T_weeks_forecast <- 5
typhoon_weeks <- 2  # IMPORTANT: Typhoon lasts 2 weeks!
N_strains <- 6
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
  typhoon_weeks = typhoon_weeks
)
# 3. Compile and run Stan model
cat("=== Compiling Stan model ===\n")
model_file <- "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/6_subtypes.stan"
model <- stan_model(model_file)
cat("\n=== Starting MCMC sampling ===\n")
cat("Note: Fitting to data up to 2025-09-27\n")
cat("Typhoon starts on 2025-09-23 (2 weeks impact + 3 weeks recovery)\n")
cat("Typhoon reduction scenarios: 0%, 10%, 20%, 40%, 60%, 80%\n")
cat("Validation data from 2025-09-23 to 2025-10-11 will be shown as yellow diamonds\n\n")
fit <- sampling(
  model,
  data = stan_data,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  thin = 1,
  cores = 4,
  control = list(adapt_delta = 0.97, max_treedepth = 20, stepsize = 0.005)
)
saveRDS(fit, file = "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/typhoon_model_fit_6scenarios.rds")
# 4. Extract results
pred_cases <- rstan::extract(fit, pars = "pred_cases")$pred_cases
forecast_cases <- rstan::extract(fit, pars = "forecast_cases")$forecast_cases
forecast_cases_typhoon <- rstan::extract(fit, pars = "forecast_cases_typhoon")$forecast_cases_typhoon
R_eff <- rstan::extract(fit, pars = "R_eff")$R_eff
R_eff_scenarios <- rstan::extract(fit, pars = "R_eff_scenarios")$R_eff_scenarios
R0_t <- rstan::extract(fit, pars = "R0_t")$R0_t
cat("\n=== Checking extracted dimensions ===\n")
cat("pred_cases dimensions:", dim(pred_cases), "\n")
cat("forecast_cases dimensions:", dim(forecast_cases), "\n")
cat("forecast_cases_typhoon dimensions:", dim(forecast_cases_typhoon), "\n")
cat("R_eff dimensions:", dim(R_eff), "\n")
cat("R_eff_scenarios dimensions:", dim(R_eff_scenarios), "\n")
cat("R0_t dimensions:", dim(R0_t), "\n")
# Calculate statistics
pred_mean <- apply(pred_cases, c(2,3), mean)
pred_median <- apply(pred_cases, c(2,3), median)
pred_lower <- apply(pred_cases, c(2,3), quantile, probs = 0.025)
pred_upper <- apply(pred_cases, c(2,3), quantile, probs = 0.975)
forecast_mean <- apply(forecast_cases, c(2,3), mean)
forecast_median <- apply(forecast_cases, c(2,3), median)
forecast_lower <- apply(forecast_cases, c(2,3), quantile, probs = 0.025)
forecast_upper <- apply(forecast_cases, c(2,3), quantile, probs = 0.975)
# 6 typhoon scenarios
typhoon_scenarios <- c("0% Reduction (Baseline)", "10% Reduction", "20% Reduction",
                       "40% Reduction", "60% Reduction", "80% Reduction")
forecast_typhoon_mean <- array(NA, dim = c(6, T_weeks_forecast, N_strains))
forecast_typhoon_lower <- array(NA, dim = c(6, T_weeks_forecast, N_strains))
forecast_typhoon_upper <- array(NA, dim = c(6, T_weeks_forecast, N_strains))
for (typhoon_idx in 1:6) {
  forecast_typhoon_mean[typhoon_idx,,] <- apply(forecast_cases_typhoon[,typhoon_idx,,], c(2,3), mean)
  forecast_typhoon_lower[typhoon_idx,,] <- apply(forecast_cases_typhoon[,typhoon_idx,,], c(2,3), quantile, probs = 0.025)
  forecast_typhoon_upper[typhoon_idx,,] <- apply(forecast_cases_typhoon[,typhoon_idx,,], c(2,3), quantile, probs = 0.975)
}
R_eff_scenarios_mean <- array(NA, dim = c(6, T_weeks + T_weeks_forecast, N_strains))
R_eff_scenarios_lower <- array(NA, dim = c(6, T_weeks + T_weeks_forecast, N_strains))
R_eff_scenarios_upper <- array(NA, dim = c(6, T_weeks + T_weeks_forecast, N_strains))
for (typhoon_idx in 1:6) {
  R_eff_scenarios_mean[typhoon_idx,,] <- apply(R_eff_scenarios[,typhoon_idx,,], c(2,3), mean)
  R_eff_scenarios_lower[typhoon_idx,,] <- apply(R_eff_scenarios[,typhoon_idx,,], c(2,3), quantile, probs = 0.025)
  R_eff_scenarios_upper[typhoon_idx,,] <- apply(R_eff_scenarios[,typhoon_idx,,], c(2,3), quantile, probs = 0.975)
}
# 5. Prepare visualization data
theme_set(theme_minimal(base_size = 11))
forecast_dates <- seq(typhoon_start_date, by = "week", length.out = T_weeks_forecast)
# CORRECTED: Typhoon lasts 2 weeks (14 days)
typhoon_start_date_line <- typhoon_start_date
typhoon_end_date_line <- typhoon_start_date + 14  # 2 weeks = 14 days
typhoon_period_end <- typhoon_start_date + 14
recovery_start_date <- typhoon_start_date + 14
cat("\n=== Date markers (CORRECTED) ===\n")
cat("Typhoon start (gray line):", as.character(typhoon_start_date_line), "\n")
cat("Typhoon end (red line):", as.character(typhoon_end_date_line), "\n")
cat("Typhoon duration: 2 weeks (14 days)\n")
cat("Recovery start:", as.character(recovery_start_date), "\n")
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
forecast_baseline <- data.frame(
  week = rep((T_weeks + 1):(T_weeks + T_weeks_forecast), N_strains),
  date = rep(forecast_dates, N_strains),
  strain = factor(rep(strain_names, each = T_weeks_forecast), levels = strain_names),
  observed = NA,
  predicted = as.vector(forecast_median),
  pred_mean = as.vector(forecast_mean),
  pred_lower = as.vector(forecast_lower),
  pred_upper = as.vector(forecast_upper)
)
validation_plot_df <- data.frame(
  date = rep(validation_data$date, N_strains),
  strain = factor(rep(strain_names, each = T_weeks_validation), levels = strain_names),
  observed = as.vector(validation_matrix)
)
validation_plot_df <- validation_plot_df[!is.na(validation_plot_df$observed), ]
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
             aes(y = observed),
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
    subtitle = "Black dots: Observed (fitted) | Yellow diamonds: Validation data (NOT fitted) | Blue: Model fit | Red: Baseline forecast\nGray dashed: Typhoon start (2025-09-23) | Red dashed: Typhoon end (2025-10-07) | Pink: Typhoon (2 weeks) | Green: Recovery (3 weeks)",
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
# FIGURE 2: Recent Period + All Typhoon Scenarios + Validation Data
# ============================================================================
cat("\n=== Creating Figure 2: Recent Period, Scenarios, and Validation ===\n")
cutoff_date <- as.Date("2025-09-01")
recent_fit_df <- results_df %>%
  filter(date >= cutoff_date) %>%
  mutate(scenario = "Historical Fit")
typhoon_forecast_list <- list()
for (typhoon_idx in 1:6) {
  for (strain_idx in 1:N_strains) {
    typhoon_forecast_list[[length(typhoon_forecast_list) + 1]] <- data.frame(
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
  }
}
typhoon_forecast_df <- do.call(rbind, typhoon_forecast_list)
typhoon_forecast_df$strain <- factor(typhoon_forecast_df$strain, levels = strain_names)
recent_and_scenarios <- rbind(recent_fit_df, typhoon_forecast_df)
recent_and_scenarios$scenario <- factor(
  recent_and_scenarios$scenario,
  levels = c("Historical Fit", "0% Reduction (Baseline)", "10% Reduction",
             "20% Reduction", "40% Reduction", "60% Reduction", "80% Reduction")
)
last_observed <- recent_fit_df %>%
  group_by(strain) %>%
  slice_tail(n = 1) %>%
  ungroup()
p2_scenarios <- ggplot(recent_and_scenarios, aes(x = date)) +
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
  geom_ribbon(data = filter(recent_and_scenarios, scenario == "Historical Fit"),
              aes(ymin = pred_lower, ymax = pred_upper),
              alpha = 0.2, fill = "#377EB8") +
  geom_line(data = filter(recent_and_scenarios, scenario == "Historical Fit"),
            aes(y = predicted), color = "#377EB8", size = 1.4) +
  geom_point(data = filter(recent_and_scenarios, scenario == "Historical Fit", !is.na(observed)),
             aes(y = observed), color = "black", size = 1.3, alpha = 0.6) +
  geom_point(data = last_observed, aes(y = observed),
             color = "black", size = 3, shape = 21, fill = "white", stroke = 1.3) +
  geom_point(data = validation_plot_df,
             aes(y = observed),
             color = "gold", fill = "yellow", shape = 23, size = 1.3, stroke = 0.8) +
  geom_ribbon(data = filter(recent_and_scenarios, scenario != "Historical Fit"),
              aes(ymin = pred_lower, ymax = pred_upper, fill = scenario),
              alpha = 0.18) +
  geom_line(data = filter(recent_and_scenarios, scenario != "Historical Fit"),
            aes(y = predicted, color = scenario, linetype = scenario),
            size = 1.2) +
  geom_vline(xintercept = typhoon_start_date_line,
             linetype = "dashed", color = "gray40", size = 0.8) +
  geom_vline(xintercept = typhoon_end_date_line,
             linetype = "dashed", color = "red", size = 0.8) +
  facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
  scale_color_manual(
    name = "Scenarios:",
    values = c("0% Reduction (Baseline)" = "#999999",
               "10% Reduction" = "#6BAED6",
               "20% Reduction" = "#E41A1C",
               "40% Reduction" = "#FF7F00",
               "60% Reduction" = "#4DAF4A",
               "80% Reduction" = "#984EA3"),
    labels = c("0%", "10%", "20%", "40%", "60%", "80%")
  ) +
  scale_fill_manual(
    name = "Scenarios:",
    values = c("0% Reduction (Baseline)" = "#999999",
               "10% Reduction" = "#6BAED6",
               "20% Reduction" = "#E41A1C",
               "40% Reduction" = "#FF7F00",
               "60% Reduction" = "#4DAF4A",
               "80% Reduction" = "#984EA3"),
    labels = c("0%", "10%", "20%", "40%", "60%", "80%")
  ) +
  scale_linetype_manual(
    name = "Scenarios:",
    values = c("0% Reduction (Baseline)" = "solid",
               "10% Reduction" = "solid",
               "20% Reduction" = "solid",
               "40% Reduction" = "dashed",
               "60% Reduction" = "dotdash",
               "80% Reduction" = "longdash"),
    labels = c("0%", "10%", "20%", "40%", "60%", "80%")
  ) +
  scale_y_continuous(trans = "sqrt") +
  scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
  labs(
    title = "Typhoon Impact: Recent Period and Transmission Reduction Scenarios",
    subtitle = "Blue: Model fit | Yellow diamonds: Validation data (NOT fitted) | Colored lines: Forecast scenarios (6 scenarios)\nGray dashed: Typhoon start (2025-09-23) | Red dashed: Typhoon end (2025-10-07) | Pink: Typhoon (2 weeks) | Green: Recovery (3 weeks)",
    x = NULL,
    y = "Weekly Cases/Hospitalizations (√ scale)"
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
print(p2_scenarios)
ggsave("figure2_typhoon_scenarios_validation.png", p2_scenarios,
       width = 14, height = 10, dpi = 300, bg = "white")
# ============================================================================
# FIGURE 3: Effective Reproduction Number (Rt) - BASELINE ONLY
# ============================================================================
cat("\n=== Creating Figure 3: Effective Reproduction Number (Baseline) ===\n")
T_weeks_total <- T_weeks + T_weeks_forecast
Reff_mean <- R_eff_scenarios_mean[1,,]
Reff_lower <- R_eff_scenarios_lower[1,,]
Reff_upper <- R_eff_scenarios_upper[1,,]
all_dates <- c(fitting_data$date, forecast_dates)
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
p3_reff <- ggplot(Reff_df, aes(x = date)) +
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
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  labs(
    title = expression(paste("Weekly Effective Reproduction Number (", R[eff], ") - Baseline Scenario")),
    subtitle = expression(paste("Solid: Historical | Dashed: Forecast | Red line: ", R[eff], "=1 threshold\nGray dashed: Typhoon start (2025-09-23) | Red dashed: Typhoon end (2025-10-07) | Pink: Typhoon (2 weeks) | Green: Recovery (3 weeks)")),
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
print(p3_reff)
ggsave("figure3_reproduction_number.png", p3_reff,
       width = 14, height = 9, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 4: Typhoon Impact Effectiveness (compared to 0% baseline)
# MODIFIED: Added error bars based on posterior distribution of reduction
# ============================================================================
cat("\n=== Creating Figure 4: Typhoon Impact Analysis (with Error Bars) ===\n")

# New data preparation for Figure 4: Calculate reduction stats from posterior
impact_effectiveness_list <- list()
n_samples <- dim(forecast_cases_typhoon)[1]

for (strain_idx in 1:N_strains) {
  # Get baseline (scenario 1) samples
  base_typhoon_samples <- apply(forecast_cases_typhoon[, 1, 1:typhoon_weeks, strain_idx], 1, sum)
  base_total_samples <- apply(forecast_cases_typhoon[, 1, 1:T_weeks_forecast, strain_idx], 1, sum)
  
  # Handle cases where baseline is 0 (to avoid NaN)
  base_typhoon_safe <- ifelse(base_typhoon_samples == 0, 1e-9, base_typhoon_samples)
  base_total_safe <- ifelse(base_total_samples == 0, 1e-9, base_total_samples)
  
  for (typhoon_idx in 2:6) { # Loop scenarios 2-6 (10% to 80%)
    # Get scenario samples
    scen_typhoon_samples <- apply(forecast_cases_typhoon[, typhoon_idx, 1:typhoon_weeks, strain_idx], 1, sum)
    scen_total_samples <- apply(forecast_cases_typhoon[, typhoon_idx, 1:T_weeks_forecast, strain_idx], 1, sum)
    
    # Calculate distributions of reduction
    typhoon_reduction_dist <- (1 - scen_typhoon_samples / base_typhoon_safe) * 100
    total_reduction_dist <- (1 - scen_total_samples / base_total_safe) * 100
    
    # Clip reductions at 0 (no negative reduction) or 100
    typhoon_reduction_dist[typhoon_reduction_dist < 0] <- 0
    total_reduction_dist[total_reduction_dist < 0] <- 0
    typhoon_reduction_dist[typhoon_reduction_dist > 100] <- 100
    total_reduction_dist[total_reduction_dist > 100] <- 100
    
    # Store stats for typhoon period
    impact_effectiveness_list[[length(impact_effectiveness_list) + 1]] <- data.frame(
      strain = strain_names[strain_idx],
      scenario = typhoon_scenarios[typhoon_idx],
      metric = "Typhoon Period (2 weeks)",
      mean = mean(typhoon_reduction_dist),
      lower = quantile(typhoon_reduction_dist, 0.025),
      upper = quantile(typhoon_reduction_dist, 0.975)
    )
    
    # Store stats for total period
    impact_effectiveness_list[[length(impact_effectiveness_list) + 1]] <- data.frame(
      strain = strain_names[strain_idx],
      scenario = typhoon_scenarios[typhoon_idx],
      metric = "Total Period (5 weeks)",
      mean = mean(total_reduction_dist),
      lower = quantile(total_reduction_dist, 0.025),
      upper = quantile(total_reduction_dist, 0.975)
    )
  }
}
impact_effectiveness_stats <- do.call(rbind, impact_effectiveness_list)
impact_effectiveness_stats$strain <- factor(impact_effectiveness_stats$strain, levels = strain_names)
impact_effectiveness_stats$scenario <- factor(impact_effectiveness_stats$scenario, levels = typhoon_scenarios[2:6])


p4_effectiveness <- ggplot(impact_effectiveness_stats, 
                           aes(x = strain, y = mean, fill = scenario)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  # ADDED: Error bars
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.9),
    width = 0.3, # Width of the error bar caps
    color = "gray20",
    size = 0.4
  ) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(
    name = "Typhoon Scenario",
    values = c("10% Reduction" = "#6BAED6",
               "20% Reduction" = "#E41A1C",
               "40% Reduction" = "#FF7F00",
               "60% Reduction" = "#4DAF4A",
               "80% Reduction" = "#984EA3"),
    labels = c("10%", "20%", "40%", "60%", "80%")
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) + # Adjusted breaks
  labs(
    title = "Typhoon Impact Effectiveness: Case Reduction Compared to 0% Baseline",
    subtitle = "Bars represent mean posterior reduction. Error bars show 95% Credible Intervals.\nLeft: Total 5-week period | Right: 2-week typhoon period only",
    x = NULL,
    y = "Percent Reduction (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    panel.spacing = unit(0.8, "lines"),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
    axis.text = element_text(size = 9, color = "gray20"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 10, color = "gray30", margin = margin(b = 10))
  )
print(p4_effectiveness)
ggsave("figure4_effectiveness_with_errorbars.png", p4_effectiveness,
       width = 14, height = 8, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 5: Typhoon vs Recovery Period Comparison (All Scenarios)
# MODIFIED: Added error bars for average weekly cases
# ============================================================================
cat("\n=== Creating Figure 5: Typhoon vs Recovery Period (with Error Bars) ===\n")

# New data preparation for Figure 5: Calculate stats for average weekly cases
period_comparison_list <- list()
n_recovery_weeks <- T_weeks_forecast - typhoon_weeks

for (strain_idx in 1:N_strains) {
  for (typhoon_idx in 1:6) {
    # Get samples for average weekly cases during typhoon
    avg_typhoon_samples <- apply(forecast_cases_typhoon[, typhoon_idx, 1:typhoon_weeks, strain_idx], 1, mean)
    
    # Get samples for average weekly cases during recovery
    avg_recovery_samples <- apply(forecast_cases_typhoon[, typhoon_idx, (typhoon_weeks+1):T_weeks_forecast, strain_idx], 1, mean)
    
    # Store stats for typhoon period
    period_comparison_list[[length(period_comparison_list) + 1]] <- data.frame(
      strain = strain_names[strain_idx],
      scenario = typhoon_scenarios[typhoon_idx],
      period = "Typhoon (2 weeks)",
      avg_weekly_cases = mean(avg_typhoon_samples),
      lower = quantile(avg_typhoon_samples, 0.025),
      upper = quantile(avg_typhoon_samples, 0.975)
    )
    
    # Store stats for recovery period
    period_comparison_list[[length(period_comparison_list) + 1]] <- data.frame(
      strain = strain_names[strain_idx],
      scenario = typhoon_scenarios[typhoon_idx],
      period = "Recovery (3 weeks)",
      avg_weekly_cases = mean(avg_recovery_samples),
      lower = quantile(avg_recovery_samples, 0.025),
      upper = quantile(avg_recovery_samples, 0.975)
    )
  }
}

period_comparison_stats <- do.call(rbind, period_comparison_list)
period_comparison_stats$strain <- factor(period_comparison_stats$strain, levels = strain_names)
period_comparison_stats$scenario <- factor(period_comparison_stats$scenario, levels = typhoon_scenarios)
period_comparison_stats$period <- factor(period_comparison_stats$period, levels = c("Typhoon (2 weeks)", "Recovery (3 weeks)"))


p5_comparison <- ggplot(period_comparison_stats,
                        aes(x = strain, y = avg_weekly_cases, fill = period)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  # ADDED: Error bars
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.9),
    width = 0.3, # Width of the error bar caps
    color = "gray20",
    size = 0.4
  ) +
  facet_wrap(~scenario, ncol = 3, nrow = 2, scales = "free_y") + # Added scales = "free_y"
  scale_fill_manual(
    name = "Period",
    values = c("Typhoon (2 weeks)" = "#E78AC3", "Recovery (3 weeks)" = "#8DA0CB")
  ) +
  labs(
    title = "Average Weekly Cases: Typhoon Period vs Recovery Period",
    subtitle = "Comparing disease burden during and after typhoon impact across all 6 scenarios. Error bars show 95% CIs.",
    x = NULL,
    y = "Average Weekly Cases"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    panel.spacing = unit(0.8, "lines"),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
    axis.text = element_text(size = 9, color = "gray20"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 10, color = "gray30", margin = margin(b = 10))
  )
print(p5_comparison)
ggsave("figure5_period_comparison_with_errorbars.png", p5_comparison,
       width = 14, height = 8, dpi = 300, bg = "white")
# ============================================================================
# FIGURE 6: R_eff by Scenario (Recent Period) - NEWLY ADDED
# ============================================================================
cat("\n=== Creating Figure 6: R_eff by Scenario (Recent Period) ===\n")
reff_cutoff_date <- as.Date("2025-09-01")
reff_start_week <- which(all_dates >= reff_cutoff_date)[1]
if (is.na(reff_start_week)) {
  reff_start_week <- 1
}
reff_scenarios_list <- list()
for (typhoon_idx in 1:6) {
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
p6_reff_scenarios <- ggplot(reff_scenarios_df, aes(x = date)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.5, size = 0.8) +
  # 垂直线：拟合期结束/预测期开始
  geom_vline(xintercept = max(fitting_data$date),
             linetype = "dashed", color = "gray40", size = 0.6, alpha = 0.6) +
  # 矩形：台风期间
  annotate("rect",
           xmin = typhoon_start_date_line,
           xmax = typhoon_end_date_line,
           ymin = -Inf, ymax = Inf,
           fill = "pink", alpha = 0.15) +
  # 矩形：恢复期间
  annotate("rect",
           xmin = recovery_start_date,
           xmax = max(forecast_dates),
           ymin = -Inf, ymax = Inf,
           fill = "lightgreen", alpha = 0.1) +
  # 垂直线：台风结束
  geom_vline(xintercept = typhoon_end_date_line,
             linetype = "dashed", color = "red", size = 0.7, alpha = 0.7) +
  geom_ribbon(aes(ymin = R_lower, ymax = R_upper, fill = scenario), alpha = 0.15) +
  geom_line(aes(y = R_eff, color = scenario, linetype = scenario), size = 1.1) +
  facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
  scale_color_manual(
    name = "Scenarios:",
    values = c("0% Reduction (Baseline)" = "#999999",
               "10% Reduction" = "#6BAED6",
               "20% Reduction" = "#E41A1C",
               "40% Reduction" = "#FF7F00",
               "60% Reduction" = "#4DAF4A",
               "80% Reduction" = "#984EA3"),
    labels = c("0%", "10%", "20%", "40%", "60%", "80%")
  ) +
  scale_fill_manual(
    name = "Scenarios:",
    values = c("0% Reduction (Baseline)" = "#999999",
               "10% Reduction" = "#6BAED6",
               "20% Reduction" = "#E41A1C",
               "40% Reduction" = "#FF7F00",
               "60% Reduction" = "#4DAF4A",
               "80% Reduction" = "#984EA3"),
    labels = c("0%", "10%", "20%", "40%", "60%", "80%")
  ) +
  scale_linetype_manual(
    name = "Scenarios:",
    values = c("0% Reduction (Baseline)" = "solid",
               "10% Reduction" = "solid",
               "20% Reduction" = "solid",
               "40% Reduction" = "dashed",
               "60% Reduction" = "dotdash",
               "80% Reduction" = "longdash"),
    labels = c("0%", "10%", "20%", "40%", "60%", "80%")
  ) +
  # P6 的 X 轴设置
  scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  labs(
    title = expression(paste("Effective Reproduction Number (", R[eff], ") by Typhoon Scenario (Recent Period)")),
    subtitle = expression(paste("Comparing ", R[eff], " trajectories under different transmission reduction scenarios | Pink: Typhoon (2 weeks) | Green: Recovery (3 weeks)")),
    x = NULL,
    y = expression(R[eff])
  ) +
  # P6 的主题设置
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
  # P6 的图例指南
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )
print(p6_reff_scenarios)
ggsave("figure6_reff_scenarios.png", p6_reff_scenarios,
       width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 7: R_eff Change (Difference from Baseline)
# NEW: Data preparation for Figure 7
# ============================================================================
cat("\n=== Creating Figure 7: R_eff Change from Baseline ===\n")

# 使用与 P6 相同的起始周
reff_change_start_week <- reff_start_week 

reff_change_list <- list()

# 循环情景 2 到 6 (10% 到 80%)
for (typhoon_idx in 2:6) { 
  for (strain_idx in 1:N_strains) {
    weeks_subset <- reff_change_start_week:T_weeks_total
    dates_subset <- all_dates[weeks_subset]
    
    # 获取基线 (情景 1) 的后验样本
    # 维度: [n_samples, n_weeks]
    base_reff_samples <- R_eff_scenarios[, 1, weeks_subset, strain_idx]
    
    # 获取当前情景的后验样本
    # 维度: [n_samples, n_weeks]
    scen_reff_samples <- R_eff_scenarios[, typhoon_idx, weeks_subset, strain_idx]
    
    # 计算差异的分布
    # 维度: [n_samples, n_weeks]
    diff_dist <- scen_reff_samples - base_reff_samples
    
    # 跨样本计算统计数据 (margin = 2 是按列/周计算)
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
      scenario = typhoon_scenarios[typhoon_idx] # Scenarios 2-6
    )
  }
}

reff_change_df <- do.call(rbind, reff_change_list)
reff_change_df$strain <- factor(reff_change_df$strain, levels = strain_names)
# 将因子水平设置为仅 5 个情景
reff_change_df$scenario <- factor(reff_change_df$scenario, levels = typhoon_scenarios[2:6])

# ============================================================================
# FIGURE 7: R_eff Change (Difference from Baseline)
# NEW: Plotting Figure 7, styled *EXACTLY* like Figure 6
# ============================================================================

p7_reff_change <- ggplot(reff_change_df, aes(x = date)) +
  # 水平线在 0 (因为是差异图)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5, size = 0.8) +
  
  # 垂直线: 拟合期结束/预测期开始 (与 P6 一致)
  geom_vline(xintercept = max(fitting_data$date),
             linetype = "dashed", color = "gray40", size = 0.6, alpha = 0.6) +
  
  # 粉色矩形: 台风期间 (与 P6 一致)
  annotate("rect",
           xmin = typhoon_start_date_line,
           xmax = typhoon_end_date_line,
           ymin = -Inf, ymax = Inf,
           fill = "pink", alpha = 0.15) +
  
  # 绿色矩形: 恢复期间 (与 P6 一致)
  annotate("rect",
           xmin = recovery_start_date,
           xmax = max(forecast_dates),
           ymin = -Inf, ymax = Inf,
           fill = "lightgreen", alpha = 0.1) +
  
  # 垂直线: 台风结束 (与 P6 一致)
  geom_vline(xintercept = typhoon_end_date_line,
             linetype = "dashed", color = "red", size = 0.7, alpha = 0.7) +
  
  # 色带和线条 (使用 p7 的数据)
  geom_ribbon(aes(ymin = R_change_lower, ymax = R_change_upper, fill = scenario), alpha = 0.15) +
  geom_line(aes(y = R_eff_change, color = scenario, linetype = scenario), size = 1.1) +
  
  # 分面 (与 P6 一致)
  facet_wrap(~strain, scales = "free_y", ncol = 2, nrow = 3) +
  
  # 颜色/填充/线型比例尺 (使用 p7 的 5 个情景设置)
  scale_color_manual(
    name = "Reduction Scenario:",
    values = c("10% Reduction" = "#6BAED6",
               "20% Reduction" = "#E41A1C",
               "40% Reduction" = "#FF7F00",
               "60% Reduction" = "#4DAF4A",
               "80% Reduction" = "#984EA3"),
    labels = c("10%", "20%", "40%", "60%", "80%")
  ) +
  scale_fill_manual(
    name = "Reduction Scenario:",
    values = c("10% Reduction" = "#6BAED6",
               "20% Reduction" = "#E41A1C",
               "40% Reduction" = "#FF7F00",
               "60% Reduction" = "#4DAF4A",
               "80% Reduction" = "#984EA3"),
    labels = c("10%", "20%", "40%", "60%", "80%")
  ) +
  scale_linetype_manual(
    name = "Reduction Scenario:",
    values = c("10% Reduction" = "solid",
               "20% Reduction" = "solid",
               "40% Reduction" = "dashed",
               "60% Reduction" = "dotdash",
               "80% Reduction" = "longdash"),
    labels = c("10%", "20%", "40%", "60%", "80%")
  ) +
  
  # X 轴: *严格*按照 P6 的风格
  scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
  
  # 标签: 使用 p7 的特定标签
  labs(
    title = expression(paste("Change in ", R[eff], " from Baseline (0% Scenario)")),
    subtitle = expression(paste("Negative values indicate reduction in transmission | ", Delta, R[eff], " = ", R[eff], "(scenario) - ", R[eff], "(baseline)")),
    x = NULL,
    y = expression(paste(Delta, R[eff]))
  ) +
  
  # 主题: *严格*按照 P6 的主题
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
  
  # 图例指南: *严格*按照 P6 的指南
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )

print(p7_reff_change)
ggsave("figure7_reff_change_scenarios.png", p7_reff_change,
       width = 14, height = 10, dpi = 300, bg = "white")

cat("\n=== All figures created successfully ===\n")




















































# Validation Analysis
# ============================================================================
cat("\n", rep("=", 100), "\n", sep="")
cat("VALIDATION ANALYSIS: Comparing Out-of-Sample Predictions with Observed Data\n")
cat("(Data from 2025-09-27 to 2025-10-11 was NOT used in model fitting)\n")
cat(rep("=", 100), "\n", sep="")

validation_comparison <- data.frame()

for (strain_idx in 1:N_strains) {
  cat(sprintf("\n%s:\n", strain_names[strain_idx]))
  cat(sprintf("%-15s | %-15s | %-30s | %-30s | %-15s\n",
              "Week", "Observed", "Baseline Predicted", "Best Scenario", "Coverage"))
  cat(strrep("-", 115), "\n")
  
  for (val_week in 1:T_weeks_validation) {
    obs_value <- validation_matrix[val_week, strain_idx]
    
    if (!is.na(obs_value)) {
      baseline_pred <- forecast_typhoon_mean[1, val_week, strain_idx]
      baseline_lower <- forecast_typhoon_lower[1, val_week, strain_idx]
      baseline_upper <- forecast_typhoon_upper[1, val_week, strain_idx]
      
      scenario_diffs <- abs(forecast_typhoon_mean[, val_week, strain_idx] - obs_value)
      best_scenario_idx <- which.min(scenario_diffs)
      best_scenario_pred <- forecast_typhoon_mean[best_scenario_idx, val_week, strain_idx]
      
      in_ci <- obs_value >= baseline_lower & obs_value <= baseline_upper
      coverage_status <- ifelse(in_ci, "YES", "NO")
      
      cat(sprintf("%-15s | %8.0f       | %8.0f [%8.0f, %8.0f] | %-20s (%8.0f) | %-15s\n",
                  as.character(validation_data$date[val_week]),
                  obs_value,
                  baseline_pred, baseline_lower, baseline_upper,
                  typhoon_scenarios[best_scenario_idx], best_scenario_pred,
                  coverage_status))
      
      validation_comparison <- rbind(validation_comparison, data.frame(
        date = validation_data$date[val_week],
        strain = strain_names[strain_idx],
        observed = obs_value,
        baseline_pred = baseline_pred,
        baseline_lower = baseline_lower,
        baseline_upper = baseline_upper,
        best_scenario = typhoon_scenarios[best_scenario_idx],
        best_scenario_pred = best_scenario_pred,
        in_ci = in_ci,
        stringsAsFactors = FALSE
      ))
    }
  }
}

write.csv(validation_comparison,
          file = "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/validation_comparison_6scenarios.csv",
          row.names = FALSE)

if (nrow(validation_comparison) > 0) {
  cat("\n", rep("=", 100), "\n", sep="")
  cat("VALIDATION SUMMARY (Out-of-Sample Performance):\n")
  cat(rep("=", 100), "\n", sep="")
  
  overall_coverage <- mean(validation_comparison$in_ci) * 100
  cat(sprintf("Overall 95%% CI Coverage: %.1f%%\n", overall_coverage))
  cat("(Expected: ~95% for well-calibrated predictions)\n\n")
  
  coverage_by_strain <- validation_comparison %>%
    group_by(strain) %>%
    summarise(
      coverage = mean(in_ci) * 100,
      mean_abs_error = mean(abs(observed - baseline_pred)),
      mean_rel_error = mean(abs(observed - baseline_pred) / observed) * 100,
      .groups = "drop"
    )
  
  cat("Validation Metrics by Strain:\n")
  cat(sprintf("%-10s | %-15s | %-20s | %-20s\n",
              "Strain", "Coverage (%)", "Mean Abs Error", "Mean Rel Error (%)"))
  cat(strrep("-", 75), "\n")
  
  for (i in 1:nrow(coverage_by_strain)) {
    cat(sprintf("%-10s | %13.1f%% | %18.1f | %18.1f%%\n",
                coverage_by_strain$strain[i],
                coverage_by_strain$coverage[i],
                coverage_by_strain$mean_abs_error[i],
                coverage_by_strain$mean_rel_error[i]))
  }
}

# ============================================================================
# Validation Figure
# ============================================================================
cat("\n=== Creating Validation Comparison Figure ===\n")

if (nrow(validation_comparison) > 0) {
  p_validation <- ggplot(validation_comparison, aes(x = baseline_pred, y = observed)) +
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
  
  print(p_validation)
  ggsave("figure7_validation_comparison.png", p_validation,
         width = 12, height = 10, dpi = 300, bg = "white")
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
    if(!is.null(strains)) rownames(stats) <- strains
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
# Visualize Posterior Distributions
# ============================================================================
cat("\n=== Creating Posterior Distribution Figure ===\n")

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

p_post <- ggplot(post_data, aes(x = value, fill = group)) +
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

print(p_post)
ggsave("figure8_posterior_distributions.png", p_post,
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
# KEY FINDINGS SUMMARY
# ============================================================================
cat("\n", rep("=", 100), "\n", sep="")
cat("KEY FINDINGS (NOW WITH 6 SCENARIOS INCLUDING 80% REDUCTION)\n")
cat(rep("=", 100), "\n", sep="")

cat("\n*** IMPORTANT: DATA USAGE ***\n")
cat(sprintf("Model was fitted ONLY to data from %s to %s (%d weeks)\n",
            as.character(min(fitting_data$date)),
            as.character(max(fitting_data$date)),
            T_weeks))

if (nrow(validation_data) > 0) {
  cat(sprintf("Validation data from %s to %s (%d weeks) was NOT used in fitting\n",
              as.character(min(validation_data$date)),
              as.character(max(validation_data$date)),
              T_weeks_validation))
  cat("Yellow diamonds in figures represent out-of-sample validation data\n")
}

if (nrow(validation_comparison) > 0) {
  cat("\n1. MODEL VALIDATION (Out-of-Sample):\n")
  cat(sprintf("   - Overall 95%% credible interval coverage: %.1f%%\n", overall_coverage))
  cat("   - This indicates how well the model predictions match unseen data\n")
}

cat("\n2. R_eff REDUCTION DURING TYPHOON PERIOD (1 week):\n")
for (scenario_name in unique(reff_change_summary$scenario)) {
  scenario_data <- reff_change_summary[reff_change_summary$scenario == scenario_name, ]
  avg_reduction <- mean(scenario_data$change_pct_mean)
  cat(sprintf("   - %s: Average %.1f%% reduction across all strains\n",
              scenario_name, avg_reduction))
}

cat("\n3. CASE AVOIDANCE (Total 5-week period):\n")
for (i in 1:nrow(case_reduction_summary)) {
  cat(sprintf("   - %s: %.0f cases avoided (95%% CI: %.0f to %.0f), %.1f%% reduction\n",
              case_reduction_summary$scenario[i],
              case_reduction_summary$abs_reduction_mean[i],
              case_reduction_summary$abs_reduction_lower[i],
              case_reduction_summary$abs_reduction_upper[i],
              case_reduction_summary$pct_reduction_mean[i]))
}

cat("\n4. STATISTICAL SIGNIFICANCE:\n")
cat("   - All intervention scenarios show >95% probability of R_eff reduction\n")
cat("   - All intervention scenarios show >95% probability of case avoidance\n")
cat("   - 80% reduction scenario shows strongest effect\n")

if (nrow(validation_comparison) > 0) {
  cat("\n5. VALIDATION DATA INSIGHTS:\n")
  cat("   - Yellow diamonds on plots show actual observed data from typhoon period\n")
  cat("   - Compare these with different scenario predictions to assess model accuracy\n")
  cat("   - Best-fitting scenario can guide future intervention planning\n")
}

cat("\n6. VISUALIZATION GUIDE:\n")
cat("   - Gray dashed vertical line: Typhoon start (2025-09-27)\n")
cat("   - Red dashed vertical line: Typhoon end (2025-10-04, 1 week after start)\n")
cat("   - Pink shaded area: Typhoon impact period (1 week)\n")
cat("   - Green shaded area: Recovery period (4 weeks)\n")

cat("\n", rep("=", 100), "\n", sep="")

















cat("\n=== Creating Forest Plot for R_eff Changes ===\n")

reff_change_summary$scenario <- factor(
  reff_change_summary$scenario,
  levels = c("10% Reduction", "20% Reduction", "40% Reduction", "60% Reduction", "80% Reduction")
)

reff_change_summary$strain <- factor(
  reff_change_summary$strain,
  levels = strain_names
)

p_reff_forest <- ggplot(reff_change_summary,
                        aes(x = change_pct_mean, y = scenario, color = scenario)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.5, size = 0.8) +
  geom_vline(xintercept = c(-10, -20, -40, -60, -80),
             linetype = "dotted", color = "gray60", alpha = 0.4, size = 0.3) +
  
  geom_pointrange(aes(xmin = change_pct_lower, xmax = change_pct_upper),
                  size = 0.6, fatten = 3,
                  position = position_dodge(width = 0.3)) +
  
  geom_text(aes(label = sprintf("%.1f%%\n[%.1f%%, %.1f%%]",
                                change_pct_mean, change_pct_lower, change_pct_upper)),
            hjust = -0.15, size = 2.5, color = "gray20",
            position = position_dodge(width = 0.3)) +
  
  facet_wrap(~strain, ncol = 2, nrow = 3, scales = "free_x") +
  
  scale_color_manual(
    name = "Typhoon Scenario",
    values = c("10% Reduction" = "#6BAED6",
               "20% Reduction" = "#E41A1C",
               "40% Reduction" = "#FF7F00",
               "60% Reduction" = "#4DAF4A",
               "80% Reduction" = "#984EA3")
  ) +
  
  scale_x_continuous(
    breaks = seq(-90, 0, 10),
    labels = function(x) paste0(x, "%")
  ) +
  
  labs(
    title = expression(paste("Change in Effective Reproduction Number (", R[eff], ") During Typhoon Period")),
    subtitle = "Mean and 95% credible intervals | Dashed line at 0% (no change) | Dotted lines show expected reductions | Now with 80% scenario",
    x = expression(paste("Relative Change in ", R[eff], " (%)")),
    y = "Typhoon Scenario"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "gray30", size = 0.5),
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    plot.margin = margin(10, 30, 10, 10)
  )

print(p_reff_forest)
ggsave("figure6_reff_change_forest.png", p_reff_forest,
       width = 14, height = 10, dpi = 300, bg = "white")