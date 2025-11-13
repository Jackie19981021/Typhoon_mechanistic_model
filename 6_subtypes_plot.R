# ============================================================================
# TYPHOON IMPACT ANALYSIS - EXTENDED WITH SUNBURST VISUALIZATION
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

# [Previous data loading code - UNCHANGED]
cat("=== Loading and preprocessing data ===\n")
data_path <- "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/HK_ILI_COVID_Sep.csv"
raw_data <- read.csv(data_path)
raw_data$date <- as.Date(raw_data$date, format = "%Y/%m/%d")

start_date <- as.Date("2023-01-07")
typhoon_start_date <- as.Date("2025-09-20")
end_date <- as.Date("2025-10-25")

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

# [Stan data preparation - UNCHANGED]
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
  typhoon_weeks = typhoon_weeks
)

# [Model compilation and fitting - UNCHANGED]
cat("=== Compiling Stan model ===\n")
model_file <- "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/6_subtypes.stan"
model <- stan_model(model_file)

cat("\n=== Starting MCMC sampling ===\n")
cat("Note: Now computing 61 scenarios (1 baseline + 60 combinations)\n\n")

# Uncomment to run new fit:
 fit <- sampling(
   model,
   data = stan_data,
   iter = 2000,
   warmup = 1000,
   chains = 4,
   thin = 1,
   cores = 4,
   control = list(adapt_delta = 0.95, max_treedepth = 20)
 )
saveRDS(fit, file = "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/typhoon_model_fit_extended.rds")

cat("=== LOADING PREVIOUS FIT ===\n")
#fit <- readRDS("/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/typhoon_model_fit_extended.rds")

# [Extract results - ORIGINAL 20 scenarios UNCHANGED]
pred_cases <- rstan::extract(fit, pars = "pred_cases")$pred_cases
forecast_cases <- rstan::extract(fit, pars = "forecast_cases")$forecast_cases
forecast_cases_typhoon <- rstan::extract(fit, pars = "forecast_cases_typhoon")$forecast_cases_typhoon
R_eff <- rstan::extract(fit, pars = "R_eff")$R_eff
R_eff_scenarios <- rstan::extract(fit, pars = "R_eff_scenarios")$R_eff_scenarios
R0_t <- rstan::extract(fit, pars = "R0_t")$R0_t

reduction_typhoon <- rstan::extract(fit, pars = "reduction_typhoon_period")$reduction_typhoon_period
avg_weekly_typhoon <- rstan::extract(fit, pars = "avg_weekly_typhoon")$avg_weekly_typhoon
avg_weekly_recovery <- rstan::extract(fit, pars = "avg_weekly_recovery")$avg_weekly_recovery

# NEW: Extract extended scenarios
# NEW: Extract extended scenarios
cases_extended <- rstan::extract(fit, pars = "cases_extended_typhoon")$cases_extended_typhoon
incidence_per10k_extended <- rstan::extract(fit, pars = "incidence_per10k_extended")$incidence_per10k_extended
reduction_extended <- rstan::extract(fit, pars = "reduction_extended_typhoon")$reduction_extended_typhoon
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

# [ALL PREVIOUS VISUALIZATIONS - COMPLETELY UNCHANGED]
# Calculate statistics
pred_mean <- apply(pred_cases, c(2,3), mean)
pred_median <- apply(pred_cases, c(2,3), median)
pred_lower <- apply(pred_cases, c(2,3), quantile, probs = 0.025)
pred_upper <- apply(pred_cases, c(2,3), quantile, probs = 0.975)

forecast_mean <- apply(forecast_cases, c(2,3), mean)
forecast_median <- apply(forecast_cases, c(2,3), median)
forecast_lower <- apply(forecast_cases, c(2,3), quantile, probs = 0.025)
forecast_upper <- apply(forecast_cases, c(2,3), quantile, probs = 0.975)

# 20 typhoon scenarios
typhoon_scenarios <- c("0% (Baseline)",
                       "3 days (10%)", "3 days (20%)", "3 days (40%)", "3 days (60%)", "3 days (80%)",
                       "5 days (10%)", "5 days (20%)", "5 days (40%)", "5 days (60%)", "5 days (80%)",
                       "7 days (10%)", "7 days (20%)", "7 days (40%)", "7 days (60%)", "7 days (80%)",
                       "7 days (60%, early 3 days)", "7 days (60%, early 5 days)",
                       "7 days (60%, delayed 3 days)", "7 days (60%, delayed 5 days)")

forecast_typhoon_mean <- array(NA, dim = c(20, T_weeks_forecast, N_strains))
forecast_typhoon_lower <- array(NA, dim = c(20, T_weeks_forecast, N_strains))
forecast_typhoon_upper <- array(NA, dim = c(20, T_weeks_forecast, N_strains))

for (typhoon_idx in 1:20) {
  forecast_typhoon_mean[typhoon_idx,,] <- apply(forecast_cases_typhoon[,typhoon_idx,,], c(2,3), mean)
  forecast_typhoon_lower[typhoon_idx,,] <- apply(forecast_cases_typhoon[,typhoon_idx,,], c(2,3), quantile, probs = 0.025)
  forecast_typhoon_upper[typhoon_idx,,] <- apply(forecast_cases_typhoon[,typhoon_idx,,], c(2,3), quantile, probs = 0.975)
}

R_eff_scenarios_mean <- array(NA, dim = c(20, T_weeks + T_weeks_forecast, N_strains))
R_eff_scenarios_lower <- array(NA, dim = c(20, T_weeks + T_weeks_forecast, N_strains))
R_eff_scenarios_upper <- array(NA, dim = c(20, T_weeks + T_weeks_forecast, N_strains))

for (typhoon_idx in 1:20) {
  R_eff_scenarios_mean[typhoon_idx,,] <- apply(R_eff_scenarios[,typhoon_idx,,], c(2,3), mean)
  R_eff_scenarios_lower[typhoon_idx,,] <- apply(R_eff_scenarios[,typhoon_idx,,], c(2,3), quantile, probs = 0.025)
  R_eff_scenarios_upper[typhoon_idx,,] <- apply(R_eff_scenarios[,typhoon_idx,,], c(2,3), quantile, probs = 0.975)
}

# 5. Prepare visualization data
theme_set(theme_minimal(base_size = 11))
forecast_dates <- seq(typhoon_start_date, by = "week", length.out = T_weeks_forecast)
all_dates <- c(fitting_data$date, forecast_dates)

typhoon_start_date_line <- typhoon_start_date
typhoon_end_date_line <- typhoon_start_date + 7
typhoon_period_end <- typhoon_start_date + 7
recovery_start_date <- typhoon_start_date + 7

cat("\n=== Date markers ===\n")
cat("Typhoon start (gray line):", as.character(typhoon_start_date_line), "\n")
cat("Typhoon end (red line):", as.character(typhoon_end_date_line), "\n")
cat("Typhoon duration: 1 week (7 days)\n")
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
# FIGURE 2: Recent Period + All Typhoon Scenarios + Validation Data - Separate for each duration
# ============================================================================
cat("\n=== Creating Figure 2: Recent Period, Scenarios, and Validation - Separate for each duration ===\n")

cutoff_date <- as.Date("2025-09-01")

recent_fit_df <- results_df %>%
  filter(date >= cutoff_date) %>%
  mutate(scenario = "Historical Fit")

typhoon_forecast_list <- list()
for (typhoon_idx in 1:20) {
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
  levels = c("Historical Fit", typhoon_scenarios)
)

last_observed <- recent_fit_df %>%
  group_by(strain) %>%
  slice_tail(n = 1) %>%
  ungroup()

# Function to create plot for selected scenarios
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
    scale_color_manual(
      name = "Scenarios:",
      values = colors
    ) +
    scale_fill_manual(
      name = "Scenarios:",
      values = colors
    ) +
    scale_linetype_manual(
      name = "Scenarios:",
      values = linetypes
    ) +
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

# Define colors and linetypes for 3 days
colors_3 <- c("0% (Baseline)" = "#999999",
              "3 days (10%)" = "#AED6F1",
              "3 days (20%)" = "#5DADE2",
              "3 days (40%)" = "#3498DB",
              "3 days (60%)" = "#21618C",
              "3 days (80%)" = "#1B4F72")

linetypes_3 <- c("0% (Baseline)" = "solid",
                 "3 days (10%)" = "solid",
                 "3 days (20%)" = "dashed",
                 "3 days (40%)" = "dotted",
                 "3 days (60%)" = "dotdash",
                 "3 days (80%)" = "longdash")

selected_3 <- c("Historical Fit", "0% (Baseline)", "3 days (10%)", "3 days (20%)", "3 days (40%)", "3 days (60%)", "3 days (80%)")

p2_3days <- create_scenario_plot(selected_3, "Typhoon Impact: Recent Period and 3 Days Duration Scenarios", colors_3, linetypes_3)
print(p2_3days)
ggsave("figure2_typhoon_scenarios_3days.png", p2_3days, width = 14, height = 10, dpi = 300, bg = "white")

# Define colors and linetypes for 5 days
colors_5 <- c("0% (Baseline)" = "#999999",
              "5 days (10%)" = "#AED6F1",
              "5 days (20%)" = "#5DADE2",
              "5 days (40%)" = "#3498DB",
              "5 days (60%)" = "#21618C",
              "5 days (80%)" = "#1B4F72")

linetypes_5 <- c("0% (Baseline)" = "solid",
                 "5 days (10%)" = "solid",
                 "5 days (20%)" = "dashed",
                 "5 days (40%)" = "dotted",
                 "5 days (60%)" = "dotdash",
                 "5 days (80%)" = "longdash")

selected_5 <- c("Historical Fit", "0% (Baseline)", "5 days (10%)", "5 days (20%)", "5 days (40%)", "5 days (60%)", "5 days (80%)")

p2_5days <- create_scenario_plot(selected_5, "Typhoon Impact: Recent Period and 5 Days Duration Scenarios", colors_5, linetypes_5)
print(p2_5days)
ggsave("figure2_typhoon_scenarios_5days.png", p2_5days, width = 14, height = 10, dpi = 300, bg = "white")

# Define colors and linetypes for 7 days
colors_7 <- c("0% (Baseline)" = "#999999",
              "7 days (10%)" = "#AED6F1",
              "7 days (20%)" = "#5DADE2",
              "7 days (40%)" = "#3498DB",
              "7 days (60%)" = "#21618C",
              "7 days (80%)" = "#1B4F72")

linetypes_7 <- c("0% (Baseline)" = "solid",
                 "7 days (10%)" = "solid",
                 "7 days (20%)" = "dashed",
                 "7 days (40%)" = "dotted",
                 "7 days (60%)" = "dotdash",
                 "7 days (80%)" = "longdash")

selected_7 <- c("Historical Fit", "0% (Baseline)", "7 days (10%)", "7 days (20%)", "7 days (40%)", "7 days (60%)", "7 days (80%)")

p2_7days <- create_scenario_plot(selected_7, "Typhoon Impact: Recent Period and 7 Days Duration Scenarios", colors_7, linetypes_7)
print(p2_7days)
ggsave("figure2_typhoon_scenarios_7days.png", p2_7days, width = 14, height = 10, dpi = 300, bg = "white")

# Define colors and linetypes for shifted scenarios
selected_shifted <- c("Historical Fit", "0% (Baseline)", "7 days (60%, early 3 days)", "7 days (60%, early 5 days)", "7 days (60%, delayed 3 days)", "7 days (60%, delayed 5 days)")

colors_shifted <- c("0% (Baseline)" = "#999999",
                    "7 days (60%, early 3 days)" = "#AED6F1",
                    "7 days (60%, early 5 days)" = "#5DADE2",
                    "7 days (60%, delayed 3 days)" = "#3498DB",
                    "7 days (60%, delayed 5 days)" = "#21618C")

linetypes_shifted <- c("0% (Baseline)" = "solid",
                       "7 days (60%, early 3 days)" = "solid",
                       "7 days (60%, early 5 days)" = "dashed",
                       "7 days (60%, delayed 3 days)" = "dotted",
                       "7 days (60%, delayed 5 days)" = "dotdash")

p2_shifted <- create_scenario_plot(selected_shifted, "Typhoon Impact: Recent Period and Shifted Timing Scenarios (7 days 60%)", colors_shifted, linetypes_shifted)
print(p2_shifted)
ggsave("figure2_typhoon_scenarios_shifted.png", p2_shifted, width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 3: Baseline vs Typhoon Forecasts - 3 days duration
# ============================================================================
cat("\n=== Creating Figure 3: Baseline vs Typhoon Forecasts - 3 days duration ===\n")

forecast_3days <- typhoon_forecast_df %>%
  filter(scenario %in% c("0% (Baseline)", "3 days (10%)", "3 days (20%)", "3 days (40%)", "3 days (60%)", "3 days (80%)"))

forecast_3days$scenario <- factor(forecast_3days$scenario, levels = c("0% (Baseline)", "3 days (10%)", "3 days (20%)", "3 days (40%)", "3 days (60%)", "3 days (80%)"))

colors_3_forecast <- c("0% (Baseline)" = "#999999",
                       "3 days (10%)" = "#AED6F1",
                       "3 days (20%)" = "#5DADE2",
                       "3 days (40%)" = "#3498DB",
                       "3 days (60%)" = "#21618C",
                       "3 days (80%)" = "#1B4F72")

linetypes_3_forecast <- c("0% (Baseline)" = "solid",
                          "3 days (10%)" = "solid",
                          "3 days (20%)" = "dashed",
                          "3 days (40%)" = "dotted",
                          "3 days (60%)" = "dotdash",
                          "3 days (80%)" = "longdash")

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
  scale_color_manual(values = colors_3_forecast) +
  scale_fill_manual(values = colors_3_forecast) +
  scale_linetype_manual(values = linetypes_3_forecast) +
  scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
  labs(
    title = "Typhoon Impact Forecasts: 3 Days Duration Scenarios",
    subtitle = "Yellow diamonds: Validation data | Gray dashed: Typhoon start | Red dashed: Typhoon end\nPink: Typhoon period | Green: Recovery period",
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

# ============================================================================
# FIGURE 4: Baseline vs Typhoon Forecasts - 5 days duration
# ============================================================================
cat("\n=== Creating Figure 4: Baseline vs Typhoon Forecasts - 5 days duration ===\n")

forecast_5days <- typhoon_forecast_df %>%
  filter(scenario %in% c("0% (Baseline)", "5 days (10%)", "5 days (20%)", "5 days (40%)", "5 days (60%)", "5 days (80%)"))

forecast_5days$scenario <- factor(forecast_5days$scenario, levels = c("0% (Baseline)", "5 days (10%)", "5 days (20%)", "5 days (40%)", "5 days (60%)", "5 days (80%)"))

colors_5_forecast <- c("0% (Baseline)" = "#999999",
                       "5 days (10%)" = "#AED6F1",
                       "5 days (20%)" = "#5DADE2",
                       "5 days (40%)" = "#3498DB",
                       "5 days (60%)" = "#21618C",
                       "5 days (80%)" = "#1B4F72")

linetypes_5_forecast <- c("0% (Baseline)" = "solid",
                          "5 days (10%)" = "solid",
                          "5 days (20%)" = "dashed",
                          "5 days (40%)" = "dotted",
                          "5 days (60%)" = "dotdash",
                          "5 days (80%)" = "longdash")

p4_5days <- ggplot(forecast_5days, aes(x = date, y = predicted, color = scenario, fill = scenario, linetype = scenario)) +
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
  scale_color_manual(values = colors_5_forecast) +
  scale_fill_manual(values = colors_5_forecast) +
  scale_linetype_manual(values = linetypes_5_forecast) +
  scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
  labs(
    title = "Typhoon Impact Forecasts: 5 Days Duration Scenarios",
    subtitle = "Yellow diamonds: Validation data | Gray dashed: Typhoon start | Red dashed: Typhoon end\nPink: Typhoon period | Green: Recovery period",
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

print(p4_5days)
ggsave("figure4_typhoon_forecasts_5days.png", p4_5days, width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 5: Baseline vs Typhoon Forecasts - 7 days duration
# ============================================================================
cat("\n=== Creating Figure 5: Baseline vs Typhoon Forecasts - 7 days duration ===\n")

forecast_7days <- typhoon_forecast_df %>%
  filter(scenario %in% c("0% (Baseline)", "7 days (10%)", "7 days (20%)", "7 days (40%)", "7 days (60%)", "7 days (80%)"))

forecast_7days$scenario <- factor(forecast_7days$scenario, levels = c("0% (Baseline)", "7 days (10%)", "7 days (20%)", "7 days (40%)", "7 days (60%)", "7 days (80%)"))

colors_7_forecast <- c("0% (Baseline)" = "#999999",
                       "7 days (10%)" = "#AED6F1",
                       "7 days (20%)" = "#5DADE2",
                       "7 days (40%)" = "#3498DB",
                       "7 days (60%)" = "#21618C",
                       "7 days (80%)" = "#1B4F72")

linetypes_7_forecast <- c("0% (Baseline)" = "solid",
                          "7 days (10%)" = "solid",
                          "7 days (20%)" = "dashed",
                          "7 days (40%)" = "dotted",
                          "7 days (60%)" = "dotdash",
                          "7 days (80%)" = "longdash")

p5_7days <- ggplot(forecast_7days, aes(x = date, y = predicted, color = scenario, fill = scenario, linetype = scenario)) +
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
  scale_color_manual(values = colors_7_forecast) +
  scale_fill_manual(values = colors_7_forecast) +
  scale_linetype_manual(values = linetypes_7_forecast) +
  scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
  labs(
    title = "Typhoon Impact Forecasts: 7 Days Duration Scenarios",
    subtitle = "Yellow diamonds: Validation data | Gray dashed: Typhoon start | Red dashed: Typhoon end\nPink: Typhoon period | Green: Recovery period",
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

print(p5_7days)
ggsave("figure5_typhoon_forecasts_7days.png", p5_7days, width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 6: Baseline vs Typhoon Forecasts - Shifted scenarios
# ============================================================================
cat("\n=== Creating Figure 6: Baseline vs Typhoon Forecasts - Shifted scenarios ===\n")

forecast_shifted <- typhoon_forecast_df %>%
  filter(scenario %in% c("0% (Baseline)", "7 days (60%, early 3 days)", "7 days (60%, early 5 days)", "7 days (60%, delayed 3 days)", "7 days (60%, delayed 5 days)"))

forecast_shifted$scenario <- factor(forecast_shifted$scenario, levels = c("0% (Baseline)", "7 days (60%, early 3 days)", "7 days (60%, early 5 days)", "7 days (60%, delayed 3 days)", "7 days (60%, delayed 5 days)"))

colors_shifted_forecast <- c("0% (Baseline)" = "#999999",
                             "7 days (60%, early 3 days)" = "#AED6F1",
                             "7 days (60%, early 5 days)" = "#5DADE2",
                             "7 days (60%, delayed 3 days)" = "#3498DB",
                             "7 days (60%, delayed 5 days)" = "#21618C")

linetypes_shifted_forecast <- c("0% (Baseline)" = "solid",
                                "7 days (60%, early 3 days)" = "solid",
                                "7 days (60%, early 5 days)" = "dashed",
                                "7 days (60%, delayed 3 days)" = "dotted",
                                "7 days (60%, delayed 5 days)" = "dotdash")

p6_shifted <- ggplot(forecast_shifted, aes(x = date, y = predicted, color = scenario, fill = scenario, linetype = scenario)) +
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
  scale_color_manual(values = colors_shifted_forecast) +
  scale_fill_manual(values = colors_shifted_forecast) +
  scale_linetype_manual(values = linetypes_shifted_forecast) +
  scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
  labs(
    title = "Typhoon Impact Forecasts: Shifted Timing Scenarios (7 days 60%)",
    subtitle = "Yellow diamonds: Validation data | Gray dashed: Typhoon start | Red dashed: Typhoon end\nPink: Typhoon period | Green: Recovery period",
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

print(p6_shifted)
ggsave("figure6_typhoon_forecasts_shifted.png", p6_shifted, width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 7: Effective Reproduction Number (Rt) - BASELINE ONLY
# ============================================================================
cat("\n=== Creating Figure 7: Effective Reproduction Number (Baseline) ===\n")

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
    subtitle = expression(paste("Solid: Historical | Dashed: Forecast | Red line: ", R[eff], "=1 threshold\nGray dashed: Typhoon start | Red dashed: Typhoon end")),
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

# ============================================================================
# FIGURE 8: Effective Reproduction Number (Rt) by Scenario (Recent Period)
# ============================================================================
cat("\n=== Creating Figure 8: R_eff by Scenario (Recent Period) - Separate for each duration ===\n")

reff_start_week <- T_weeks - 3

reff_scenarios_list <- list()
for (typhoon_idx in 1:20) {
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
    scale_color_manual(
      name = "Scenarios:",
      values = colors
    ) +
    scale_fill_manual(
      name = "Scenarios:",
      values = colors
    ) +
    scale_linetype_manual(
      name = "Scenarios:",
      values = linetypes
    ) +
    scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
    scale_y_continuous(breaks = seq(0, 10, 1)) +
    coord_cartesian(ylim = c(0, 5)) +
    labs(
      title = plot_title,
      subtitle = "Comparing R_eff trajectories under different transmission reduction scenarios\nPink: Typhoon | Green: Recovery",
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

# Define colors and linetypes for 3 days reff
colors_3_reff <- colors_3
linetypes_3_reff <- linetypes_3
selected_3_reff <- c("0% (Baseline)", "3 days (10%)", "3 days (20%)", "3 days (40%)", "3 days (60%)", "3 days (80%)")

p8_3days <- create_reff_plot(selected_3_reff, "Effective Reproduction Number (R_eff) by Typhoon Scenario - 3 Days Duration", colors_3_reff, linetypes_3_reff)
print(p8_3days)
ggsave("figure8_reff_scenarios_3days.png", p8_3days, width = 14, height = 10, dpi = 300, bg = "white")

# Define colors and linetypes for 5 days reff
colors_5_reff <- colors_5
linetypes_5_reff <- linetypes_5
selected_5_reff <- c("0% (Baseline)", "5 days (10%)", "5 days (20%)", "5 days (40%)", "5 days (60%)", "5 days (80%)")

p8_5days <- create_reff_plot(selected_5_reff, "Effective Reproduction Number (R_eff) by Typhoon Scenario - 5 Days Duration", colors_5_reff, linetypes_5_reff)
print(p8_5days)
ggsave("figure8_reff_scenarios_5days.png", p8_5days, width = 14, height = 10, dpi = 300, bg = "white")

# Define colors and linetypes for 7 days reff
colors_7_reff <- colors_7
linetypes_7_reff <- linetypes_7
selected_7_reff <- c("0% (Baseline)", "7 days (10%)", "7 days (20%)", "7 days (40%)", "7 days (60%)", "7 days (80%)")

p8_7days <- create_reff_plot(selected_7_reff, "Effective Reproduction Number (R_eff) by Typhoon Scenario - 7 Days Duration", colors_7_reff, linetypes_7_reff)
print(p8_7days)
ggsave("figure8_reff_scenarios_7days.png", p8_7days, width = 14, height = 10, dpi = 300, bg = "white")

# Define colors and linetypes for shifted reff
colors_shifted_reff <- colors_shifted
linetypes_shifted_reff <- linetypes_shifted
selected_shifted_reff <- c("0% (Baseline)", "7 days (60%, early 3 days)", "7 days (60%, early 5 days)", "7 days (60%, delayed 3 days)", "7 days (60%, delayed 5 days)")

p8_shifted <- create_reff_plot(selected_shifted_reff, "Effective Reproduction Number (R_eff) by Typhoon Scenario - Shifted Timing Scenarios", colors_shifted_reff, linetypes_shifted_reff)
print(p8_shifted)
ggsave("figure8_reff_scenarios_shifted.png", p8_shifted, width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 9: Change in R_eff from Baseline - Separate for each duration
# ============================================================================
cat("\n=== Creating Figure 9: Change in R_eff from Baseline ===\n")

reff_change_start_week <- T_weeks - 4

reff_change_list <- list()
for (typhoon_idx in 2:20) {
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
reff_change_df$scenario <- factor(reff_change_df$scenario, levels = typhoon_scenarios[2:20])

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
    scale_color_manual(
      name = "Scenarios:",
      values = colors
    ) +
    scale_fill_manual(
      name = "Scenarios:",
      values = colors
    ) +
    scale_linetype_manual(
      name = "Scenarios:",
      values = linetypes
    ) +
    scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") +
    labs(
      title = plot_title,
      subtitle = "Negative values indicate reduction in transmission | ΔR_eff = R_eff(scenario) - R_eff(baseline)",
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

# Define colors and linetypes for 3 days change
colors_3_change <- c("3 days (10%)" = "#AED6F1",
                     "3 days (20%)" = "#5DADE2",
                     "3 days (40%)" = "#3498DB",
                     "3 days (60%)" = "#21618C",
                     "3 days (80%)" = "#1B4F72")

linetypes_3_change <- c("3 days (10%)" = "solid",
                        "3 days (20%)" = "dashed",
                        "3 days (40%)" = "dotted",
                        "3 days (60%)" = "dotdash",
                        "3 days (80%)" = "longdash")

selected_3_change <- c("3 days (10%)", "3 days (20%)", "3 days (40%)", "3 days (60%)", "3 days (80%)")

p9_3days <- create_reff_change_plot(selected_3_change, "Change in R_eff from Baseline - 3 Days Duration", colors_3_change, linetypes_3_change)
print(p9_3days)
ggsave("figure9_reff_change_3days.png", p9_3days, width = 14, height = 10, dpi = 300, bg = "white")

# Define colors and linetypes for 5 days change
colors_5_change <- c("5 days (10%)" = "#AED6F1",
                     "5 days (20%)" = "#5DADE2",
                     "5 days (40%)" = "#3498DB",
                     "5 days (60%)" = "#21618C",
                     "5 days (80%)" = "#1B4F72")

linetypes_5_change <- c("5 days (10%)" = "solid",
                        "5 days (20%)" = "dashed",
                        "5 days (40%)" = "dotted",
                        "5 days (60%)" = "dotdash",
                        "5 days (80%)" = "longdash")

selected_5_change <- c("5 days (10%)", "5 days (20%)", "5 days (40%)", "5 days (60%)", "5 days (80%)")

p9_5days <- create_reff_change_plot(selected_5_change, "Change in R_eff from Baseline - 5 Days Duration", colors_5_change, linetypes_5_change)
print(p9_5days)
ggsave("figure9_reff_change_5days.png", p9_5days, width = 14, height = 10, dpi = 300, bg = "white")

# Define colors and linetypes for 7 days change
colors_7_change <- c("7 days (10%)" = "#AED6F1",
                     "7 days (20%)" = "#5DADE2",
                     "7 days (40%)" = "#3498DB",
                     "7 days (60%)" = "#21618C",
                     "7 days (80%)" = "#1B4F72")

linetypes_7_change <- c("7 days (10%)" = "solid",
                        "7 days (20%)" = "dashed",
                        "7 days (40%)" = "dotted",
                        "7 days (60%)" = "dotdash",
                        "7 days (80%)" = "longdash")

selected_7_change <- c("7 days (10%)", "7 days (20%)", "7 days (40%)", "7 days (60%)", "7 days (80%)")

p9_7days <- create_reff_change_plot(selected_7_change, "Change in R_eff from Baseline - 7 Days Duration", colors_7_change, linetypes_7_change)
print(p9_7days)
ggsave("figure9_reff_change_7days.png", p9_7days, width = 14, height = 10, dpi = 300, bg = "white")

# Define colors and linetypes for shifted change
colors_shifted_change <- c("7 days (60%, early 3 days)" = "#AED6F1",
                           "7 days (60%, early 5 days)" = "#5DADE2",
                           "7 days (60%, delayed 3 days)" = "#3498DB",
                           "7 days (60%, delayed 5 days)" = "#21618C")

linetypes_shifted_change <- c("7 days (60%, early 3 days)" = "solid",
                              "7 days (60%, early 5 days)" = "dashed",
                              "7 days (60%, delayed 3 days)" = "dotted",
                              "7 days (60%, delayed 5 days)" = "dotdash")

selected_shifted_change <- c("7 days (60%, early 3 days)", "7 days (60%, early 5 days)", "7 days (60%, delayed 3 days)", "7 days (60%, delayed 5 days)")

p9_shifted <- create_reff_change_plot(selected_shifted_change, "Change in R_eff from Baseline - Shifted Timing Scenarios", colors_shifted_change, linetypes_shifted_change)
print(p9_shifted)
ggsave("figure9_reff_change_shifted.png", p9_shifted, width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 10 (MODIFIED): Typhoon Impact Effectiveness - TYPHOON PERIOD ONLY
# ============================================================================
cat("\n=== Creating Figure 10: Typhoon Impact Effectiveness (Typhoon Period Only) ===\n")

# MODIFIED: Extract reduction for TYPHOON PERIOD ONLY
reduction_typhoon_mean <- apply(reduction_typhoon, c(2,3), mean)
reduction_typhoon_lower <- apply(reduction_typhoon, c(2,3), quantile, probs = 0.025)
reduction_typhoon_upper <- apply(reduction_typhoon, c(2,3), quantile, probs = 0.975)

# Function to create reduction plot for a specific duration (MODIFIED SUBTITLE)
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
    facet_wrap(~ strain, scales = "free_y") +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title = plot_title,
      subtitle = "Error bars: 95% credible intervals | Based on Typhoon Period ONLY (1 week) vs Baseline",
      x = "Typhoon Reduction Scenario",
      y = "Reduction During Typhoon Period (%)",
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

# Define colors for 3 days reduction
colors_3_reduction <- c("0% (Baseline)" = "#999999", 
                        "3 days (10%)" = "#AED6F1", 
                        "3 days (20%)" = "#5DADE2", 
                        "3 days (40%)" = "#3498DB", 
                        "3 days (60%)" = "#21618C", 
                        "3 days (80%)" = "#1B4F72")

p10_3days <- create_reduction_plot(c(1, 2:6), 
                                   "Percentage Reduction in Cases (Typhoon Period) - 3 Days Duration", 
                                   colors_3_reduction)
print(p10_3days)
ggsave("figure10_typhoon_reduction_3days_typhoon_only.png", p10_3days, 
       width = 12, height = 12, dpi = 300, bg = "white")

# Define colors for 5 days reduction
colors_5_reduction <- c("0% (Baseline)" = "#999999", 
                        "5 days (10%)" = "#AED6F1", 
                        "5 days (20%)" = "#5DADE2", 
                        "5 days (40%)" = "#3498DB", 
                        "5 days (60%)" = "#21618C", 
                        "5 days (80%)" = "#1B4F72")

p10_5days <- create_reduction_plot(c(1, 7:11), 
                                   "Percentage Reduction in Cases (Typhoon Period) - 5 Days Duration", 
                                   colors_5_reduction)
print(p10_5days)
ggsave("figure10_typhoon_reduction_5days_typhoon_only.png", p10_5days, 
       width = 12, height = 12, dpi = 300, bg = "white")

# Define colors for 7 days reduction
colors_7_reduction <- c("0% (Baseline)" = "#999999", 
                        "7 days (10%)" = "#AED6F1", 
                        "7 days (20%)" = "#5DADE2", 
                        "7 days (40%)" = "#3498DB", 
                        "7 days (60%)" = "#21618C", 
                        "7 days (80%)" = "#1B4F72")

p10_7days <- create_reduction_plot(c(1, 12:16), 
                                   "Percentage Reduction in Cases (Typhoon Period) - 7 Days Duration", 
                                   colors_7_reduction)
print(p10_7days)
ggsave("figure10_typhoon_reduction_7days_typhoon_only.png", p10_7days, 
       width = 12, height = 12, dpi = 300, bg = "white")

# Define colors for shifted reduction
colors_shifted_reduction <- c("0% (Baseline)" = "#999999", 
                              "7 days (60%, early 3 days)" = "#AED6F1", 
                              "7 days (60%, early 5 days)" = "#5DADE2", 
                              "7 days (60%, delayed 3 days)" = "#3498DB", 
                              "7 days (60%, delayed 5 days)" = "#21618C")

p10_shifted <- create_reduction_plot(c(1, 17:20), 
                                     "Percentage Reduction in Cases (Typhoon Period) - Shifted Timing", 
                                     colors_shifted_reduction)
print(p10_shifted)
ggsave("figure10_typhoon_reduction_shifted_typhoon_only.png", p10_shifted, 
       width = 12, height = 12, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 11: Average Weekly Cases - Typhoon vs Recovery Period
# ============================================================================
cat("\n=== Creating Figure 11: Average Weekly Cases by Period - Separate for each duration ===\n")

avg_typhoon_mean <- apply(avg_weekly_typhoon, c(2,3), mean)
avg_typhoon_lower <- apply(avg_weekly_typhoon, c(2,3), quantile, probs = 0.025)
avg_typhoon_upper <- apply(avg_weekly_typhoon, c(2,3), quantile, probs = 0.975)

avg_recovery_mean <- apply(avg_weekly_recovery, c(2,3), mean)
avg_recovery_lower <- apply(avg_weekly_recovery, c(2,3), quantile, probs = 0.025)
avg_recovery_upper <- apply(avg_weekly_recovery, c(2,3), quantile, probs = 0.975)

# Function to create avg cases plot for selected scenarios
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
      subtitle = "Error bars: 95% credible intervals | Compared across scenarios",
      x = "Typhoon Reduction Scenario",
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

selected_3_avg_indices <- c(1, 2:6)
p11_3days <- create_avg_cases_plot(selected_3_avg_indices, "Average Weekly Cases: Typhoon vs Recovery Period - 3 Days Duration")
print(p11_3days)
ggsave("figure11_avg_weekly_cases_3days.png", p11_3days, width = 14, height = 10, dpi = 300, bg = "white")

selected_5_avg_indices <- c(1, 7:11)
p11_5days <- create_avg_cases_plot(selected_5_avg_indices, "Average Weekly Cases: Typhoon vs Recovery Period - 5 Days Duration")
print(p11_5days)
ggsave("figure11_avg_weekly_cases_5days.png", p11_5days, width = 14, height = 10, dpi = 300, bg = "white")

selected_7_avg_indices <- c(1, 12:16)
p11_7days <- create_avg_cases_plot(selected_7_avg_indices, "Average Weekly Cases: Typhoon vs Recovery Period - 7 Days Duration")
print(p11_7days)
ggsave("figure11_avg_weekly_cases_7days.png", p11_7days, width = 14, height = 10, dpi = 300, bg = "white")

selected_shifted_avg_indices <- c(1, 17:20)
p11_shifted <- create_avg_cases_plot(selected_shifted_avg_indices, "Average Weekly Cases: Typhoon vs Recovery Period - Shifted Timing Scenarios")
print(p11_shifted)
ggsave("figure11_avg_weekly_cases_shifted.png", p11_shifted, width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================================
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
          file = "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/validation_comparison_typhoon_only.csv",
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
# FIGURE 12: Validation Comparison
# ============================================================================
cat("\n=== Creating Figure 12: Validation Comparison ===\n")

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
# FIGURE 13: Visualize Posterior Distributions
# ============================================================================
cat("\n=== Creating Figure 13: Posterior Distribution Figure ===\n")

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

p13_post <- ggplot(post_data, aes(x = value, fill = group)) +
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

print(p13_post)
ggsave("figure13_posterior_distributions.png", p13_post,
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
cat("\n=== Typhoon Period Reductions (1 week only) ===\n")

for (typhoon_idx in 2:20) {
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
# FIGURE 14A: COMBINED SUNBURST CHART - All 6 strains in one figure
# Using manual approach for independent color scales per strain
# ============================================================================
cat("\n=== Creating Figure 14A: Combined Sunburst Chart (6 strains, 2x3 layout) ===\n")

# Build scenario labels
durations <- c(3, 5, 7)
intensities <- c(0.2, 0.4, 0.6, 0.8)
shifts <- c(-7, -5, -3, 0, 3, 5, 7)

cat("Total scenarios:", nrow(incidence_per10k_mean), "\n")

# ============================================================================
# Create sunburst data with variable outer ring extensions
# ============================================================================
create_extension_sunburst <- function(strain_idx) {
  sunburst_data <- data.frame()
  
  # Build hierarchical structure
  for (dur_idx in 1:3) {
    dur_days <- durations[dur_idx]
    
    for (int_idx in 1:4) {
      int_val <- intensities[int_idx]
      
      # Extract all 7 timing incidences for this group
      temp_incidences <- numeric(7)
      for (shift_idx in 1:7) {
        scenario_idx <- 1 + (dur_idx-1)*28 + (int_idx-1)*7 + shift_idx
        temp_incidences[shift_idx] <- incidence_per10k_mean[scenario_idx, strain_idx]
      }
      
      # Calculate group-wise min/max for normalization
      max_inc_group <- max(temp_incidences)
      min_inc_group <- min(temp_incidences)
      
      # Now process each of the 7 timing scenarios
      for (shift_idx in 1:7) {
        shift_val <- shifts[shift_idx]
        scenario_idx <- 1 + (dur_idx-1)*28 + (int_idx-1)*7 + shift_idx
        incidence_val <- incidence_per10k_mean[scenario_idx, strain_idx]
        
        # Timing labels
        if (shift_val < 0) {
          shift_label <- paste0("E", abs(shift_val))
        } else if (shift_val > 0) {
          shift_label <- paste0("D", shift_val)
        } else {
          shift_label <- "On"
        }
        
        segment_idx <- (dur_idx-1)*28 + (int_idx-1)*7 + shift_idx
        
        # Group-wise normalization for EXTENSION LENGTH
        if (max_inc_group > min_inc_group) {
          extension_length <- 0.3 + 0.9 * (incidence_val - min_inc_group) /
            (max_inc_group - min_inc_group)
        } else {
          extension_length <- 0.75
        }
        
        # Calculate angular positions (84 equal segments)
        theta_start <- (segment_idx - 1) * (360 / 84)
        theta_end <- segment_idx * (360 / 84)
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
          group_id = paste(dur_idx, int_idx, sep = "_"),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(sunburst_data)
}

# Create individual plots for each strain, then combine with patchwork
library(patchwork)

create_individual_sunburst_subplot <- function(strain_idx) {
  strain_name <- strain_names[strain_idx]
  
  # Get data for this strain
  sunburst_data <- create_extension_sunburst(strain_idx)
  
  # Calculate aggregated data for inner and middle rings
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
  
  # Create plot with independent color scale
  p <- ggplot() +
    # Outer ring: variable xmax for extension
    geom_rect(
      data = sunburst_data,
      aes(xmin = 2.5,
          xmax = 2.5 + extension_length * 1.5,
          ymin = theta_start,
          ymax = theta_end,
          fill = incidence),
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
    
    # Labels for duration
    geom_text(
      data = duration_data,
      aes(x = 1, y = theta_mid,
          label = paste0(duration, "d"),
          angle = ifelse(theta_mid > 180, theta_mid - 90, theta_mid + 90)),
      size = 2.5, fontface = "bold", color = "black"
    ) +
    
    # Labels for intensity
    geom_text(
      data = intensity_data,
      aes(x = 2, y = theta_mid,
          label = paste0(intensity, "%"),
          angle = ifelse(theta_mid > 180, theta_mid - 90, theta_mid + 90)),
      size = 1.8, color = "black"
    ) +
    
    # Scenario names INSIDE the outer ring
    geom_text(
      data = sunburst_data,
      aes(x = 2.5 + (extension_length * 1.5) * 0.5,
          y = theta_mid,
          label = shift_label,
          angle = ifelse(theta_mid > 180, theta_mid - 90, theta_mid + 90)),
      size = 1.0, color = "white", fontface = "bold"
    ) +
    
    # Incidence values OUTSIDE the extended outer ring
    geom_text(
      data = sunburst_data,
      aes(x = 2.5 + extension_length * 1.5 + 0.15,
          y = theta_mid,
          label = sprintf("%.1f", incidence),
          angle = ifelse(theta_mid > 180, theta_mid - 90, theta_mid + 90)),
      size = 0.8, color = "gray20", fontface = "plain"
    ) +
    
    # Color scale: independent for each strain
    scale_fill_gradientn(
      colors = c("#1A9641", "#A6D96A", "#FFFFBF", "#FDAE61", "#D7191C"),
      values = scales::rescale(c(0, 25, 50, 75, 100)),
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
    
    labs(title = strain_name) +
    
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

cat("\n=== Creating individual sunburst subplots ===\n")

# Create all 6 subplots
p1 <- create_individual_sunburst_subplot(1)  # B
p2 <- create_individual_sunburst_subplot(2)  # H3
p3 <- create_individual_sunburst_subplot(3)  # H1
p4 <- create_individual_sunburst_subplot(4)  # COVID
p5 <- create_individual_sunburst_subplot(5)  # RSV
p6 <- create_individual_sunburst_subplot(6)  # HFMD

# Combine using patchwork
p14a_combined <- (p1 | p2 | p3) / (p4 | p5 | p6) +
  plot_annotation(
    title = "Typhoon Impact Scenarios: Case Incidence per 10,000 Population",
    subtitle = "Inner: Duration (3d/5d/7d) | Middle: Intensity (20%/40%/60%/80%) | Outer: Timing (E=Early, D=Delayed, On=On time)\nExtension length = relative incidence within intensity group | Color: green = lower, red = higher | Each strain has independent color scale",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30", lineheight = 1.2)
    )
  )

print(p14a_combined)

ggsave(
  filename = "figure14a_sunburst_combined_6strains.png",
  plot = p14a_combined,
  width = 20,
  height = 13,
  dpi = 300,
  bg = "white"
)

cat("\nSaved: figure14a_sunburst_combined_6strains.png\n")

# ============================================================================
# Print verification
# ============================================================================
cat("\n=== Verification of Variable Outer Ring Extension Design ===\n")
cat("Key design:\n")
cat("1. All 6 strains in one figure (2 rows × 3 columns)\n")
cat("2. Each strain has independent color scale and color bar\n")
cat("3. Outer ring: FIXED xmin=2.5, VARIABLE xmax based on incidence\n")
cat("4. Longer extension = higher incidence within each intensity group\n")
cat("5. Scenario names inside extended ring\n")
cat("6. Incidence values outside the extended ring\n\n")

# Get sample data for verification
sample_data <- create_extension_sunburst(1)  # B strain
sample_group <- sample_data %>%
  filter(duration == 3, intensity == 20) %>%
  arrange(shift) %>%
  select(shift_label, incidence, extension_length, theta_mid)

cat("Sample verification for strain B, 3 days, 20% intensity:\n")
print(sample_group)

cat("\nExtension length verification:\n")
cat("Shortest extension:", sample_group$shift_label[which.min(sample_group$extension_length)],
    "=", min(sample_group$extension_length), "\n")
cat("Longest extension:", sample_group$shift_label[which.max(sample_group$extension_length)],
    "=", max(sample_group$extension_length), "\n")
cat("Extension range: 0.3 to 1.2 units (multiplied by 1.5 in plot)\n")

# ============================================================================
# FIGURE 14B: Heatmap Grid for all scenarios
# ============================================================================
cat("\n=== Creating Figure 14B: Reduction Heatmap ===\n")

grid_data_list <- list()
for (dur_idx in 1:3) {
  for (int_idx in 1:4) {
    for (shift_idx in 1:7) {
      scenario_idx <- 1 + (dur_idx-1)*28 + (int_idx-1)*7 + shift_idx
      
      for (strain_idx in 1:N_strains) {
        grid_data_list[[length(grid_data_list) + 1]] <- data.frame(
          duration = factor(durations[dur_idx], levels = durations),
          intensity = factor(paste0(as.integer(intensities[int_idx]*100), "%")),
          shift = factor(shifts[shift_idx]),
          reduction = reduction_extended_mean[scenario_idx, strain_idx],
          strain = strain_names[strain_idx],
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

grid_df <- do.call(rbind, grid_data_list)
grid_df$strain <- factor(grid_df$strain, levels = strain_names)

p14b_grid <- ggplot(grid_df, aes(x = shift, y = intensity, fill = reduction)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%.0f", reduction)), size = 2.5, color = "black") +
  facet_grid(duration ~ strain,
             labeller = labeller(duration = function(x) paste(x, "days"))) +
  scale_fill_viridis_c(option = "plasma", begin = 0.1, end = 0.9,
                       name = "Reduction\n(%)",
                       limits = c(0, max(grid_df$reduction))) +
  scale_x_discrete(labels = c("-7" = "E7", "-5" = "E5", "-3" = "E3", "0" = "On", 
                              "3" = "D3", "5" = "D5", "7" = "D7")) +
  labs(
    title = "Typhoon Impact: Case Reduction by Duration, Intensity, and Timing",
    subtitle = "Values show percentage reduction in cases during typhoon period (compared to baseline)\nColumns: Strains | Rows: Duration | X-axis: Timing (E=Early, D=Delayed, On=On time) | Y-axis: Intensity",
    x = "Timing Shift (days)",
    y = "Transmission Reduction Intensity"
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

print(p14b_grid)
ggsave("figure14b_reduction_heatmap.png", p14b_grid,
       width = 16, height = 10, dpi = 300, bg = "white")

cat("\n=== Figure 14 Complete ===\n")
cat("Generated files:\n")
cat("  Combined sunburst (1 file):\n")
cat("    - figure14a_sunburst_combined_6strains.png (2×3 layout, each with independent color bar)\n")
cat("  Heatmap grid (1 file):\n")
cat("    - figure14b_reduction_heatmap.png\n")

cat("\n=== All Figure 14 visualizations complete ===\n")
cat("IMPORTANT: Each strain has its own independent color scale in the combined figure.\n")
cat("The outer ring extension length represents relative incidence within each intensity group.\n")