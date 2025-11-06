# ============================================================================
# TYPHOON IMPACT ANALYSIS - CORRECTED ERROR BAR CALCULATION - ALL FIGURES
# ============================================================================
#Sys.setenv('R_MAX_VSIZE' = 32 * 1024^3)
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

# 2. Build Stan data
T_weeks <- nrow(fitting_data)
T_weeks_validation <- nrow(validation_data)
T_weeks_forecast <- 8
typhoon_weeks <- 1
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
cat("Typhoon duration and reduction scenarios: baseline, and combinations of 3,7,10 days with 10,20,40,60,80% reductions\n")
cat("Additional scenarios: early 3 days, early 5 days, delayed 3 days, delayed 5 days (7 days duration, 60% reduction)\n")
cat("Validation data from 2025-09-23 to 2025-10-11 will be shown as yellow diamonds\n\n")

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

saveRDS(fit, file = "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/typhoon_model_fit_corrected.rds")


#cat("=== LOADING PREVIOUS FIT (typhoon_model_fit_corrected.rds) ===\n")
#fit <- readRDS("/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/typhoon_model_fit_corrected.rds")

# 4. Extract results - INCLUDING THE NEW QUANTITIES
pred_cases <- rstan::extract(fit, pars = "pred_cases")$pred_cases
forecast_cases <- rstan::extract(fit, pars = "forecast_cases")$forecast_cases
forecast_cases_typhoon <- rstan::extract(fit, pars = "forecast_cases_typhoon")$forecast_cases_typhoon
R_eff <- rstan::extract(fit, pars = "R_eff")$R_eff
R_eff_scenarios <- rstan::extract(fit, pars = "R_eff_scenarios")$R_eff_scenarios
R0_t <- rstan::extract(fit, pars = "R0_t")$R0_t

# NEW: Extract the properly calculated reduction percentages
reduction_typhoon <- rstan::extract(fit, pars = "reduction_typhoon_period")$reduction_typhoon_period
reduction_total <- rstan::extract(fit, pars = "reduction_total_period")$reduction_total_period
avg_weekly_typhoon <- rstan::extract(fit, pars = "avg_weekly_typhoon")$avg_weekly_typhoon
avg_weekly_recovery <- rstan::extract(fit, pars = "avg_weekly_recovery")$avg_weekly_recovery

cat("\n=== Checking extracted dimensions ===\n")
cat("pred_cases dimensions:", dim(pred_cases), "\n")
cat("forecast_cases dimensions:", dim(forecast_cases), "\n")
cat("forecast_cases_typhoon dimensions:", dim(forecast_cases_typhoon), "\n")
cat("R_eff dimensions:", dim(R_eff), "\n")
cat("R_eff_scenarios dimensions:", dim(R_eff_scenarios), "\n")
cat("R0_t dimensions:", dim(R0_t), "\n")
cat("reduction_typhoon dimensions:", dim(reduction_typhoon), "\n")
cat("reduction_total dimensions:", dim(reduction_total), "\n")
cat("avg_weekly_typhoon dimensions:", dim(avg_weekly_typhoon), "\n")
cat("avg_weekly_recovery dimensions:", dim(avg_weekly_recovery), "\n")

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
                       "7 days (10%)", "7 days (20%)", "7 days (40%)", "7 days (60%)", "7 days (80%)",
                       "10 days (10%)", "10 days (20%)", "10 days (40%)", "10 days (60%)", "10 days (80%)",
                       "7 days (60%, early 3 days)", "7 days (60%, early 5 days)", "7 days (60%, delayed 3 days)", "7 days (60%, delayed 5 days)")

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
    subtitle = "Black dots: Observed (fitted) | Yellow diamonds: Validation data (NOT fitted) | Blue: Model fit | Red: Baseline forecast\nGray dashed: Typhoon start | Red dashed: Typhoon end",
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
               aes(y = observed),
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
      subtitle = "Blue: Model fit | Yellow diamonds: Validation data (NOT fitted) | Colored lines: Forecast scenarios\nGray dashed: Typhoon start (2025-09-23) | Red dashed: Typhoon end (2025-10-07) | Pink: Typhoon (2 weeks) | Green: Recovery (3 weeks)",
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

# Define colors and linetypes for 10 days
colors_10 <- c("0% (Baseline)" = "#999999",
               "10 days (10%)" = "#AED6F1",
               "10 days (20%)" = "#5DADE2",
               "10 days (40%)" = "#3498DB",
               "10 days (60%)" = "#21618C",
               "10 days (80%)" = "#1B4F72")
linetypes_10 <- c("0% (Baseline)" = "solid",
                  "10 days (10%)" = "solid",
                  "10 days (20%)" = "dashed",
                  "10 days (40%)" = "dotted",
                  "10 days (60%)" = "dotdash",
                  "10 days (80%)" = "longdash")
selected_10 <- c("Historical Fit", "0% (Baseline)", "10 days (10%)", "10 days (20%)", "10 days (40%)", "10 days (60%)", "10 days (80%)")
p2_10days <- create_scenario_plot(selected_10, "Typhoon Impact: Recent Period and 10 Days Duration Scenarios", colors_10, linetypes_10)
print(p2_10days)
ggsave("figure2_typhoon_scenarios_10days.png", p2_10days, width = 14, height = 10, dpi = 300, bg = "white")

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
p2_shifted <- create_scenario_plot(selected_shifted, "Typhoon Impact: Recent Period and Shifted Timing Scenarios (7 days, 60% reduction)", colors_shifted, linetypes_shifted)
print(p2_shifted)
ggsave("figure2_typhoon_scenarios_shifted.png", p2_shifted, width = 14, height = 10, dpi = 300, bg = "white")

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
# FIGURE 4: Typhoon Impact Effectiveness - CORRECTED WITH PROPER ERROR BARS - Separate for each duration
# ============================================================================
cat("\n=== Creating Figure 4: Typhoon Impact Effectiveness - Separate for each duration ===\n")

# Prepare reduction data for typhoon period and total period
reduction_typhoon_mean <- apply(reduction_typhoon, c(2,3), mean)
reduction_typhoon_lower <- apply(reduction_typhoon, c(2,3), quantile, probs = 0.025)
reduction_typhoon_upper <- apply(reduction_typhoon, c(2,3), quantile, probs = 0.975)

reduction_total_mean <- apply(reduction_total, c(2,3), mean)
reduction_total_lower <- apply(reduction_total, c(2,3), quantile, probs = 0.025)
reduction_total_upper <- apply(reduction_total, c(2,3), quantile, probs = 0.975)

# Function to create reduction plot for a specific duration
create_reduction_plot <- function(scenario_indices, plot_title, colors) {
  selected_scenarios <- typhoon_scenarios[scenario_indices]
  
  reduction_df_typhoon <- data.frame(
    scenario = rep(selected_scenarios, N_strains),
    strain = rep(strain_names, each = length(selected_scenarios)),
    mean = as.vector(reduction_typhoon_mean[scenario_indices,]),
    lower = as.vector(reduction_typhoon_lower[scenario_indices,]),
    upper = as.vector(reduction_typhoon_upper[scenario_indices,]),
    period = "Typhoon Period (2 weeks)"
  )
  
  reduction_df_total <- data.frame(
    scenario = rep(selected_scenarios, N_strains),
    strain = rep(strain_names, each = length(selected_scenarios)),
    mean = as.vector(reduction_total_mean[scenario_indices,]),
    lower = as.vector(reduction_total_lower[scenario_indices,]),
    upper = as.vector(reduction_total_upper[scenario_indices,]),
    period = "Total Period (5 weeks)"
  )
  
  reduction_df <- rbind(reduction_df_typhoon, reduction_df_total)
  reduction_df$scenario <- factor(reduction_df$scenario, levels = selected_scenarios)
  reduction_df$strain <- factor(reduction_df$strain, levels = strain_names)
  reduction_df$period <- factor(reduction_df$period, levels = c("Typhoon Period (2 weeks)", "Total Period (5 weeks)"))
  
  p <- ggplot(reduction_df, aes(x = scenario, y = mean, fill = scenario)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.8), width = 0.25, size = 0.6) +
    facet_grid(strain ~ period, scales = "free_y") +
    scale_fill_manual(
      values = colors
    ) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title = plot_title,
      subtitle = "Error bars: 95% credible intervals | Calculated per MCMC draw for accurate uncertainty",
      x = "Typhoon Reduction Scenario",
      y = "Reduction (%)",
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

# *** MODIFICATION: Added baseline ("#999999") to color maps ***
# Define colors for 3 days reduction
colors_3_reduction <- c("0% (Baseline)" = "#999999", "3 days (10%)" = "#AED6F1", "3 days (20%)" = "#5DADE2", "3 days (40%)" = "#3498DB", "3 days (60%)" = "#21618C", "3 days (80%)" = "#1B4F72")
# *** MODIFICATION: Changed 2:6 to c(1, 2:6) to include baseline ***
p4_3days <- create_reduction_plot(c(1, 2:6), "Percentage Reduction in Cases Compared to Baseline - 3 Days Duration", colors_3_reduction)
print(p4_3days)
ggsave("figure4_typhoon_reduction_3days.png", p4_3days, width = 12, height = 12, dpi = 300, bg = "white")

# Define colors for 7 days reduction
colors_7_reduction <- c("0% (Baseline)" = "#999999", "7 days (10%)" = "#AED6F1", "7 days (20%)" = "#5DADE2", "7 days (40%)" = "#3498DB", "7 days (60%)" = "#21618C", "7 days (80%)" = "#1B4F72")
# *** MODIFICATION: Changed 7:11 to c(1, 7:11) to include baseline ***
p4_7days <- create_reduction_plot(c(1, 7:11), "Percentage Reduction in Cases Compared to Baseline - 7 Days Duration", colors_7_reduction)
print(p4_7days)
ggsave("figure4_typhoon_reduction_7days.png", p4_7days, width = 12, height = 12, dpi = 300, bg = "white")

# Define colors for 10 days reduction
colors_10_reduction <- c("0% (Baseline)" = "#999999", "10 days (10%)" = "#AED6F1", "10 days (20%)" = "#5DADE2", "10 days (40%)" = "#3498DB", "10 days (60%)" = "#21618C", "10 days (80%)" = "#1B4F72")
# *** MODIFICATION: Changed 12:16 to c(1, 12:16) to include baseline ***
p4_10days <- create_reduction_plot(c(1, 12:16), "Percentage Reduction in Cases Compared to Baseline - 10 Days Duration", colors_10_reduction)
print(p4_10days)
ggsave("figure4_typhoon_reduction_10days.png", p4_10days, width = 12, height = 12, dpi = 300, bg = "white")

# Define colors for shifted reduction
colors_shifted_reduction <- c("0% (Baseline)" = "#999999", "7 days (60%, early 3 days)" = "#AED6F1", "7 days (60%, early 5 days)" = "#5DADE2", "7 days (60%, delayed 3 days)" = "#3498DB", "7 days (60%, delayed 5 days)" = "#21618C")
# *** MODIFICATION: Changed 17:20 to c(1, 17:20) to include baseline ***
p4_shifted <- create_reduction_plot(c(1, 17:20), "Percentage Reduction in Cases Compared to Baseline - Shifted Timing Scenarios", colors_shifted_reduction)
print(p4_shifted)
ggsave("figure4_typhoon_reduction_shifted.png", p4_shifted, width = 12, height = 12, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 5: Average Weekly Cases - Typhoon vs Recovery Period - Separate for each duration
# ============================================================================
cat("\n=== Creating Figure 5: Average Weekly Cases by Period - Separate for each duration ===\n")

avg_typhoon_mean <- apply(avg_weekly_typhoon, c(2,3), mean)
avg_typhoon_lower <- apply(avg_weekly_typhoon, c(2,3), quantile, probs = 0.025)
avg_typhoon_upper <- apply(avg_weekly_typhoon, c(2,3), quantile, probs = 0.975)

avg_recovery_mean <- apply(avg_weekly_recovery, c(2,3), mean)
avg_recovery_lower <- apply(avg_weekly_recovery, c(2,3), quantile, probs = 0.025)
avg_recovery_upper <- apply(avg_weekly_recovery, c(2,3), quantile, probs = 0.975)

# Function to create avg cases plot for selected scenarios
# (This function was corrected in the previous step and is correct)
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

# *** MODIFICATION: Changed 2:6 to c(1, 2:6) to include baseline ***
selected_3_avg_indices <- c(1, 2:6)
p5_3days <- create_avg_cases_plot(selected_3_avg_indices, "Average Weekly Cases: Typhoon vs Recovery Period - 3 Days Duration")
print(p5_3days)
ggsave("figure5_avg_weekly_cases_3days.png", p5_3days, width = 14, height = 10, dpi = 300, bg = "white")

# *** MODIFICATION: Changed 7:11 to c(1, 7:11) to include baseline ***
selected_7_avg_indices <- c(1, 7:11)
p5_7days <- create_avg_cases_plot(selected_7_avg_indices, "Average Weekly Cases: Typhoon vs Recovery Period - 7 Days Duration")
print(p5_7days)
ggsave("figure5_avg_weekly_cases_7days.png", p5_7days, width = 14, height = 10, dpi = 300, bg = "white")

# *** MODIFICATION: Changed 12:16 to c(1, 12:16) to include baseline ***
selected_10_avg_indices <- c(1, 12:16)
p5_10days <- create_avg_cases_plot(selected_10_avg_indices, "Average Weekly Cases: Typhoon vs Recovery Period - 10 Days Duration")
print(p5_10days)
ggsave("figure5_avg_weekly_cases_10days.png", p5_10days, width = 14, height = 10, dpi = 300, bg = "white")

# *** MODIFICATION: Changed 17:20 to c(1, 17:20) to include baseline ***
selected_shifted_avg_indices <- c(1, 17:20)
p5_shifted <- create_avg_cases_plot(selected_shifted_avg_indices, "Average Weekly Cases: Typhoon vs Recovery Period - Shifted Timing Scenarios")
print(p5_shifted)
ggsave("figure5_avg_weekly_cases_shifted.png", p5_shifted, width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 6: Effective Reproduction Number (Rt) by Scenario (Recent Period) - Separate for each duration
# ============================================================================
cat("\n=== Creating Figure 6: R_eff by Scenario (Recent Period) - Separate for each duration ===\n")

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
    scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
    labs(
      title = plot_title,
      subtitle = "Comparing R_eff trajectories under different transmission reduction scenarios | Pink: Typhoon (2 weeks) | Green: Recovery (3 weeks)",
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
colors_3_reff <- c("0% (Baseline)" = "#999999",
                   "3 days (10%)" = "#AED6F1",
                   "3 days (20%)" = "#5DADE2",
                   "3 days (40%)" = "#3498DB",
                   "3 days (60%)" = "#21618C",
                   "3 days (80%)" = "#1B4F72")
linetypes_3_reff <- c("0% (Baseline)" = "solid",
                      "3 days (10%)" = "solid",
                      "3 days (20%)" = "dashed",
                      "3 days (40%)" = "dotted",
                      "3 days (60%)" = "dotdash",
                      "3 days (80%)" = "longdash")
selected_3_reff <- c("0% (Baseline)", "3 days (10%)", "3 days (20%)", "3 days (40%)", "3 days (60%)", "3 days (80%)")
p6_3days <- create_reff_plot(selected_3_reff, "Effective Reproduction Number (R_eff) by Typhoon Scenario (Recent Period) - 3 Days Duration", colors_3_reff, linetypes_3_reff)
print(p6_3days)
ggsave("figure6_reff_scenarios_3days.png", p6_3days, width = 14, height = 10, dpi = 300, bg = "white")

# Define colors and linetypes for 7 days reff
colors_7_reff <- c("0% (Baseline)" = "#999999",
                   "7 days (10%)" = "#AED6F1",
                   "7 days (20%)" = "#5DADE2",
                   "7 days (40%)" = "#3498DB",
                   "7 days (60%)" = "#21618C",
                   "7 days (80%)" = "#1B4F72")
linetypes_7_reff <- c("0% (Baseline)" = "solid",
                      "7 days (10%)" = "solid",
                      "7 days (20%)" = "dashed",
                      "7 days (40%)" = "dotted",
                      "7 days (60%)" = "dotdash",
                      "7 days (80%)" = "longdash")
selected_7_reff <- c("0% (Baseline)", "7 days (10%)", "7 days (20%)", "7 days (40%)", "7 days (60%)", "7 days (80%)")
p6_7days <- create_reff_plot(selected_7_reff, "Effective Reproduction Number (R_eff) by Typhoon Scenario (Recent Period) - 7 Days Duration", colors_7_reff, linetypes_7_reff)
print(p6_7days)
ggsave("figure6_reff_scenarios_7days.png", p6_7days, width = 14, height = 10, dpi = 300, bg = "white")

# Define colors and linetypes for 10 days reff
colors_10_reff <- c("0% (Baseline)" = "#999999",
                    "10 days (10%)" = "#AED6F1",
                    "10 days (20%)" = "#5DADE2",
                    "10 days (40%)" = "#3498DB",
                    "10 days (60%)" = "#21618C",
                    "10 days (80%)" = "#1B4F72")
linetypes_10_reff <- c("0% (Baseline)" = "solid",
                       "10 days (10%)" = "solid",
                       "10 days (20%)" = "dashed",
                       "10 days (40%)" = "dotted",
                       "10 days (60%)" = "dotdash",
                       "10 days (80%)" = "longdash")
selected_10_reff <- c("0% (Baseline)", "10 days (10%)", "10 days (20%)", "10 days (40%)", "10 days (60%)", "10 days (80%)")
p6_10days <- create_reff_plot(selected_10_reff, "Effective Reproduction Number (R_eff) by Typhoon Scenario (Recent Period) - 10 Days Duration", colors_10_reff, linetypes_10_reff)
print(p6_10days)
ggsave("figure6_reff_scenarios_10days.png", p6_10days, width = 14, height = 10, dpi = 300, bg = "white")

# Define colors and linetypes for shifted reff
colors_shifted_reff <- c("0% (Baseline)" = "#999999",
                         "7 days (60%, early 3 days)" = "#AED6F1",
                         "7 days (60%, early 5 days)" = "#5DADE2",
                         "7 days (60%, delayed 3 days)" = "#3498DB",
                         "7 days (60%, delayed 5 days)" = "#21618C")
linetypes_shifted_reff <- c("0% (Baseline)" = "solid",
                            "7 days (60%, early 3 days)" = "solid",
                            "7 days (60%, early 5 days)" = "dashed",
                            "7 days (60%, delayed 3 days)" = "dotted",
                            "7 days (60%, delayed 5 days)" = "dotdash")
selected_shifted_reff <- c("0% (Baseline)", "7 days (60%, early 3 days)", "7 days (60%, early 5 days)", "7 days (60%, delayed 3 days)", "7 days (60%, delayed 5 days)")
p6_shifted <- create_reff_plot(selected_shifted_reff, "Effective Reproduction Number (R_eff) by Typhoon Scenario (Recent Period) - Shifted Timing Scenarios", colors_shifted_reff, linetypes_shifted_reff)
print(p6_shifted)
ggsave("figure6_reff_scenarios_shifted.png", p6_shifted, width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================================
# FIGURE 7: R_eff Change (Difference from Baseline) - Separate for each duration
# ============================================================================
cat("\n=== Creating Figure 7: R_eff Change from Baseline - Separate for each duration ===\n")

reff_change_start_week <- reff_start_week 
reff_change_list <- list()

# Loop scenarios 2 to 20
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
    
    # Calculate statistics across samples (margin = 2 is by column/week)
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
p7_3days <- create_reff_change_plot(selected_3_change, "Change in R_eff from Baseline - 3 Days Duration", colors_3_change, linetypes_3_change)
print(p7_3days)
ggsave("figure7_reff_change_3days.png", p7_3days, width = 14, height = 10, dpi = 300, bg = "white")

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
p7_7days <- create_reff_change_plot(selected_7_change, "Change in R_eff from Baseline - 7 Days Duration", colors_7_change, linetypes_7_change)
print(p7_7days)
ggsave("figure7_reff_change_7days.png", p7_7days, width = 14, height = 10, dpi = 300, bg = "white")

# Define colors and linetypes for 10 days change
colors_10_change <- c("10 days (10%)" = "#AED6F1",
                      "10 days (20%)" = "#5DADE2",
                      "10 days (40%)" = "#3498DB",
                      "10 days (60%)" = "#21618C",
                      "10 days (80%)" = "#1B4F72")
linetypes_10_change <- c("10 days (10%)" = "solid",
                         "10 days (20%)" = "dashed",
                         "10 days (40%)" = "dotted",
                         "10 days (60%)" = "dotdash",
                         "10 days (80%)" = "longdash")
selected_10_change <- c("10 days (10%)", "10 days (20%)", "10 days (40%)", "10 days (60%)", "10 days (80%)")
p7_10days <- create_reff_change_plot(selected_10_change, "Change in R_eff from Baseline - 10 Days Duration", colors_10_change, linetypes_10_change)
print(p7_10days)
ggsave("figure7_reff_change_10days.png", p7_10days, width = 14, height = 10, dpi = 300, bg = "white")

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
p7_shifted <- create_reff_change_plot(selected_shifted_change, "Change in R_eff from Baseline - Shifted Timing Scenarios", colors_shifted_change, linetypes_shifted_change)
print(p7_shifted)
ggsave("figure7_reff_change_shifted.png", p7_shifted, width = 14, height = 10, dpi = 300, bg = "white")

cat("\n=== All figures created successfully ===\n")

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
          file = "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/validation_comparison_corrected.csv",
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
# FIGURE 8: Validation Comparison
# ============================================================================
cat("\n=== Creating Figure 8: Validation Comparison ===\n")

if (nrow(validation_comparison) > 0) {
  p8_validation <- ggplot(validation_comparison, aes(x = baseline_pred, y = observed)) +
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
  
  print(p8_validation)
  ggsave("figure8_validation_comparison.png", p8_validation,
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
# FIGURE 9: Visualize Posterior Distributions
# ============================================================================
cat("\n=== Creating Figure 9: Posterior Distribution Figure ===\n")

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

p9_post <- ggplot(post_data, aes(x = value, fill = group)) +
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

print(p9_post)
ggsave("figure9_posterior_distributions.png", p9_post,
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
cat("\n=== VERIFICATION: Reduction Statistics (CORRECTED) ===\n")
cat("Typhoon Period Reductions (2 weeks):\n")
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

cat("\n=== Total Period Reductions (5 weeks) ===\n")
for (typhoon_idx in 2:20) {
  cat(sprintf("\n%s:\n", typhoon_scenarios[typhoon_idx]))
  for (strain_idx in 1:N_strains) {
    red_dist <- reduction_total[, typhoon_idx, strain_idx]
    cat(sprintf("  %s: %.1f%% [%.1f%%, %.1f%%]\n",
                strain_names[strain_idx],
                mean(red_dist),
                quantile(red_dist, 0.025),
                quantile(red_dist, 0.975)))
  }
}

cat("\n=== Analysis Complete ===\n")
cat("All 9 figures have been generated successfully.\n")
cat("The error bars now correctly represent uncertainty in the reduction percentages,\n")
cat("calculated from the posterior distribution for each MCMC draw.\n")