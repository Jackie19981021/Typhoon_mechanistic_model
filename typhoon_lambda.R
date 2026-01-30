# ============================================================================
# TYPHOON vs NON-TYPHOON COMPARISON ANALYSIS (CORRECTED)
# Full fitting period: 2025-06-01 to 2025-12-13
# Key changes: 
# 1. Typhoon knot placed at TYPHOON END date
# 2. Beta ratio calculated in Stan for correct CI
# 3. Visualization ensures smooth connection at typhoon start date
# 4. Added Cori Rt estimation with pathogen-specific SI parameters
# 5. Using spline interpolation for weekly to daily conversion
# 6. Added weekly incidence visualization below R_eff plots
# ============================================================================
library(rstan)
library(ggplot2)
library(dplyr)
library(lubridate)
library(tidyr)
library(patchwork)
library(EpiEstim)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# ============================================================================
# 1. DATA LOADING
# ============================================================================
cat("=== Loading data ===\n")
data_path <- "/Users/chenjiaqi/Desktop/COVID-19_HK/typhoon/HK_ILI_COVID_Sep.csv"
raw_data <- read.csv(data_path)
raw_data$date <- as.Date(raw_data$date, format = "%Y/%m/%d")

# Full fitting period
start_date <- as.Date("2025-06-01")
end_date <- as.Date("2025-12-13")
typhoon_start_date <- as.Date("2025-09-20")
typhoon_duration_days <- 7
typhoon_end_date <- typhoon_start_date + typhoon_duration_days

fitting_data <- raw_data[raw_data$date >= start_date & raw_data$date <= end_date, ]
fitting_data <- fitting_data[complete.cases(fitting_data), ]

strain_names <- c("B", "H3", "H1", "COVID", "RSV", "HFMD")
fitting_data <- fitting_data[, c("date", strain_names)]

cat("Full fitting period:", as.character(range(fitting_data$date)), "\n")
cat("Total weeks:", nrow(fitting_data), "\n")
cat("Typhoon starts:", as.character(typhoon_start_date), "\n")
cat("Typhoon ends:", as.character(typhoon_end_date), "\n")
cat("Typhoon knot will be placed at END date\n\n")

# ============================================================================
# 2. PREPARE DATA FOR BOTH SCENARIOS
# ============================================================================
T_weeks <- nrow(fitting_data)
N_strains <- 6

# Calculate typhoon start week
typhoon_start_week <- which(fitting_data$date == typhoon_start_date)
cat("Typhoon start week index:", typhoon_start_week, "\n\n")

cases_matrix <- as.matrix(fitting_data[, strain_names])
cases_matrix[cases_matrix < 0] <- 0
cases_matrix <- round(cases_matrix)

# ============================================================================
# 3. DEFINE PATHOGEN-SPECIFIC SERIAL INTERVAL PARAMETERS
# ============================================================================
cat("=== Serial Interval Parameters (Days) ===\n")

# Serial Interval 参数（日尺度）
si_params <- list(
  COVID = list(mean_si = 5.2, std_si = 1.72),
  RSV = list(mean_si = 6.3, std_si = 2.5),
  HFMD = list(mean_si = 4.0, std_si = 2.0),
  H1 = list(mean_si = 2.6, std_si = 1.5),
  H3 = list(mean_si = 2.6, std_si = 1.5),
  B = list(mean_si = 2.6, std_si = 1.5)
)

# Print SI parameters
for (strain in strain_names) {
  cat(sprintf("  %s: Mean SI = %.1f days, SD = %.2f days\n", 
              strain, si_params[[strain]]$mean_si, si_params[[strain]]$std_si))
}
cat("\n")

# ============================================================================
# 4. CALCULATE CORI Rt FOR ALL STRAINS
# ============================================================================
cat("=== Calculating Cori Rt with pathogen-specific SI ===\n")

# Hong Kong total population
hk_population <- 7524000

# Function to convert weekly cumulative cases to daily per 10,000 population using spline interpolation
weekly_to_daily_per10k_spline <- function(weekly_cases, weekly_dates, population) {
  # Convert to per 10,000 population
  weekly_per10k <- (weekly_cases / population) * 12000
  
  # Create time points for weekly data (in days from start)
  time_weekly <- as.numeric(weekly_dates - weekly_dates[1])
  
  # Create daily time points
  time_daily <- seq(from = 0, 
                    to = time_weekly[length(time_weekly)] + 6, 
                    by = 1)
  
  # Use splinefun with monoH.FC for shape-preserving interpolation
  # This preserves monotonicity in each segment without requiring overall monotonicity
  spline_func <- splinefun(x = time_weekly, 
                           y = weekly_per10k, 
                           method = "monoH.FC")
  
  # Apply the spline function to get daily values
  daily_per10k <- spline_func(time_daily)
  
  # Ensure non-negative values
  daily_per10k <- pmax(daily_per10k, 0)
  
  return(daily_per10k)
}

# Function to calculate Cori Rt with pathogen-specific SI
calculate_cori_rt <- function(weekly_cases, weekly_dates, population, strain_name, si_mean, si_std, window = 6) {
  cat(sprintf("  Processing %s (SI: mean=%.1f, sd=%.2f, window=%d)...\n", 
              strain_name, si_mean, si_std, window))
  
  # Convert to daily per 10,000 using spline interpolation
  daily_cases <- weekly_to_daily_per10k_spline(weekly_cases, weekly_dates, population)
  
  # Create incidence object
  dates_daily <- seq(from = weekly_dates[1], 
                     to = weekly_dates[1] + days(length(daily_cases) - 1), 
                     by = "day")
  
  incidence_data <- data.frame(
    dates = dates_daily,
    I = pmax(daily_cases, 0.1)  # Avoid zeros
  )
  
  # Define configuration with pathogen-specific SI and window size
  window_size <- window
  config <- make_config(
    list(
      mean_si = si_mean,
      std_si = si_std,
      t_start = seq(2, length(daily_cases) - window_size),
      t_end = seq(2 + window_size, length(daily_cases))
      
      # 左对齐（前向窗口）：时间点t的Rt基于 [t, t+6] 的数据（window=6时）
    )
  )
  
  # Estimate Rt
  tryCatch({
    res <- estimate_R(
      incid = incidence_data,
      method = "parametric_si",
      config = config
    )
    
    # Extract Rt estimates (weekly aggregation)
    rt_daily <- data.frame(
      date = res$dates[res$R$t_start],  # 使用 t_start 作为时间点（左对齐）
      R_mean = res$R$`Mean(R)`,
      R_lower = res$R$`Quantile.0.025(R)`,
      R_upper = res$R$`Quantile.0.975(R)`
    )
    
    # Aggregate to weekly (matching original weekly dates)
    rt_daily$week <- ceiling(as.numeric(rt_daily$date - weekly_dates[1] + 1) / 7)
    
    rt_weekly <- rt_daily %>%
      filter(week <= length(weekly_dates)) %>%
      group_by(week) %>%
      summarise(
        R_mean = mean(R_mean, na.rm = TRUE),
        R_lower = mean(R_lower, na.rm = TRUE),
        R_upper = mean(R_upper, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Match to original weekly dates
    rt_weekly$date <- weekly_dates[rt_weekly$week]
    
    cat(sprintf("    Success! Computed %d weekly Rt estimates\n", nrow(rt_weekly)))
    return(rt_weekly)
    
  }, error = function(e) {
    cat(sprintf("    Warning: Cori Rt calculation failed: %s\n", e$message))
    return(data.frame(
      date = weekly_dates,
      R_mean = NA,
      R_lower = NA,
      R_upper = NA
    ))
  })
}

# Calculate Cori Rt for all strains with strain-specific SI
cori_rt_list <- list()
for (i in 1:N_strains) {
  strain <- strain_names[i]
  cori_rt_list[[strain]] <- calculate_cori_rt(
    weekly_cases = cases_matrix[, i],
    weekly_dates = fitting_data$date,
    population = hk_population,
    strain_name = strain,
    si_mean = si_params[[strain]]$mean_si,
    si_std = si_params[[strain]]$std_si,
    window = 6
  )
}

# Combine into single dataframe
df_cori_rt <- bind_rows(
  lapply(names(cori_rt_list), function(strain) {
    df <- cori_rt_list[[strain]]
    df$strain <- strain
    return(df)
  })
)

df_cori_rt$strain <- factor(df_cori_rt$strain, levels = strain_names)
df_cori_rt$scenario <- "Cori Rt"

cat("\nCori Rt calculation complete (using spline interpolation with LEFT-aligned window)!\n\n")

# ============================================================================
# 5. PREPARE WEEKLY INCIDENCE DATA (PER 10,000 POPULATION)
# ============================================================================
cat("=== Preparing weekly incidence data (per 10,000 population) ===\n")

# Calculate weekly incidence per 10,000 population
weekly_incidence_per10k <- (cases_matrix / hk_population) * 120000

# Create dataframe for visualization
df_incidence <- data.frame(
  date = rep(fitting_data$date, N_strains),
  strain = factor(rep(strain_names, each = T_weeks), levels = strain_names),
  incidence_per10k = as.vector(weekly_incidence_per10k)
)

cat("Weekly incidence data prepared.\n\n")

# ============================================================================
# 6. STAN MODEL COMPILATION
# ============================================================================
cat("=== Compiling Stan model ===\n")
model_file <- "/Users/chenjiaqi/SPH Dropbox/Jackie Chen/Shared RESV_HK_Typhoon/Fig_rds/6_subtypes_last_epi.stan"
model <- stan_model(model_file)

# ============================================================================
# 7. FIT MODEL WITH TYPHOON KNOT (at typhoon END)
# ============================================================================
cat("\n=== Fitting model WITH typhoon knot (at typhoon END) ===\n")
stan_data_typhoon <- list(
  T_weeks = T_weeks,
  N_strains = N_strains,
  cases = cases_matrix,
  num_knots = 15,
  spline_degree = 3,
  population = 7524000,
  child_ratio = 0.029,
  typhoon_start_week = typhoon_start_week,
  typhoon_duration_days = typhoon_duration_days,
  use_typhoon_knot = 1,
  typhoon_knot_multiplicity = 3
)

fit_typhoon <- sampling(
  model,
  data = stan_data_typhoon,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  thin = 1,
  cores = 4,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

saveRDS(fit_typhoon, 
        file = "/Users/chenjiaqi/SPH Dropbox/Jackie Chen/rds/fitting_result/fit_with_typhoon_knot.rds")

# ============================================================================
# 8. FIT MODEL WITHOUT TYPHOON KNOT
# ============================================================================
cat("\n=== Fitting model WITHOUT typhoon knot ===\n")
stan_data_normal <- list(
  T_weeks = T_weeks,
  N_strains = N_strains,
  cases = cases_matrix,
  num_knots = 15,
  spline_degree = 3,
  population = 7524000,
  child_ratio = 0.029,
  typhoon_start_week = typhoon_start_week,
  typhoon_duration_days = typhoon_duration_days,
  use_typhoon_knot = 0,
  typhoon_knot_multiplicity = 3
)

fit_normal <- sampling(
  model,
  data = stan_data_normal,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  thin = 1,
  cores = 4,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

saveRDS(fit_normal, 
        file = "/Users/chenjiaqi/SPH Dropbox/Jackie Chen/rds/fitting_result/fit_without_typhoon_knot.rds")

# ============================================================================
# 9. EXTRACT SAMPLES & CALCULATE STATISTICS
# ============================================================================
cat("\n=== Extracting samples and calculating statistics ===\n")

# Extract predictions
pred_typhoon <- rstan::extract(fit_typhoon, pars = "pred_cases")$pred_cases
pred_normal  <- rstan::extract(fit_normal, pars = "pred_cases")$pred_cases

# Extract R_eff
R_eff_typhoon <- rstan::extract(fit_typhoon, pars = "R_eff")$R_eff
R_eff_normal  <- rstan::extract(fit_normal, pars = "R_eff")$R_eff

# Extract typhoon knot information
typhoon_knot_coeff <- rstan::extract(fit_typhoon, pars = "typhoon_knot_coeff")$typhoon_knot_coeff
typhoon_knot_index <- rstan::extract(fit_typhoon, pars = "typhoon_knot_index")$typhoon_knot_index
typhoon_knot_position <- rstan::extract(fit_typhoon, pars = "typhoon_knot_position")$typhoon_knot_position

# Extract beta values calculated in Stan
beta_at_start_typhoon <- rstan::extract(fit_typhoon, pars = "beta_at_typhoon_start")$beta_at_typhoon_start
beta_at_end_typhoon <- rstan::extract(fit_typhoon, pars = "beta_at_typhoon_end")$beta_at_typhoon_end
ln_beta_ratio_typhoon <- rstan::extract(fit_typhoon, pars = "ln_beta_ratio")$ln_beta_ratio

beta_at_start_normal <- rstan::extract(fit_normal, pars = "beta_at_typhoon_start")$beta_at_typhoon_start
beta_at_end_normal <- rstan::extract(fit_normal, pars = "beta_at_typhoon_end")$beta_at_typhoon_end
ln_beta_ratio_normal <- rstan::extract(fit_normal, pars = "ln_beta_ratio")$ln_beta_ratio

# Calculate statistics
calc_stats <- function(x) {
  list(
    median = apply(x, c(2,3), median),
    lower  = apply(x, c(2,3), quantile, probs = 0.025),
    upper  = apply(x, c(2,3), quantile, probs = 0.975)
  )
}

stats_pred_typhoon <- calc_stats(pred_typhoon)
stats_pred_normal  <- calc_stats(pred_normal)
stats_reff_typhoon <- calc_stats(R_eff_typhoon)
stats_reff_normal  <- calc_stats(R_eff_normal)

# Typhoon knot coefficient summary
typhoon_coeff_summary <- data.frame(
  Strain   = strain_names,
  Mean     = apply(typhoon_knot_coeff, 2, mean),
  Lower_95 = apply(typhoon_knot_coeff, 2, quantile, probs = 0.025),
  Upper_95 = apply(typhoon_knot_coeff, 2, quantile, probs = 0.975)
)

cat("\n=== Typhoon Knot Coefficient Summary (at typhoon END) ===\n")
print(typhoon_coeff_summary)
cat("\nTyphoon knot position (normalized):", mean(typhoon_knot_position), "\n")
cat("Typhoon knot basis index:", round(mean(typhoon_knot_index)), "\n\n")

# ============================================================================
# 10. PREPARE VISUALIZATION DATA WITH SMOOTH CONNECTION
# ============================================================================
cat("=== Preparing visualization data with smooth connection ===\n")

# --- Cases Data Preparation ---
# Normal model: full period
df_cases_normal <- data.frame(
  date      = rep(fitting_data$date, N_strains),
  strain    = factor(rep(strain_names, each = T_weeks), levels = strain_names),
  observed  = as.vector(cases_matrix),
  predicted = as.vector(stats_pred_normal$median),
  lower     = as.vector(stats_pred_normal$lower),
  upper     = as.vector(stats_pred_normal$upper),
  scenario  = "Without Typhoon Knot"
)

# Typhoon model: only from typhoon start onwards
df_cases_typhoon <- data.frame(
  date      = rep(fitting_data$date, N_strains),
  strain    = factor(rep(strain_names, each = T_weeks), levels = strain_names),
  observed  = as.vector(cases_matrix),
  predicted = as.vector(stats_pred_typhoon$median),
  lower     = as.vector(stats_pred_typhoon$lower),
  upper     = as.vector(stats_pred_typhoon$upper),
  scenario  = "With Typhoon Knot"
)

# Filter typhoon model to start from typhoon_start_date
df_cases_typhoon_display <- df_cases_typhoon %>%
  filter(date >= typhoon_start_date)

# Create connection point at typhoon start date
connection_point_cases <- df_cases_normal %>%
  filter(date == typhoon_start_date) %>%
  mutate(scenario = "With Typhoon Knot")

# Combine: Normal (full) + Connection point + Typhoon (from start)
df_cases_combined <- bind_rows(
  df_cases_normal,
  connection_point_cases,
  df_cases_typhoon_display %>% filter(date > typhoon_start_date)
) %>%
  mutate(across(c(predicted, lower, upper), ~pmax(.x, 0)))

# --- R_eff Data Preparation ---
df_reff_normal <- data.frame(
  date     = rep(fitting_data$date, N_strains),
  strain   = factor(rep(strain_names, each = T_weeks), levels = strain_names),
  R_eff    = as.vector(stats_reff_normal$median),
  lower    = as.vector(stats_reff_normal$lower),
  upper    = as.vector(stats_reff_normal$upper),
  scenario = "Without Typhoon Knot"
)

df_reff_typhoon <- data.frame(
  date     = rep(fitting_data$date, N_strains),
  strain   = factor(rep(strain_names, each = T_weeks), levels = strain_names),
  R_eff    = as.vector(stats_reff_typhoon$median),
  lower    = as.vector(stats_reff_typhoon$lower),
  upper    = as.vector(stats_reff_typhoon$upper),
  scenario = "With Typhoon Knot"
)

# Filter typhoon model to start from typhoon_start_date
df_reff_typhoon_display <- df_reff_typhoon %>%
  filter(date >= typhoon_start_date)

# Create connection point for R_eff at typhoon start date
connection_point_reff <- df_reff_normal %>%
  filter(date == typhoon_start_date) %>%
  mutate(scenario = "With Typhoon Knot")

# Combine: Normal (full) + Connection point + Typhoon (from start) + Cori Rt
df_reff_combined <- bind_rows(
  df_reff_normal,
  connection_point_reff,
  df_reff_typhoon_display %>% filter(date > typhoon_start_date),
  df_cori_rt %>% 
    rename(R_eff = R_mean, lower = R_lower, upper = R_upper) %>%
    select(date, strain, R_eff, lower, upper, scenario)
) %>%
  mutate(lower = pmax(lower, 0))

# Set factor levels (3 scenarios now)
df_cases_combined$scenario <- factor(df_cases_combined$scenario, 
                                     levels = c("Without Typhoon Knot", "With Typhoon Knot"))
df_reff_combined$scenario  <- factor(df_reff_combined$scenario, 
                                     levels = c("Without Typhoon Knot", "With Typhoon Knot", "Cori Rt"))

cat("Connection point added at:", as.character(typhoon_start_date), "\n")
cat("Typhoon model now smoothly connects with Normal model\n")
cat("Cori Rt added as third comparison (spline interpolation, pathogen-specific SI, window=6)\n\n")

# ============================================================================
# 11. VISUALIZATION: CASES COMPARISON - INDIVIDUAL PLOTS FOR EACH STRAIN
# ============================================================================
cat("\n=== Creating Individual Figures for Each Strain ===\n")

for (strain_name in strain_names) {
  # Filter data for current strain
  df_strain <- df_cases_combined %>%
    filter(strain == strain_name)
  
  df_observed <- df_cases_normal %>%
    filter(strain == strain_name)
  
  # Create plot for this strain
  p <- ggplot(df_strain, aes(x = date)) +
    # Typhoon period shading
    annotate("rect", 
             xmin = typhoon_start_date, 
             xmax = typhoon_end_date,
             ymin = -Inf, ymax = Inf,
             fill = "gray60", alpha = 0.15) +
    
    # Typhoon START line (connection point)
    geom_vline(xintercept = as.numeric(typhoon_start_date),
               linetype = "dotted", color = "darkgreen", size = 0.5, alpha = 0.6) +
    
    # Typhoon END line (where knot is placed)
    geom_vline(xintercept = as.numeric(typhoon_end_date),
               linetype = "dashed", color = "red", size = 0.5, alpha = 0.6) +
    
    # Confidence intervals
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = scenario), alpha = 0.2) +
    
    # Prediction lines
    geom_line(aes(y = predicted, color = scenario), size = 1.1) +
    
    # Observed points
    geom_point(data = df_observed, aes(y = observed),
               color = "black", size = 0.8, alpha = 0.4) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_x_date(date_labels = "%m-%d", date_breaks = "1 month") +
    scale_color_manual(values = c("Without Typhoon Knot" = "#377EB8",
                                  "With Typhoon Knot" = "#E41A1C")) +
    scale_fill_manual(values = c("Without Typhoon Knot" = "#377EB8",
                                 "With Typhoon Knot" = "#E41A1C")) +
    labs(
      title = paste0("Epidemic Fitting: ", strain_name),
      subtitle = paste0("Baseline (blue): entire period | Typhoon model (red): from typhoon START with smooth connection\n",
                        "Green dotted line: connection point (", format(typhoon_start_date, "%m-%d"),
                        ") | Red dashed line: typhoon knot (", format(typhoon_end_date, "%m-%d"),
                        ") | Gray: typhoon period"),
      x = NULL,
      y = "Weekly Cases",
      color = "Model",
      fill = "Model"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      panel.border = element_rect(color = "gray30", fill = NA, size = 0.6),
      legend.position = "bottom",
      plot.subtitle = element_text(size = 9)
    )
  
  print(p)
}

cat("\nAll individual plots displayed successfully!\n")

# ============================================================================
# 12. VISUALIZATION: R_eff + INCIDENCE COMBINED PLOTS (ALIGNED TYPHOON SHADING)
# ============================================================================
cat("\n=== Creating Combined Figure: R_eff + Incidence for each strain ===\n")

reff_start_date <- as.Date("2025-06-01")

# Filter data for the visualization period
df_reff_filtered <- df_reff_combined %>%
  filter(date >= reff_start_date)

df_incidence_filtered <- df_incidence %>%
  filter(date >= reff_start_date)

# Define common x-axis limits to ensure alignment
x_axis_limits <- c(reff_start_date, max(df_reff_filtered$date))

# Create combined plots for each strain
for (strain_name in strain_names) {
  cat(sprintf("  Creating combined plot for %s...\n", strain_name))
  
  # Filter data for current strain
  df_reff_strain <- df_reff_filtered %>%
    filter(strain == strain_name)
  
  df_incidence_strain <- df_incidence_filtered %>%
    filter(strain == strain_name)
  
  # === TOP PANEL: R_eff ===
  p_reff <- ggplot(df_reff_strain, aes(x = date)) +
    # Typhoon period shading (ALIGNED)
    annotate("rect",
             xmin = typhoon_start_date,
             xmax = typhoon_end_date,
             ymin = -Inf, ymax = Inf,
             fill = "gray60", alpha = 0.15) +
    
    # Typhoon START line
    geom_vline(xintercept = as.numeric(typhoon_start_date),
               linetype = "dotted", color = "darkgreen", linewidth = 0.5, alpha = 0.6) +
    
    # Typhoon END line
    geom_vline(xintercept = as.numeric(typhoon_end_date),
               linetype = "dashed", color = "red", linewidth = 0.5, alpha = 0.6) +
    
    # R_eff = 1 reference line
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "gray50", linewidth = 0.5) +
    
    # Confidence intervals and lines
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = scenario), alpha = 0.15) +
    geom_line(aes(y = R_eff, color = scenario), linewidth = 1.1) +
    
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_x_date(
      limits = x_axis_limits,
      date_labels = "%m-%d", 
      date_breaks = "1 month",
      expand = c(0, 0)  # Remove padding
    ) +
    scale_color_manual(
      values = c(
        "Without Typhoon Knot" = "#377EB8",
        "With Typhoon Knot" = "#E41A1C",
        "Cori Rt" = "#4DAF4A"
      )
    ) +
    scale_fill_manual(
      values = c(
        "Without Typhoon Knot" = "#377EB8",
        "With Typhoon Knot" = "#E41A1C",
        "Cori Rt" = "#4DAF4A"
      )
    ) +
    labs(
      title = paste0(strain_name, ": Effective Reproduction Number (R_eff)"),
      x = NULL,
      y = expression(R[eff]),
      color = "Model",
      fill = "Model"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      panel.border = element_rect(color = "gray30", fill = NA, linewidth = 0.6),
      legend.position = "top",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(5, 5, 0, 5)
    )
  
  # === BOTTOM PANEL: Weekly Incidence (displayed as ILI+proxy) ===
  p_incidence <- ggplot(df_incidence_strain, aes(x = date, y = incidence_per10k)) +
    # Typhoon period shading (ALIGNED WITH TOP PANEL)
    annotate("rect",
             xmin = typhoon_start_date,
             xmax = typhoon_end_date,
             ymin = -Inf, ymax = Inf,
             fill = "gray60", alpha = 0.15) +
    
    # Typhoon START line
    geom_vline(xintercept = as.numeric(typhoon_start_date),
               linetype = "dotted", color = "darkgreen", linewidth = 0.5, alpha = 0.6) +
    
    # Typhoon END line
    geom_vline(xintercept = as.numeric(typhoon_end_date),
               linetype = "dashed", color = "red", linewidth = 0.5, alpha = 0.6) +
    
    # Bar chart
    geom_col(fill = "#2E86AB", alpha = 0.7, width = 6) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_x_date(
      limits = x_axis_limits,  # Same limits as top panel
      date_labels = "%m-%d", 
      date_breaks = "1 month",
      expand = c(0, 0)  # Remove padding
    ) +
    labs(
      title = NULL,
      x = "Date",
      y = "Weekly ILI+proxy"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      panel.border = element_rect(color = "gray30", fill = NA, linewidth = 0.6),
      plot.margin = margin(0, 5, 5, 5)
    )
  
  # Combine plots vertically with aligned x-axes
  p_combined <- p_reff / p_incidence + 
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
      subtitle = paste0("Starting from ", format(reff_start_date, "%Y-%m-%d"),
                        " | Blue: Normal | Red: Typhoon | Green: Cori Rt\n",
                        "Green dotted: typhoon start | Red dashed: typhoon end | Gray: typhoon period")
    )
  
  print(p_combined)
}


# ============================================================================
# 13. BETA RATIO ANALYSIS (CALCULATED IN STAN - CORRECT METHOD)
# ============================================================================
cat("\n\n")
cat("============================================================================\n")
cat("BETA RATIO ANALYSIS: ln(Beta_End / Beta_Start)\n")
cat("Calculated in Stan's generated quantities for CORRECT confidence intervals\n")
cat("============================================================================\n\n")

# Function to format output
format_result <- function(samples) {
  mean_val <- mean(samples)
  lower_val <- quantile(samples, probs = 0.025)
  upper_val <- quantile(samples, probs = 0.975)
  sprintf("%.4f (%.4f, %.4f)", mean_val, lower_val, upper_val)
}

# A. TYPHOON MODEL: Beta change from start to end
cat("=== A. TYPHOON MODEL: Beta change from start to end ===\n\n")

typhoon_beta_summary <- data.frame(
  Strain = strain_names,
  Beta_Start = character(N_strains),
  Beta_End = character(N_strains),
  Ln_Ratio = character(N_strains),
  Exp_Ratio = character(N_strains),
  Percent_Change = character(N_strains),
  stringsAsFactors = FALSE
)

for (i in 1:N_strains) {
  typhoon_beta_summary$Beta_Start[i] <- format_result(beta_at_start_typhoon[, i])
  typhoon_beta_summary$Beta_End[i] <- format_result(beta_at_end_typhoon[, i])
  typhoon_beta_summary$Ln_Ratio[i] <- format_result(ln_beta_ratio_typhoon[, i])
  
  exp_ratio <- exp(ln_beta_ratio_typhoon[, i])
  typhoon_beta_summary$Exp_Ratio[i] <- format_result(exp_ratio)
  
  pct_change <- (exp_ratio - 1) * 100
  typhoon_beta_summary$Percent_Change[i] <- format_result(pct_change)
}

print(typhoon_beta_summary)

# B. NORMAL MODEL: Beta change from start to end
cat("\n\n=== B. NORMAL MODEL: Beta change from start to end ===\n\n")

normal_beta_summary <- data.frame(
  Strain = strain_names,
  Beta_Start = character(N_strains),
  Beta_End = character(N_strains),
  Ln_Ratio = character(N_strains),
  Exp_Ratio = character(N_strains),
  Percent_Change = character(N_strains),
  stringsAsFactors = FALSE
)

for (i in 1:N_strains) {
  normal_beta_summary$Beta_Start[i] <- format_result(beta_at_start_normal[, i])
  normal_beta_summary$Beta_End[i] <- format_result(beta_at_end_normal[, i])
  normal_beta_summary$Ln_Ratio[i] <- format_result(ln_beta_ratio_normal[, i])
  
  exp_ratio <- exp(ln_beta_ratio_normal[, i])
  normal_beta_summary$Exp_Ratio[i] <- format_result(exp_ratio)
  
  pct_change <- (exp_ratio - 1) * 100
  normal_beta_summary$Percent_Change[i] <- format_result(pct_change)
}

print(normal_beta_summary)

# C. COMPARISON: Difference in beta change between models
cat("\n\n=== C. COMPARISON: Difference in beta change between models ===\n\n")
cat("ln(Beta_Ratio_Typhoon) - ln(Beta_Ratio_Normal)\n")
cat("Positive values indicate GREATER increase in typhoon model\n\n")

comparison_summary <- data.frame(
  Strain = strain_names,
  Ln_Ratio_Typhoon = character(N_strains),
  Ln_Ratio_Normal = character(N_strains),
  Difference = character(N_strains),
  stringsAsFactors = FALSE
)

for (i in 1:N_strains) {
  comparison_summary$Ln_Ratio_Typhoon[i] <- format_result(ln_beta_ratio_typhoon[, i])
  comparison_summary$Ln_Ratio_Normal[i] <- format_result(ln_beta_ratio_normal[, i])
  
  difference <- ln_beta_ratio_typhoon[, i] - ln_beta_ratio_normal[, i]
  comparison_summary$Difference[i] <- format_result(difference)
}

print(comparison_summary)

# D. VISUALIZATION: Beta ratio comparison
cat("\n\n=== Creating visualization: Beta ratio comparison ===\n")

beta_ratio_df <- data.frame(
  Strain = rep(strain_names, 2),
  Model = rep(c("Typhoon", "Normal"), each = N_strains),
  Mean = c(apply(ln_beta_ratio_typhoon, 2, mean),
           apply(ln_beta_ratio_normal, 2, mean)),
  Lower = c(apply(ln_beta_ratio_typhoon, 2, quantile, probs = 0.025),
            apply(ln_beta_ratio_normal, 2, quantile, probs = 0.025)),
  Upper = c(apply(ln_beta_ratio_typhoon, 2, quantile, probs = 0.975),
            apply(ln_beta_ratio_normal, 2, quantile, probs = 0.975))
)

beta_ratio_df$Strain <- factor(beta_ratio_df$Strain, levels = strain_names)

p_beta_ratio <- ggplot(beta_ratio_df, aes(x = Strain, y = Mean, color = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(ymin = Lower, ymax = Upper),
                  position = position_dodge(width = 0.3),
                  size = 0.8) +
  scale_color_manual(values = c("Typhoon" = "#E41A1C", "Normal" = "#377EB8")) +
  labs(
    title = "Beta Change: ln(Beta_End / Beta_Start)",
    subtitle = "Comparison between Typhoon and Normal models",
    x = "Pathogen",
    y = "ln(Beta Ratio)",
    color = "Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.6),
    legend.position = "bottom"
  )

print(p_beta_ratio)

# ============================================================================
# 14. ADDITIONAL ANALYSIS: CROSS-MODEL RATIO
# ============================================================================
cat("\n\n")
cat("============================================================================\n")
cat("ADDITIONAL ANALYSIS: ln(Beta_End_Typhoon / Beta_Start_Normal)\n")
cat("============================================================================\n\n")
cat("WARNING: This ratio mixes different models and time points:\n")
cat("  - Numerator:   Beta from TYPHOON model at typhoon END (", 
    as.character(typhoon_end_date), ")\n", sep="")
cat("  - Denominator: Beta from NORMAL model at typhoon START (", 
    as.character(typhoon_start_date), ")\n\n", sep="")

# Calculate cross-model ratio
cross_model_summary <- data.frame(
  Strain = strain_names,
  Beta_End_Typhoon = character(N_strains),
  Beta_Start_Normal = character(N_strains),
  Ln_Ratio = character(N_strains),
  Multiplicative_Change = character(N_strains),
  stringsAsFactors = FALSE
)

for (i in 1:N_strains) {
  # Extract samples
  beta_end_t <- beta_at_end_typhoon[, i]
  beta_start_n <- beta_at_start_normal[, i]
  
  # Format individual betas
  cross_model_summary$Beta_End_Typhoon[i] <- format_result(beta_end_t)
  cross_model_summary$Beta_Start_Normal[i] <- format_result(beta_start_n)
  
  # Calculate ln(ratio) - pairing samples from independent chains
  ln_ratio <- log(beta_end_t / beta_start_n)
  cross_model_summary$Ln_Ratio[i] <- format_result(ln_ratio)
  
  # Multiplicative change
  mult_change <- exp(ln_ratio)
  cross_model_summary$Multiplicative_Change[i] <- sprintf("%.2fx (%.2fx, %.2fx)",
                                                          mean(mult_change),
                                                          quantile(mult_change, 0.025),
                                                          quantile(mult_change, 0.975))
}

cat("=== CROSS-MODEL RATIO RESULTS ===\n\n")
print(cross_model_summary)

cat("\n\nINTERPRETATION:\n")
cat("This ratio quantifies the combined effect of:\n")
cat("  1. Model difference (Typhoon vs Normal structure)\n")
cat("  2. Temporal change (7-day typhoon period)\n")
cat("Positive Ln_Ratio: Beta increased from Normal-Start to Typhoon-End\n")
cat("Multiplicative_Change: Fold-change in beta\n\n")
cat("CAUTION: Confidence intervals combine uncertainty from:\n")
cat("  - Two independent MCMC chains (no pairing)\n")
cat("  - Model structural differences\n")
cat("  - Temporal dynamics\n")

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n\n")
cat("============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================================\n")
cat("\nKey outputs:\n")
cat("1. Individual epidemic curves for each strain (with smooth connection)\n")
cat("2. Combined R_eff + Weekly Incidence plots for each strain:\n")
cat("   - Top panel: R_eff comparison (3 curves)\n")
cat("     * Without Typhoon Knot (Blue)\n")
cat("     * With Typhoon Knot (Red)\n")
cat("     * Cori Rt (Green)\n")
cat("   - Bottom panel: Weekly incidence per 10,000 population (bar chart)\n")
cat("3. Pathogen-specific SI parameters:\n")
for (strain in strain_names) {
  cat(sprintf("   - %s: SI mean=%.1f days, SD=%.2f days\n", 
              strain, si_params[[strain]]$mean_si, si_params[[strain]]$std_si))
}
cat("4. Within-model beta change analysis (Typhoon vs Normal)\n")
cat("5. Cross-model beta ratio: ln(Beta_End_Typhoon / Beta_Start_Normal)\n")
cat("6. Beta ratio comparison visualization\n")
cat("\nVisualization features:\n")
cat("  - Green dotted line: Connection point at typhoon START (", 
    as.character(typhoon_start_date), ")\n", sep="")
cat("  - Red dashed line: Typhoon knot at typhoon END (", 
    as.character(typhoon_end_date), ")\n", sep="")
cat("  - Gray shading: Typhoon period\n")
cat("  - Spline interpolation: monoH.FC (shape-preserving)\n")
cat("  - Window alignment: LEFT-aligned (forward window)\n")
cat("============================================================================\n")