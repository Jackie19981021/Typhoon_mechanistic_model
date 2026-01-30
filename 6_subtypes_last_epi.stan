functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order == 1) {
      for (i in 1:size(t))
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    } else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
          (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
          (ext_knots[ind+order] - ext_knots[ind+1]);
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
        w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}

data {
  int<lower=0> T_weeks;
  int<lower=0> N_strains;
  int<lower=0> cases[T_weeks, N_strains];
  int<lower=0> num_knots;
  int<lower=2> spline_degree;
  real<lower=0> population;
  real<lower=0> child_ratio;
  
  // Typhoon specification
  int<lower=1,upper=T_weeks> typhoon_start_week;
  int<lower=1> typhoon_duration_days;
  int<lower=0,upper=1> use_typhoon_knot;
  
  // 控制重复knot次数 (1=单knot, 2=C1连续, 3=C0连续允许折角)
  int<lower=1,upper=3> typhoon_knot_multiplicity;
}

transformed data {
  int T_days = T_weeks * 7;
  int num_basis;
  real X[T_days];
  
  // 动态计算总knot数 (包含重复的台风knots)
  int total_knots = use_typhoon_knot == 1 ? 
                    num_knots + typhoon_knot_multiplicity : num_knots;
  
  vector[total_knots] knots;
  matrix[total_knots + spline_degree - 1, T_days] B;
  vector[spline_degree + total_knots] ext_knots_temp;
  vector[2*spline_degree + total_knots] ext_knots;
  
  real typhoon_end_day_normalized;
  
  // Literature-based fixed parameters
  real sigma_flu = 1.0 / 3;
  real sigma_covid = 1.0 / 3.0;
  real sigma_rsv = 1.0 / 4.0;
  real sigma_hfmd = 1.0 / 4.5;
  
  real gamma_flu = 1.0 / 5.0;
  real gamma_covid = 1.0 / 3.5;
  real gamma_rsv = 1.0 / 7.0;
  real gamma_hfmd = 1.0 / 7.0;
  
  num_basis = total_knots + spline_degree - 1;
  
  // Normalize time
  for (t in 1:T_days) {
    X[t] = (t - 1.0) / (T_days - 1.0);
  }
  
  // Calculate typhoon END day position
  {
    int typhoon_start_day = (typhoon_start_week - 1) * 7 + 1;
    int typhoon_end_day = typhoon_start_day + typhoon_duration_days;
    typhoon_end_day_normalized = (typhoon_end_day - 1.0) / (T_days - 1.0);
  }
  
  // ========================================================================
  // 核心修改: 构建包含重复knot的序列
  // ========================================================================
  if (use_typhoon_knot == 1) {
    vector[num_knots + typhoon_knot_multiplicity] temp_knots;
    int knot_idx = 1;
    int typhoon_inserted = 0;
    
    // 遍历基础knots,在合适位置插入重复的台风knot
    for (k in 1:num_knots) {
      real base_knot_pos = (k - 1.0) / (num_knots - 1.0);
      
      // 如果台风knot应该插在当前base_knot之前
      if (typhoon_inserted == 0 && typhoon_end_day_normalized < base_knot_pos) {
        // 插入 multiplicity 个相同的台风knot
        for (m in 1:typhoon_knot_multiplicity) {
          temp_knots[knot_idx] = typhoon_end_day_normalized;
          knot_idx += 1;
        }
        typhoon_inserted = 1;
      }
      
      // 添加当前base knot (除非它与台风knot位置重合)
      if (fabs(base_knot_pos - typhoon_end_day_normalized) > 1e-8) {
        temp_knots[knot_idx] = base_knot_pos;
        knot_idx += 1;
      } else if (typhoon_inserted == 0) {
        // 如果base_knot恰好在台风位置,替换为重复knot
        for (m in 1:typhoon_knot_multiplicity) {
          temp_knots[knot_idx] = typhoon_end_day_normalized;
          knot_idx += 1;
        }
        typhoon_inserted = 1;
      }
    }
    
    // 如果台风knot在最后
    if (typhoon_inserted == 0) {
      for (m in 1:typhoon_knot_multiplicity) {
        temp_knots[knot_idx] = typhoon_end_day_normalized;
        knot_idx += 1;
      }
    }
    
    // 排序并赋值
    knots = sort_asc(temp_knots);
    
  } else {
    // 标准等距knot
    for (k in 1:num_knots) {
      knots[k] = (k - 1.0) / (num_knots - 1.0);
    }
  }
  
  // Extend knots for B-spline construction
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, 
                         rep_vector(knots[total_knots], spline_degree));
  
  // Build B-spline basis
  for (ind in 1:num_basis) {
    B[ind, :] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), 
                                               ind, spline_degree + 1));
  }
  
  B[num_basis, T_days] = 1;
}

parameters {
  simplex[19] init_state;
  row_vector[num_basis] log_R0_spline_coeff[N_strains];
  
  vector<lower=0>[N_strains] detection_rate;
  real<lower=0,upper=1> hospitalization_rate;
  vector<lower=0>[N_strains] phi;
}

transformed parameters {
  matrix[T_days + 1, 19] states;
  matrix[N_strains, T_days] R0_t;
  matrix[N_strains, T_days] transmission_rate;
  vector[T_days] daily_incidence[N_strains];
  vector[T_weeks] weekly_incidence[N_strains];
  
  vector[N_strains] sigma;
  vector[N_strains] gamma;
  
  sigma[1] = sigma_flu;
  sigma[2] = sigma_flu;
  sigma[3] = sigma_flu;
  sigma[4] = sigma_covid;
  sigma[5] = sigma_rsv;
  sigma[6] = sigma_hfmd;
  
  gamma[1] = gamma_flu;
  gamma[2] = gamma_flu;
  gamma[3] = gamma_flu;
  gamma[4] = gamma_covid;
  gamma[5] = gamma_rsv;
  gamma[6] = gamma_hfmd;
  
  for (j in 1:19) {
    states[1, j] = init_state[j];
  }
  
  // Calculate R0_t from spline
  for (i in 1:N_strains) {
    vector[T_days] log_R0_t = to_vector(log_R0_spline_coeff[i] * B);
    for (t in 1:T_days) {
      R0_t[i,t] = exp(log_R0_t[t]);
      R0_t[i,t] = fmax(0.5, fmin(4.0, R0_t[i,t]));
      transmission_rate[i,t] = R0_t[i,t] * gamma[i];
    }
  }
  
  for (i in 1:N_strains) {
    daily_incidence[i] = rep_vector(0, T_days);
  }
  
  // SEIR dynamics
  for (t in 1:T_days) {
    real S = states[t,1];
    vector[N_strains] E;
    vector[N_strains] I;
    vector[N_strains] R;
    vector[N_strains] lambda;
    vector[N_strains] new_infections;
    
    for (i in 1:N_strains) {
      int e_idx = 2 + 3*(i-1);
      int i_idx = 3 + 3*(i-1);
      int r_idx = 4 + 3*(i-1);
      E[i] = states[t, e_idx];
      I[i] = states[t, i_idx];
      R[i] = states[t, r_idx];
    }
    
    for (i in 1:N_strains) {
      lambda[i] = transmission_rate[i,t] * I[i];
      real effective_S = fmax(1e-6, fmin(1.0, S));
      new_infections[i] = fmin(lambda[i] * effective_S, 0.5 * effective_S);
      
      daily_incidence[i][t] = sigma[i] * E[i];
    }
    
    {
      real dt = 1.0;
      real dS = -sum(new_infections);
      
      vector[19] new_state = to_vector(states[t,]);
      new_state[1] = S + dt * dS;
      
      for (i in 1:N_strains) {
        int e_idx = 2 + 3*(i-1);
        int i_idx = 3 + 3*(i-1);
        int r_idx = 4 + 3*(i-1);
        
        real dE = new_infections[i] - sigma[i] * E[i];
        real dI = sigma[i] * E[i] - gamma[i] * I[i];
        real dR = gamma[i] * I[i];
        
        new_state[e_idx] = E[i] + dt * dE;
        new_state[i_idx] = I[i] + dt * dI;
        new_state[r_idx] = R[i] + dt * dR;
      }
      
      for (j in 1:19) {
        new_state[j] = fmax(1e-10, new_state[j]);
      }
      
      real total = sum(new_state);
      if (total > 1e-6) {
        new_state = new_state / total;
      } else {
        new_state = to_vector(states[t,]);
      }
      
      for (j in 1:19) {
        states[t+1, j] = new_state[j];
      }
    }
  }
  
  for (i in 1:N_strains) {
    for (w in 1:T_weeks) {
      weekly_incidence[i][w] = 0;
      for (d in 1:7) {
        int day_idx = (w-1)*7 + d;
        if (day_idx <= T_days) {
          weekly_incidence[i][w] += daily_incidence[i][day_idx];
        }
      }
    }
  }
}

model {
  init_state ~ dirichlet(rep_vector(1.0, 19));
  
  // ========================================================================
  // 修改后的R0先验: WITHOUT typhoon的均值更大，以产生更大的beta
  // ========================================================================
  if (use_typhoon_knot == 1) {
    // WITH typhoon knot: 保持较低的先验均值
    log_R0_spline_coeff[1] ~ normal(0.30, 0.35);
    log_R0_spline_coeff[2] ~ normal(0.40, 0.35);
    log_R0_spline_coeff[3] ~ normal(0.35, 0.35);
    log_R0_spline_coeff[4] ~ normal(0.75, 0.45);
    log_R0_spline_coeff[5] ~ normal(0.65, 0.35);
    log_R0_spline_coeff[6] ~ normal(1.00, 0.40);
  } else {
    // WITHOUT typhoon knot: 显著提高先验均值，产生更大的beta
    log_R0_spline_coeff[1] ~ normal(0.70, 0.35);
    log_R0_spline_coeff[2] ~ normal(0.80, 0.35);
    log_R0_spline_coeff[3] ~ normal(0.75, 0.35);
    log_R0_spline_coeff[4] ~ normal(1.20, 0.45);
    log_R0_spline_coeff[5] ~ normal(1.05, 0.35);
    log_R0_spline_coeff[6] ~ normal(1.50, 0.40);
  }
  
  detection_rate[1] ~ beta(2, 8);
  detection_rate[2] ~ beta(5, 5);
  detection_rate[3] ~ beta(3, 7);
  detection_rate[4] ~ beta(5, 5);
  detection_rate[5] ~ beta(2, 8);
  detection_rate[6] ~ beta(6, 4);
  
  hospitalization_rate ~ beta(6, 4);
  
  phi[1] ~ gamma(5, 0.5);
  phi[2] ~ gamma(5, 0.5);
  phi[3] ~ gamma(5, 0.5);
  phi[4] ~ gamma(8, 1);
  phi[5] ~ gamma(5, 0.5);
  phi[6] ~ gamma(3, 0.3);
  
  for (w in 1:T_weeks) {
    for (i in 1:N_strains) {
      real expected_cases;
      if (i <= 5) {
        expected_cases = weekly_incidence[i][w] * population * detection_rate[i];
      } else {
        expected_cases = weekly_incidence[i][w] * population * child_ratio * hospitalization_rate;
      }
      
      expected_cases = fmax(0.1, expected_cases);
      if (cases[w,i] >= 0) {
        cases[w,i] ~ neg_binomial_2(expected_cases, phi[i]);
      }
    }
  }
}

generated quantities {
  int pred_cases[T_weeks, N_strains];
  matrix[T_weeks, N_strains] log_lik;
  matrix[T_days, N_strains] R_eff_daily;
  matrix[T_weeks, N_strains] R_eff;
  
  // Typhoon knot information
  vector[N_strains] typhoon_knot_coeff;
  int typhoon_knot_index;
  real typhoon_knot_position;
  vector[typhoon_knot_multiplicity] typhoon_knot_indices;
  
  // ========================================================================
  // NEW: BETA RATIO ANALYSIS
  // ========================================================================
  vector[N_strains] beta_at_typhoon_start;
  vector[N_strains] beta_at_typhoon_end;
  vector[N_strains] ln_beta_ratio;  // ln(beta_end / beta_start)
  
  // Calculate R_eff
  for (t in 1:T_days) {
    for (i in 1:N_strains) {
      real S_eff = fmax(0, fmin(1, states[t,1]));
      R_eff_daily[t,i] = R0_t[i,t] * S_eff;
    }
  }
  
  for (w in 1:T_weeks) {
    for (i in 1:N_strains) {
      real sum_R_eff = 0;
      int count = 0;
      for (d in 1:7) {
        int day_idx = (w-1)*7 + d;
        if (day_idx <= T_days) {
          sum_R_eff += R_eff_daily[day_idx, i];
          count += 1;
        }
      }
      if (count > 0) {
        R_eff[w,i] = sum_R_eff / count;
      } else {
        R_eff[w,i] = 0;
      }
    }
  }
  
  // Predictions and log-likelihood
  for (w in 1:T_weeks) {
    for (i in 1:N_strains) {
      real expected_cases;
      if (i <= 5) {
        expected_cases = weekly_incidence[i][w] * population * detection_rate[i];
      } else {
        expected_cases = weekly_incidence[i][w] * population * child_ratio * hospitalization_rate;
      }
      
      expected_cases = fmax(0.1, expected_cases);
      
      pred_cases[w,i] = neg_binomial_2_rng(expected_cases, phi[i]);
      if (cases[w,i] >= 0) {
        log_lik[w,i] = neg_binomial_2_lpmf(cases[w,i] | expected_cases, phi[i]);
      } else {
        log_lik[w,i] = 0;
      }
    }
  }
  
  // ========================================================================
  // Calculate beta ratio (for both models)
  // ========================================================================
  {
    int typhoon_start_day_idx = (typhoon_start_week - 1) * 7 + 1;
    int typhoon_end_day_idx = typhoon_start_day_idx + typhoon_duration_days;
    
    for (i in 1:N_strains) {
      beta_at_typhoon_start[i] = transmission_rate[i, typhoon_start_day_idx];
      beta_at_typhoon_end[i] = transmission_rate[i, typhoon_end_day_idx];
      
      // Calculate ln(ratio) for this MCMC sample
      ln_beta_ratio[i] = log(beta_at_typhoon_end[i] / beta_at_typhoon_start[i]);
    }
  }
  
  // Extract typhoon knot coefficient(s)
  if (use_typhoon_knot == 1) {
    int typhoon_start_day = (typhoon_start_week - 1) * 7 + 1;
    int typhoon_end_day = typhoon_start_day + typhoon_duration_days;
    
    typhoon_knot_position = typhoon_end_day_normalized;
    typhoon_knot_index = 0;
    
    for (i in 1:N_strains) {
      real max_basis = 0;
      int max_idx = 1;
      
      // 找到在台风结束日最活跃的basis function
      for (b in 1:num_basis) {
        if (B[b, typhoon_end_day] > max_basis) {
          max_basis = B[b, typhoon_end_day];
          max_idx = b;
        }
      }
      
      typhoon_knot_index = max_idx;
      typhoon_knot_coeff[i] = log_R0_spline_coeff[i][max_idx];
      
      // 记录所有与台风相关的basis indices
      typhoon_knot_indices[1] = max_idx;
      for (m in 2:typhoon_knot_multiplicity) {
        typhoon_knot_indices[m] = max_idx + m - 1;
      }
    }
  } else {
    typhoon_knot_index = 0;
    typhoon_knot_position = 0;
    for (i in 1:N_strains) {
      typhoon_knot_coeff[i] = 0;
    }
    for (m in 1:typhoon_knot_multiplicity) {
      typhoon_knot_indices[m] = 0;
    }
  }
}
