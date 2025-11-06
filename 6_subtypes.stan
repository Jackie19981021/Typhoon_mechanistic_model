functions {
  // B-spline construction function
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
  int<lower=0> T_weeks_forecast;
  int<lower=0> N_strains;
  int<lower=0> cases[T_weeks, N_strains];
  int<lower=0> num_knots;
  int<lower=2> spline_degree;
  real<lower=0> population;
  int<lower=0> typhoon_weeks;
}

transformed data {
  int T_weeks_total = T_weeks + T_weeks_forecast;
  int T_days = T_weeks * 7;
  int T_days_forecast = T_weeks_forecast * 7;
  int T_days_total = T_weeks_total * 7;
  int num_basis = num_knots + spline_degree - 1;
  real X[T_days];
  vector[num_knots] knots;
  matrix[num_basis, T_days] B;
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots;
  
  for (t in 1:T_days) {
    X[t] = (t - 1.0) / (T_days - 1.0);
  }
  
  for (k in 1:num_knots) {
    knots[k] = (k - 1.0) / (num_knots - 1.0);
  }
  
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  
  for (ind in 1:num_basis) {
    B[ind, :] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));
  }
  
  B[num_basis, T_days] = 1;
}

parameters {
  simplex[19] init_state;
  row_vector[num_basis] log_R0_spline_coeff[N_strains];
  real<lower=0,upper=1> cross_immunity_flu;
  real<lower=0,upper=1> cross_immunity_flu_covid;
  real<lower=0,upper=1> cross_immunity_flu_rsv;
  real<lower=0,upper=1> cross_immunity_covid_rsv;
  vector<lower=0.05,upper=1>[N_strains] sigma;
  vector<lower=0,upper=1>[N_strains] gamma;
  vector<lower=0,upper=0.05>[N_strains] mu;
  vector<lower=0,upper=1>[5] detection_rate;
  real<lower=0,upper=0.01> hospitalization_rate;
  vector<lower=0>[N_strains] phi;
}

transformed parameters {
  matrix[T_days + 1, 19] states;
  matrix[N_strains, T_days] R0_t;
  matrix[N_strains, T_days] transmission_rate;
  vector[T_days] daily_incidence[N_strains];
  vector[T_weeks] weekly_incidence[N_strains];
  matrix[N_strains, N_strains] cross_immunity;
  
  for (i in 1:N_strains) {
    for (j in 1:N_strains) {
      if (i == j) {
        cross_immunity[i,j] = 1.0;
      } else if (i <= 3 && j <= 3) {
        cross_immunity[i,j] = cross_immunity_flu;
      } else if ((i <= 3 && j == 4) || (i == 4 && j <= 3)) {
        cross_immunity[i,j] = cross_immunity_flu_covid;
      } else if ((i <= 3 && j == 5) || (i == 5 && j <= 3)) {
        cross_immunity[i,j] = cross_immunity_flu_rsv;
      } else if ((i == 4 && j == 5) || (i == 5 && j == 4)) {
        cross_immunity[i,j] = cross_immunity_covid_rsv;
      } else {
        cross_immunity[i,j] = 0.0;
      }
    }
  }
  
  for (j in 1:19) {
    states[1, j] = init_state[j];
  }
  
  for (i in 1:N_strains) {
    vector[T_days] log_R0_t = to_vector(log_R0_spline_coeff[i] * B);
    
    for (t in 1:T_days) {
      R0_t[i,t] = exp(log_R0_t[t]);
      R0_t[i,t] = fmax(0.5, fmin(10.0, R0_t[i,t]));
      transmission_rate[i,t] = R0_t[i,t] * gamma[i];
    }
  }
  
  for (i in 1:N_strains) {
    daily_incidence[i] = rep_vector(0, T_days);
  }
  
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
      
      real effective_S = S;
      
      for (j in 1:N_strains) {
        if (j != i) {
          real immunity_effect = R[j] * (1 - cross_immunity[j,i]);
          effective_S += immunity_effect;
        }
      }
      
      effective_S = fmax(1e-6, fmin(1.0, effective_S));
      new_infections[i] = fmin(lambda[i] * effective_S, 0.5 * effective_S);
      
      daily_incidence[i][t] = sigma[i] * E[i];
    }
    
    {
      real dt = 1.0;
      real dS = -sum(new_infections) + sum(mu .* R);
      
      vector[19] new_state = to_vector(states[t,]);
      new_state[1] = S + dt * dS;
      
      for (i in 1:N_strains) {
        int e_idx = 2 + 3*(i-1);
        int i_idx = 3 + 3*(i-1);
        int r_idx = 4 + 3*(i-1);
        
        real dE = new_infections[i] - sigma[i] * E[i];
        real dI = sigma[i] * E[i] - gamma[i] * I[i];
        real dR = gamma[i] * I[i] - mu[i] * R[i];
        
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
  
  log_R0_spline_coeff[1] ~ normal(0.55, 0.4);
  log_R0_spline_coeff[2] ~ normal(0.65, 0.4);
  log_R0_spline_coeff[3] ~ normal(0.57, 0.4);
  log_R0_spline_coeff[4] ~ normal(0.83, 0.5);
  log_R0_spline_coeff[5] ~ normal(0.75, 0.4);
  log_R0_spline_coeff[6] ~ normal(1.04, 0.45);
  
  cross_immunity_flu ~ beta(2, 6);
  cross_immunity_flu_covid ~ beta(1, 5);
  cross_immunity_flu_rsv ~ beta(1, 5);
  cross_immunity_covid_rsv ~ beta(1, 6);
  
  for (i in 1:N_strains) {
    if (i <= 3) {
      sigma[i] ~ beta(2, 2);
      gamma[i] ~ beta(2, 10);
    } else if (i == 4) {
      sigma[i] ~ beta(2, 5);
      gamma[i] ~ beta(2, 15);
    } else if (i == 5) {
      sigma[i] ~ beta(2, 3);
      gamma[i] ~ beta(2, 12);
    } else {
      sigma[i] ~ beta(2, 3);
      gamma[i] ~ beta(2, 8);
    }
  }
  
  mu[1] ~ exponential(500);
  mu[2] ~ exponential(300);
  mu[3] ~ exponential(400);
  mu[4] ~ exponential(800);
  mu[5] ~ exponential(400);
  mu[6] ~ exponential(200);
  
  detection_rate[1] ~ beta(2, 8);
  detection_rate[2] ~ beta(3, 7);
  detection_rate[3] ~ beta(3, 7);
  detection_rate[4] ~ beta(5, 5);
  detection_rate[5] ~ beta(2, 8);
  hospitalization_rate ~ beta(1, 99);
  
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
        expected_cases = weekly_incidence[i][w] * population * hospitalization_rate;
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
  
  matrix[T_days_total, N_strains] R_eff_daily;
  matrix[T_weeks_total, N_strains] R_eff;
  
  matrix[T_weeks_total, N_strains] R_eff_scenarios[20];
  
  matrix[T_days_total + 1, 19] states_forecast;
  vector[T_days_forecast] daily_incidence_forecast[N_strains];
  vector[T_weeks_forecast] weekly_incidence_forecast[N_strains];
  int forecast_cases[T_weeks_forecast, N_strains];
  
  real typhoon_effects[20] = {0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.1, 0.2, 0.4, 0.6, 0.8, 0.1, 0.2, 0.4, 0.6, 0.8, 0.6, 0.6, 0.6, 0.6};
  int typhoon_days_arr[20] = {0, 3, 3, 3, 3, 3, 7, 7, 7, 7, 7, 10, 10, 10, 10, 10, 7, 7, 7, 7};
  int shift_arr[20] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, -5, 3, 5};
  int forecast_cases_typhoon[20, T_weeks_forecast, N_strains];
  
  // NEW: Store raw case counts for proper error bar calculation
  real cases_typhoon_period[20, N_strains];
  real cases_total_period[20, N_strains];
  real cases_baseline_typhoon_period[N_strains];
  real cases_baseline_total_period[N_strains];
  
  // NEW: Store reductions as percentages for each draw
  real reduction_typhoon_period[20, N_strains];
  real reduction_total_period[20, N_strains];
  
  // NEW: Store average weekly cases for period comparison
  real avg_weekly_typhoon[20, N_strains];
  real avg_weekly_recovery[20, N_strains];
  
  // 1. Calculate historical R_eff
  for (t in 1:T_days) {
    for (i in 1:N_strains) {
      real S_eff = states[t,1];
      
      for (j in 1:N_strains) {
        if (j != i) {
          int r_idx = 4 + 3*(j-1);
          real immunity_contribution = states[t, r_idx] * (1 - cross_immunity[j,i]);
          S_eff += immunity_contribution;
        }
      }
      
      S_eff = fmax(0, fmin(1, S_eff));
      R_eff_daily[t,i] = R0_t[i,t] * S_eff;
    }
  }
  
  // 2. Forecast state simulation (baseline)
  for (t in 1:(T_days + 1)) {
    for (j in 1:19) {
      states_forecast[t, j] = states[t, j];
    }
  }
  
  matrix[N_strains, T_days_forecast] R0_forecast;
  matrix[N_strains, T_days_forecast] transmission_rate_forecast;
  for (i in 1:N_strains) {
    real last_R0 = R0_t[i, T_days];
    for (t in 1:T_days_forecast) {
      R0_forecast[i, t] = last_R0;
      transmission_rate_forecast[i, t] = last_R0 * gamma[i];
    }
  }
  
  for (i in 1:N_strains) {
    daily_incidence_forecast[i] = rep_vector(0, T_days_forecast);
  }
  
  for (t_rel in 1:T_days_forecast) {
    int t = T_days + t_rel;
    real S = states_forecast[t, 1];
    vector[N_strains] E;
    vector[N_strains] I;
    vector[N_strains] R;
    vector[N_strains] lambda;
    vector[N_strains] new_infections;
    
    for (i in 1:N_strains) {
      int e_idx = 2 + 3*(i-1);
      int i_idx = 3 + 3*(i-1);
      int r_idx = 4 + 3*(i-1);
      E[i] = states_forecast[t, e_idx];
      I[i] = states_forecast[t, i_idx];
      R[i] = states_forecast[t, r_idx];
    }
    
    for (i in 1:N_strains) {
      lambda[i] = transmission_rate_forecast[i, t_rel] * I[i];
      
      real effective_S = S;
      for (j in 1:N_strains) {
        if (j != i) {
          real immunity_effect = R[j] * (1 - cross_immunity[j,i]);
          effective_S += immunity_effect;
        }
      }
      
      effective_S = fmax(1e-6, fmin(1.0, effective_S));
      new_infections[i] = fmin(lambda[i] * effective_S, 0.5 * effective_S);
      
      daily_incidence_forecast[i][t_rel] = sigma[i] * E[i];
    }
    
    for (i in 1:N_strains) {
      real S_eff = S;
      for (j in 1:N_strains) {
        if (j != i) {
          real immunity_contribution = R[j] * (1 - cross_immunity[j,i]);
          S_eff += immunity_contribution;
        }
      }
      S_eff = fmax(0, fmin(1, S_eff));
      R_eff_daily[t, i] = R0_forecast[i, t_rel] * S_eff;
    }
    
    {
      real dt = 1.0;
      real dS = -sum(new_infections) + sum(mu .* R);
      
      vector[19] new_state = to_vector(states_forecast[t,]);
      new_state[1] = S + dt * dS;
      
      for (i in 1:N_strains) {
        int e_idx = 2 + 3*(i-1);
        int i_idx = 3 + 3*(i-1);
        int r_idx = 4 + 3*(i-1);
        
        real dE = new_infections[i] - sigma[i] * E[i];
        real dI = sigma[i] * E[i] - gamma[i] * I[i];
        real dR = gamma[i] * I[i] - mu[i] * R[i];
        
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
        new_state = to_vector(states_forecast[t,]);
      }
      
      for (j in 1:19) {
        states_forecast[t+1, j] = new_state[j];
      }
    }
  }
  
  // 3. Calculate weekly R_eff (baseline)
  for (w in 1:T_weeks_total) {
    for (i in 1:N_strains) {
      real sum_R_eff = 0;
      int count = 0;
      for (d in 1:7) {
        int day_idx = (w-1)*7 + d;
        if (day_idx <= T_days_total) {
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
  
  // 4. Historical predictions and log-likelihood
  for (w in 1:T_weeks) {
    for (i in 1:N_strains) {
      real expected_cases;
      
      if (i <= 5) {
        expected_cases = weekly_incidence[i][w] * population * detection_rate[i];
      } else {
        expected_cases = weekly_incidence[i][w] * population * hospitalization_rate;
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
  
  // 5. Baseline forecast
  for (i in 1:N_strains) {
    for (w in 1:T_weeks_forecast) {
      weekly_incidence_forecast[i][w] = 0;
      for (d in 1:7) {
        int day_idx = (w-1)*7 + d;
        if (day_idx <= T_days_forecast) {
          weekly_incidence_forecast[i][w] += daily_incidence_forecast[i][day_idx];
        }
      }
    }
  }
  
  for (w in 1:T_weeks_forecast) {
    for (i in 1:N_strains) {
      real expected_cases;
      if (i <= 5) {
        expected_cases = weekly_incidence_forecast[i][w] * population * detection_rate[i];
      } else {
        expected_cases = weekly_incidence_forecast[i][w] * population * hospitalization_rate;
      }
      expected_cases = fmax(0.1, expected_cases);
      forecast_cases[w, i] = neg_binomial_2_rng(expected_cases, phi[i]);
    }
  }
  
  // NEW: Calculate baseline totals first (scenario 1 = 0% reduction)
  for (i in 1:N_strains) {
    cases_baseline_typhoon_period[i] = 0;
    cases_baseline_total_period[i] = 0;
  }
  
  // 6. Typhoon scenario simulations (WITH 20 SCENARIOS)
  for (typhoon_idx in 1:20) {
    real typhoon_reduction = typhoon_effects[typhoon_idx];
    int typhoon_days = typhoon_days_arr[typhoon_idx];
    int shift = shift_arr[typhoon_idx];
    matrix[T_days_total + 1, 19] states_typhoon;
    vector[T_days_forecast] daily_incidence_typhoon[N_strains];
    matrix[T_days_forecast, N_strains] R_eff_daily_scenario;
    
    for (i in 1:N_strains) {
      daily_incidence_typhoon[i] = rep_vector(0, T_days_forecast);
    }
    
    if (shift >= 0) {
      // Normal or delayed
      for (t in 1:(T_days + 1)) {
        for (j in 1:19) {
          states_typhoon[t, j] = states[t, j];
        }
      }
      
      for (t_rel in 1:T_days_forecast) {
        int t = T_days + t_rel;
        real S = states_typhoon[t, 1];
        vector[N_strains] E;
        vector[N_strains] I;
        vector[N_strains] R;
        vector[N_strains] lambda;
        vector[N_strains] new_infections;
        
        real reduction_factor = ((t_rel > shift) && (t_rel <= shift + typhoon_days)) ? 1.0 - typhoon_reduction : 1.0;
        
        for (i in 1:N_strains) {
          int e_idx = 2 + 3*(i-1);
          int i_idx = 3 + 3*(i-1);
          int r_idx = 4 + 3*(i-1);
          E[i] = states_typhoon[t, e_idx];
          I[i] = states_typhoon[t, i_idx];
          R[i] = states_typhoon[t, r_idx];
        }
        
        for (i in 1:N_strains) {
          lambda[i] = transmission_rate_forecast[i, t_rel] * reduction_factor * I[i];
          
          real effective_S = S;
          for (j in 1:N_strains) {
            if (j != i) {
              real immunity_effect = R[j] * (1 - cross_immunity[j,i]);
              effective_S += immunity_effect;
            }
          }
          
          effective_S = fmax(1e-6, fmin(1.0, effective_S));
          new_infections[i] = fmin(lambda[i] * effective_S, 0.5 * effective_S);
          
          daily_incidence_typhoon[i][t_rel] = sigma[i] * E[i];
        }
        
        for (i in 1:N_strains) {
          real S_eff = S;
          for (j in 1:N_strains) {
            if (j != i) {
              real immunity_contribution = R[j] * (1 - cross_immunity[j,i]);
              S_eff += immunity_contribution;
            }
          }
          S_eff = fmax(0, fmin(1, S_eff));
          R_eff_daily_scenario[t_rel, i] = R0_forecast[i, t_rel] * reduction_factor * S_eff;
        }
        
        {
          real dt = 1.0;
          real dS = -sum(new_infections) + sum(mu .* R);
          
          vector[19] new_state = to_vector(states_typhoon[t,]);
          new_state[1] = S + dt * dS;
          
          for (i in 1:N_strains) {
            int e_idx = 2 + 3*(i-1);
            int i_idx = 3 + 3*(i-1);
            int r_idx = 4 + 3*(i-1);
            
            real dE = new_infections[i] - sigma[i] * E[i];
            real dI = sigma[i] * E[i] - gamma[i] * I[i];
            real dR = gamma[i] * I[i] - mu[i] * R[i];
            
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
            new_state = to_vector(states_typhoon[t,]);
          }
          
          for (j in 1:19) {
            states_typhoon[t+1, j] = new_state[j];
          }
        }
      }
    } else {
      // Early
      int pre_days = -shift;
      int start_historical = T_days - pre_days + 1;
      vector[19] current_state = to_vector(states[start_historical, ]);
      
      // Simulate pre_days with reduction
      for (pre_t in 1:pre_days) {
        int historical_t = start_historical + pre_t - 1;
        real S = current_state[1];
        vector[N_strains] E;
        vector[N_strains] I;
        vector[N_strains] R;
        vector[N_strains] lambda;
        vector[N_strains] new_infections;
        
        real reduction_factor = 1.0 - typhoon_reduction;
        
        for (i in 1:N_strains) {
          int e_idx = 2 + 3*(i-1);
          int i_idx = 3 + 3*(i-1);
          int r_idx = 4 + 3*(i-1);
          E[i] = current_state[e_idx];
          I[i] = current_state[i_idx];
          R[i] = current_state[r_idx];
        }
        
        for (i in 1:N_strains) {
          lambda[i] = transmission_rate[i, historical_t] * reduction_factor * I[i];
          
          real effective_S = S;
          for (j in 1:N_strains) {
            if (j != i) {
              real immunity_effect = R[j] * (1 - cross_immunity[j,i]);
              effective_S += immunity_effect;
            }
          }
          
          effective_S = fmax(1e-6, fmin(1.0, effective_S));
          new_infections[i] = fmin(lambda[i] * effective_S, 0.5 * effective_S);
        }
        
        {
          real dt = 1.0;
          real dS = -sum(new_infections) + sum(mu .* R);
          
          vector[19] new_state = current_state;
          new_state[1] = S + dt * dS;
          
          for (i in 1:N_strains) {
            int e_idx = 2 + 3*(i-1);
            int i_idx = 3 + 3*(i-1);
            int r_idx = 4 + 3*(i-1);
            
            real dE = new_infections[i] - sigma[i] * E[i];
            real dI = sigma[i] * E[i] - gamma[i] * I[i];
            real dR = gamma[i] * I[i] - mu[i] * R[i];
            
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
            new_state = current_state;
          }
          
          current_state = new_state;
        }
      }
      
      // Set initial forecast state
      for (j in 1:19) {
        states_typhoon[T_days + 1, j] = current_state[j];
      }
      
      int remaining_days = typhoon_days - pre_days;
      
      for (t_rel in 1:T_days_forecast) {
        int t = T_days + t_rel;
        real S = states_typhoon[t, 1];
        vector[N_strains] E;
        vector[N_strains] I;
        vector[N_strains] R;
        vector[N_strains] lambda;
        vector[N_strains] new_infections;
        
        real reduction_factor = (t_rel <= remaining_days) ? 1.0 - typhoon_reduction : 1.0;
        
        for (i in 1:N_strains) {
          int e_idx = 2 + 3*(i-1);
          int i_idx = 3 + 3*(i-1);
          int r_idx = 4 + 3*(i-1);
          E[i] = states_typhoon[t, e_idx];
          I[i] = states_typhoon[t, i_idx];
          R[i] = states_typhoon[t, r_idx];
        }
        
        for (i in 1:N_strains) {
          lambda[i] = transmission_rate_forecast[i, t_rel] * reduction_factor * I[i];
          
          real effective_S = S;
          for (j in 1:N_strains) {
            if (j != i) {
              real immunity_effect = R[j] * (1 - cross_immunity[j,i]);
              effective_S += immunity_effect;
            }
          }
          
          effective_S = fmax(1e-6, fmin(1.0, effective_S));
          new_infections[i] = fmin(lambda[i] * effective_S, 0.5 * effective_S);
          
          daily_incidence_typhoon[i][t_rel] = sigma[i] * E[i];
        }
        
        for (i in 1:N_strains) {
          real S_eff = S;
          for (j in 1:N_strains) {
            if (j != i) {
              real immunity_contribution = R[j] * (1 - cross_immunity[j,i]);
              S_eff += immunity_contribution;
            }
          }
          S_eff = fmax(0, fmin(1, S_eff));
          R_eff_daily_scenario[t_rel, i] = R0_forecast[i, t_rel] * reduction_factor * S_eff;
        }
        
        {
          real dt = 1.0;
          real dS = -sum(new_infections) + sum(mu .* R);
          
          vector[19] new_state = to_vector(states_typhoon[t,]);
          new_state[1] = S + dt * dS;
          
          for (i in 1:N_strains) {
            int e_idx = 2 + 3*(i-1);
            int i_idx = 3 + 3*(i-1);
            int r_idx = 4 + 3*(i-1);
            
            real dE = new_infections[i] - sigma[i] * E[i];
            real dI = sigma[i] * E[i] - gamma[i] * I[i];
            real dR = gamma[i] * I[i] - mu[i] * R[i];
            
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
            new_state = to_vector(states_typhoon[t,]);
          }
          
          for (j in 1:19) {
            states_typhoon[t+1, j] = new_state[j];
          }
        }
      }
    }
    
    for (w in 1:T_weeks_total) {
      for (i in 1:N_strains) {
        if (w <= T_weeks) {
          R_eff_scenarios[typhoon_idx][w,i] = R_eff[w,i];
        } else {
          int w_forecast = w - T_weeks;
          real sum_R_eff = 0;
          int count = 0;
          for (d in 1:7) {
            int day_idx = (w_forecast-1)*7 + d;
            if (day_idx <= T_days_forecast) {
              sum_R_eff += R_eff_daily_scenario[day_idx, i];
              count += 1;
            }
          }
          if (count > 0) {
            R_eff_scenarios[typhoon_idx][w,i] = sum_R_eff / count;
          } else {
            R_eff_scenarios[typhoon_idx][w,i] = 0;
          }
        }
      }
    }
    
    // NEW: Calculate metrics for this scenario for each strain
    for (i in 1:N_strains) {
      real typhoon_total = 0;
      real total_period = 0;
      real recovery_total = 0;
      
      for (w in 1:T_weeks_forecast) {
        real weekly_inc = 0;
        for (d in 1:7) {
          int day_idx = (w-1)*7 + d;
          if (day_idx <= T_days_forecast) {
            weekly_inc += daily_incidence_typhoon[i][day_idx];
          }
        }
        
        real expected_cases;
        if (i <= 5) {
          expected_cases = weekly_inc * population * detection_rate[i];
        } else {
          expected_cases = weekly_inc * population * hospitalization_rate;
        }
        
        expected_cases = fmax(0.1, expected_cases);
        forecast_cases_typhoon[typhoon_idx, w, i] = neg_binomial_2_rng(expected_cases, phi[i]);
        
        // Accumulate for metrics
        total_period += expected_cases;
        if (w <= typhoon_weeks) {
          typhoon_total += expected_cases;
        } else {
          recovery_total += expected_cases;
        }
      }
      
      // Store raw case counts
      cases_typhoon_period[typhoon_idx, i] = typhoon_total;
      cases_total_period[typhoon_idx, i] = total_period;
      
      // Store baseline (scenario 1)
      if (typhoon_idx == 1) {
        cases_baseline_typhoon_period[i] = typhoon_total;
        cases_baseline_total_period[i] = total_period;
      }
      
      // Calculate percentage reductions for this draw
      if (typhoon_idx > 1) {
        // Avoid division by zero
        real base_typhoon_safe = fmax(cases_baseline_typhoon_period[i], 0.01);
        real base_total_safe = fmax(cases_baseline_total_period[i], 0.01);
        
        reduction_typhoon_period[typhoon_idx, i] = 
          fmax(0, fmin(100, (1 - typhoon_total / base_typhoon_safe) * 100));
        reduction_total_period[typhoon_idx, i] = 
          fmax(0, fmin(100, (1 - total_period / base_total_safe) * 100));
      } else {
        // Baseline scenario has 0% reduction
        reduction_typhoon_period[typhoon_idx, i] = 0;
        reduction_total_period[typhoon_idx, i] = 0;
      }
      
      // Calculate average weekly cases
      int n_recovery = T_weeks_forecast - typhoon_weeks;
      avg_weekly_typhoon[typhoon_idx, i] = typhoon_total / typhoon_weeks;
      avg_weekly_recovery[typhoon_idx, i] = recovery_total / n_recovery;
    }
  }
}
