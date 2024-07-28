library(haven)
library(tidyverse)
library(Synth)
library(devtools)
library(fect)

canadian_provinces = c('alberta', 'british columbia', 'manitoba',
                       'New brunswick', 'Newfoundland and Labrador', 'nova scotia',
                       'ontario', 'Prince Edward Island', 'quebec', 'saskatchewan')

# Read data
# energy_panel = read_dta("energy_panel_stata_order.dta", col_select=-c(16:19))
energy_panel = read_csv("energy_panel_new.csv")

# Encode region as factor
region_levels = unique(energy_panel$region)
region_factors = factor(energy_panel$region, levels=region_levels)
energy_panel$region_id = as.numeric(region_factors)
treatment_region_id = which(region_levels == 'british columbia')

# Find which index of region_levels is in canadian_provinces
canadian_provinces_ids = which(region_levels %in% canadian_provinces)

# Derive secondary metrics
energy_panel = energy_panel %>%
  mutate(
    lpop = log(pop),
    empinenergy = 100 * employment_mogas / emp_tot
  )

energy_panel = energy_panel %>%
  group_by(region) %>%
  mutate(emindex = 100 * emissions / emissions[year == 1998])

# Prep data for synth
energy_panel = as.data.frame(energy_panel)
dataprep_out = dataprep(
  foo = energy_panel,
  special.predictors = list(
    list("emindex", c(2000), "mean"),
    list("emindex", c(2004), "mean"),
    list("emindex", c(2007), "mean"),
    list("rgdp", 1998:2007, "mean"),
    list("unemprate", 1998:2007, "mean"),
    list("lpop", 2007, "mean"),
    list("empinenergy", 2007, "mean")
    ),
  dependent = "emindex",
  unit.variable = "region_id",
  unit.names.variable = "region",
  treatment.identifier = treatment_region_id,
  controls.identifier = c(1:(treatment_region_id-1), (treatment_region_id+1):61),
  time.variable = "year",
  time.predictors.prior = 1998:2007,
  time.optimize.ssr = 1998:2007,
  time.plot = 1998:2017
)

# Run synthetic control generation
synth_out = synth(data.prep.obj = dataprep_out)
sqrt(synth_out$loss.v)

# Get predictor & region weights
predictor_weights = synth_out$solution.v
predictor_weights
region_weights = data.frame(region=region_levels[-treatment_region_id], weight=synth_out$solution.w)
region_weights

# Confirm pre-treatment RMSPE from weights
ycontrol_pre_r = dataprep_out$Z0
ytreat_pre_r = dataprep_out$Z1
ypred_pre_r = ycontrol_pre_r %*% synth_out$solution.w
sqrt(mean((ypred_pre_r - ytreat_pre_r)^2))
cat('ytreat_pre_r:', mean(ytreat_pre_r), '\nypred_pre_r:', mean(ypred_pre_r))

# Double-check pre-treatment RMSPE from Stata
region_weights_stata = as.data.frame(read_csv("stata_region_weights.csv"))
ytreat_pre_stata = dataprep_out$Z1
ypred_pre_stata = dataprep_out$Z0 %*% region_weights_stata$weights
sqrt(mean((ypred_pre_stata - ytreat_pre_stata)^2))
cat('ytreat_pre_stata:', mean(ytreat_pre_stata), '\nypred_pre_stata:', mean(ypred_pre_stata))

# Calculate post-treatment RMSPE
ycontrol_post_r = dataprep_out$Y0[as.character(2008:2017), ]
ytreat_post_r = dataprep_out$Y1[as.character(2008:2017), ]
ypred_post_r = ycontrol_post_r %*% synth_out$solution.w
sqrt(mean((ypred_post_r - ytreat_post_r)^2))
cat('ytreat_post_r:', mean(ytreat_post_r), '\nypred_post_r:', mean(ypred_post_r))

# Double-check post-treatment RMSPE from Stata
ys_stata = as.data.frame(read_csv("ys_stata.csv")) %>% filter(Period == 1)
# Average Y_treated & Y_synthetic from ytreat_post_stata
ytreat_post_stata = as.numeric(ys_stata$Y_treated)
ypred_post_stata = as.numeric(ys_stata$Y_synthetic)
cat('ytreat_post_stata:', mean(ytreat_post_stata), '\nypred_post_stata:', mean(ypred_post_stata))

# Calculate ATT
pre_diff = mean(ytreat_pre_r) - mean(ypred_pre_r)
post_diff = mean(ytreat_post_r) - mean(ypred_post_r)
att = post_diff - pre_diff
cat('Pre-treatment difference:', pre_diff, '\nPost-treatment difference:', post_diff, '\nATT:', att)

# Plot treated vs synthetic CO2 emissions index
png('output/synth_path_plot.png')
path.plot(synth_out, dataprep_out, Ylim=c(80,120), Ylab='Emission index')
abline(v=2008, lty=1)
abline(v=2001, lty=2)
abline(v=2004, lty=2)
dev.off()
gaps.plot(synth_out, dataprep_out)

### Cross validation
# Define function that takes into the list of years of emindex to include in the synthetic control, 
# then prep the data, train the model, and evalute its RMSPE on training & validation set
synth_cross_val = function(outcome_type, covariate_type, region_type, verbose=FALSE) {
  # Create list of outcome covariates to include in the synthetic control
  
  if (outcome_type == 'all') {
  outcome_covariates = lapply(1999:2003, function(year) {list("emindex", year, "mean")})
  }
  else if (outcome_type == 'last') {
    outcome_covariates = list(list("emindex", 2003, "mean"))
  }
  else if (outcome_type == 'mean') {
    outcome_covariates = list(list("emindex", 1998:2003, "mean"))
  }
  else if (outcome_type == 'first_mid_last') {
    outcome_covariates = lapply(c(1999, 2001, 2003), function(year) {list("emindex", year, "mean")})
  }
  
  # Create list of other covariates to include in the synthetic control
  if (covariate_type == 'none') {
    other_covariates = list()
  }
  else if (covariate_type == 'econ') {
    other_covariates = list(
      list("rgdp", 1998:2003, "mean"),
      list("unemprate", 1998:2003, "mean"),
      list("lpop", 2003, "mean")
    )
  }
  else if (covariate_type == 'all') {
    other_covariates = list(
      list("rgdp", 1998:2003, "mean"),
      list("unemprate", 1998:2003, "mean"),
      list("lpop", 2003, "mean"),
      list("empinenergy", 2003, "mean")
    )
  }
  
  special_predictors = c(outcome_covariates, other_covariates)
  if (verbose) {
    print('Special predictors:')
    print(special_predictors)}
  
  if (region_type == 'all') {
    input_data = energy_panel
    control_region_ids = c(1:(treatment_region_id-1), (treatment_region_id+1):61)
  }
  else if (region_type == 'canadian') {
    input_data = energy_panel[energy_panel$region %in% canadian_provinces, ]
    control_region_ids = canadian_provinces_ids[-which(canadian_provinces_ids == treatment_region_id)]
  }
  
  # Prep data
  dataprep_out = dataprep(
    foo = input_data,
    special.predictors = special_predictors,
    dependent = "emindex",
    unit.variable = "region_id",
    unit.names.variable = "region",
    treatment.identifier = treatment_region_id,
    controls.identifier = control_region_ids,
    time.variable = "year",
    time.predictors.prior = 1998:2003,
    time.optimize.ssr = 1998:2003,
    time.plot = 1998:2007
  )
  if (verbose) {
    cat('input_data dim', dim(input_data))
    cat('X0 dim:', dim(dataprep_out$X0))
    cat('Z0 dim:', dim(dataprep_out$Z0))
    }
  
  # Run synthetic control generation
  synth_out = synth(data.prep.obj=dataprep_out, verbose=verbose)
  
  # Calculate pre-treatment RMSPE
  ycontrol_pre = dataprep_out$Z0
  ytreat_pre = dataprep_out$Z1
  ypred_pre = ycontrol_pre %*% synth_out$solution.w
  pre_rmspe = sqrt(mean((ypred_pre - ytreat_pre)^2))
  
  # Calculate post-treatment RMSPE
  ycontrol_post = dataprep_out$Y0[as.character(2004:2007), ]
  ytreat_post = dataprep_out$Y1[as.character(2004:2007), ]
  ypred_post = ycontrol_post %*% synth_out$solution.w
  post_rmspe = sqrt(mean((ypred_post - ytreat_post)^2))
  
  if (verbose) {
    cat('Pre RMSPE:', pre_rmspe, '\nPost RMSPE:', post_rmspe)
    }
  return(c(pre_rmspe, post_rmspe))
}

# Try different cross-validation scenarios with different outcome_type, covariate_type, region_type, and store the pre & post RMSPE for each configuration alongside the configuration values
cross_val_results = data.frame()
for (outcome_type in c('all', 'last', 'mean', 'first_mid_last')) {
  for (covariate_type in c('none', 'econ', 'all')) {
    for (region_type in c('all', 'canadian')) {
      cat('outcome_type:', outcome_type, ', covariate_type:', covariate_type, ', region_type:', region_type)
      rmspes = synth_cross_val(outcome_type, covariate_type, region_type, verbose=FALSE)
      cross_val_results = rbind(cross_val_results, c(outcome_type, covariate_type, region_type, rmspes))
    }
  }
}
colnames(cross_val_results) = c('outcome_type', 'covariate_type', 'region_type', 'pre_rmspe', 'post_rmspe')
write.csv(cross_val_results, 'output/cross_val_results.csv', row.names=FALSE)

### Inference
# Run synthetic control using last year of pre-treatment data as only predictor
synth_last = function(pre_years, post_years, dependent="emindex", verbose=FALSE) {
  dataprep_out = dataprep(
    foo = as.data.frame(energy_panel),
    special.predictors = list(
      list(dependent, c(max(pre_years)), "mean")
    ),
    dependent = dependent,
    unit.variable = "region_id",
    unit.names.variable = "region",
    treatment.identifier = treatment_region_id,
    controls.identifier = c(1:(treatment_region_id-1), (treatment_region_id+1):61),
    time.variable = "year",
    time.predictors.prior = pre_years,
    time.optimize.ssr = pre_years,
    time.plot = c(pre_years, post_years)
  )
  
  # Run synthetic control generation
  synth_out = synth(data.prep.obj=dataprep_out, verbose=verbose)
  
  # Calculate pre-treatment RMSPE
  ycontrol_pre = dataprep_out$Z0
  ytreat_pre = dataprep_out$Z1
  ypred_pre = ycontrol_pre %*% synth_out$solution.w
  pre_rmspe = sqrt(mean((ypred_pre - ytreat_pre)^2))
  
  # Calculate post-treatment RMSPE
  ycontrol_post = dataprep_out$Y0[as.character(post_years), ]
  ytreat_post = dataprep_out$Y1[as.character(post_years), ]
  ypred_post = ycontrol_post %*% synth_out$solution.w
  post_rmspe = sqrt(mean((ypred_post - ytreat_post)^2))
  
  if (verbose) {
  cat('Pre RMSPE:', pre_rmspe, '\nPost RMSPE:', post_rmspe)
  }
  
  return(list(
    synth_out=synth_out,
    dataprep_out=dataprep_out,
    ys=data.frame(
      year=c(pre_years, post_years),
      treated=rep(c(0, 1), c(length(pre_years), length(post_years))),
      ytreat=c(dataprep_out$Y1),
      ypred=c(ypred_pre, ypred_post))))
}

calculate_att = function(ys) {
  ys_avg = ys %>% group_by(treated) %>% summarise(across(c(ytreat, ypred), mean))
  ytreat_pre_avg = ys_avg %>% filter(treated==0) %>% select(ytreat) %>% as.numeric()
  ypred_pre_avg = ys_avg %>% filter(treated==0) %>% select(ypred) %>% as.numeric()
  ytreat_post_avg = ys_avg %>% filter(treated==1) %>% select(ytreat) %>% as.numeric()
  ypred_post_avg = ys_avg %>% filter(treated==1) %>% select(ypred) %>% as.numeric()
  
  pre_diff = ytreat_pre_avg - ypred_pre_avg
  post_diff = ytreat_post_avg - ypred_post_avg
  att = post_diff - pre_diff
  return(list(
    ytreat_pre_avg=ytreat_pre_avg,
    ypred_pre_avg=ypred_pre_avg,
    ytreat_post_avg=ytreat_post_avg,
    ypred_post_avg=ypred_post_avg,
    pre_diff=pre_diff,
    post_diff=post_diff,
    att=att))
}

synth_last_result = synth_last(1998:2007, 2008:2017)
synth_last_att = calculate_att(synth_last_result$ys)
synth_last_att

# Save path.plot as png and set y limit from 80 to 120
png('output/synth_last_path_plot.png')
path.plot(synth_last_result$synth_out, synth_last_result$dataprep_out, 
          Ylim=c(80, 120), Ylab='emindex')
abline(v=2008, lty=2)
dev.off()

# Save gaps plot as png
png('output/synth_last_gaps_plot.png')
gaps.plot(synth_last_result$synth_out, synth_last_result$dataprep_out,
          Ylab='Treatment Effect')
abline(v=2008, lty=2)
dev.off()

synth_out_last = synth_last_result$synth_out
dataprep_out_last = synth_last_result$dataprep_out
predictor_weights_last = synth_out_last$solution.v
region_weights_last = data.frame(region=region_levels[-treatment_region_id], weight=synth_out_last$solution.w)
region_weights_last

# 3-fold t-test
fold_years = list(1998:2000, 2001:2003, 2004:2007)

scm_t_test = function(fold_years, dependent="emindex") {
  fold_results = data.frame()
  
  for (pre_years in fold_years) {
    synth_last_result = synth_last(pre_years, 2008:2017, dependent=dependent, verbose=FALSE)
    synth_last_att = calculate_att(synth_last_result$ys)
    fold_results = rbind(fold_results, c(
      min(pre_years),
      max(pre_years),
      synth_last_att$ytreat_pre_avg,
      synth_last_att$ypred_pre_avg,
      synth_last_att$ytreat_post_avg,
      synth_last_att$ypred_post_avg,
      synth_last_att$pre_diff,
      synth_last_att$post_diff,
      synth_last_att$att))
  }
  
  colnames(fold_results) = c('start_year', 'end_year', 'ytreat_pre_avg', 'ypred_pre_avg', 'ytreat_post_avg', 'ypred_post_avg', 'pre_diff', 'post_diff', 'att')
  
  K = length(fold_years)
  T0 = 10
  T = 20
  att_t = mean(fold_results$att)
  sigma = sqrt((1 + (K*min(floor(T0/K), T-T0-1)/(T-T0-1))) * 1/(K-1) * sum((fold_results$att - mean(fold_results$att))^2))
  t = sqrt(K) * att_t / sigma
  
  # Find p-value from t-value with 3 degrees of freedom
  p_value = (1 - pt(abs(t), df=K)) * 2
  
  # Find lower & uper CI of t-test
  se = sigma / sqrt(K)
  lower_ci_t = att_t - qt(0.975, df=K) * se
  upper_ci_t = att_t + qt(0.975, df=K) * se
  cat('att_t', att_t, 'p_value', p_value, 'lower_ci_t', lower_ci_t, 'upper_ci_t', upper_ci_t)
  return (list(fold_results=fold_results, att_t=att_t, p_value=p_value, lower_ci_t=lower_ci_t, upper_ci_t=upper_ci_t))
}

scm_t_test_result = scm_t_test(fold_years, dependent="emindex")
scm_t_test_result

### Generalized SCM
# Prep data
dataprep_genscm = energy_panel[c('region', 'year', 'emindex')]
dataprep_genscm$treated = ifelse(dataprep_genscm$year %in% 2008:2017 & dataprep_genscm$region == 'british columbia', 1, 0)

train_gen_scm = function(dataprep_genscm, method, dependent = 'emindex', threshold='upper') {
  gen_scm = fect(data = dataprep_genscm,
       Y = dependent,
       D = 'treated',
       index = c('region', 'year'),
       se = TRUE,
       nboots = 100,
       method = method)
  #Obtain 95% CI bounds
  lower_ci_genscm = gen_scm$att.avg - qnorm(0.975)*sd(gen_scm$att.avg.boot)
  upper_ci_genscm = gen_scm$att.avg + qnorm(0.975)*sd(gen_scm$att.avg.boot)
  
  threshold = ifelse(threshold == 'upper', upper_ci_genscm, lower_ci_genscm)

  # Save plot
  png(paste("output/gen_scm_", method, "_", dependent, ".png", sep=""))
  print(plot(gen_scm, main = paste("Estimated ATT:", method), ylab = "Effect of D on Y",
             cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8,
             bound = "both", tost.threshold = threshold))
  dev.off()
  
  return(list(
    gen_scm=gen_scm,
    att=gen_scm$att.avg,
    lower_ci=lower_ci_genscm,
    upper_ci=upper_ci_genscm)
  )
}

set.seed(12)
fe = train_gen_scm(dataprep_genscm, "fe", threshold='lower')
fe
ife = train_gen_scm(dataprep_genscm, "ife", threshold='lower')
ife



# Compile ATT, lower CI & upper CI of fe, ife and t-test together in one table
att_table = data.frame(
  method=c('FE', 'IFE', 'T-test'),
  att=c(fe$att, ife$att, scm_t_test_result$att_t),
  lower_ci=c(fe$lower_ci, ife$lower_ci, scm_t_test_result$lower_ci_t),
  upper_ci=c(fe$upper_ci, ife$upper_ci, scm_t_test_result$upper_ci_t),
  p_value=c(fe$gen_scm$est.avg[, 'p.value'], ife$gen_scm$est.avg[, 'p.value'], scm_t_test_result$p_value)
)
write.csv(att_table, "output/att_comparison.csv", row.names=FALSE)

# Plot att and associated confidence interval lower & upper bounds for different methods on the same graph\
ggplot(att_table, aes(x = method, y = att)) +
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
  labs(title = "Point Estimates with Confidence Intervals",
       x = "Method",
       y = "ATT")
ggsave("output/att_comparison.png")


### Diagnose parallel trend violation
# Plot realgdp and emissions of british columbia over years
energy_panel = energy_panel %>%
  group_by(region) %>%
  mutate(gdpindex = 100 * gdp_new / gdp_new[year == 1998])

dataprep_gdpindex = dataprep(
  foo = as.data.frame(energy_panel),
  special.predictors = list(
    list("gdpindex", c(2007), "mean")
  ),
  dependent = "gdpindex",
  unit.variable = "region_id",
  unit.names.variable = "region",
  treatment.identifier = treatment_region_id,
  controls.identifier = c(1:(treatment_region_id-1), (treatment_region_id+1):61),
  time.variable = "year",
  time.predictors.prior = 1998:2007,
  time.optimize.ssr = 1998:2007,
  time.plot = c(1998:2017)
)


gdp_control = dataprep_gdpindex$Y0
gdp_treat = dataprep_gdpindex$Y1
gdp_pred = gdp_control %*% synth_out_last$solution.w


# Plot gdp_treat and gdp_pred over year
gdp_data = data.frame(rownames(gdp_treat), gdp_treat, gdp_pred)
colnames(gdp_data) = c('year', 'gdp_treat', 'gdp_pred')
gdp_data$year = as.numeric(gdp_data$year)
ggplot(gdp_data, aes(x=year)) +
  geom_line(aes(y=gdp_treat, color='British Columbia')) +
  geom_line(aes(y=gdp_pred, color='Synthetic Control')) +
  labs(title='Actual vs Predicted GDP of British Columbia',
       x='Year',
       y='Value') +
  scale_color_manual(name='Region', values=c('British Columbia'='blue', 'Synthetic Control'='red')) +
  geom_vline(xintercept = 2008, linetype = "dotted", color = "black")

# Plot gdp_diff against year in gdp_data, with gdp_diff as GDP Difference in chart legends
gdp_data$gdp_diff = gdp_data$gdp_treat - gdp_data$gdp_pred
ggplot(gdp_data, aes(x=year, y=gdp_diff)) +
  geom_line(color='black', size=1) +
  labs(title='GDP Difference between Actual and Predicted',
       x='Year',
       y='Value') +
  geom_vline(xintercept = 2008, linetype = "dotted", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")


### Rerun the synthetic control method with emission per GDP as dependent variable
energy_panel$em_per_gdp = energy_panel$emissions / energy_panel$gdp_new
em_per_gdp_last_result = synth_last(1998:2007, 2008:2017, dependent = 'em_per_gdp')
em_per_gdp_last_att = calculate_att(em_per_gdp_last_result$ys)
em_per_gdp_last_att

# Save path.plot as png and set y limit from 80 to 120
png('output/em_per_gdp_last_path_plot.png')
path.plot(em_per_gdp_last_result$synth_out, em_per_gdp_last_result$dataprep_out, 
          Ylab='Emission per GDP')
abline(v=2008, lty=2)
dev.off()

# Save gaps plot as png
png('output/em_per_gdp_last_gaps_plot.png')
gaps.plot(em_per_gdp_last_result$synth_out, em_per_gdp_last_result$dataprep_out,
          Ylab='Treatment Effect')
abline(v=2008, lty=2)
dev.off()

# Generalized SCM for emission per GDP
dataprep_genscm_em_per_gdp = energy_panel[c('region', 'year', 'em_per_gdp')]
dataprep_genscm_em_per_gdp$treated = ifelse(dataprep_genscm_em_per_gdp$year %in% 2008:2017 & dataprep_genscm_em_per_gdp$region == 'british columbia', 1, 0)

set.seed(12)
fe_em_per_gdp = train_gen_scm(dataprep_genscm_em_per_gdp, "fe", dependent="em_per_gdp", threshold='upper')
fe_em_per_gdp
ife_em_per_gdp = train_gen_scm(dataprep_genscm_em_per_gdp, "ife", dependent="em_per_gdp", threshold='upper')
ife_em_per_gdp

scm_t_test_result_em_per_gdp = scm_t_test(fold_years, dependent="em_per_gdp")
scm_t_test_result_em_per_gdp
att_table_em_per_gdp = data.frame(
  method=c('FE', 'IFE', 'T-test'),
  att=c(fe_em_per_gdp$att, ife_em_per_gdp$att, scm_t_test_result_em_per_gdp$att_t),
  lower_ci=c(fe_em_per_gdp$lower_ci, ife_em_per_gdp$lower_ci, scm_t_test_result_em_per_gdp$lower_ci_t),
  upper_ci=c(fe_em_per_gdp$upper_ci, ife_em_per_gdp$upper_ci, scm_t_test_result_em_per_gdp$upper_ci_t),
  p_value=c(fe_em_per_gdp$gen_scm$est.avg[, 'p.value'], ife_em_per_gdp$gen_scm$est.avg[, 'p.value'], scm_t_test_result_em_per_gdp$p_value)
)
write.csv(att_table_em_per_gdp, "output/att_table_em_per_gdp.csv", row.names=FALSE)
ggplot(att_table_em_per_gdp, aes(x = method, y = att)) +
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
  labs(title = "Point Estimates with Confidence Intervals",
       x = "Method",
       y = "ATT")
  