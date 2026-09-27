##############################################################################
# ANCHORED INDIRECT TREATMENT COMPARISON ON A BODY-WEIGHT ENDPOINT
# MAIC (two weighting methods) vs STC vs naive comparison
#
# Scenario: Trial 1 (Drug A vs Placebo, IPD available) and 
#          Trial 2 (Drug B vs Placebo, AgD only).
# Placebo anchors both. Estimand: d_AB = d_AC - d_BC (Bucher method).
# 
# Key feature: Analysis uses ONLY what a real analyst would receive from
# Trial 2's publication (arm sizes, outcome means/SEs, covariate means/SDs).
# Trial 2's patient-level data is used only to generate a realistic
# aggregate summary and to validate the final result at the end.
##############################################################################

set.seed(123)

## ---------------------------------------------------------------------------
## STEP 1: Simulate Trial 1 IPD (Drug A vs Placebo)
## ---------------------------------------------------------------------------
# 300 patients, balanced allocation
n1 <- 300
trial1_ipd <- data.frame(
  USUBJID = 1:n1,
  ARM = rep(c("Drug A", "Placebo"), c(150, 150)),
  age = rnorm(n1, mean = 45, sd = 10),
  sex = rbinom(n1, 1, 0.50),  # 1 = Female
  bmi = rnorm(n1, mean = 32, sd = 4)
)

# Generate outcomes
# Placebo response: baseline + age effect + BMI effect + sex effect
placebo_response <- 5 - 0.05 * trial1_ipd$age - 0.30 * (trial1_ipd$bmi - 32) - 2 * trial1_ipd$sex

# Drug A adds an extra benefit that DEPENDS on BMI (effect modifier)
drug_a_benefit <- ifelse(trial1_ipd$ARM == "Drug A", 
                         -5 - 0.4 * (trial1_ipd$bmi - 32), 
                         0)

trial1_ipd$chg_bwt <- placebo_response + drug_a_benefit + rnorm(n1, 0, 3)

head(trial1_ipd)


## ---------------------------------------------------------------------------
## STEP 2: Simulate Trial 2, then extract ONLY aggregate summary
## ---------------------------------------------------------------------------
# Trial 2 has a different population (older, more female, heavier)
n2 <- 280  # Different sample size
trial2_ipd_hidden <- data.frame(
  USUBJID = (n1+1):(n1+n2),
  ARM = rep(c("Drug B", "Placebo"), c(140, 140)),
  age = rnorm(n2, mean = 55, sd = 10),
  sex = rbinom(n2, 1, 0.65),  # higher proportion female
  bmi = rnorm(n2, mean = 36, sd = 4)
)

# Generate outcomes
placebo_response2 <- 5 - 0.05 * trial2_ipd_hidden$age - 0.30 * (trial2_ipd_hidden$bmi - 32) - 2 * trial2_ipd_hidden$sex
drug_b_benefit <- ifelse(trial2_ipd_hidden$ARM == "Drug B", 
                         -7 - 0.4 * (trial2_ipd_hidden$bmi - 32), 
                         0)
trial2_ipd_hidden$chg_bwt <- placebo_response2 + drug_b_benefit + rnorm(n2, 0, 3)

# ===== EXTRACT AGGREGATE DATA (what a real analyst receives) =====
# This is the ONLY Trial 2 information used in the analysis below
trial2_agd <- list(
  n_B = sum(trial2_ipd_hidden$ARM == "Drug B"),
  n_C = sum(trial2_ipd_hidden$ARM == "Placebo"),
  age_mean = mean(trial2_ipd_hidden$age),
  age_sd = sd(trial2_ipd_hidden$age),
  sex_prop = mean(trial2_ipd_hidden$sex),  # proportion female
  bmi_mean = mean(trial2_ipd_hidden$bmi),
  bmi_sd = sd(trial2_ipd_hidden$bmi),
  mean_y_B = mean(trial2_ipd_hidden$chg_bwt[trial2_ipd_hidden$ARM == "Drug B"]),
  se_y_B = sd(trial2_ipd_hidden$chg_bwt[trial2_ipd_hidden$ARM == "Drug B"]) / sqrt(sum(trial2_ipd_hidden$ARM == "Drug B")),
  mean_y_C = mean(trial2_ipd_hidden$chg_bwt[trial2_ipd_hidden$ARM == "Placebo"]),
  se_y_C = sd(trial2_ipd_hidden$chg_bwt[trial2_ipd_hidden$ARM == "Placebo"]) / sqrt(sum(trial2_ipd_hidden$ARM == "Placebo"))
)

# Trial 2 treatment effect (observed difference)
d_BC <- trial2_agd$mean_y_B - trial2_agd$mean_y_C
se_d_BC <- sqrt(trial2_agd$se_y_B^2 + trial2_agd$se_y_C^2)

# Target population characteristics (Trial 2's published means)
target_means <- c(age = trial2_agd$age_mean, 
                  sex = trial2_agd$sex_prop, 
                  bmi = trial2_agd$bmi_mean)

# True effect (for final validation only - never used in analysis)
true_dAB <- -5 - (-7)  


## ---------------------------------------------------------------------------
## STEP 3: Naive (unadjusted) comparison - Trial 1 only
## ---------------------------------------------------------------------------
# Fit an unadjusted model to Trial 1
naive_mod <- lm(chg_bwt ~ ARM, data = trial1_ipd)
d_AC_naive <- coef(naive_mod)["ARMPlacebo"] * -1  # flip sign for Drug A effect


## ---------------------------------------------------------------------------
## STEP 4: MAIC Method 1 - Method of Moments
## ---------------------------------------------------------------------------
# Goal: Weight Trial 1 patients so their WEIGHTED covariate distribution
# matches Trial 2's reported covariate means exactly.
# 
# Method: Find weights w_i = exp(x_i' * beta) that minimize the entropy
# of the weight distribution subject to the moment constraint.
# This is solved via convex optimization (minimizing sum of exponential weights).

get_maic_mom_weights <- function(ipd, target) {
  # Center covariates on target means
  X <- cbind(
    age = ipd$age - target["age"],
    sex = ipd$sex - target["sex"],
    bmi = ipd$bmi - target["bmi"]
  )
  
  # Objective: minimize sum(exp(X %*% beta))
  # (convex, solved via BFGS)
  obj <- function(beta) sum(exp(X %*% beta))
  
  opt <- optim(par = c(0, 0, 0), fn = obj, method = "BFGS")
  
  # Return exponential weights
  as.numeric(exp(X %*% opt$par))
}

w_mom <- get_maic_mom_weights(trial1_ipd, target_means)
ess_mom <- sum(w_mom)^2 / sum(w_mom^2)  # Effective Sample Size

# Fit weighted regression
mom_mod <- lm(chg_bwt ~ ARM, data = trial1_ipd, weights = w_mom)
d_AC_mom <- coef(mom_mod)["ARMPlacebo"] * -1


## ---------------------------------------------------------------------------
## STEP 5: Build synthetic Trial 2 population from aggregate data
## ---------------------------------------------------------------------------
# Since we only have Trial 2's reported means/SDs/proportions (not patient-level
# data), simulate a large pseudo-population to stand in for Trial 2's 
# covariate distribution. Use this in both IPW weighting and STC below.
# (Making it large avoids Monte Carlo noise in estimation.)

simulate_synthetic_population <- function(agd, n_synthetic = 5000) {
  data.frame(
    age = rnorm(n_synthetic, agd$age_mean, agd$age_sd),
    sex = rbinom(n_synthetic, 1, agd$sex_prop),
    bmi = rnorm(n_synthetic, agd$bmi_mean, agd$bmi_sd)
  )
}

set.seed(999)
synth_trial2 <- simulate_synthetic_population(trial2_agd, n_synthetic = 5000)

# Quick check: synthetic sample means should be close to reported Trial 2 means
cat("\nSynthetic population baseline balance check:\n")
print(rbind(synthetic = colMeans(synth_trial2), published = target_means))


## ---------------------------------------------------------------------------
## STEP 6: MAIC Method 2 - Inverse Probability Weighting (IPW)
## ---------------------------------------------------------------------------
# Goal: Re-weight Trial 1 patients to look like Trial 2's population.
#
# Method: Fit a propensity score P(Trial 1 | covariates) by pooling Trial 1's
# real IPD with the synthetic Trial 2 population and fitting logistic regression.
# Then use odds (1-p)/p to weight Trial 1 patients.
#
# Note: This uses the SYNTHETIC Trial 2, not hidden IPD. A real analyst would
# replicate Trial 2's AgD arm-wise to mimic individual data, then pool with
# Trial 1 for the propensity model (alternative: use the synthetic population
# as we do here for numerical stability).

get_ipw_weights <- function(trial1_ipd, synthetic_trial2) {
  # Pool Trial 1 (source=1) and synthetic Trial 2 (source=0)
  pooled <- rbind(
    data.frame(trial1_ipd[, c("age", "sex", "bmi")], source = 1),
    data.frame(synthetic_trial2[, c("age", "sex", "bmi")], source = 0)
  )
  
  # Fit propensity score model
  ps_model <- glm(source ~ age + sex + bmi, data = pooled, family = binomial)
  
  # Extract propensity scores for Trial 1 patients only
  ps <- predict(ps_model, newdata = trial1_ipd[, c("age", "sex", "bmi")], 
                type = "response")
  
  # IPW: weight by (1-p)/p (down-weight high-p patients, up-weight low-p patients)
  (1 - ps) / ps
}

w_ipw <- get_ipw_weights(trial1_ipd, synth_trial2)
ess_ipw <- sum(w_ipw)^2 / sum(w_ipw^2)

# Fit weighted regression
ipw_mod <- lm(chg_bwt ~ ARM, data = trial1_ipd, weights = w_ipw)
d_AC_ipw <- coef(ipw_mod)["ARMPlacebo"] * -1

# Check correlation between the two weighting schemes
cat("\nCorrelation between MoM and IPW weights:", round(cor(w_mom, w_ipw), 3), "\n")


## ---------------------------------------------------------------------------
## STEP 7: Outcome Regression & Standardization (STC)
## ---------------------------------------------------------------------------
# Goal: Fit an outcome model on Trial 1 IPD with treatment-covariate
# interactions, then predict what the effect would be at Trial 2's
# population characteristics.
#
# Steps:
#   1. FIT: outcome model with treatment interactions (allows effect modification)
#   2. STANDARDIZE: predict for each synthetic Trial 2 patient under both
#      treatment and control, then average (integrating over target distribution)
#   3. CONTRAST: difference of the two averages = population-averaged effect

stc_model <- lm(chg_bwt ~ ARM * (age + sex + bmi), data = trial1_ipd)

# Predict for all synthetic Trial 2 patients under Drug A
synth_drugA <- synth_trial2; synth_drugA$ARM <- "Drug A"
pred_drugA <- predict(stc_model, newdata = synth_drugA)

# Predict for all synthetic Trial 2 patients under Placebo
synth_placebo <- synth_trial2; synth_placebo$ARM <- "Placebo"
pred_placebo <- predict(stc_model, newdata = synth_placebo)

# STC effect: marginal difference averaged over target population
d_AC_stc <- mean(pred_drugA) - mean(pred_placebo)

# Shortcut check (mean-centering): For linear models with identity link,
# E[f(X)] = f(E[X]), so centering covariates and reading the coefficient
# should give the same answer (Monte Carlo noise only)
trial1_centered <- trial1_ipd
trial1_centered$age_c <- trial1_ipd$age - target_means["age"]
trial1_centered$sex_c <- trial1_ipd$sex - target_means["sex"]
trial1_centered$bmi_c <- trial1_ipd$bmi - target_means["bmi"]
stc_shortcut_mod <- lm(chg_bwt ~ ARM * (age_c + sex_c + bmi_c), data = trial1_centered)
d_AC_stc_shortcut <- coef(stc_shortcut_mod)["ARMPlacebo"] * -1

cat("\nSTC full standardization:", round(d_AC_stc, 4), "\n")
cat("STC mean-centering shortcut:", round(d_AC_stc_shortcut, 4), 
    "(should closely match for linear models)\n")


## ---------------------------------------------------------------------------
## STEP 8: Bootstrap standard errors
## ---------------------------------------------------------------------------
# Resample Trial 1 IPD with replacement; keep synthetic Trial 2 fixed
# (it represents a known published population, not a random sample).

compute_all_estimates <- function(ipd_sample, synth_pop, target_cov) {
  # Naive
  mod_naive <- lm(chg_bwt ~ ARM, data = ipd_sample)
  d_naive <- coef(mod_naive)["ARMPlacebo"] * -1
  
  # MAIC MoM
  w1 <- get_maic_mom_weights(ipd_sample, target_cov)
  mod_mom <- lm(chg_bwt ~ ARM, data = ipd_sample, weights = w1)
  d_mom <- coef(mod_mom)["ARMPlacebo"] * -1
  
  # MAIC IPW
  w2 <- get_ipw_weights(ipd_sample, synth_pop)
  mod_ipw <- lm(chg_bwt ~ ARM, data = ipd_sample, weights = w2)
  d_ipw <- coef(mod_ipw)["ARMPlacebo"] * -1
  
  # STC (using mean-centering shortcut for speed)
  ipd_centered <- ipd_sample
  ipd_centered$age_c <- ipd_sample$age - target_cov["age"]
  ipd_centered$sex_c <- ipd_sample$sex - target_cov["sex"]
  ipd_centered$bmi_c <- ipd_sample$bmi - target_cov["bmi"]
  mod_stc <- lm(chg_bwt ~ ARM * (age_c + sex_c + bmi_c), data = ipd_centered)
  d_stc <- coef(mod_stc)["ARMPlacebo"] * -1
  
  c(naive = d_naive, mom = d_mom, ipw = d_ipw, stc = d_stc)
}

# Point estimates
point_ests <- compute_all_estimates(trial1_ipd, synth_trial2, target_means)

# Bootstrap
n_boot <- 500
boot_results <- matrix(NA_real_, nrow = n_boot, ncol = 4,
                       dimnames = list(NULL, names(point_ests)))

set.seed(2024)
for (b in seq_len(n_boot)) {
  idx <- sample(seq_len(nrow(trial1_ipd)), replace = TRUE)
  boot_results[b, ] <- compute_all_estimates(trial1_ipd[idx, ], synth_trial2, target_means)
}

se_estimates <- apply(boot_results, 2, sd)


## ---------------------------------------------------------------------------
## STEP 9: Summary table (Bucher: d_AB = d_AC - d_BC)
## ---------------------------------------------------------------------------
summary_table <- data.frame(
  Method = c("Naive", "MAIC (MoM)", "MAIC (IPW)", "STC"),
  d_AC = point_ests,
  SE_AC = se_estimates,
  d_BC = d_BC,
  SE_BC = se_d_BC,
  d_AB = point_ests - d_BC,
  SE_AB = sqrt(se_estimates^2 + se_d_BC^2)
)

summary_table$CI_lower <- summary_table$d_AB - 1.96 * summary_table$SE_AB
summary_table$CI_upper <- summary_table$d_AB + 1.96 * summary_table$SE_AB

cat("\n=== ANCHORED INDIRECT COMPARISON: DRUG A vs DRUG B ===\n")
cat("Endpoint: Change in body weight (kg, negative = loss)\n")
cat("Comparator: Placebo (anchor)\n\n")
print(summary_table, row.names = FALSE, digits = 3)

cat("\n=== EFFECTIVE SAMPLE SIZES ===\n")
cat("MAIC (MoM)  ESS:", round(ess_mom, 1), "out of", n1, "patients\n")
cat("MAIC (IPW)  ESS:", round(ess_ipw, 1), "out of", n1, "patients\n")


## ---------------------------------------------------------------------------
## STEP 10: Diagnostic plots
## ---------------------------------------------------------------------------
par(mfrow = c(1, 2))

# Weight distributions (rescaled to mean=1 for comparison)
hist(w_mom / mean(w_mom), breaks = 20, 
     main = "MAIC Weights\n(Method of Moments)",
     xlab = "Weight (rescaled, mean=1)", col = "skyblue", xlim = c(0, 15))

hist(w_ipw / mean(w_ipw), breaks = 30, 
     main = "MAIC Weights\n(IPW from synthetic AgD)",
     xlab = "Weight (rescaled, mean=1)", col = "lightcoral", xlim = c(0, 15))

par(mfrow = c(1, 1))

# Forest plot: point estimates and CIs
par(mar = c(4, 16, 3, 2))
y_pos <- 1:nrow(summary_table)
x_lim <- range(c(summary_table$CI_lower, summary_table$CI_upper, true_dAB))

plot(summary_table$d_AB, y_pos, 
     xlim = x_lim, yaxt = "n", ylab = "", 
     xlab = "d_AB: Drug A vs Drug B (kg)",
     pch = 19, cex = 1.2, ylim = c(0.5, nrow(summary_table) + 0.7),
     main = "Anchored Indirect Comparison: All Methods")

axis(2, at = y_pos, labels = summary_table$Method, las = 1, cex.axis = 0.9)
segments(summary_table$CI_lower, y_pos, summary_table$CI_upper, y_pos, lwd = 2)

abline(v = true_dAB, lty = 2, col = "red", lwd = 2)
text(true_dAB, nrow(summary_table) + 0.55, "Simulation truth", col = "red", cex = 0.85)

par(mar = c(5, 4, 4, 2))  # reset margins


## ---------------------------------------------------------------------------
## FINAL VALIDATION (for this simulation only)
## ---------------------------------------------------------------------------
cat("\n=== SIMULATION GROUND-TRUTH CHECK ===\n")
cat("True simulated d_AB:", true_dAB, "kg\n")
cat("Compare to d_AB column in table above.\n")
cat("\nNote: The hidden Trial 2 IPD was used ONLY to generate the aggregate\n")
cat("summary (means, SDs, arm sizes) and this final validation.\n")
cat("The analysis itself used only published aggregate data.\n")
