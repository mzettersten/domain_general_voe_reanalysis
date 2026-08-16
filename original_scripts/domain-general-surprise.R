########################################################################
# Domain-General Surprise Analysis Script
#
# Description:
# This script analyzes infant looking time data from a violation-of-expectation 
# (VoE) paradigm battery to test whether surprise responses reflect a 
# domain-general latent factor.
#
# Analyses include:
# - Data cleaning and transformation
# - Confirmatory Factor Analysis (CFA)
# - Multivariate regressions predicting later curiosity traits
########################################################################


### Load Required Libraries ###
library(tidyverse)     # For data-wrangling and plotting
library(lavaan)        # For confirmatory factor analysis (CFA)
library(psych)         # For factor analysis diagnostics
library(semTools)      # For reliability estimates
library(brms)          # For Bayesian regression modeling
library(lm.beta)       # For standardized coefficients
library(ggpubr)
library(interactions)  # For Johnson-Neyman analysis


# ----------------------------------------------------------------------
# Load and Clean Raw VoE Data
# ----------------------------------------------------------------------

# a. Load raw frame-level dataset
raw_voe_data <- read.csv("~/path_to/raw_voe_data.csv")

# b. Aggregate to trial-level means for each subject/scene/outcome
primary_voe_data <- raw_voe_data %>%
  group_by(subj, age_months, gender, eu_order, domain, scene, outcome, test) %>%
  summarize(look = mean(look, na.rm = TRUE)) %>%
  ungroup() %>%
  spread(test, look) %>%
  dplyr::select(
    subject = subj,
    age_months,
    gender,
    outcome_order = eu_order,
    domain,
    task = scene,
    outcome,
    familiarization = fam,
    test_window = test
  )

# c. Retain only trials where infants looked at least 50% of the familiarization phase
primary_voe_data <- primary_voe_data %>%
  filter(familiarization > 0.5) %>%
  dplyr::select(-familiarization)  # Drop column after filter

# d. Compute log-looking time to normalize skew and reduce influence of long looks
# - Trials must exceed 5s minimum looking time
# - Compute difference score: unexpected - expected
window_size_s <- 30                      # Freeze frame duration per test trial (in seconds)
minimum_looking_time <- 5                # Minimum threshold (in seconds)


primary_voe_data <- primary_voe_data %>%
  mutate(look = log(window_size_s * test_window)) %>%
  dplyr::select(-test_window) %>%
  filter(look > log(minimum_looking_time))  # Exclude very low-looking trials

# e. Compute difference between looking to unexpected and expected outcomes
primary_voe_data <- primary_voe_data %>%
  spread(outcome, look) %>%
  mutate(diff_log = unexpected - expected) %>%
  na.omit() %>%
  ungroup()

# write.csv(primary_voe_data, "primary_voe_data.csv")


# ----------------------------------------------------------------------
# Demographics
# ----------------------------------------------------------------------

# ---- Initial VoE Sample ----
voe_demog <- primary_voe_data %>%
  group_by(subject) %>%
  summarize(age_months = first(age_months),
            gender = first(gender)) %>%
  mutate(gender = tolower(gender))

nrow(voe_demog)
min(voe_demog$age_months, na.rm = TRUE)
max(voe_demog$age_months, na.rm = TRUE)
mean(voe_demog$age_months, na.rm = TRUE)
sd(voe_demog$age_months, na.rm = TRUE)
sum(voe_demog$gender == "female", na.rm = TRUE)

# ---- Follow-up Sample (EMCS questionnaire) ----
# emcs_raw <- read.csv("~emcs_raw.csv")

emcs_joined <- emcs_raw %>%
  rename(subject = subj) %>%
  left_join(voe_demog, by = "subject") %>%
  filter(!is.na(age_months), !is.na(age_surv)) %>%
  mutate(age_surv = age_surv * 12) # conversion to months

nrow(emcs_joined)
min(emcs_joined$age_surv - emcs_joined$age_months, na.rm = TRUE)
max(emcs_joined$age_surv - emcs_joined$age_months, na.rm = TRUE)
mean(emcs_joined$age_surv - emcs_joined$age_months, na.rm = TRUE)
mean(emcs_joined$age_surv, na.rm = TRUE)
sd(emcs_joined$age_surv, na.rm = TRUE)
sum(emcs_joined$gender == "female", na.rm = TRUE)

# ----------------------------------------------------------------------
# T-Tests — Is surprise (unexpected > expected) present in each task?
# ----------------------------------------------------------------------

# Note: 'voe_data' should be 'primary_voe_data'
# Run one-sample one-sided t-tests on diff_log per task

t.test(primary_voe_data$diff_log[primary_voe_data$task == "solidity"], alternative = "greater")
t.test(primary_voe_data$diff_log[primary_voe_data$task == "support"], alternative = "greater")
t.test(primary_voe_data$diff_log[primary_voe_data$task == "continuity"], alternative = "greater")
t.test(primary_voe_data$diff_log[primary_voe_data$task == "containment"], alternative = "greater")
t.test(primary_voe_data$diff_log[primary_voe_data$task == "reach"], alternative = "greater")
t.test(primary_voe_data$diff_log[primary_voe_data$task == "snack"], alternative = "greater")
t.test(primary_voe_data$diff_log[primary_voe_data$task == "drink"], alternative = "greater")
t.test(primary_voe_data$diff_log[primary_voe_data$task == "fairness"], alternative = "greater")

# ----------------------------------------------------------------------
# Prepare Data for Confirmatory Factor Analysis (CFA)
# ----------------------------------------------------------------------

# Create wide-format dataset with one row per subject and columns for each task
# Values are z-scored difference scores (unexpected - expected) per task
dat.lav <- primary_voe_data %>%
  dplyr::select(subject, task, age_months, diff_log) %>%
  group_by(task) %>%
  mutate(diff_log = scale(diff_log)) %>%       # z-score within task
  spread(task, diff_log)

# ----------------------------------------------------------------------
# CFA — Compare 1-Factor vs. 2-Factor Models
# ----------------------------------------------------------------------

# a. Define CFA models
model_gen <- '
gen_factor =~ continuity + support + solidity + containment + reach + snack + drink + fairness
'

model_two <- '
physical =~ continuity + solidity + containment + support
agent    =~ reach + snack + drink + fairness
'

uncorrelated <- 'physical ~~ 0 * agent'

# b. Test factorability
cortest.bartlett(dat.lav[, 3:10])   # Bartlett’s test of sphericity
KMO(dat.lav[, 3:10])                # Kaiser-Meyer-Olkin measure of sampling adequacy

# c. Fit models using full-information maximum likelihood (FIML)
fit1.lav <- cfa(model_gen, missing = "FIML", estimator = "MLR", data = dat.lav, std.lv = TRUE)
fit2.lav <- cfa(model_two, missing = "FIML", estimator = "MLR", data = dat.lav, std.lv = TRUE)
fit3.lav <- cfa(c(model_two, uncorrelated), missing = "FIML", estimator = "MLR", data = dat.lav, std.lv = TRUE)


# d. Summarize model fits
summary(fit1.lav)
fitMeasures(fit1.lav, c("cfi.robust", "tli.robust", "rmsea.robust", "srmr", "aic", "bic"))

summary(fit2.lav)
fitMeasures(fit2.lav, c("cfi.robust", "tli.robust", "rmsea.robust", "srmr", "aic", "bic"))

summary(fit3.lav)
fitMeasures(fit3.lav, c("cfi.robust", "tli.robust", "rmsea.robust", "srmr", "aic", "bic"))

# e. Compare models
anova(fit1.lav, fit2.lav, fit3.lav)

# f. Estimate reliability (internal consistency of latent factor)
semTools::reliability(fit1.lav, what = c("alpha", "omega", "ave"))

# g. Bootstrapped single factor model
fit1.boot <- lavaan::sem(
  model_gen, data = dat.lav, se = "bootstrap", bootstrap = 2000,
  estimator = "ML", missing = "FIML", std.lv = TRUE
)
summary(fit1.boot, fit.measures = TRUE, standardized = TRUE, ci = TRUE)

# ----------------------------------------------------------------------
# Create Composite Surprise Scores by Domain
# ----------------------------------------------------------------------

dat.cross <- primary_voe_data %>%
  group_by(subject, age_months, domain) %>%
  summarize(diff.m = mean(diff_log, na.rm = TRUE)) %>%
  spread(domain, diff.m) %>%
  na.omit() %>%
  ungroup() %>%
  mutate(
    age.z    = scale(age_months),
    physical = scale(physical),
    social   = scale(social)
  )

# ----------------------------------------------------------------------
# Test Interaction Between Domains and Age
# ----------------------------------------------------------------------

mod_linear <- lm(physical ~ social + age.z, data = dat.cross)
mod_int    <- lm(physical ~ social * age.z, data = dat.cross)

# Likelihood ratio test: is the interaction model better?
anova(mod_linear, mod_int)

# Extract standardized coefficients and confidence intervals
lm.beta(mod_int)
confint(mod_int)

# ----------------------------------------------------------------------
# Visualize Domain Integration Across Age
# ----------------------------------------------------------------------

# A. Johnson-Neyman interaction plot
mod_int <- lm(physical ~ social * age_months, data = dat.cross)
jn <- johnson_neyman(model = mod_int, pred = social, modx = age_months)

cb <- jn$cbands %>%
  rename(age       = `age_months`,
         slope     = `Slope of social`,
         conf.low  = Lower,
         conf.high = Upper)

jn_plot <- ggplot(cb, aes(x = age, y = slope)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .15, fill = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = jn$bounds, colour = "firebrick",
             linetype = "dotted", size = 1) +
  labs(x = "Age (months)",
       y = "Between-Domain Correlation",
       title = "A") +
  xlim(12,30) +
  theme_classic() +
  theme(aspect.ratio = 1)

median_age <- median(dat.cross$age_months)

# B. Younger infants (< median age)
younger <- dat.cross %>%
  filter(age_months <= median_age) %>%
  ggplot(aes(x = scale(physical), y = social)) +
  geom_point(color = "blue2") +
  geom_smooth(method = 'lm', se = FALSE, color = "blue2") +
  labs(title = "B", x = "Physical Surprise Scores", y = "Social Surprise Scores") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  xlim(-2.5, 2.5) + ylim(-3.5, 3.5)

# C. Older infants (> median age)
older <- dat.cross %>%
  filter(age_months > median_age) %>%
  ggplot(aes(x = scale(physical), y = scale(social))) +
  geom_point(color = "blue3") +
  geom_smooth(method = 'lm', se = FALSE, color = "blue3") +
  labs(title = "C", x = "Physical Surprise Scores", y = "Social Surprise Scores") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  xlim(-2.5, 2.5) + ylim(-3.5, 3.5)

plot3 <- ggarrange(jn_plot, younger, older, ncol = 3, nrow = 1, align = 'h')
# install.packages("svglite")
ggsave(file = "plot3.svg", plot = plot3, width = 9, height = 3)

# ----------------------------------------------------------------------
# Merge VoE Data with Follow-Up Curiosity Questionnaire
# ----------------------------------------------------------------------

emcs_raw <- read.csv("~/path_to/emcs_raw.csv")

dat.brms <- left_join(primary_voe_data, emcs_raw, by = c("subject" = "subj")) %>%
  ungroup() %>%
  mutate(
    sociality         = scale(emcs_soc),
    info_seeking      = scale(emcs_seek),
    broad_exploration = scale(emcs_exp),
    persistence       = scale(emcs_per),
    vocab_tot         = scale(vocab_tot)
  ) %>%
  group_by(subject, age_months, vocab_tot,
           broad_exploration, info_seeking,
           persistence, sociality, age_surv) %>%
  summarize(
    gen_factor = mean(diff_log, na.rm = TRUE)  # average surprise score per subject
  ) %>%
  na.omit() %>%
  ungroup() %>%
  mutate(
    gen_factor   = scale(gen_factor),
    age_months.z = scale(age_months),
    age_surv.z   = scale(age_surv),
    age_split    = if_else(age_months.z > 0, "older", "younger")
  )

# ----------------------------------------------------------------------
# Multivariate Bayesian Regression — Does Early Surprise Predict Later Curiosity?
# ----------------------------------------------------------------------

set.seed(20250711)  # for reproducibility

mod_curiosity <- brm(
  bf(info_seeking      ~ gen_factor + gen_factor:age_months.z + age_surv.z) +
    bf(broad_exploration ~ gen_factor + gen_factor:age_months.z + age_surv.z) +
    bf(persistence       ~ gen_factor + gen_factor:age_months.z + age_surv.z) +
    bf(sociality         ~ gen_factor + gen_factor:age_months.z + age_surv.z) +
    bf(vocab_tot         ~ gen_factor + gen_factor:age_months.z + age_surv.z) +
    set_rescor(TRUE),
  data = dat.brms,
  family = gaussian(),
  chains = 4,
  cores = 4,
  iter = 10000,
  control = list(adapt_delta = 0.98)
)

# ----------------------------------------------------------------------
# Summarize Posterior Results
# ----------------------------------------------------------------------

# Convergence
brms::rhat(mod_curiosity) |> range()
brms::neff_ratio(mod_curiosity) |> range()

# Posterior summaries for all fixed effects, including vocabulary
summary(mod_curiosity)
posterior_summary(mod_curiosity)
conditional_effects(mod_curiosity, effects = "gen_factor:age_months.z")

write.csv(dat.brms, "dat.brms.csv")





