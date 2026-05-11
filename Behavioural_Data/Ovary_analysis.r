## =========================================================
## Ovary analysis: status, score, and length
## =========================================================

library(lme4)
library(ordinal)
library(censReg)
library(car)
library(emmeans)
library(dplyr)
library(ggplot2)
library(patchwork)

## ---- Data prep ----

dat <- read.csv(
  "~/Dropbox/Data/kristi_ovaries/reproductive_data_int.csv",
  stringsAsFactors = TRUE
)

lvl <- c("Control", "6-azacytidine", "Decitabine")
dat$Drug   <- factor(dat$Drug, levels = lvl)
dat$Colony <- factor(dat$Colony)
dat$Box    <- factor(dat$Box)

dat$y <- as.integer(dat$Reproductive == "Reproductive")
dat$Score_ord <- factor(dat$Score, levels = 0:4, ordered = TRUE)

treat_cols <- c(
  "Control"       = "#333333",
  "6-azacytidine" = "#1b9e77",
  "Decitabine"    = "#d95f02"
)

## =========================================================
## 1. Reproductive status (binary) -- GLMM
## =========================================================

m_status <- glmer(
  y ~ Drug + (1 | Colony/Box),
  family  = binomial,
  data    = dat,
  control = glmerControl(optimizer = "bobyqa")
)

summary(m_status)
car::Anova(m_status, type = "II")

emm_status <- emmeans(m_status, ~ Drug, type = "response")
pairs(emm_status, adjust = "tukey")

status_df <- as.data.frame(emm_status) %>%
  transmute(Drug, estimate = prob, lower = asymp.LCL, upper = asymp.UCL) %>%
  mutate(Drug = factor(Drug, levels = lvl))

## =========================================================
## 2. Reproductive score (0-4 ordinal) -- CLMM
## =========================================================

m_score <- clmm(
  Score_ord ~ Drug + (1 | Colony/Box),
  data = dat,
  link = "logit"
)

summary(m_score)

m_score_null <- clmm(
  Score_ord ~ 1 + (1 | Colony/Box),
  data = dat,
  link = "logit"
)
anova(m_score_null, m_score)

emm_score <- emmeans(m_score, ~ Drug, mode = "latent")
pairs(emm_score, adjust = "tukey")

score_df <- as.data.frame(emm_score) %>%
  transmute(Drug, estimate = emmean, lower = asymp.LCL, upper = asymp.UCL) %>%
  mutate(Drug = factor(Drug, levels = lvl))


emm_cat  <- emmeans(m_score, ~ Drug, mode = "mean.class")
score_df <- as.data.frame(emm_cat) %>%
  transmute(Drug,
            estimate = mean.class,
            lower    = asymp.LCL,
            upper    = asymp.UCL) %>%
  mutate(Drug = factor(Drug, levels = lvl))

## =========================================================
## 3. Ovary length -- Tobit (left-censored at 0)
## =========================================================

m_length <- censReg(
  Length ~ Drug + Colony,
  left = 0,
  data = dat
)

summary(m_length)
margEff(m_length)

## Extract Drug estimates with 95% CIs on the latent scale.
## censReg doesn't play with emmeans, so we pull coefficients directly
## and compute marginal means relative to the Control reference.
coefs  <- coef(summary(m_length))
vc     <- vcov(m_length)

# Intercept = Control mean (at reference Colony level)
mu_control  <- coefs["(Intercept)", "Estimate"]
se_control  <- coefs["(Intercept)", "Std. error"]

# Other treatments = intercept + Drug coefficient
b_aza  <- coefs["Drug6-azacytidine", "Estimate"]
b_dec  <- coefs["DrugDecitabine", "Estimate"]

mu_aza <- mu_control + b_aza
mu_dec <- mu_control + b_dec

# SEs using the covariance matrix
se_aza <- sqrt(vc["(Intercept)", "(Intercept)"] +
               vc["Drug6-azacytidine", "Drug6-azacytidine"] +
               2 * vc["(Intercept)", "Drug6-azacytidine"])
se_dec <- sqrt(vc["(Intercept)", "(Intercept)"] +
               vc["DrugDecitabine", "DrugDecitabine"] +
               2 * vc["(Intercept)", "DrugDecitabine"])

length_df <- data.frame(
  Drug     = factor(lvl, levels = lvl),
  estimate = c(mu_control, mu_aza, mu_dec),
  se       = c(se_control, se_aza, se_dec)
) %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se
  )




library(lmtest)

m_length      <- censReg(Length ~ Drug + Colony, left = 0, data = dat)
m_length_null <- censReg(Length ~ Colony,        left = 0, data = dat)

lrtest(m_length_null, m_length)

## =========================================================
## Unified figure: three panels, one measure of ovary per panel
## =========================================================

base_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(size = 10, face = "plain")
  )

p_status <- ggplot(status_df, aes(Drug, estimate, colour = Drug)) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  scale_colour_manual(values = treat_cols) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "P(reproductive)",
       title = "Reproductive status (GLMM)") +
  base_theme



p_score <- ggplot(score_df, aes(Drug, estimate, colour = Drug)) +
  geom_jitter(
    data = dat,
    aes(x = Drug, y = as.numeric(as.character(Score)), colour = Drug),
    width = 0.08, height = 0.08, alpha = 0.25, inherit.aes = FALSE
  ) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  scale_colour_manual(values = treat_cols) +
  scale_y_continuous(
    breaks = 0:4,
    limits = c(-0.3, 4.3)
  ) +
  labs(x = NULL, y = "Reproductive score (0–4)",
       title = "Reproductive score (CLMM)") +
  base_theme





p_length <- ggplot(length_df, aes(Drug, estimate, colour = Drug)) +
  geom_jitter(
    data = dat,
    aes(x = Drug, y = Length, colour = Drug),
    width = 0.08, height = 0, alpha = 0.25, inherit.aes = FALSE
  ) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  scale_colour_manual(values = treat_cols) +
  labs(x = NULL, y = "Ovary length (mm)",
       title = "Ovary length (Tobit)") +
  base_theme

combined <- p_status + p_score + p_length +
  plot_annotation(tag_levels = "A")

ggsave(
  "~/Dropbox/Data/kristi_ovaries/ovary_three_panel.pdf",
  combined,
  width  = 10,
  height = 4
)
