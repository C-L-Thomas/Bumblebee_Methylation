library(ggplot2)
library(glmmTMB)
library(dplyr)
library(splines)
library(emmeans)




setwd("/Users/emb3/Dropbox/Projects/kristi_behaviour")
comp_data <- read.csv(file = "complied_aggression.csv")

comp_data$Treatment   <- factor(comp_data$Treatment)
comp_data$Colony <- factor(comp_data$Colony)
comp_data$Box.Name    <- factor(comp_data$Box.Name)
comp_data$Batch <- factor(comp_data$Batch)
comp_data$Date <- as.Date(comp_data$Date, format = "%d/%m/%Y")

comp_data$DayNum <- as.integer(sub("D", "", comp_data$Day))

comp_data$Day <- as.integer(comp_data$Date - min(comp_data$Date, na.rm = TRUE))
summary(comp_data$Day)

comp_data$Treatment <- relevel(comp_data$Treatment, ref = "Control")

# m_total <- glmmTMB(
#   Total ~ Treatment + Time +
#     (1 | DayNum) +
#     (1 | Batch/Colony/Box.Name),
#   data = comp_data,
#   family = nbinom2
# )
# summary(m_total)
# performance::check_overdispersion(m_total)
# 
# 
# m_total_clean2 <- glmmTMB(
#   Total ~ Treatment + Time + poly(DayNum, 2) +
#     (1 | Colony/Box.Name),
#   data = comp_data,
#   family = nbinom2
# )
# summary(m_total_clean2)
# performance::check_overdispersion(m_total_clean2)



m_total_spline <- glmmTMB(
  Total ~ Treatment + Time + ns(DayNum, df = 3) +
    (1 | Colony/Box.Name),
  data = comp_data,
  family = nbinom2
)
summary(m_total_spline) #Aggression changes strongly and non-linearly across experimental days.
performance::check_overdispersion(m_total_spline)

#Aggressive behaviour counts were analysed using negative binomial mixed-effects models with treatment, time of day, and a natural spline for experimental day as fixed effects, and box nested within colony as random intercepts to account for repeated measures.

#Aggressive behaviour differed significantly among treatments. Relative to control boxes, aggression increased by approximately 36% under 6-azacytidine and more than doubled under Decitabine (negative binomial mixed-effects model, both p < 0.05). Aggression also varied strongly and non-linearly across experimental days, consistent with habituation effects, but did not differ systematically across times of day.




emm_treat <- emmeans(
  m_total_spline,
  ~ Treatment,
  type = "response"
)

emm_treat

pairs(emm_treat)

#Model-estimated aggressive behaviour increased from a mean of 1.48 events per observation in control boxes to 2.02 under 6-azacytidine and 3.28 under Decitabine. Pairwise comparisons indicated a strong increase under Decitabine relative to both control (rate ratio ≈ 2.21, p < 0.001) and 6-azacytidine (rate ratio ≈ 1.62, p < 0.001). The increase under 6-azacytidine relative to control was more modest (rate ratio ≈ 1.36) and marginal after Tukey correction (p = 0.067).




####################################
# --- Raw data dots + model lines/ribbon (population-level) ---

# Model predictions as you already built them
day_seq   <- sort(unique(comp_data$DayNum))
time_vals <- sort(unique(comp_data$Time))

newdat2 <- expand.grid(
  Treatment = levels(comp_data$Treatment),
  DayNum    = day_seq,
  Time      = time_vals
)

pred_link2 <- predict(
  m_total_spline,
  newdata = newdat2,
  type    = "link",
  se.fit  = TRUE,
  re.form = NA
)

newdat2 <- newdat2 %>%
  mutate(
    eta = pred_link2$fit,
    se  = pred_link2$se.fit,
    mu  = exp(eta),
    lcl = exp(eta - 1.96 * se),
    ucl = exp(eta + 1.96 * se)
  )

newdat_avg <- newdat2 %>%
  group_by(Treatment, DayNum) %>%
  summarise(
    mu  = mean(mu),
    lcl = mean(lcl),
    ucl = mean(ucl),
    .groups = "drop"
  )



# --- Marginal treatment means (right-hand inset) ---
emm_treat <- emmeans(m_total_spline, ~ Treatment, type = "response")
emm_df <- as.data.frame(emm_treat) %>%
  transmute(
    Treatment,
    mu  = response,
    lcl = asymp.LCL,
    ucl = asymp.UCL
  )

# --- x placement for the inset panel ---
x_min <- min(comp_data$DayNum, na.rm = TRUE)
x_max <- max(comp_data$DayNum, na.rm = TRUE)
x_rng <- x_max - x_min

x_inset_left  <- x_max + 0.03 * x_rng
x_inset_right <- x_max + 0.18 * x_rng
x_inset_mid   <- (x_inset_left + x_inset_right) / 2

emm_df$x <- x_inset_mid

# --- y helper for label placement ---
y_top <- max(c(comp_data$Total, newdat_avg$ucl, emm_df$ucl), na.rm = TRUE)

# --- colours (edit if you want different ones) ---
treat_cols <- c(
  "Control"       = "#333333",
  "6-azacytidine"  = "#1b9e77",
  "Decitabine"     = "#d95f02"
)

p <- ggplot() +
  
  # Inset background (white box)
  annotate(
    "rect",
    xmin = x_inset_left, xmax = x_inset_right,
    ymin = -Inf, ymax = Inf,
    fill = "white", colour = NA
  ) +
  
  # Separator line to mark inset boundary
  geom_vline(xintercept = x_inset_left, linewidth = 0.4, colour = "grey70") +
  
  # Model CI ribbon (behind)
  geom_ribbon(
    data = newdat_avg,
    aes(x = DayNum, ymin = lcl, ymax = ucl, fill = Treatment),
    alpha = 0.18,
    colour = NA
  ) +
  
  # Model mean line
  geom_line(
    data = newdat_avg,
    aes(x = DayNum, y = mu, colour = Treatment),
    linewidth = 1.0
  ) +
  
  # Raw observations
  geom_point(
    data = comp_data,
    aes(x = DayNum, y = Total, colour = Treatment),
    alpha = 0.18,
    size = 1.2,
    position = position_jitter(width = 0.10, height = 0)
  ) +
  
  # Marginal treatment means (95% CI) in inset
  geom_pointrange(
    data = emm_df,
    aes(x = x, y = mu, ymin = lcl, ymax = ucl, colour = Treatment),
    linewidth = 0.8
  ) +
  
  # Label for inset: left-aligned inside inset (won't overlap separator)
  annotate(
    "text",
    x = x_inset_left + 0.01 * x_rng,
    y = y_top * 1.04,
    label = "Marginal mean\n(95% CI)",
    hjust = 0,
    size = 3.2
  ) +
  
  scale_colour_manual(values = treat_cols) +
  scale_fill_manual(values = treat_cols) +
  
  labs(
    x = "Experimental day (within batch)",
    y = "Aggressive interactions (count)"
  ) +
  
  coord_cartesian(
    xlim = c(x_min, x_inset_right),
    ylim = c(0, NA),     # start y-axis at 0; remove if you don't want this
    clip = "off"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.margin = margin(5.5, 10, 5.5, 5.5)
  )

p

ggsave("aggression.pdf", p, width = 6.8, height = 4.6, units = "in")

