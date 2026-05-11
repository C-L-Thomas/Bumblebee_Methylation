library(data.table)
library(lme4)
library(vegan)
library(survival)
library(splines)
library(uwot)
library(ggplot2)
library(dplyr)
library(rrcov)
library(ggbreak)
library(patchwork)

# ============================================================
# Bumblebee methylation structure and reproductive status
# - Compute CpG methylation proportions from count matrices
# - Reduce genome-wide CpG methylation to robust principal components
# - Visualise sample structure with UMAP on robust PCA space
# - Test association between methylation PCs and reproductive status
#   using conditional logistic regression on matched pairs
# - Test for global treatment shifts in methylation structure (PERMANOVA)
# ============================================================

# ------------------- SETTINGS -------------------
k_extract <- 5   # how many robust PCs to carry into dat

# ------------------- INPUT ----------------------
#C <- read.table(
#  "/PATH/methylated_counts_matrix.tsv",
#  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
#)# home

C <- read.table(
  "/PATH/methylated_counts_matrix.tsv",
  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
)# work

#N <- read.table(
#  "/PATH/coverage_counts_matrix.tsv",
#  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
#)

N <- read.table(
  "/PATH/coverage_counts_matrix.tsv",
  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
)#work

#covar <- fread(
#  "/PATH/covariates.tsv"
#)

covar <- fread(
  "/PATH/covariates.tsv"
)#work

# ------------------- ALIGN SAMPLES --------------
common <- intersect(colnames(C), covar$sample)

C <- C[, common, drop = FALSE]
N <- N[, common, drop = FALSE]
covar <- covar[match(common, covar$sample)]

stopifnot(identical(colnames(C), covar$sample))
stopifnot(identical(colnames(N), covar$sample))

# ------------------- METHYLATION PROPORTION -----
C_mat <- as.matrix(C)
storage.mode(C_mat) <- "double"

N_mat <- as.matrix(N)
storage.mode(N_mat) <- "double"

meth_prop <- C_mat / N_mat
meth_prop[!is.finite(meth_prop)] <- NA

# Impute missing per CpG with row mean
if (anyNA(meth_prop)) {
  rm <- rowMeans(meth_prop, na.rm = TRUE)
  na_rows <- which(rowSums(is.na(meth_prop)) > 0)
  for (i in na_rows) {
    meth_prop[i, is.na(meth_prop[i, ])] <- rm[i]
  }
}

# ------------------- MATRIX FOR PCA -------------
# samples in rows, CpGs in columns
pca_mat <- t(meth_prop)
rownames(pca_mat) <- covar$sample

# Remove CpGs with zero variance across samples
sd_cols <- apply(pca_mat, 2, sd)
keep <- is.finite(sd_cols) & sd_cols > 0
pca_mat <- pca_mat[, keep, drop = FALSE]

# ------------------- ROBUST PCA ------------------
# k must not exceed number of components available
k_extract <- min(k_extract, nrow(pca_mat) - 1, ncol(pca_mat))
pca_rob <- PcaHubert(pca_mat, k = k_extract)

# Extract robust PC scores
rob_scores <- as.data.frame(pca_rob@scores)
colnames(rob_scores) <- paste0("PC", seq_len(ncol(rob_scores)))
rob_scores$sample <- rownames(rob_scores)

# ------------------- MERGE METADATA -------------
dat <- merge(as.data.frame(covar), rob_scores, by = "sample", sort = FALSE)

if (!("treatment" %in% names(dat))) {
  if ("group" %in% names(dat)) {
    dat$treatment <- dat$group
  } else {
    stop(
      "Neither 'treatment' nor 'group' found in covariates.tsv. Columns are: ",
      paste(names(dat), collapse = ", ")
    )
  }
}

dat$status <- factor(dat$status)
dat$treatment <- factor(dat$treatment)
dat$colony <- factor(dat$colony)

###################
###################

# ------------------- ROBUST PCA PLOT ------------
ggplot(dat, aes(PC1, PC2, group = interaction(colony, treatment))) +
  geom_line(
    aes(linetype = treatment),
    colour = "grey30",
    linewidth = 0.7,
    alpha = 0.6
  ) +
  geom_point(
    aes(colour = colony, fill = status),
    shape = 21,
    size = 2.2,, alpha = 0.9,
    stroke = 1
  ) +
  scale_fill_manual(
    values = c("repro" = "black", "sterile" = "white")
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Robust PC1",
    y = "Robust PC2",
    colour = "Colony",
    fill = "Reproductive status",
    linetype = "Treatment"
  )


ggsave(
  "/PATH/robust_pca.pdf"
)

# ------------------- UMAP VISUALISATION ---------
# Use robust PC scores instead of standard PCA scores
n_umap_pcs <- min(10, ncol(pca_rob@scores))
pc_umap <- pca_rob@scores[, 1:n_umap_pcs, drop = FALSE]
pc_umap <- scale(pc_umap)

set.seed(1)

um <- umap(
  pc_umap,
  n_neighbors = 6,
  min_dist = 0.3,
  metric = "euclidean"
)

dat$UMAP1 <- um[, 1]
dat$UMAP2 <- um[, 2]

dat$fill_col <- ifelse(dat$status == "repro",
                       as.character(dat$colony),
                       NA)
dat$fill_col <- factor(dat$fill_col)

ggplot(dat, aes(UMAP1, UMAP2, group = interaction(colony, treatment))) +
  geom_line(
    aes(linetype = treatment),
    colour = "grey30",
    linewidth = 0.7,
    alpha = 0.7
  ) +
  geom_point(
    aes(colour = colony, fill = fill_col),
    shape = 21,
    size = 1.2, alpha = 0.9,
    stroke = 1.0
  ) +
  #scale_x_break(c(-8.5,-4.5,0, 16.5)) +
  scale_fill_discrete(na.value = "transparent") +
  theme_classic(base_size = 11) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    linetype = "Treatment"
  ) +
  guides(
    colour = "none",
    fill   = "none"
  )

ggsave(
  "/PATH/umap_robust_pca.pdf"
)






# ------------------- MATCHED-PAIR SETUP ---------
dat$colony <- factor(dat$colony)
dat$treatment <- factor(dat$treatment)

dat$status01 <- as.integer(dat$status == "repro")
dat$stratum <- interaction(dat$colony, dat$treatment, drop = TRUE)

with(dat, table(stratum, status01)) |> head()

# ------------------- CONDITIONAL LOGISTIC MODELS -----------------
m_cl <- clogit(status01 ~ PC1 + PC2 + strata(stratum), data = dat)
summary(m_cl)

m_spline2 <- clogit(
  status01 ~ ns(PC1, df = 2) + PC2 + strata(stratum),
  data = dat
)
summary(m_spline2)

anova(m_cl, m_spline2, test = "Chisq")

m_no_pc2 <- clogit(
  status01 ~ ns(PC1, df = 2) + strata(stratum),
  data = dat
)
anova(m_no_pc2, m_spline2, test = "Chisq")

m_spline3 <- clogit(
  status01 ~ ns(PC1, df = 2) + ns(PC2, df = 2) + strata(stratum),
  data = dat
)
anova(m_spline2, m_spline3)

m_null <- clogit(status01 ~ strata(stratum), data = dat)

data.frame(
  Model = c("Null", "Linear", "Spline", "Spline+PC2"),
  logLik = c(logLik(m_null), logLik(m_cl), logLik(m_no_pc2), logLik(m_spline2)),
  AIC = c(AIC(m_null), AIC(m_cl), AIC(m_no_pc2), AIC(m_spline2))
)

# ------------------- GLOBAL TREATMENT EFFECT ON METHYLATION STRUCTURE ----------
# PERMANOVA on robust PC space, blocked by colony
K <- min(10, ncol(pca_rob@scores))
D <- dist(scale(pca_rob@scores[, 1:K, drop = FALSE]))

adonis2(D ~ treatment, data = dat, permutations = 999, strata = dat$colony)




m_spline3df <- clogit(
  status01 ~ ns(PC1, df = 3) + PC2 + strata(stratum),
  data = dat
)

anova(m_spline2, m_spline3df)

m_quad <- clogit(
  status01 ~ PC1 + I(PC1^2) + PC2 + strata(stratum),
  data = dat
)

anova(m_cl, m_quad)
#AIC(m_spline2, m_quad)
summary(m_quad)






## ------------------------------------------------------------
## Robustness summary for m_quad
## Model: status01 ~ PC1 + I(PC1^2) + PC2 + strata(stratum)
## Produces:
##   1) leave-one-stratum-out coefficients
##   2) stratum bootstrap coefficients
##   3) compact summary table
## ------------------------------------------------------------


orig_coef <- coef(m_quad)

# ------------------------------------------------------------
# 1. Leave-one-stratum-out
# ------------------------------------------------------------
strata_levels <- unique(dat$stratum)

loso_list <- lapply(strata_levels, function(s) {
  fit <- try(
    clogit(
      status01 ~ PC1 + I(PC1^2) + PC2 + strata(stratum),
      data = dat[dat$stratum != s, ]
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(c(PC1 = NA_real_, `I(PC1^2)` = NA_real_, PC2 = NA_real_))
  }
  
  coef(fit)
})

loso_coefs <- do.call(rbind, loso_list)
rownames(loso_coefs) <- strata_levels

# Summaries for leave-one-stratum-out
loso_summary <- data.frame(
  term = colnames(loso_coefs),
  orig = orig_coef[colnames(loso_coefs)],
  mean = colMeans(loso_coefs, na.rm = TRUE),
  sd = apply(loso_coefs, 2, sd, na.rm = TRUE),
  median = apply(loso_coefs, 2, median, na.rm = TRUE),
  q025 = apply(loso_coefs, 2, quantile, probs = 0.025, na.rm = TRUE),
  q975 = apply(loso_coefs, 2, quantile, probs = 0.975, na.rm = TRUE),
  prop_positive = sapply(colnames(loso_coefs), function(x) {
    mean(loso_coefs[, x] > 0, na.rm = TRUE)
  }),
  n_success = sapply(colnames(loso_coefs), function(x) {
    sum(!is.na(loso_coefs[, x]))
  })
)

# ------------------------------------------------------------
# 2. Stratum bootstrap
# ------------------------------------------------------------
set.seed(1)

n_boot <- 1000

boot_list <- lapply(seq_len(n_boot), function(b) {
  samp_strata <- sample(strata_levels, length(strata_levels), replace = TRUE)
  
  boot_dat <- do.call(
    rbind,
    lapply(seq_along(samp_strata), function(i) {
      tmp <- dat[dat$stratum == samp_strata[i], ]
      tmp$stratum_boot <- paste0("boot_", i)
      tmp
    })
  )
  
  fit <- try(
    clogit(
      status01 ~ PC1 + I(PC1^2) + PC2 + strata(stratum_boot),
      data = boot_dat
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(c(PC1 = NA_real_, `I(PC1^2)` = NA_real_, PC2 = NA_real_))
  }
  
  coef(fit)
})

boot_coefs <- do.call(rbind, boot_list)

# Summaries for bootstrap
boot_summary <- data.frame(
  term = colnames(boot_coefs),
  orig = orig_coef[colnames(boot_coefs)],
  mean = colMeans(boot_coefs, na.rm = TRUE),
  sd = apply(boot_coefs, 2, sd, na.rm = TRUE),
  median = apply(boot_coefs, 2, median, na.rm = TRUE),
  q025 = apply(boot_coefs, 2, quantile, probs = 0.025, na.rm = TRUE),
  q975 = apply(boot_coefs, 2, quantile, probs = 0.975, na.rm = TRUE),
  prop_positive = sapply(colnames(boot_coefs), function(x) {
    mean(boot_coefs[, x] > 0, na.rm = TRUE)
  }),
  n_success = sapply(colnames(boot_coefs), function(x) {
    sum(!is.na(boot_coefs[, x]))
  })
)

# ------------------------------------------------------------
# 3. Combined compact table
# ------------------------------------------------------------
robustness_table <- bind_rows(
  mutate(loso_summary, method = "leave_one_stratum_out"),
  mutate(boot_summary, method = "stratum_bootstrap")
) %>%
  select(method, term, orig, mean, sd, median, q025, q975, prop_positive, n_success)

robustness_table






p_global <- ggplot(dat, aes(UMAP1, UMAP2, group = interaction(colony, treatment))) +
  geom_line(
    aes(linetype = treatment),
    colour = "grey30",
    linewidth = 0.7,
    alpha = 0.7
  ) +
  geom_point(
    aes(colour = colony, fill = fill_col),
    shape = 21,
    size = 1.8,
    alpha = 0.9,
    stroke = 0.8
  ) +
  scale_colour_discrete(guide = "none") +
  scale_fill_discrete(na.value = "transparent", guide = "none") +
  theme_classic(base_size = 11) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    linetype = "Treatment"
  )

p_facet <- ggplot(dat, aes(UMAP1, UMAP2, group = interaction(colony, treatment))) +
  geom_line(
    aes(linetype = treatment),
    colour = "grey30",
    linewidth = 0.7,
    alpha = 0.7
  ) +
  geom_point(
    aes(colour = colony, fill = fill_col),
    shape = 21,
    size = 1.8,
    alpha = 0.9,
    stroke = 0.8
  ) +
  facet_wrap(~ colony, scales = "free") +
  scale_colour_discrete(guide = "none") +
  scale_fill_discrete(na.value = "transparent", guide = "none") +
  theme_classic(base_size = 11) +
  theme(
    strip.text = element_blank()   # remove colony labels
  ) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    linetype = "Treatment"
  ) +
  guides(linetype = "none")

p_facet <- p_facet +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

combined_plot <- (p_global + p_facet) +
  plot_layout(widths = c(1.2, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") &   # adds A, B
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.tag = element_text(size = 14, face = "bold")  # style A/B
  )

combined_plot

ggsave(
  "/PATH/facet_umap_robust_pca.pdf",
  combined_plot,
  width = 11.4,
  height = 4.9
)




##################################
# ------------------- PREDICTED PROBABILITY CURVE FOR m_quad -------------------

b <- coef(m_quad)

# Vertex of the quadratic
vertex <- -b["PC1"] / (2 * b["I(PC1^2)"])

# Jitter the raw status points slightly for visibility
set.seed(42)
dat$jitter_status <- jitter(dat$status01, amount = 0.02)

# Global prediction grid
pc1_seq <- seq(min(dat$PC1), max(dat$PC1), length.out = 200)
pred_grid <- data.frame(PC1 = pc1_seq, PC2 = median(dat$PC2))
lp <- b["PC1"] * pred_grid$PC1 +
      b["I(PC1^2)"] * pred_grid$PC1^2 +
      b["PC2"] * pred_grid$PC2
pred_grid$prob <- plogis(lp - median(lp))

# Stratum-level prediction grid
pred_grid2 <- do.call(rbind, lapply(levels(dat$stratum), function(s) {
  sub <- dat[dat$stratum == s, ]
  pc1_seq <- seq(min(sub$PC1), max(sub$PC1), length.out = 100)
  lp <- b["PC1"]*pc1_seq + b["I(PC1^2)"]*pc1_seq^2 + b["PC2"]*median(dat$PC2)
  data.frame(stratum = s, PC1 = pc1_seq, prob = plogis(lp - median(lp)))
}))

# Colony-level prediction grid
pred_grid_colony <- do.call(rbind, lapply(levels(dat$colony), function(col) {
  sub <- dat[dat$colony == col, ]
  pc1_seq <- seq(min(sub$PC1), max(sub$PC1), length.out = 200)
  lp <- b["PC1"]*pc1_seq + b["I(PC1^2)"]*pc1_seq^2 + b["PC2"]*median(dat$PC2)
  data.frame(colony = col, PC1 = pc1_seq, prob = plogis(lp - median(lp)))
}))

# Segment data with arrow direction depending on vertex
segment_dat <- dat %>%
  group_by(stratum, colony) %>%
  summarise(
    x_sterile = PC1[status01 == 0],
    y_sterile = jitter_status[status01 == 0],
    x_repro   = PC1[status01 == 1],
    y_repro   = jitter_status[status01 == 1],
    .groups = "drop"
  ) %>%
  mutate(
    pc1_mid       = (x_sterile + x_repro) / 2,
    before_vertex = pc1_mid < vertex,
    x_start = ifelse(before_vertex, x_repro,   x_sterile),
    y_start = ifelse(before_vertex, y_repro,   y_sterile),
    x_end   = ifelse(before_vertex, x_sterile, x_repro),
    y_end   = ifelse(before_vertex, y_sterile, y_repro),
    x_mid   = (x_start + x_end) / 2,
    y_mid   = (y_start + y_end) / 2
  )

# ------------------- PLOT -------------------
ggplot() +
  # First half: start -> midpoint (carries the arrowhead)
  geom_segment(
    data = segment_dat,
    aes(x = x_start, y = y_start,
        xend = x_mid, yend = y_mid,
        colour = colony),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 0.6,
    alpha = 0.6
  ) +
  # Second half: midpoint -> end (no arrowhead)
  geom_segment(
    data = segment_dat,
    aes(x = x_mid, y = y_mid,
        xend = x_end, yend = y_end,
        colour = colony),
    linewidth = 0.6,
    alpha = 0.6
  ) +
  # Observed points
  geom_point(
    data = dat,
    aes(x = PC1, y = jitter_status, colour = colony, fill = colony),
    shape = 21,
    size = 2.2,
    alpha = 0.7,
    stroke = 0.8
  ) +
  # Predicted probability curve
  geom_line(
    data = pred_grid_colony,
    aes(x = PC1, y = prob),
    linewidth = 1.1,
    colour = "black"
  ) +
  geom_vline(xintercept = vertex, linetype = "dashed", colour = "grey50") +
  facet_wrap(~ colony, scales = "free_x") +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.5", "0.75", "1"),
    limits = c(-0.05, 1.05)
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text       = element_blank(),
    strip.background = element_blank(),
    legend.position  = "none"
  ) +
  labs(
    x = "Robust PC1",
    y = "Reproductive status / predicted probability"
  )

ggsave(
  "/PATH/pc1_quad_predicted_prob_by_colony.pdf",
  width = 12,
  height = 9
)








## ============================================================
## PREP: define the fitted minimum from the quadratic model
## ============================================================

# make sure m_quad already exists
# m_quad <- clogit(status01 ~ PC1 + I(PC1^2) + PC2 + strata(stratum), data = dat)

b <- coef(m_quad)
pc1_min <- -b["PC1"] / (2 * b["I(PC1^2)"])
pc1_min


## ============================================================
## 1. Compare residual multivariate dispersion of sterile vs reproductive workers
##    after accounting for the U-shaped PC1 relationship
## ============================================================

# align PC1 to fitted minimum
dat$PC1_aligned <- dat$PC1 - pc1_min
dat$PC1_dev <- abs(dat$PC1_aligned)

# residualise higher PCs against distance from optimum
resid_mat <- sapply(c("PC2", "PC3", "PC4", "PC5"), function(pc) {
  fit <- lm(dat[[pc]] ~ dat$PC1_dev)
  resid(fit)
})

resid_mat <- as.data.frame(resid_mat)
colnames(resid_mat) <- c("rPC2", "rPC3", "rPC4", "rPC5")

D_resid <- dist(scale(resid_mat))
bd_resid <- betadisper(D_resid, group = dat$status)

anova(bd_resid)
permutest(bd_resid, permutations = 999)

# optional plot
disp_dat <- dat
disp_dat$dist_to_centroid_resid <- bd_resid$distances

ggplot(disp_dat, aes(x = status, y = dist_to_centroid_resid, fill = status)) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.08, height = 0, alpha = 0.8, size = 2) +
  theme_classic(base_size = 12) +
  labs(
    x = "Reproductive status",
    y = "Residual distance to status centroid"
  ) +
  guides(fill = "none")


## ============================================================
## 3. Do the position of the minimum differ by colony?
## ============================================================

dat$colony <- factor(dat$colony)
dat$stratum <- interaction(dat$colony, dat$treatment, drop = TRUE)

m_quad_common <- clogit(
  status01 ~ PC1 + I(PC1^2) + PC2 + strata(stratum),
  data = dat
)

m_quad_colony_min <- clogit(
  status01 ~ PC1 + I(PC1^2) + PC2 +
    PC1:colony +
    strata(stratum),
  data = dat
)

anova(m_quad_common, m_quad_colony_min, test = "Chisq")
AIC(m_quad_common, m_quad_colony_min)
