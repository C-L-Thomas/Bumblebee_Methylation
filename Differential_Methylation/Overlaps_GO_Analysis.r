   # =========================================================
  # GO term overlap analysis with colony excluded
  # - reads GO enrichment tables
  # - extracts unique GO.ID sets
  # - computes pairwise overlaps
  # - stores actual shared GO.IDs and Terms
  # - makes pairwise overlap heatmap
  # - makes UpSet plot
  # - returns Decitabine- and 6-aza-related reproductive overlaps
  # =========================================================
  
  library(readr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(UpSetR)
  
  # -----------------------------
  # Read files
  # -----------------------------
    
    
    colony <- read_delim(
      "DSS_topGO_Colony_effect_BP.tsv",
      delim = "\t", escape_double = FALSE, trim_ws = TRUE
    )
    
  reprovssterile <- read_delim(
    "DSS_topGO_repro_vs_sterile_BP.tsv",
    delim = "\t", escape_double = FALSE, trim_ws = TRUE
  )
  
  sixvsdeci <- read_delim(
    "DSS_topGO_5aza_vs_6aza_BP.tsv",
    delim = "\t", escape_double = FALSE, trim_ws = TRUE
  )
  
  decivscontrol <- read_delim(
    "DSS_topGO_Control_vs_5aza_BP.tsv",
    delim = "\t", escape_double = FALSE, trim_ws = TRUE
  )
  
  sixvscontrol <- read_delim(
    "DSS_topGO_Control_vs_6aza_BP.tsv",
    delim = "\t", escape_double = FALSE, trim_ws = TRUE
  )
  
  # -----------------------------
  # GO.ID lookup table
  # -----------------------------
  go_lookup <- bind_rows(
    reprovssterile %>% dplyr::select(GO.ID, Term),
    sixvsdeci %>% dplyr::select(GO.ID, Term),
    decivscontrol %>% dplyr::select(GO.ID, Term),
    sixvscontrol %>% dplyr::select(GO.ID, Term)
  ) %>%
    dplyr::filter(!is.na(GO.ID)) %>%
    distinct(GO.ID, .keep_all = TRUE)
  
  # -----------------------------
  # Extract unique GO.ID vectors
  # -----------------------------
  go_lists <- list(
    repro_vs_sterile = unique(reprovssterile$GO.ID),
    six_vs_deci      = unique(sixvsdeci$GO.ID),
    deci_vs_control  = unique(decivscontrol$GO.ID),
    six_vs_control   = unique(sixvscontrol$GO.ID)
  )
  
  go_lists <- lapply(go_lists, function(x) x[!is.na(x)])
  
  # -----------------------------
  # Pairwise overlaps: actual GO.IDs
  # -----------------------------
  pair_names <- combn(names(go_lists), 2, simplify = FALSE)
  
  pairwise_overlap_ids <- purrr::map(
    pair_names,
    ~ intersect(go_lists[[.x[1]]], go_lists[[.x[2]]])
  ) %>%
    set_names(purrr::map_chr(pair_names, ~ paste(.x, collapse = "___")))
  
  # -----------------------------
  # Pairwise overlaps as tidy dataframe with Term
  # -----------------------------
  pairwise_overlap_df <- purrr::imap_dfr(pairwise_overlap_ids, function(ids, nm) {
    parts <- strsplit(nm, "___")[[1]]
    tibble(
      comparison_1 = parts[1],
      comparison_2 = parts[2],
      GO.ID = ids
    )
  }) %>%
    dplyr::left_join(go_lookup, by = "GO.ID") %>%
    distinct()
  
  # -----------------------------
  # Helper function: get overlap with GO.ID + Term
  # -----------------------------
  get_overlap_terms <- function(df, x, y) {
    df %>%
      dplyr::filter(
        (comparison_1 == x & comparison_2 == y) |
          (comparison_1 == y & comparison_2 == x)
      ) %>%
      dplyr::select(GO.ID, Term) %>%
      distinct()
  }
  
  # -----------------------------
  # Pairwise overlap counts and Jaccard index
  # -----------------------------
  pairwise_count_df <- expand.grid(
    comparison_1 = names(go_lists),
    comparison_2 = names(go_lists),
    stringsAsFactors = FALSE
  ) %>%
    rowwise() %>%
    mutate(
      n_overlap = length(intersect(
        go_lists[[comparison_1]],
        go_lists[[comparison_2]]
      )),
      n_1 = length(go_lists[[comparison_1]]),
      n_2 = length(go_lists[[comparison_2]]),
      jaccard = n_overlap / (n_1 + n_2 - n_overlap)
    ) %>%
    ungroup()
  
  pairwise_counts <- pairwise_count_df %>%
    dplyr::select(comparison_1, comparison_2, n_overlap) %>%
    tidyr::pivot_wider(names_from = comparison_2, values_from = n_overlap)
  
  # -----------------------------
  # Save overlap tables
  # -----------------------------
  write_csv(
    pairwise_overlap_df,
    "pairwise_overlap_GO_IDs_and_terms_no_colony.csv"
  )
  
  write_csv(
    pairwise_count_df,
    "pairwise_overlap_counts_no_colony.csv"
  )
  
  write_csv(
    pairwise_counts,
    "pairwise_overlap_matrix_no_colony.csv"
  )
  
  # -----------------------------
  # Heatmap of pairwise overlap counts
  # -----------------------------
  p_heatmap <- ggplot(pairwise_count_df, aes(comparison_1, comparison_2, fill = n_overlap)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_overlap), size = 4) +
    theme_minimal(base_size = 11) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Shared\nGO IDs"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  print(p_heatmap)
  
  ggsave(
    "pairwise_overlap_heatmap_no_colony.pdf",
    p_heatmap,
    width = 5.8,
    height = 4.8
  )
  
  # -----------------------------
  # Build binary membership table for UpSetR
  # -----------------------------
  all_go_ids <- sort(unique(unlist(go_lists)))
  
  upset_df <- tibble(GO.ID = all_go_ids)
  
  for (nm in names(go_lists)) {
    upset_df[[nm]] <- as.integer(all_go_ids %in% go_lists[[nm]])
  }
  
  write_csv(
    upset_df,
    "upset_input_GO_IDs_no_colony.csv"
  )
  
  # -----------------------------
  # UpSet plot
  # -----------------------------
  pdf(
    "GO_overlap_upset_no_colony.pdf",
    width = 8.5,
    height = 5.2
  )
  
  suppressWarnings({
    upset(
      as.data.frame(upset_df[, -1]),
      sets = names(go_lists),
      order.by = "freq",
      keep.order = TRUE,
      mb.ratio = c(0.65, 0.35),
      text.scale = 1.2
    )
  })
  
  dev.off()
  
  # -----------------------------
  # Extract biologically useful overlaps
  # -----------------------------
  six_repro_df  <- get_overlap_terms(pairwise_overlap_df, "six_vs_control", "repro_vs_sterile")
  deci_repro_df <- get_overlap_terms(pairwise_overlap_df, "deci_vs_control", "repro_vs_sterile")
  six_deci_df   <- get_overlap_terms(pairwise_overlap_df, "six_vs_deci", "repro_vs_sterile")
  
  # Save them
  write_csv(
    six_repro_df,
    "six_control_repro_overlap_GO_IDs_and_terms_no_colony.csv"
  )
  
  write_csv(
    deci_repro_df,
    "deci_control_repro_overlap_GO_IDs_and_terms_no_colony.csv"
  )
  
  write_csv(
    six_deci_df,
    "six_deci_repro_overlap_GO_IDs_and_terms_no_colony.csv"
  )
  
  # -----------------------------
  # Decitabine/reproductive overlap unique versus 6-aza/reproductive overlap
  # -----------------------------
  deci_unique_df <- deci_repro_df %>%
    dplyr::filter(!GO.ID %in% six_repro_df$GO.ID) %>%
    distinct()
  
  write_csv(
    deci_unique_df,
    "deci_unique_repro_overlap_GO_IDs_and_terms_no_colony.csv"
  )
  
  # -----------------------------
  # 6-aza/reproductive overlap unique versus Decitabine/reproductive overlap
  # -----------------------------
  six_unique_df <- six_repro_df %>%
    dplyr::filter(!GO.ID %in% deci_repro_df$GO.ID) %>%
    distinct()
  
  write_csv(
    six_unique_df,
    "six_unique_repro_overlap_GO_IDs_and_terms_no_colony.csv"
  )
  
  # -----------------------------
  # Exact two-way only intersections (excluding all other sets)
  # -----------------------------
  deci_repro_only_df <- upset_df %>%
    dplyr::filter(
      deci_vs_control == 1,
      repro_vs_sterile == 1,
      six_vs_control == 0,
      six_vs_deci == 0
    ) %>%
    dplyr::left_join(go_lookup, by = "GO.ID") %>%
    dplyr::select(GO.ID, Term) %>%
    distinct()
  
  six_repro_only_df <- upset_df %>%
    dplyr::filter(
      six_vs_control == 1,
      repro_vs_sterile == 1,
      deci_vs_control == 0,
      six_vs_deci == 0
    ) %>%
    dplyr::left_join(go_lookup, by = "GO.ID") %>%
    dplyr::select(GO.ID, Term) %>%
    distinct()
  
  write_csv(
    deci_repro_only_df,
    "deci_repro_exact_two_way_only_GO_IDs_and_terms_no_colony.csv"
  )
  
  write_csv(
    six_repro_only_df,
    "six_repro_exact_two_way_only_GO_IDs_and_terms_no_colony.csv"
  )
  
  # -----------------------------
  # Print key objects
  # -----------------------------
  print(pairwise_overlap_df)
  print(pairwise_count_df)
  print(pairwise_counts)
  print(six_repro_df)
  print(deci_repro_df)
  print(deci_unique_df)
  print(six_unique_df)
  print(deci_repro_only_df)
  print(six_repro_only_df)
  
  
  
  
  # ============================================================
  # Hypergeometric test for GO term overlap
  #
  # Tests whether the overlap between drug-associated and
  # reproductive-status-associated enriched GO terms exceeds
  # chance, and whether the Decitabine-reproductive overlap is
  # significantly larger than the 6-azacytidine-reproductive overlap.
  #
  # Approach:
  #   - Universe = union of all GO terms appearing in any of the
  #     four enrichment tables (conservative; true universe would
  #     be all GO terms tested in topGO, including non-significant
  #     ones — if you have the full topGO output, swap in those
  #     counts for N below)
  #   - phyper() for each drug-repro overlap
  #   - Fisher's exact test (2x2) for each drug-repro overlap
  #   - Direct comparison of Decitabine vs 6-aza overlap sizes
  #     via Fisher's exact test on a 2x2 contingency table
  # ============================================================
  
  library(readr)
  library(dplyr)
  
  # ------------------- PATHS ----------------------------------
  # Adjust base_path if running from a different working directory
  base_path <- "GO_analysis"
  
  reprovssterile <- read_delim(
    file.path(base_path, "DSS_topGO_repro_vs_sterile_BP.tsv"),
    delim = "\t", escape_double = FALSE, trim_ws = TRUE
  )
  
  decivscontrol <- read_delim(
    file.path(base_path, "DSS_topGO_Control_vs_5aza_BP.tsv"),
    delim = "\t", escape_double = FALSE, trim_ws = TRUE
  )
  
  sixvscontrol <- read_delim(
    file.path(base_path, "DSS_topGO_Control_vs_6aza_BP.tsv"),
    delim = "\t", escape_double = FALSE, trim_ws = TRUE
  )
  
  sixvsdeci <- read_delim(
    file.path(base_path, "DSS_topGO_5aza_vs_6aza_BP.tsv"),
    delim = "\t", escape_double = FALSE, trim_ws = TRUE
  )
  
  # ------------------- EXTRACT ENRICHED SETS ------------------
  go_repro <- unique(na.omit(reprovssterile$GO.ID))
  go_deci  <- unique(na.omit(decivscontrol$GO.ID))
  go_six   <- unique(na.omit(sixvscontrol$GO.ID))
  go_sd    <- unique(na.omit(sixvsdeci$GO.ID))
  
  # ------------------- UNIVERSE -------------------------------
  # Union of all enriched terms across all four contrasts.
  # NOTE: this is a lower bound on the true universe. If you have
  # the full topGO result tables (all tested terms, not just
  # significant ones), use nrow() of those tables instead of
  # length(go_universe) below for a more conservative test.
  go_universe <- unique(c(go_repro, go_deci, go_six, go_sd))
  N <- length(go_universe)
  
  cat("Universe size (union of all enriched GO terms):", N, "\n")
  cat("Enriched in repro_vs_sterile:", length(go_repro), "\n")
  cat("Enriched in deci_vs_control: ", length(go_deci), "\n")
  cat("Enriched in six_vs_control:  ", length(go_six), "\n\n")
  
  # ------------------- OVERLAPS --------------------------------
  overlap_deci <- intersect(go_repro, go_deci)
  overlap_six  <- intersect(go_repro, go_six)
  
  cat("Observed overlap — Decitabine x Reproductive:", length(overlap_deci), "\n")
  cat("Observed overlap — 6-azacytidine x Reproductive:", length(overlap_six), "\n\n")
  
  # ------------------- HYPERGEOMETRIC TESTS -------------------
  # phyper(q, m, n, k, lower.tail = FALSE)
  #   q = observed overlap - 1 (number of white balls drawn - 1)
  #   m = K  = size of "success" set in universe (repro GO terms)
  #   n = N - K = remainder of universe
  #   k = number drawn = size of drug GO term set
  
  hg_deci <- phyper(
    q            = length(overlap_deci) - 1,
    m            = length(go_repro),
    n            = N - length(go_repro),
    k            = length(go_deci),
    lower.tail   = FALSE
  )
  
  hg_six <- phyper(
    q            = length(overlap_six) - 1,
    m            = length(go_repro),
    n            = N - length(go_repro),
    k            = length(go_six),
    lower.tail   = FALSE
  )
  
  cat("--- Hypergeometric test ---\n")
  cat("Decitabine x Reproductive  p =", formatC(hg_deci, format = "e", digits = 3), "\n")
  cat("6-azacytidine x Reproductive p =", formatC(hg_six, format = "e", digits = 3), "\n\n")
  
  # ------------------- FISHER'S EXACT TESTS -------------------
  # 2x2 table for each drug:
  #
  #                  | in repro | not in repro
  #  in drug         |    a     |      b
  #  not in drug     |    c     |      d
  #
  # Where a + b + c + d = N (universe)
  
  fisher_2x2 <- function(go_drug, go_repro, N) {
    a <- length(intersect(go_drug, go_repro))
    b <- length(setdiff(go_drug, go_repro))
    c <- length(setdiff(go_repro, go_drug))
    d <- N - a - b - c
    mat <- matrix(c(a, b, c, d), nrow = 2,
                  dimnames = list(c("in_drug", "not_in_drug"),
                                  c("in_repro", "not_in_repro")))
    list(table = mat, test = fisher.test(mat, alternative = "greater"))
  }
  
  ft_deci <- fisher_2x2(go_deci, go_repro, N)
  ft_six  <- fisher_2x2(go_six,  go_repro, N)
  
  cat("--- Fisher's exact test (one-sided: overlap > chance) ---\n")
  cat("Decitabine x Reproductive\n")
  print(ft_deci$table)
  cat("  p =", formatC(ft_deci$test$p.value, format = "e", digits = 3),
      " OR =", round(ft_deci$test$estimate, 2), "\n\n")
  
  cat("6-azacytidine x Reproductive\n")
  print(ft_six$table)
  cat("  p =", formatC(ft_six$test$p.value, format = "e", digits = 3),
      " OR =", round(ft_six$test$estimate, 2), "\n\n")
  
  # ------------------- DIRECT COMPARISON: DECI vs 6-AZA -------
  # Tests whether the Decitabine-reproductive overlap is significantly
  # larger than the 6-azacytidine-reproductive overlap.
  #
  # Each repro GO term is paired: it either overlaps with Decitabine
  # or not, AND either overlaps with 6-aza or not. These are paired
  # binary observations from the same repro set, so McNemar's test
  # is the appropriate approach — it tests whether the marginal
  # proportions differ (i.e., does deci capture more repro terms
  # than 6-aza does?).
  #
  # McNemar table (rows = in deci, cols = in 6-aza):
  #
  #                  | in 6-aza | not in 6-aza
  #  in deci         |    a     |      b
  #  not in deci     |    c     |      d
  #
  # b = in deci AND not in 6-aza (deci-unique repro overlap)
  # c = not in deci AND in 6-aza (6-aza-unique repro overlap)
  # McNemar chi-sq = (b - c)^2 / (b + c)
  
  a_mc <- length(intersect(intersect(go_deci, go_six), go_repro))  # both
  b_mc <- length(intersect(setdiff(go_deci, go_six), go_repro))    # deci only
  c_mc <- length(intersect(setdiff(go_six, go_deci), go_repro))    # six only
  d_mc <- length(setdiff(go_repro, union(go_deci, go_six)))         # neither
  
  mcnemar_mat <- matrix(c(a_mc, c_mc, b_mc, d_mc), nrow = 2,
                        dimnames = list(c("in_deci", "not_in_deci"),
                                        c("in_6aza", "not_in_6aza")))
  
  cat("--- Direct comparison: Decitabine vs 6-aza overlap with Reproductive ---\n")
  cat("(McNemar's test on paired repro GO terms)\n")
  print(mcnemar_mat)
  cat(sprintf("  b (deci-unique repro overlap) = %d\n", b_mc))
  cat(sprintf("  c (6-aza-unique repro overlap) = %d\n", c_mc))
  
  mc_result <- mcnemar.test(mcnemar_mat)
  cat("  McNemar chi-sq =", round(mc_result$statistic, 3),
      " df =", mc_result$parameter,
      " p =", formatC(mc_result$p.value, format = "e", digits = 3), "\n\n")
  
  # ------------------- SUMMARY TABLE --------------------------
  results <- data.frame(
    comparison       = c("Decitabine x Reproductive",
                         "6-azacytidine x Reproductive",
                         "Deci vs 6-aza (McNemar, within repro)"),
    n_drug_set       = c(length(go_deci), length(go_six), NA),
    n_repro_set      = c(length(go_repro), length(go_repro), length(go_repro)),
    n_overlap        = c(length(overlap_deci), length(overlap_six), NA),
    universe_N       = c(N, N, NA),
    hypergeometric_p = c(hg_deci, hg_six, NA),
    fisher_p         = c(ft_deci$test$p.value, ft_six$test$p.value, NA),
    fisher_OR        = c(ft_deci$test$estimate, ft_six$test$estimate, NA),
    mcnemar_p        = c(NA, NA, mc_result$p.value),
    mcnemar_chisq    = c(NA, NA, mc_result$statistic)
  )
  
  print(results)
  
  write_csv(
    results,
    file.path(base_path, "go_overlap_hypergeometric_results.csv")
  )
  
  cat("\nDone. Results written to go_overlap_hypergeometric_results.csv\n")
