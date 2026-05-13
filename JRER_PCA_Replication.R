## JRER Revision: SLD PCA Replication
## Reads from original shapefile DBF (attribute table only — no geometry needed for PCA)
## Key addition: per-CBSA loading exports to verify Factor 1/2 consistency (addresses Reviewer 1 Comment 1)

require(foreign)
require(tidyverse)
require(FactoMineR)
require(factoextra)

# ── Paths ─────────────────────────────────────────────────────────────────────
DBF_PATH <- "C:/Users/jerem/OneDrive/Documents/JRER2026_GKRS/data/SLD_sf/SLD_sf/SmartLocationDb.dbf"
RCA_PATH <- "C:/Users/jerem/OneDrive/Documents/JRER2026_GKRS/data/RCAdata.csv"
OUT_DIR  <- "C:/Users/jerem/OneDrive/Documents/JRER2026_GKRS/results"

cat("=== Step 1: Load SLD from original shapefile DBF ===\n")
sld_raw <- read.dbf(DBF_PATH, as.is = TRUE)

# Drop metadata columns (matching original Urban_Form_PCA.R exactly)
sld_meta <- as_tibble(sld_raw) %>%
  select(-TRFIPS, -CFIPS, -SFIPS, -CBSA_EMP, -CBSA_POP, -CBSA_WRK,
         -Shape_Leng, -Shape_Area, -D5br_Flag, -D5be_Flag, -D1_flag)

# Recode missing-data placeholders
sld_meta[sld_meta == -99999] <- NA
sld_meta$CBSA[sld_meta$CBSA == 0] <- NA
sld_meta$CSA[sld_meta$CSA == 0]   <- NA

cat("SLD loaded:", nrow(sld_meta), "CBGs,", ncol(sld_meta), "columns\n")

# ── Step 2: Identify paper CBSAs (for loading consistency analysis only) ──────
cat("\n=== Step 2: Identify paper CBSAs from RCA ===\n")
rca <- read_csv(RCA_PATH, show_col_types = FALSE)
paper_cbsa_codes <- rca %>%
  filter(!is.na(CBSA_cd)) %>%
  distinct(CBSA_cd) %>%
  pull(CBSA_cd)
cat("Found", length(paper_cbsa_codes), "unique CBSAs in RCA data\n")

# ── Step 3: Prepare subsets ────────────────────────────────────────────────────
cat("\n=== Step 3: Prepare data subsets ===\n")

# Variables with no missing data across all CBGs
sld_all <- sld_meta %>%
  select(which(colSums(is.na(.)) == 0))
cat("sld_all:", nrow(sld_all), "rows,", ncol(sld_all), "cols\n")

# ALL CBGs in ANY CBSA, complete cases — matches original Urban_Form_PCA.R sld_cbsa_all
# (original ran PCA on all 900+ CBSAs so every US CBG gets a CBSA factor score)
sld_cbsa_all <- sld_meta %>%
  drop_na(CBSA) %>%
  select(which(colSums(is.na(.)) == 0))
cat("sld_cbsa_all:", nrow(sld_cbsa_all), "rows,", ncol(sld_cbsa_all), "cols\n")

# Paper-CBSAs-only subset (used only for loading consistency analysis in Step 7)
sld_paper_cbsa <- sld_meta %>%
  filter(CBSA %in% paper_cbsa_codes) %>%
  select(which(colSums(is.na(.)) == 0))
cat("sld_paper_cbsa:", nrow(sld_paper_cbsa), "rows,", ncol(sld_paper_cbsa), "cols\n")

# ── Step 4: National-scale PCA ─────────────────────────────────────────────────
cat("\n=== Step 4: National-scale PCA ===\n")
sld_pca_national <- sld_all

drop_cols <- c("GEOID10", "CSA", "CSA_Name", "CBSA", "CBSA_Name")
national_pca_out <- PCA(select(sld_pca_national, -any_of(drop_cols)),
                         scale.unit = TRUE, ncp = 5, graph = FALSE)

# Extract and save national loadings (correlations)
national_contrib <- get_pca_var(national_pca_out)
write.csv(national_contrib$coord,   file.path(OUT_DIR, "National_PCA_Correlations.csv"))
write.csv(national_contrib$contrib, file.path(OUT_DIR, "National_PCA_Contributions.csv"))

# Top 10 variables for Factor 1 and Factor 2 nationally
top_national <- as.data.frame(national_contrib$coord[, 1:2]) %>%
  rownames_to_column("variable") %>%
  rename(Factor1_loading = Dim.1, Factor2_loading = Dim.2) %>%
  arrange(desc(abs(Factor1_loading)))
cat("\nTop 10 variables by |loading| on national Factor 1:\n")
print(head(top_national[order(abs(top_national$Factor1_loading), decreasing=TRUE), ], 10))
cat("\nTop 10 variables by |loading| on national Factor 2:\n")
print(head(top_national[order(abs(top_national$Factor2_loading), decreasing=TRUE), ], 10))

# National factor scores
national_scores <- get_pca_ind(national_pca_out)$coord %>%
  as.data.frame() %>%
  bind_cols(select(sld_pca_national, GEOID10)) %>%
  select(GEOID10, everything()) %>%
  setNames(c("GEOID10", paste0("Factor_", 1:5, "_National")))

cat("National PCA complete.\n")

# ── Step 5: CBSA-scale PCA loop (paper CBSAs only) ────────────────────────────
# All CBSAs in the SLD (matches original Urban_Form_PCA.R scope)
all_cbsa_names <- sld_cbsa_all %>%
  distinct(CBSA, CBSA_Name) %>%
  drop_na() %>%
  arrange(CBSA_Name)
cat("\n=== Step 5: CBSA-scale PCA loop (", nrow(all_cbsa_names), "CBSAs) ===\n")

all_cbsa_scores    <- tibble()
all_cbsa_loadings  <- tibble()   # per-CBSA loadings (paper CBSAs only, for Step 7)

for (i in seq_len(nrow(all_cbsa_names))) {
  cbsa_code <- all_cbsa_names$CBSA[i]
  cbsa_name <- all_cbsa_names$CBSA_Name[i]

  cbsa_data <- sld_cbsa_all %>%
    filter(CBSA == cbsa_code) %>%
    select(-any_of(c("GEOID10", "CSA", "CSA_Name", "CBSA", "CBSA_Name")))

  if (nrow(cbsa_data) < 10) {
    cat("  Skipping", cbsa_name, "(too few CBGs:", nrow(cbsa_data), ")\n")
    next
  }

  if (i %% 50 == 0 || cbsa_code %in% paper_cbsa_codes) {
    cat(sprintf("  [%d/%d] %s (n=%d)\n", i, nrow(all_cbsa_names), cbsa_name, nrow(cbsa_data)))
  }

  pca_out <- PCA(cbsa_data, scale.unit = TRUE, ncp = 5, graph = FALSE)

  # Factor scores (use sld_cbsa_all for GEOIDs)
  scores <- pca_out$ind$coord %>%
    as.data.frame() %>%
    bind_cols(sld_cbsa_all %>% filter(CBSA == cbsa_code) %>% select(GEOID10)) %>%
    mutate(CBSA = cbsa_code, CBSA_Name = cbsa_name) %>%
    select(GEOID10, CBSA, CBSA_Name, everything()) %>%
    setNames(c("GEOID10", "CBSA", "CBSA_Name", paste0("Factor_", 1:5, "_CBSA_raw")))

  all_cbsa_scores <- bind_rows(all_cbsa_scores, scores)

  # Loadings — only capture for paper CBSAs (for consistency analysis in Step 7)
  if (cbsa_code %in% paper_cbsa_codes) {
    contrib <- get_pca_var(pca_out)
    loadings <- as.data.frame(contrib$coord[, 1:2]) %>%
      rownames_to_column("variable") %>%
      rename(Factor1_loading = Dim.1, Factor2_loading = Dim.2) %>%
      mutate(CBSA = cbsa_code, CBSA_Name = cbsa_name,
             abs_F1 = abs(Factor1_loading),
             abs_F2 = abs(Factor2_loading))
    all_cbsa_loadings <- bind_rows(all_cbsa_loadings, loadings)
  }
}

cat("CBSA PCA loop complete.\n")

# ── Step 6: Normalize CBSA factor scores to N(0,1) within each CBSA ───────────
cat("\n=== Step 6: Normalize CBSA factor scores ===\n")
normalise <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

cbsa_factors_norm <- all_cbsa_scores %>%
  group_by(CBSA_Name) %>%
  mutate(
    Factor_1_CBSA_Norm = normalise(Factor_1_CBSA_raw),
    Factor_2_CBSA_Norm = normalise(Factor_2_CBSA_raw),
    Factor_3_CBSA_Norm = normalise(Factor_3_CBSA_raw),
    Factor_4_CBSA_Norm = normalise(Factor_4_CBSA_raw),
    Factor_5_CBSA_Norm = normalise(Factor_5_CBSA_raw)
  ) %>%
  ungroup()

# ── Step 7: Consistency analysis (for Reviewer 1 Comment 1) ───────────────────
cat("\n=== Step 7: CBSA loading consistency analysis ===\n")

# For each CBSA, get rank of each variable on Factor 1 and Factor 2
cbsa_ranks <- all_cbsa_loadings %>%
  group_by(CBSA_Name) %>%
  mutate(
    rank_F1 = rank(-abs_F1),   # rank 1 = highest loading
    rank_F2 = rank(-abs_F2)
  ) %>%
  ungroup()

# Summary: median rank of key variables across CBSAs
key_f1_vars <- c("EMPTOT", "D1C", "D1C5_Svc10", "D1C8_Svc10", "E5_OFF10", "E8_OFF10",
                 "D1C5_Off10", "E5_SVC10", "E8_SVC10")
key_f2_vars <- c("PCT_AO0", "AUTOOWN0", "PCT_AO2P", "AUTOOWN2P", "D1A", "D1B",
                 "D2B_E8MIX", "D2B_E5MIX")

f1_consistency <- cbsa_ranks %>%
  filter(variable %in% key_f1_vars) %>%
  group_by(variable) %>%
  summarise(
    median_rank_F1 = median(rank_F1),
    mean_abs_loading_F1 = mean(abs_F1),
    pct_in_top10_F1 = mean(rank_F1 <= 10) * 100,
    n_cbsa = n()
  ) %>%
  arrange(median_rank_F1)

f2_consistency <- cbsa_ranks %>%
  filter(variable %in% key_f2_vars) %>%
  group_by(variable) %>%
  summarise(
    median_rank_F2 = median(rank_F2),
    mean_abs_loading_F2 = mean(abs_F2),
    pct_in_top10_F2 = mean(rank_F2 <= 10) * 100,
    n_cbsa = n()
  ) %>%
  arrange(median_rank_F2)

cat("\nFactor 1 key variable consistency across CBSAs:\n")
print(f1_consistency)
cat("\nFactor 2 key variable consistency across CBSAs:\n")
print(f2_consistency)

# ── Step 8: Save outputs ───────────────────────────────────────────────────────
cat("\n=== Step 8: Saving outputs ===\n")

write_csv(national_scores,
          file.path(OUT_DIR, "National_Factor_Scores_AllCBG.csv"))
write_csv(cbsa_factors_norm,
          file.path(OUT_DIR, paste0("SLD_CBSA_Factors_Normalized_", Sys.Date(), ".csv")))
write_csv(all_cbsa_loadings,
          file.path(OUT_DIR, "CBSA_PCA_Loadings_AllMarkets.csv"))
write_csv(f1_consistency,
          file.path(OUT_DIR, "CBSA_Factor1_Consistency.csv"))
write_csv(f2_consistency,
          file.path(OUT_DIR, "CBSA_Factor2_Consistency.csv"))

# Also merge with national scores for a complete factor file
all_factors <- left_join(
  cbsa_factors_norm %>% select(GEOID10, CBSA, CBSA_Name,
                                Factor_1_CBSA_Norm, Factor_2_CBSA_Norm,
                                Factor_3_CBSA_Norm, Factor_4_CBSA_Norm, Factor_5_CBSA_Norm),
  national_scores,
  by = "GEOID10"
)
write_csv(all_factors,
          file.path(OUT_DIR, paste0("SLD_Factors_Complete_", Sys.Date(), ".csv")))

cat("\n=== Done! Files saved to:", OUT_DIR, "===\n")
cat("Outputs:\n")
cat("  SLD_CBSA_Factors_Normalized_[date].csv  — normalized CBSA factor scores per CBG\n")
cat("  SLD_Factors_Complete_[date].csv          — CBSA + national factors merged\n")
cat("  CBSA_PCA_Loadings_AllMarkets.csv         — per-CBSA loadings for all variables\n")
cat("  CBSA_Factor1_Consistency.csv             — Factor 1 consistency summary\n")
cat("  CBSA_Factor2_Consistency.csv             — Factor 2 consistency summary\n")
cat("  National_PCA_Correlations.csv            — national-scale loadings\n")
cat("  National_PCA_Contributions.csv           — national-scale contributions\n")
