## JRER_PCA_Replication_V3.R
##
## Replicates JRER_PCA_Replication.R using the EPA Smart Location Database
## Version 3.0 (May 2021), which uses 2020 Census block group boundaries.
##
## Key V3 differences from V2 (handled below):
##   - Source format: File Geodatabase (GDB), read with sf::read_sf()
##   - V3 carries both GEOID10 (2010 BG) and GEOID20 (2020 BG); 220,739 records
##   - Variable renames: employment/density vars dropped the "10" suffix
##     (e.g. EMPTOT -> TotEmp, E8_SVC10 -> E8_Svc, D1C8_Svc10 -> D1C8_SVC)
##   - New derived variables (VMT, GHG, model coefficients, walkability index)
##     are excluded from the PCA to preserve comparability with V2
##   - Removed: TRFIPS/CFIPS/SFIPS (now STATEFP/COUNTYFP/TRACTCE/BLKGRPCE),
##     all E_FED* federal employment variables, D5br_Flag, D5be_Flag
##
## Prerequisite: SmartLocationDatabaseV3.gdb must be extracted to the data
## folder. If starting from the zip, extract SmartLocationDatabaseV3.zip
## using Windows Explorer or:
##   Add-Type -AssemblyName System.IO.Compression.FileSystem
##   [System.IO.Compression.ZipFile]::ExtractToDirectory("...V3.zip", "...dst")
## then copy the resulting SmartLocationDatabase.gdb folder to:
##   data/SmartLocationDatabaseV3.gdb
##
## Output files are written with "_V3" suffixes so they never overwrite V2 results.

library(tidyverse)
library(sf)
library(FactoMineR)
library(factoextra)

# ── Paths ──────────────────────────────────────────────────────────────────────
GDB_PATH <- "C:/Users/jerem/OneDrive/Documents/JRER2026_GKRS/data/SmartLocationDatabaseV3.gdb"
RCA_PATH <- "C:/Users/jerem/OneDrive/Documents/JRER2026_GKRS/data/RCAdata.csv"
OUT_DIR  <- "C:/Users/jerem/OneDrive/Documents/JRER2026_GKRS/results"

normalise <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

# ── Step 1: Load V3 SLD from GDB ──────────────────────────────────────────────
cat("=== Step 1: Load V3 SLD from GDB ===\n")
cat("Reading layer EPA_SLD_Database_V3 (this may take a few minutes)...\n")
sld_raw <- read_sf(GDB_PATH, layer = "EPA_SLD_Database_V3")
sld_raw <- st_drop_geometry(sld_raw)   # attributes only; no geometry needed for PCA
cat("Loaded:", nrow(sld_raw), "rows,", ncol(sld_raw), "columns\n")

# Drop V3-specific derived/model variables that have no V2 equivalent,
# plus the standard metadata columns (updated names for V3).
# Keeps all core D1-D5 SLD variables for PCA comparability with V2.
v3_drop <- c(
  # Geographic identifiers / metadata
  "GEOID10", "GEOID20", "STATEFP", "COUNTYFP", "TRACTCE", "BLKGRPCE",
  "CSA", "CSA_Name", "CBSA", "CBSA_Name", "CBSA_EMP", "CBSA_POP", "CBSA_WRK",
  "Region",
  # Geometry
  "Shape_Length", "Shape_Area",
  # Flag
  "D1_FLAG",
  # Ranked index variables (derived from core SLD vars)
  "D2A_Ranked", "D2B_Ranked", "D3B_Ranked", "D4A_Ranked", "NatWalkInd",
  # Behavioral model inputs (demographics/prices not in V2 PCA)
  "Households", "Workers_1", "Residents", "Drivers", "Vehicles",
  "White", "Male", "Lowwage", "Medwage", "Highwage",
  "W_P_Lowwage", "W_P_Medwage", "W_P_Highwage", "GasPrice",
  # Log-transformed and per-capita intermediates
  "logd1a", "logd1c", "logd3aao", "logd3apo", "d4bo25",
  "d5dei_1", "logd4d", "UPTpercap",
  # Behavioral model coefficients
  "B_C_constant", "B_C_male", "B_C_ld1c", "B_C_drvmveh", "B_C_ld1a",
  "B_C_ld3apo", "B_C_inc1", "B_C_gasp",
  "B_N_constant", "B_N_inc2", "B_N_inc3", "B_N_white", "B_N_male",
  "B_N_drvmveh", "B_N_gasp", "B_N_ld1a", "B_N_ld1c", "B_N_ld3aao",
  "B_N_ld3apo", "B_N_d4bo25", "B_N_d5dei", "B_N_UPTpc",
  # Calibrated rates
  "C_R_Households", "C_R_Pop", "C_R_Workers", "C_R_Drivers", "C_R_Vehicles",
  "C_R_White", "C_R_Male", "C_R_Lowwage", "C_R_Medwage", "C_R_Highwage", "C_R_DrmV",
  # VMT / GHG outputs
  "NonCom_VMT_Per_Worker", "Com_VMT_Per_Worker", "VMT_per_worker",
  "VMT_tot_min", "VMT_tot_max", "VMT_tot_avg",
  "GHG_per_worker", "Annual_GHG", "SLC_score"
)

sld_meta <- as_tibble(sld_raw) %>%
  select(-any_of(v3_drop))

# Recode missing-data placeholder (-99999) to NA
sld_meta[sld_meta == -99999] <- NA

# Add GEOID10/CBSA/CSA back as join keys (not included in PCA).
# Convert to numeric first so character "0" is handled the same as integer 0.
sld_meta$GEOID10   <- sld_raw$GEOID10
sld_meta$CBSA      <- { x <- suppressWarnings(as.numeric(sld_raw$CBSA));
                         if_else(is.na(x) | x == 0, NA_real_, x) }
sld_meta$CBSA_Name <- sld_raw$CBSA_Name
sld_meta$CSA       <- { x <- suppressWarnings(as.numeric(sld_raw$CSA));
                         if_else(is.na(x) | x == 0, NA_real_, x) }
sld_meta$CSA_Name  <- sld_raw$CSA_Name

cat("After dropping derived/metadata columns:", ncol(sld_meta), "columns remain\n")

# ── Step 2: Identify paper CBSAs ──────────────────────────────────────────────
cat("\n=== Step 2: Identify paper CBSAs from RCA ===\n")
rca <- read_csv(RCA_PATH, show_col_types = FALSE)
paper_cbsa_codes <- rca %>%
  filter(!is.na(CBSA_cd)) %>%
  distinct(CBSA_cd) %>%
  pull(CBSA_cd)
cat("Found", length(paper_cbsa_codes), "unique CBSAs in RCA data\n")

# ── Step 3: Prepare subsets ────────────────────────────────────────────────────
cat("\n=== Step 3: Prepare data subsets ===\n")
meta_cols <- c("GEOID10", "CSA", "CSA_Name", "CBSA", "CBSA_Name")
drop_cols  <- meta_cols   # excluded from PCA inputs

# Keep meta columns always; drop only data columns that have any NA
keep_complete <- function(df) {
  data_cols <- setdiff(names(df), meta_cols)
  good      <- data_cols[colSums(is.na(df[data_cols])) == 0]
  select(df, any_of(meta_cols), all_of(good))
}

sld_all      <- keep_complete(sld_meta)
sld_cbsa_all <- keep_complete(filter(sld_meta, !is.na(CBSA)))
sld_csa_all  <- keep_complete(filter(sld_meta, !is.na(CSA)))
cat("sld_all:", nrow(sld_all), "rows,", ncol(sld_all), "cols\n")
cat("sld_cbsa_all:", nrow(sld_cbsa_all), "rows,", ncol(sld_cbsa_all), "cols\n")
cat("sld_csa_all:", nrow(sld_csa_all), "rows,", ncol(sld_csa_all), "cols\n")

# ── Step 4: National-scale PCA ─────────────────────────────────────────────────
cat("\n=== Step 4: National-scale PCA ===\n")
national_pca_out <- PCA(select(sld_all, -any_of(drop_cols)),
                        scale.unit = TRUE, ncp = 5, graph = FALSE)

national_contrib <- get_pca_var(national_pca_out)
write.csv(national_contrib$coord,   file.path(OUT_DIR, "National_PCA_Correlations_V3.csv"))
write.csv(national_contrib$contrib, file.path(OUT_DIR, "National_PCA_Contributions_V3.csv"))

top_national <- as.data.frame(national_contrib$coord[, 1:2]) %>%
  rownames_to_column("variable") %>%
  rename(Factor1_loading = Dim.1, Factor2_loading = Dim.2)
cat("\nTop 10 variables by |loading| on national Factor 1:\n")
print(head(top_national[order(abs(top_national$Factor1_loading), decreasing = TRUE), ], 10))
cat("\nTop 10 variables by |loading| on national Factor 2:\n")
print(head(top_national[order(abs(top_national$Factor2_loading), decreasing = TRUE), ], 10))

national_scores <- get_pca_ind(national_pca_out)$coord %>%
  as.data.frame() %>%
  bind_cols(select(sld_all, GEOID10)) %>%
  select(GEOID10, everything()) %>%
  setNames(c("GEOID10", paste0("Factor_", 1:5, "_National"))) %>%
  mutate(
    Factor_1_National_Norm = normalise(Factor_1_National),
    Factor_2_National_Norm = normalise(Factor_2_National)
  )

write_csv(national_scores, file.path(OUT_DIR, "National_Factor_Scores_AllCBG_V3.csv"))
cat("National PCA complete. Rows:", nrow(national_scores), "\n")

# ── Step 5: CBSA-scale PCA loop ────────────────────────────────────────────────
all_cbsa_names <- sld_cbsa_all %>%
  distinct(CBSA, CBSA_Name) %>%
  drop_na() %>%
  arrange(CBSA_Name)
cat("\n=== Step 5: CBSA PCA loop (", nrow(all_cbsa_names), "CBSAs) ===\n")

all_cbsa_scores   <- tibble()
all_cbsa_loadings <- tibble()

for (i in seq_len(nrow(all_cbsa_names))) {
  cbsa_code <- all_cbsa_names$CBSA[i]
  cbsa_name <- all_cbsa_names$CBSA_Name[i]

  cbsa_data <- sld_cbsa_all %>%
    filter(CBSA == cbsa_code) %>%
    select(-any_of(drop_cols))

  if (nrow(cbsa_data) < 10) {
    cat("  Skipping", cbsa_name, "(n =", nrow(cbsa_data), ")\n")
    next
  }
  if (i %% 50 == 0 || cbsa_code %in% paper_cbsa_codes)
    cat(sprintf("  [%d/%d] %s (n=%d)\n", i, nrow(all_cbsa_names), cbsa_name, nrow(cbsa_data)))

  pca_out <- PCA(cbsa_data, scale.unit = TRUE, ncp = 5, graph = FALSE)

  scores <- pca_out$ind$coord %>%
    as.data.frame() %>%
    bind_cols(sld_cbsa_all %>% filter(CBSA == cbsa_code) %>% select(GEOID10)) %>%
    mutate(CBSA = cbsa_code, CBSA_Name = cbsa_name) %>%
    select(GEOID10, CBSA, CBSA_Name, everything()) %>%
    setNames(c("GEOID10", "CBSA", "CBSA_Name", paste0("Factor_", 1:5, "_CBSA_raw")))

  all_cbsa_scores <- bind_rows(all_cbsa_scores, scores)

  if (cbsa_code %in% paper_cbsa_codes) {
    contrib <- get_pca_var(pca_out)
    loadings <- as.data.frame(contrib$coord[, 1:2]) %>%
      rownames_to_column("variable") %>%
      rename(Factor1_loading = Dim.1, Factor2_loading = Dim.2) %>%
      mutate(CBSA = cbsa_code, CBSA_Name = cbsa_name,
             abs_F1 = abs(Factor1_loading), abs_F2 = abs(Factor2_loading))
    all_cbsa_loadings <- bind_rows(all_cbsa_loadings, loadings)
  }
}
cat("CBSA PCA loop complete.\n")

# ── Step 6: Normalize CBSA scores ─────────────────────────────────────────────
cat("\n=== Step 6: Normalize CBSA factor scores ===\n")
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

write_csv(cbsa_factors_norm,
          file.path(OUT_DIR, paste0("SLD_CBSA_Factors_Normalized_V3_", Sys.Date(), ".csv")))
cat("CBSA normalized scores saved. Rows:", nrow(cbsa_factors_norm), "\n")

# ── Step 7: CSA-scale PCA loop ─────────────────────────────────────────────────
all_csa_names <- sld_csa_all %>%
  distinct(CSA, CSA_Name) %>%
  drop_na() %>%
  arrange(CSA_Name)
cat("\n=== Step 7: CSA PCA loop (", nrow(all_csa_names), "CSAs) ===\n")

all_csa_scores <- tibble()

for (i in seq_len(nrow(all_csa_names))) {
  csa_code <- all_csa_names$CSA[i]
  csa_name <- all_csa_names$CSA_Name[i]

  csa_data <- sld_csa_all %>%
    filter(CSA == csa_code) %>%
    select(-any_of(drop_cols))

  if (nrow(csa_data) < 10) {
    cat("  Skipping", csa_name, "(n =", nrow(csa_data), ")\n")
    next
  }
  if (i %% 25 == 0)
    cat(sprintf("  [%d/%d] %s (n=%d)\n", i, nrow(all_csa_names), csa_name, nrow(csa_data)))

  pca_out <- PCA(csa_data, scale.unit = TRUE, ncp = 5, graph = FALSE)

  scores <- pca_out$ind$coord %>%
    as.data.frame() %>%
    bind_cols(sld_csa_all %>% filter(CSA == csa_code) %>% select(GEOID10)) %>%
    mutate(CSA = csa_code, CSA_Name = csa_name) %>%
    select(GEOID10, CSA, CSA_Name, everything()) %>%
    setNames(c("GEOID10", "CSA", "CSA_Name", paste0("Factor_", 1:5, "_CSA_raw")))

  all_csa_scores <- bind_rows(all_csa_scores, scores)
}
cat("CSA PCA loop complete.\n")

# ── Step 8: Normalize CSA scores ──────────────────────────────────────────────
cat("\n=== Step 8: Normalize CSA factor scores ===\n")
csa_factors_norm <- all_csa_scores %>%
  group_by(CSA_Name) %>%
  mutate(
    Factor_1_CSA_Norm = normalise(Factor_1_CSA_raw),
    Factor_2_CSA_Norm = normalise(Factor_2_CSA_raw),
    Factor_3_CSA_Norm = normalise(Factor_3_CSA_raw),
    Factor_4_CSA_Norm = normalise(Factor_4_CSA_raw),
    Factor_5_CSA_Norm = normalise(Factor_5_CSA_raw)
  ) %>%
  ungroup()

write_csv(csa_factors_norm,
          file.path(OUT_DIR, paste0("SLD_CSA_Factors_Normalized_V3_", Sys.Date(), ".csv")))
cat("CSA normalized scores saved. Rows:", nrow(csa_factors_norm), "\n")

# ── Step 9: CBSA loading consistency analysis ──────────────────────────────────
cat("\n=== Step 9: CBSA loading consistency analysis ===\n")

# V3 variable names for consistency analysis (same concepts as V2, updated names)
key_f1_vars <- c("TotEmp", "D1C", "D1C5_SVC", "D1C8_SVC", "E5_Off", "E8_off",
                 "D1C5_OFF", "E5_Svc", "E8_Svc")
key_f2_vars <- c("Pct_AO0", "AutoOwn0", "Pct_AO2p", "AutoOwn2p", "D1A", "D1B",
                 "D2B_E8MIX", "D2B_E5MIX")

cbsa_ranks <- all_cbsa_loadings %>%
  group_by(CBSA_Name) %>%
  mutate(rank_F1 = rank(-abs_F1), rank_F2 = rank(-abs_F2)) %>%
  ungroup()

f1_consistency <- cbsa_ranks %>%
  filter(variable %in% key_f1_vars) %>%
  group_by(variable) %>%
  summarise(
    median_rank_F1      = median(rank_F1),
    mean_abs_loading_F1 = mean(abs_F1),
    pct_in_top10_F1     = mean(rank_F1 <= 10) * 100,
    n_cbsa = n()
  ) %>%
  arrange(median_rank_F1)

f2_consistency <- cbsa_ranks %>%
  filter(variable %in% key_f2_vars) %>%
  group_by(variable) %>%
  summarise(
    median_rank_F2      = median(rank_F2),
    mean_abs_loading_F2 = mean(abs_F2),
    pct_in_top10_F2     = mean(rank_F2 <= 10) * 100,
    n_cbsa = n()
  ) %>%
  arrange(median_rank_F2)

cat("\nFactor 1 key variable consistency (V3 SLD):\n"); print(f1_consistency)
cat("\nFactor 2 key variable consistency (V3 SLD):\n"); print(f2_consistency)

# ── Step 10: Assemble supplementary file ──────────────────────────────────────
cat("\n=== Step 10: Assemble supplementary factor score file ===\n")

nat_slim  <- national_scores %>%
  select(GEOID10, Factor_1_National_Norm, Factor_2_National_Norm)
cbsa_slim <- cbsa_factors_norm %>%
  select(GEOID10, CBSA, CBSA_Name, Factor_1_CBSA_Norm, Factor_2_CBSA_Norm)
csa_slim  <- csa_factors_norm %>%
  select(GEOID10, CSA, CSA_Name, Factor_1_CSA_Norm, Factor_2_CSA_Norm)

supplementary <- nat_slim %>%
  left_join(cbsa_slim, by = "GEOID10") %>%
  left_join(csa_slim,  by = "GEOID10") %>%
  select(GEOID10,
    Factor_1_National_Norm, Factor_2_National_Norm,
    CBSA, CBSA_Name, Factor_1_CBSA_Norm, Factor_2_CBSA_Norm,
    CSA,  CSA_Name,  Factor_1_CSA_Norm,  Factor_2_CSA_Norm)

cat("Output rows:", nrow(supplementary), "\n")
cat("CBGs with CBSA score:", sum(!is.na(supplementary$Factor_1_CBSA_Norm)), "\n")
cat("CBGs with CSA score: ", sum(!is.na(supplementary$Factor_1_CSA_Norm)),  "\n")
cat("CBGs national only:  ",
    sum(is.na(supplementary$Factor_1_CBSA_Norm) & is.na(supplementary$Factor_1_CSA_Norm)), "\n")

# ── Step 11: Save all outputs ──────────────────────────────────────────────────
cat("\n=== Step 11: Saving outputs ===\n")

supp_path <- file.path(OUT_DIR,
  paste0("SLD_Factor_Scores_Supplementary_V3_", Sys.Date(), ".csv"))
write_csv(supplementary, supp_path)

write_csv(all_cbsa_loadings, file.path(OUT_DIR, "CBSA_PCA_Loadings_AllMarkets_V3.csv"))
write_csv(f1_consistency,    file.path(OUT_DIR, "CBSA_Factor1_Consistency_V3.csv"))
write_csv(f2_consistency,    file.path(OUT_DIR, "CBSA_Factor2_Consistency_V3.csv"))

all_factors <- left_join(
  cbsa_factors_norm %>% select(GEOID10, CBSA, CBSA_Name,
    Factor_1_CBSA_Norm, Factor_2_CBSA_Norm,
    Factor_3_CBSA_Norm, Factor_4_CBSA_Norm, Factor_5_CBSA_Norm),
  national_scores,
  by = "GEOID10"
)
write_csv(all_factors,
  file.path(OUT_DIR, paste0("SLD_Factors_Complete_V3_", Sys.Date(), ".csv")))

cat("\n=== Done! Files saved to:", OUT_DIR, "===\n")
cat("Key output:\n")
cat(" ", supp_path, "\n")
cat("Supporting outputs:\n")
cat("  National_Factor_Scores_AllCBG_V3.csv\n")
cat("  National_PCA_Correlations_V3.csv\n")
cat("  SLD_CBSA_Factors_Normalized_V3_[date].csv\n")
cat("  SLD_CSA_Factors_Normalized_V3_[date].csv\n")
cat("  CBSA_PCA_Loadings_AllMarkets_V3.csv\n")
cat("  CBSA_Factor1/2_Consistency_V3.csv\n")
cat("  SLD_Factors_Complete_V3_[date].csv\n")
