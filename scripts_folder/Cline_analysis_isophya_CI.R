############################# META-CLINE + HZAR #############################

## ====== CONFIG ======================================================================
setwd("/home/ismail/Research/isophya/cline_analysis/")

# A) Meta-cline inputs (minor-allele freq across altitudes)
freq_path <- "isophya.pops.assoc.loci.mafs.tsv"
n_path    <- "isophya.pops.assoc.loci.sample.number.tsv"
n_is_individuals <- FALSE     # TRUE if counts are individuals; converts to chromosomes

# B) HZAR inputs built from genotypes + sample altitudes
geno_path <- "isophya71.assoc.geno.tsv"   # wide: Chr Pos Maj Min Sample_1 …
info_path <- "isophya71.info"             # columns: sample_id, altitude[, site_id]
recode_to_global_minor <- TRUE            # keep focal allele = global minor
min_sites_per_locus     <- 3              # require ≥3 altitude sites for fitting

# Outputs
out_meta_pdf   <- "Figure_4_meta_cline.pdf"
out_hzar_dir   <- "hzar_inputs_all"
out_param_tsv  <- "hzar_ML_summary.tsv"             # HZAR params + LL/AIC (no CI)
out_join_tsv   <- "hzar_summary_with_trend.tsv"     # add trend (Spearman rho) + Δp
out_group_tsv  <- "cline_summary_by_trend.tsv"      # Increasing vs Decreasing comparison
out_violin_pdf <- "Figure_S4_center_width_violin_filtered_linear.pdf"

# HZAR model options (kept simple: free asymptotes, no tails)
fit_single_tail_too <- FALSE

# Parallel + batching (keeps RAM flat on typical workstations)
library(parallel)
total_cores <- parallel::detectCores(logical = TRUE)
workers     <- max(1L, floor(total_cores / 2))   # half cores → smoother memory use
chunk_size  <- 20L                                # tune 10–30

# Reproducibility
set.seed(1)

## ====== LIBS ========================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(fs)
  library(hzar)
  library(patchwork)
  library(scales)
})

## ====== PREVENT THREAD OVER-SUBSCRIPTION ===========================================
Sys.setenv("OMP_NUM_THREADS" = "1", "OPENBLAS_NUM_THREADS" = "1", "MKL_NUM_THREADS" = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)

## ====== HELPERS =====================================================================

# Convert letter genotypes to minor-allele counts (0/1/2; NA for missing).
convert_to_minor_counts <- function(gt_chr, maj, min) {
  gt_chr <- toupper(trimws(gt_chr))
  out <- rep(NA_integer_, length(gt_chr))
  mm <- paste0(maj, maj); mn <- paste0(maj, min)
  nm <- paste0(min, maj); nn <- paste0(min, min)
  out[gt_chr == mm] <- 0L
  out[gt_chr == mn] <- 1L
  out[gt_chr == nm] <- 1L
  out[gt_chr == nn] <- 2L
  out[gt_chr %in% c("NN","NA",".","")] <- NA_integer_
  out
}

# Plain TSV writer (stable for journal supplements).
write_tsv_plain <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

# Sigmoidal cline (consistent with HZAR parameterization used here).
sigmoid_cline <- function(x, center, width, pMin, pMax) {
  pMin + (pMax - pMin) * (1 / (1 + exp(-4 * (x - center) / width)))
}

## ------------------------------------------------------------------------------------
## PART 1. META-CLINE (two-slope regression; main Figure 4)
## ------------------------------------------------------------------------------------
freq <- read.table(freq_path, header=TRUE, sep="\t", check.names=FALSE)
ns   <- read.table(n_path,   header=TRUE, sep="\t", check.names=FALSE)
stopifnot(all(c("chromo","position") %in% names(freq)),
          all(c("chromo","position") %in% names(ns)))

alt_map <- function(df) {
  tibble(colname = names(df), altitude = readr::parse_number(names(df))) %>%
    filter(!is.na(altitude))
}
freq_cols <- alt_map(freq)
ns_cols   <- alt_map(ns)
common_alts <- intersect(freq_cols$altitude, ns_cols$altitude)
stopifnot(length(common_alts) > 0)
freq_cols <- freq_cols %>% filter(altitude %in% common_alts)
ns_cols   <- ns_cols   %>% filter(altitude %in% common_alts)

freq_long <- freq %>%
  pivot_longer(all_of(freq_cols$colname), names_to="alt_col", values_to="p") %>%
  left_join(freq_cols, by=c("alt_col"="colname")) %>%
  transmute(chromo, position, altitude = as.numeric(altitude), p = as.numeric(p))

n_long <- ns %>%
  pivot_longer(all_of(ns_cols$colname), names_to="alt_col", values_to="n") %>%
  left_join(ns_cols, by=c("alt_col"="colname")) %>%
  transmute(chromo, position, altitude = as.numeric(altitude), n = as.numeric(n))

if (n_is_individuals) n_long <- n_long %>% mutate(n = 2*n)

df_meta <- inner_join(freq_long, n_long, by=c("chromo","position","altitude")) %>%
  filter(is.finite(p), is.finite(n), n > 0) %>%
  mutate(p = pmin(1, pmax(0, p)))

# Classify loci by endpoint direction (robust; avoids overfitting).
trend_by_locus <- df_meta %>%
  group_by(chromo, position) %>%
  summarise(
    p_lo = p[which.min(altitude)],
    p_hi = p[which.max(altitude)],
    trend = if_else(p_hi > p_lo, "Increasing", "Decreasing"),
    .groups = "drop"
  )

df_meta <- df_meta %>% inner_join(trend_by_locus, by=c("chromo","position"))

# Collapse to locus means per altitude (avoid pseudo-replication).
locus_means <- df_meta %>%
  group_by(trend, chromo, position, altitude) %>%
  summarise(p_locus = weighted.mean(p, w = n), .groups = "drop")

# Slope-difference test (interaction).
fit_meta <- lm(p_locus ~ altitude * trend, data = locus_means)
print(summary(fit_meta))

# Figure 4 (two lines with SE ribbons).
summary_by_alt <- locus_means %>%
  group_by(trend, altitude) %>%
  summarise(mean_p = mean(p_locus), se = sd(p_locus)/sqrt(n()), .groups = "drop")

p_two_slope <- ggplot(summary_by_alt,
                      aes(x = altitude, y = mean_p, color = trend, fill = trend)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_ribbon(aes(ymin = mean_p - 1.96*se, ymax = mean_p + 1.96*se),
              alpha = 0.20, color = NA) +
  scale_color_manual(values = c(Increasing = "#00441B", Decreasing = "#66C2A5")) +
  scale_fill_manual(values  = c(Increasing = "#00441B", Decreasing = "#66C2A5")) +
  labs(x = "Altitude (m)", y = "Mean frequency of globally defined minor allele",
       color = "Locus group", fill = "Locus group") +
  theme_minimal(base_size = 11)
ggsave(out_meta_pdf, p_two_slope, width = 140, height = 90, units = "mm")
message("Wrote: ", out_meta_pdf)

## ------------------------------------------------------------------------------------
## PART 2. BUILD HZAR INPUTS (from genotypes) FOR ALL LOCI
## ------------------------------------------------------------------------------------
geno_raw <- read_tsv(
  geno_path,
  col_types = cols(
    Chr = col_character(),
    Pos = col_character(),
    Maj = col_character(),
    Min = col_character(),
    .default = col_character()
  ),
  trim_ws = TRUE,
  na = c("NA","Na","na",".","")
)
info <- read_tsv(info_path, col_types = cols(.default = col_character()))
stopifnot(all(c("Chr","Pos","Maj","Min") %in% names(geno_raw)))
stopifnot(all(c("sample_id","altitude") %in% names(info)))
info <- info %>% mutate(altitude = as.numeric(altitude))
if (!"site_id" %in% names(info)) info <- info %>% mutate(site_id = paste0("ALT_", altitude))

sample_cols <- setdiff(names(geno_raw), c("Chr","Pos","Maj","Min"))
stopifnot(length(sample_cols) > 0)

geno_long <- geno_raw %>%
  pivot_longer(all_of(sample_cols), names_to = "sample_id", values_to = "GT") %>%
  mutate(
    Chr = as.character(Chr),
    Pos = as.character(Pos),
    locus_id = paste(Chr, Pos, sep = "_")
  ) %>%
  inner_join(info[, c("sample_id","altitude","site_id")], by = "sample_id") %>%
  mutate(Maj = toupper(Maj), Min = toupper(Min))

# Sanity: one Maj/Min per locus.
chk <- geno_long %>%
  group_by(locus_id) %>%
  summarise(nMaj = n_distinct(Maj), nMin = n_distinct(Min), .groups = "drop") %>%
  filter(nMaj != 1 | nMin != 1)
if (nrow(chk) > 0) {
  warning("Loci with inconsistent Maj/Min: ",
          paste(head(chk$locus_id, 10), collapse = ", "),
          if (nrow(chk) > 10) " …")
}

dir_create(out_hzar_dir)
n_written <- 0L
for (loc in unique(geno_long$locus_id)) {
  sub <- geno_long %>% filter(locus_id == loc)
  
  maj_allele <- sub %>% distinct(Maj) %>% pull()
  min_allele <- sub %>% distinct(Min) %>% pull()
  if (length(maj_allele) != 1 || length(min_allele) != 1) next
  maj_allele <- maj_allele[1]; min_allele <- min_allele[1]
  
  g_counts <- convert_to_minor_counts(sub$GT, maj = maj_allele, min = min_allele)
  
  # Optionally recode to global minor (stable focal allele across altitudes).
  if (recode_to_global_minor) {
    ALT_k_global <- sum(g_counts, na.rm = TRUE)
    n_chr_global <- 2L * sum(!is.na(g_counts))
    ALT_p_global <- ifelse(n_chr_global > 0, ALT_k_global / n_chr_global, NA_real_)
    if (!is.na(ALT_p_global) && ALT_p_global > 0.5) {
      g_counts <- ifelse(is.na(g_counts), NA_integer_, 2L - g_counts)
    }
  }
  
  sub$k_minor <- g_counts
  sub$n_chrom <- ifelse(is.na(g_counts), 0L, 2L)
  
  per_site <- sub %>%
    filter(!is.na(altitude)) %>%
    group_by(site_id, altitude) %>%
    summarise(k_minor = sum(k_minor, na.rm = TRUE),
              n_chrom = sum(n_chrom, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(altitude)
  
  if (nrow(per_site) < min_sites_per_locus) next
  
  out_file <- fs::path(out_hzar_dir, paste0(gsub("[^A-Za-z0-9_.:-]+","_", loc), ".tsv"))
  write_tsv_plain(per_site, out_file)
  n_written <- n_written + 1L
}
message("Built ", n_written, " per-locus TSVs in ", out_hzar_dir)

## ------------------------------------------------------------------------------------
## PART 3. FIT HZAR CLINES (ML, TRIMMED MCMC) + INLINE LL/AIC (NO CI)
## ------------------------------------------------------------------------------------
read_one_locus <- function(path_tsv) {
  df <- tryCatch(read.table(path_tsv, header = TRUE, sep = "\t", check.names = FALSE),
                 error = function(e) NULL)
  if (is.null(df)) return(NULL)
  req <- c("site_id","altitude","k_minor","n_chrom")
  if (!all(req %in% names(df))) return(NULL)
  df <- df %>%
    mutate(altitude = as.numeric(altitude),
           k_minor  = as.numeric(k_minor),
           n_chrom  = as.numeric(n_chrom)) %>%
    filter(is.finite(altitude), is.finite(k_minor), is.finite(n_chrom), n_chrom > 0) %>%
    arrange(altitude)
  if (nrow(df) < 3) return(NULL)
  eps <- 1e-6
  df <- df %>% mutate(p = pmin(1 - eps, pmax(eps, k_minor / n_chrom)))
  df
}

fit_hzar_basic <- function(x, p, n,
                           chain_length = 2e5L,
                           burnin       = 2e4L,
                           thin         = 100L) {
  obs   <- hzar.doMolecularData1DPops(distance = x, pObs = p, n)  # 3rd arg positional
  model <- hzar.makeCline1DFreq(obs, scaling = "free", tails = "none")
  req   <- hzar.first.fitRequest.old.ML(model, obs)
  
  # Trim MCMC (biggest runtime win; keeps ML stable enough for ranks/summary).
  if (!is.null(req$mcmcParam)) {
    req$mcmcParam$chainLength <- as.integer(chain_length)
    req$mcmcParam$burnin      <- as.integer(burnin)
    req$mcmcParam$thin        <- as.integer(thin)
  }
  
  fit <- try(hzar.doFit(req), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  ml  <- try(hzar.get.ML.cline(fit), silent = TRUE)
  if (inherits(ml, "try-error")) return(NULL)
  list(ml = ml)
}

extract_params <- function(ml) {
  p_all <- try(ml$param.all, silent = TRUE)
  if (inherits(p_all, "try-error") || is.null(p_all)) return(NULL)
  list(center = suppressWarnings(as.numeric(p_all$center)),
       width  = suppressWarnings(as.numeric(p_all$width)),
       pMin   = suppressWarnings(as.numeric(p_all$pMin)),
       pMax   = suppressWarnings(as.numeric(p_all$pMax)))
}

fit_one_file <- function(path_tsv) {
  df <- read_one_locus(path_tsv); if (is.null(df)) return(NULL)
  x <- df$altitude; p <- df$p; n <- df$n_chrom; k <- df$k_minor
  
  res_basic <- fit_hzar_basic(x, p, n)
  if (is.null(res_basic)) return(NULL)
  params <- extract_params(res_basic$ml)
  if (is.null(params) || any(!is.finite(unlist(params)))) return(NULL)
  
  # LL/AIC at observed sites (computed here to avoid post-processing).
  p_hat <- sigmoid_cline(x, params$center, params$width, params$pMin, params$pMax)
  eps <- 1e-12; p_hat <- pmin(1 - eps, pmax(eps, p_hat))
  ll  <- sum(k * log(p_hat) + (n - k) * log(1 - p_hat))
  aic <- 2*4 - 2*ll
  
  tibble(locus_id = path_file(path_tsv) %>% path_ext_remove(),
         model    = "free_noTails",
         center   = params$center,
         width    = params$width,
         pMin     = params$pMin,
         pMax     = params$pMax,
         logLik   = ll,
         AIC      = aic)
}

# Discover files
in_dir <- out_hzar_dir
tsv_files <- dir_ls(in_dir, glob = "*.tsv")
stopifnot(length(tsv_files) > 0)

# Batched parallel execution (steady memory profile).
chunks <- split(tsv_files, ceiling(seq_along(tsv_files) / chunk_size))
res_list <- vector("list", length(tsv_files))
idx <- 1L
for (ci in seq_along(chunks)) {
  ch <- chunks[[ci]]
  cat("Batch", ci, "of", length(chunks), " | n =", length(ch), " | workers =", workers, "\n")
  r <- mclapply(ch, fit_one_file, mc.cores = workers, mc.preschedule = TRUE)
  for (i in seq_along(r)) { res_list[[idx]] <- r[[i]]; idx <- idx + 1L }
  rm(r); gc()
}
res_list <- res_list[!vapply(res_list, is.null, logical(1))]
res <- bind_rows(res_list)
write_tsv_plain(res %>% arrange(locus_id), out_param_tsv)
message("Wrote: ", out_param_tsv)

## ------------------------------------------------------------------------------------
## PART 4. TREND (Spearman), QC FILTER, FIGURE (LINEAR), GROUP COMPARISON
## ------------------------------------------------------------------------------------

# Spearman trend from per-locus TSV (robust to outliers & nonlinearity).
trend_from_file <- function(path_tsv, rho_threshold = 0.10) {
  df <- tryCatch(read.table(path_tsv, header = TRUE, sep = "\t", check.names = FALSE),
                 error = function(e) NULL)
  if (is.null(df)) return(NULL)
  df <- df %>%
    mutate(altitude = as.numeric(altitude),
           k_minor  = as.numeric(k_minor),
           n_chrom  = as.numeric(n_chrom),
           p        = k_minor / n_chrom) %>%
    filter(is.finite(altitude), is.finite(p), is.finite(n_chrom), n_chrom > 0)
  if (nrow(df) < 3) return(tibble(trend = "Flat", rho = NA_real_))
  rho <- suppressWarnings(cor(df$altitude, df$p, method = "spearman"))
  trend <- dplyr::case_when(
    is.finite(rho) & rho >  0.10 ~ "Increasing",
    is.finite(rho) & rho < -0.10 ~ "Decreasing",
    TRUE                         ~ "Flat"
  )
  tibble(trend = trend, rho = rho)
}

# Compute trend for all fitted loci.
paths <- fs::path(in_dir, paste0(gsub("[^A-Za-z0-9_.:-]+","_", res$locus_id), ".tsv"))
names(paths) <- res$locus_id
trend_tbl <- bind_rows(lapply(names(paths), function(lid) {
  out <- trend_from_file(paths[[lid]]); if (is.null(out)) return(NULL)
  mutate(out, locus_id = lid)
}))

# Join trend + Δp and save.
res_join <- res %>%
  mutate(delta_p = pMax - pMin) %>%
  left_join(trend_tbl, by = "locus_id") %>%
  arrange(locus_id)
write_tsv_plain(res_join, out_join_tsv)
message("Wrote: ", out_join_tsv)

# ---- QC thresholds (tune here if needed) ----
QC_RANGE_BUFFER_PROP <- 0.10     # allow center within ±10% of sampled range
QC_WIDTH_MAX_FACTOR  <- 5        # max width allowed as multiple of sampled range
QC_WIDTH_MIN_ABS     <- 10       # widths < 10 m are too sharp for site spacing
QC_DP_MIN            <- 0.02     # minimal fitted amplitude
QC_DP_OBS_MIN        <- 0.02     # minimal observed amplitude

# Build per-locus altitude bounds and observed Δp directly from TSVs.
get_locus_meta <- function(locus_id) {
  f <- fs::path(in_dir, paste0(gsub("[^A-Za-z0-9_.:-]+","_", locus_id), ".tsv"))
  if (!file.exists(f)) return(NULL)
  df <- tryCatch(read.table(f, header=TRUE, sep="\t", check.names=FALSE), error=function(e) NULL)
  if (is.null(df)) return(NULL)
  df <- df %>% mutate(
    altitude = as.numeric(altitude),
    k_minor  = as.numeric(k_minor),
    n_chrom  = as.numeric(n_chrom),
    p        = k_minor / n_chrom
  ) %>% filter(is.finite(altitude), is.finite(p), is.finite(n_chrom), n_chrom>0)
  if (nrow(df) < 3) return(NULL)
  tibble(
    locus_id = locus_id,
    alt_min  = min(df$altitude, na.rm=TRUE),
    alt_max  = max(df$altitude, na.rm=TRUE),
    alt_rng  = max(df$altitude, na.rm=TRUE) - min(df$altitude, na.rm=TRUE),
    dp_obs   = diff(range(df$p, na.rm=TRUE))
  )
}

meta_list <- lapply(res_join$locus_id, get_locus_meta)
meta_list <- meta_list[!vapply(meta_list, is.null, logical(1))]
meta_tbl  <- bind_rows(meta_list)

# QC flags and filtered table.
res_qc <- res_join %>%
  left_join(meta_tbl, by="locus_id") %>%
  mutate(
    alt_rng = ifelse(is.na(alt_rng), NA_real_, alt_rng),
    flag_center_extrap = ifelse(is.finite(center) & is.finite(alt_min) & is.finite(alt_max) & is.finite(alt_rng),
                                center < (alt_min - QC_RANGE_BUFFER_PROP*alt_rng) |
                                  center > (alt_max + QC_RANGE_BUFFER_PROP*alt_rng),
                                TRUE),
    flag_width_inflate = ifelse(is.finite(width) & is.finite(alt_rng),
                                width > QC_WIDTH_MAX_FACTOR*alt_rng | width < QC_WIDTH_MIN_ABS,
                                TRUE),
    flag_delta_p_small = !is.finite(delta_p) | delta_p < QC_DP_MIN,
    flag_dp_obs_small  = !is.finite(dp_obs)  | dp_obs  < QC_DP_OBS_MIN,
    flagged = flag_center_extrap | flag_width_inflate | flag_delta_p_small | flag_dp_obs_small
  )

write_tsv_plain(res_qc %>% arrange(locus_id), "hzar_ML_summary_QC.tsv")

res_qc_filt <- res_qc %>%
  filter(!flagged, is.finite(center), is.finite(width), width > 0)

write_tsv_plain(res_qc_filt %>% arrange(locus_id), "hzar_ML_summary_QC_filtered.tsv")
message("QC kept ", nrow(res_qc_filt), " / ", nrow(res_qc), " loci (",
        round(100*nrow(res_qc_filt)/max(1,nrow(res_qc)),1), "%).")

# Figure S4 (linear axis, more ticks, filtered).
res_plot_filt <- res_qc_filt %>%
  filter(trend %in% c("Increasing","Decreasing","Flat"))

p_center_box_lin <- ggplot(res_plot_filt, aes(trend, center, fill = trend)) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  scale_y_continuous(
    name = "Cline center (m)",
    breaks = pretty_breaks(n = 8),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_fill_manual(values = c(Increasing = "#1b9e77", Decreasing = "#d95f02", Flat = "grey60")) +
  theme_minimal(base_size = 12) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 11),
        axis.text.y  = element_text(size = 10),
        legend.position = "none")

p_width_box_lin <- ggplot(res_plot_filt, aes(trend, width, fill = trend)) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  scale_y_continuous(
    name = "Cline width (m)",
    breaks = pretty_breaks(n = 8),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_fill_manual(values = c(Increasing = "#1b9e77", Decreasing = "#d95f02", Flat = "grey60")) +
  theme_minimal(base_size = 12) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 11),
        axis.text.y  = element_text(size = 10),
        legend.position = "none")

ggsave(out_violin_pdf, (p_center_box_lin | p_width_box_lin), width = 180, height = 95, units = "mm")
message("Wrote: ", out_violin_pdf)

# Comparative summary by trend (filtered set only).
group_cmp <- res_qc_filt %>%
  filter(trend %in% c("Increasing","Decreasing")) %>%
  group_by(trend) %>%
  summarise(
    n_loci      = n(),
    center_med  = median(center, na.rm = TRUE),
    center_IQR  = IQR(center, na.rm = TRUE),
    width_med   = median(width, na.rm = TRUE),
    width_IQR   = IQR(width, na.rm = TRUE),
    dP_med      = median(delta_p, na.rm = TRUE),
    dP_IQR      = IQR(delta_p, na.rm = TRUE),
    AIC_med     = median(AIC, na.rm = TRUE),
    .groups     = "drop"
  )
write_tsv_plain(group_cmp, out_group_tsv)
message("Wrote: ", out_group_tsv)

# Done.
####################################################################################################


# Load filtered HZAR summary with trend + QC
res_qc_filt <- read.table("hzar_ML_summary_QC_filtered.tsv", header=TRUE, sep="\t", check.names=FALSE)

# If trend is not present (should be, from your script), merge it in:
if (!"trend" %in% names(res_qc_filt)) {
  trend_tbl <- read.table("hzar_summary_with_trend.tsv", header=TRUE, sep="\t", check.names=FALSE)[, c("locus_id","trend")]
  res_qc_filt <- merge(res_qc_filt, trend_tbl, by="locus_id", all.x=TRUE)
}

# Counts
N_keep <- nrow(res_qc_filt)
tab_trend <- table(res_qc_filt$trend)
N_inc <- as.integer(tab_trend[["Increasing"]])
N_dec <- as.integer(tab_trend[["Decreasing"]])

# Summaries
center_med <- median(res_qc_filt$center, na.rm=TRUE)
center_q   <- quantile(res_qc_filt$center, probs=c(0.25, 0.75), na.rm=TRUE)
width_med  <- median(res_qc_filt$width,  na.rm=TRUE)
width_q    <- quantile(res_qc_filt$width,  probs=c(0.25, 0.75), na.rm=TRUE)
dP_med     <- median(res_qc_filt$pMax - res_qc_filt$pMin, na.rm=TRUE)
dP_q       <- quantile(res_qc_filt$pMax - res_qc_filt$pMin, probs=c(0.25, 0.75), na.rm=TRUE)

list(
  N_keep = N_keep,
  N_inc  = N_inc,
  N_dec  = N_dec,
  center_med = center_med,
  center_q25 = center_q[[1]],
  center_q75 = center_q[[2]],
  width_med  = width_med,
  width_q25  = width_q[[1]],
  width_q75  = width_q[[2]],
  dP_med     = dP_med,
  dP_q25     = dP_q[[1]],
  dP_q75     = dP_q[[2]]
)
