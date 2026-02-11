############################################
# 4_NMDS_Turnover_and_CommunitySpace_with_Uncertainty_Outhwaite
# Standalone script (run top-to-bottom)
#
# INPUT:
#   - Combined BUGS summary CSV containing psi.fs.r_{Zone}[YearIndex]
#     with columns at least: Parameter, Species, mean, sd
#
# UNCERTAINTY PROPAGATION (Cooke/Outhwaite-style approximation):
#   - For each Species × Zone × Year cell, reconstruct a Beta distribution
#     from BUGS posterior mean + sd
#   - Monte Carlo draw one occupancy per cell per draw
#   - For each draw:
#       * build Zone-Year × Species community matrix
#       * NMDS (Bray–Curtis, k=2)
#       * Procrustes-align to a reference NMDS (computed on posterior means)
#       * compute interannual turnover per zone = Euclidean step distance in aligned NMDS space
#       * save aligned NMDS site coordinates per Zone-Year for "community space" visualisation
#   - Summaries:
#       * Turnover (Zone×Year): median + 95% CrI
#       * NMDS site coords (Zone×Year): median + 95% CrI for NMDS1 and NMDS2
#
# OUTPUTS (all prefixed "4_"):
#   - Reference ordination (means):
#       4_NMDS_Ordination_REFERENCE_MEANS_2000_2023.png/.jpg
#   - Uncertainty-aware ordination:
#       4_NMDS_Ordination_WITH_UNCERTAINTY_2000_2023.png/.jpg
#       4_NMDS_Ordination_DrawCloud_2000_2023.png/.jpg
#   - Turnover:
#       4_NMDS_Interannual_Turnover_WITH_UNCERTAINTY_2000_2023.png/.jpg
#   - Tables:
#       4_NMDS_Turnover_Draws_Long_2000_2023.csv
#       4_NMDS_Turnover_Summary_ZoneYear_2000_2023.csv
#       4_NMDS_SiteScores_Draws_Long_2000_2023.csv
#       4_NMDS_SiteScores_Summary_2000_2023.csv
#       4_NMDS_Summary_Tables_2000_2023.xlsx
#
# NOTES:
#   - Running NMDS hundreds/thousands of times is computationally expensive.
#     Start with n_draws_nmds = 200–400 and increase if needed.
#   - Procrustes alignment is essential: NMDS solutions are only defined up to
#     rotation/reflection/translation. We align each draw to the reference.
#   - PERMANOVA is run on the reference (posterior mean) matrix only.
############################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(vegan)
  library(ggplot2)
  library(readr)
  library(writexl)
  library(scales)  # alpha()
})

# ----------------------------
# 0) USER SETTINGS
# ----------------------------
setwd("C:/Users/georg/OneDrive - University of Reading/George Allen - PhD Master Folder/Year Three/Chapter 2 - Occupancy Modelling/Analysis Nov 2025/Nov_Outhwaite_Outputs")

bugs_file <- "Combined_BUGS_Data_Outhwaite_finalclustv5_32000iterations_FINAL_WITH_TAXCORR.csv"

start_year <- 2000
end_year   <- 2023
year_index_offset <- 1999

zones_keep <- c("Core","Buffer","Outside")
zone_colors <- c("Core"="#b2182b","Buffer"="black","Outside"="grey70")
zone_fill_colors <- c("Core"=alpha("#b2182b",0.20),
                      "Buffer"=alpha("black",0.20),
                      "Outside"=alpha("grey40",0.30))

# Regex for parameters
psi_zone_regex <- "^psi\\.fs\\.r_(Core|Buffer|Outside)\\[(\\d+)\\]$"

# Monte Carlo draws for uncertainty-aware NMDS turnover + space
set.seed(123)
n_draws_nmds <- 300        # increase to 500–1000 if runtime allows
eps <- 1e-6

# Optional: save aligned site coordinates only every k draws (1 = keep all)
save_every_k <- 1

# NMDS settings
k_dim <- 2
trymax_nmds <- 40          # higher = more stable but slower
trace_nmds <- 0

# Plot cloud thinning (fraction of draw-points to keep)
cloud_frac <- 0.25

# Outputs (prefix "4_")
out_ref_ord_png <- paste0("4_NMDS_Ordination_REFERENCE_MEANS_", start_year, "_", end_year, ".png")
out_ref_ord_jpg <- paste0("4_NMDS_Ordination_REFERENCE_MEANS_", start_year, "_", end_year, ".jpg")

out_unc_ord_png <- paste0("4_NMDS_Ordination_WITH_UNCERTAINTY_", start_year, "_", end_year, ".png")
out_unc_ord_jpg <- paste0("4_NMDS_Ordination_WITH_UNCERTAINTY_", start_year, "_", end_year, ".jpg")

out_cloud_png   <- paste0("4_NMDS_Ordination_DrawCloud_", start_year, "_", end_year, ".png")
out_cloud_jpg   <- paste0("4_NMDS_Ordination_DrawCloud_", start_year, "_", end_year, ".jpg")

out_turn_png <- paste0("4_NMDS_Interannual_Turnover_WITH_UNCERTAINTY_", start_year, "_", end_year, ".png")
out_turn_jpg <- paste0("4_NMDS_Interannual_Turnover_WITH_UNCERTAINTY_", start_year, "_", end_year, ".jpg")

out_turn_draws_csv <- paste0("4_NMDS_Turnover_Draws_Long_", start_year, "_", end_year, ".csv")
out_turn_sum_csv   <- paste0("4_NMDS_Turnover_Summary_ZoneYear_", start_year, "_", end_year, ".csv")

out_sites_draws_csv <- paste0("4_NMDS_SiteScores_Draws_Long_", start_year, "_", end_year, ".csv")
out_sites_sum_csv   <- paste0("4_NMDS_SiteScores_Summary_", start_year, "_", end_year, ".csv")

out_xlsx <- paste0("4_NMDS_Summary_Tables_", start_year, "_", end_year, ".xlsx")

# ----------------------------
# 1) HELPERS (same philosophy as your other scripts)
# ----------------------------
beta_ab_from_mean_sd <- function(mu, sd) {
  # returns list(a,b) or NA if invalid
  if (!is.finite(mu) || !is.finite(sd)) return(list(a=NA_real_, b=NA_real_))
  if (mu <= 0 || mu >= 1) return(list(a=NA_real_, b=NA_real_))
  if (sd <= 0) return(list(a=NA_real_, b=NA_real_))
  v <- sd^2
  vmax <- mu * (1 - mu)
  if (!is.finite(v) || !is.finite(vmax) || vmax <= 0) return(list(a=NA_real_, b=NA_real_))
  if (v >= vmax) v <- 0.95 * vmax
  k <- mu * (1 - mu) / v - 1
  if (!is.finite(k) || k <= 0) return(list(a=NA_real_, b=NA_real_))
  a <- mu * k
  b <- (1 - mu) * k
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(list(a=NA_real_, b=NA_real_))
  list(a=a, b=b)
}

draw_row <- function(mu, sd, a, b, n = 1) {
  # Prefer Beta if feasible; fallback to truncated normal if Beta invalid
  if (is.finite(a) && is.finite(b) && a > 0 && b > 0) {
    d <- rbeta(n, a, b)
  } else {
    sduse <- ifelse(is.finite(sd) && sd > 0, sd, 0.05)
    d <- rnorm(n, mean = mu, sd = sduse)
  }
  pmin(1 - eps, pmax(eps, d))
}

euclid_step <- function(x1, y1, x0, y0) sqrt((x1 - x0)^2 + (y1 - y0)^2)

# ----------------------------
# 2) LOAD BUGS + EXTRACT psi.fs.r_{Zone}[YearIndex]
# ----------------------------
bugs <- read.csv(bugs_file, stringsAsFactors = FALSE)
stopifnot(all(c("Parameter","Species","mean","sd") %in% names(bugs)))

psi_long <- bugs %>%
  filter(str_detect(Parameter, "^psi\\.fs\\.r_")) %>%
  mutate(
    Zone = str_match(Parameter, psi_zone_regex)[,2],
    YearIndex = suppressWarnings(as.integer(str_match(Parameter, psi_zone_regex)[,3]))
  ) %>%
  filter(!is.na(Zone), !is.na(YearIndex)) %>%
  mutate(
    Year = YearIndex + year_index_offset,
    Zone = factor(Zone, levels = zones_keep),
    Species = as.character(Species),
    psi_mean = suppressWarnings(as.numeric(mean)),
    psi_sd   = suppressWarnings(as.numeric(sd))
  ) %>%
  filter(Year >= start_year, Year <= end_year) %>%
  filter(!is.na(Species), Species != "", !is.na(Zone),
         is.finite(psi_mean), is.finite(psi_sd)) %>%
  mutate(
    psi_mean = pmin(pmax(psi_mean, eps), 1 - eps),
    psi_sd   = pmax(psi_sd, 1e-12)
  ) %>%
  select(Species, Zone, Year, psi_mean, psi_sd)

cat("\n========================================================\n")
cat("Extracted psi.fs.r_{Zone}[YearIndex]\n")
cat("========================================================\n")
cat("Rows:", nrow(psi_long), "\n")
cat("Species:", n_distinct(psi_long$Species), "\n")
cat("Years:", paste(range(psi_long$Year), collapse="–"), "\n")
cat("Zones:\n"); print(table(psi_long$Zone, useNA="ifany"))

# Add Beta parameters
psi_par <- psi_long %>%
  rowwise() %>%
  mutate(
    ab = list(beta_ab_from_mean_sd(psi_mean, psi_sd)),
    a = ab$a,
    b = ab$b
  ) %>%
  ungroup() %>%
  select(-ab)

# ----------------------------
# 3) BUILD REFERENCE COMMUNITY MATRIX (posterior MEANS)
#    - rows: Zone-Year, cols: Species
# ----------------------------
comm_ref_df <- psi_par %>%
  group_by(Zone, Year, Species) %>%
  summarise(val = mean(psi_mean, na.rm = TRUE), .groups = "drop") %>%
  unite("Zone_Year", Zone, Year, sep = "_") %>%
  pivot_wider(names_from = Species, values_from = val, values_fill = 0)

metadata_ref <- comm_ref_df %>%
  separate(Zone_Year, into = c("Zone","Year"), sep = "_") %>%
  mutate(
    Zone = factor(Zone, levels = zones_keep),
    Year = as.numeric(Year)
  )

comm_ref <- comm_ref_df %>% select(-Zone_Year) %>% as.matrix()
comm_ref[is.na(comm_ref)] <- 0

# Remove all-zero rows/cols (and lock them for all draws)
keep_rows <- rowSums(comm_ref) > 0
comm_ref <- comm_ref[keep_rows, , drop = FALSE]
metadata_ref <- metadata_ref[keep_rows, , drop = FALSE]

keep_cols <- colSums(comm_ref) > 0
comm_ref <- comm_ref[, keep_cols, drop = FALSE]

stopifnot(nrow(comm_ref) == nrow(metadata_ref))

# Locked ordering for all draws
species_levels <- colnames(comm_ref)
sample_levels  <- paste(metadata_ref$Zone, metadata_ref$Year, sep = "_")

# Prepare indexing for fast fill-in per draw
psi_indexed <- psi_par %>%
  mutate(
    Zone_Year = paste(Zone, Year, sep = "_"),
    sample_id = match(Zone_Year, sample_levels),
    species_id = match(Species, species_levels)
  ) %>%
  filter(!is.na(sample_id), !is.na(species_id)) %>%
  select(sample_id, species_id, psi_mean, psi_sd, a, b)

cat("\nLocked community matrix size:\n")
cat("Samples (Zone-Year):", length(sample_levels), "\n")
cat("Species columns:", length(species_levels), "\n")
cat("Cells used:", nrow(psi_indexed), "\n")

# ----------------------------
# 4) REFERENCE NMDS (on posterior means) + plot
# ----------------------------
set.seed(123)
nmds_ref <- metaMDS(comm_ref, distance = "bray", k = k_dim, trymax = trymax_nmds, trace = trace_nmds)

ref_sites <- as.data.frame(scores(nmds_ref, display = "sites"))
ref_sites <- bind_cols(ref_sites, metadata_ref)

cat("\nReference NMDS stress:", round(nmds_ref$stress, 4), "\n")

g_ref <- ggplot(ref_sites, aes(x = NMDS1, y = NMDS2, color = Zone)) +
  geom_point(size = 3, alpha = 0.75) +
  stat_ellipse(type = "t", linetype = 2, linewidth = 1) +
  scale_color_manual(values = zone_colors) +
  labs(
    title = "Reference NMDS (posterior mean occupancy)",
    subtitle = paste0("Bray–Curtis, Zone-Year samples; stress = ", round(nmds_ref$stress, 3)),
    x = "NMDS1", y = "NMDS2"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(color = "gray30"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.position = "none"
  )

print(g_ref)
ggsave(out_ref_ord_png, g_ref, width = 7.5, height = 6, dpi = 900)
ggsave(out_ref_ord_jpg, g_ref, width = 7.5, height = 6, dpi = 900)

# Reference XY for Procrustes target
ref_xy <- as.matrix(ref_sites[, c("NMDS1","NMDS2")])
rownames(ref_xy) <- sample_levels

# ----------------------------
# 5) UNCERTAINTY PROPAGATION:
#    Repeat NMDS across draws, Procrustes-align, compute turnover + save site coords
# ----------------------------
turnover_draws <- vector("list", n_draws_nmds)
site_draws     <- vector("list", n_draws_nmds)

set.seed(999)
for (d in seq_len(n_draws_nmds)) {
  
  # 5a) Draw one occupancy for each indexed cell
  vals <- numeric(nrow(psi_indexed))
  for (i in seq_len(nrow(psi_indexed))) {
    vals[i] <- draw_row(
      mu = psi_indexed$psi_mean[i],
      sd = psi_indexed$psi_sd[i],
      a  = psi_indexed$a[i],
      b  = psi_indexed$b[i],
      n  = 1
    )
  }
  
  # 5b) Fill community matrix for this draw
  comm_d <- matrix(0, nrow = nrow(comm_ref), ncol = ncol(comm_ref))
  comm_d[cbind(psi_indexed$sample_id, psi_indexed$species_id)] <- vals
  
  # 5c) NMDS for this draw
  nmds_d <- try(
    metaMDS(comm_d, distance = "bray", k = k_dim, trymax = trymax_nmds, trace = trace_nmds),
    silent = TRUE
  )
  
  if (inherits(nmds_d, "try-error")) {
    turnover_draws[[d]] <- tibble(Zone=character(0), Year=numeric(0), draw_id=integer(0), delta=numeric(0))
    site_draws[[d]]     <- tibble(draw_id=integer(0), Zone=character(0), Year=numeric(0), NMDS1=numeric(0), NMDS2=numeric(0))
    next
  }
  
  cur_sites <- as.data.frame(scores(nmds_d, display = "sites"))
  cur_xy <- as.matrix(cur_sites[, 1:2, drop = FALSE])
  
  # 5d) Procrustes-align current sites to reference sites
  proc <- procrustes(ref_xy, cur_xy, symmetric = TRUE)
  cur_rot <- proc$Yrot
  
  # 5e) Save aligned site coordinates
  if (d %% save_every_k == 0) {
    site_draws[[d]] <- metadata_ref %>%
      mutate(
        draw_id = d,
        NMDS1 = cur_rot[,1],
        NMDS2 = cur_rot[,2]
      ) %>%
      select(draw_id, Zone, Year, NMDS1, NMDS2)
  } else {
    site_draws[[d]] <- tibble(draw_id=integer(0), Zone=character(0), Year=numeric(0), NMDS1=numeric(0), NMDS2=numeric(0))
  }
  
  # 5f) Compute year-to-year step distances per zone (aligned coordinates)
  tmp_turn <- metadata_ref %>%
    mutate(
      NMDS1 = cur_rot[,1],
      NMDS2 = cur_rot[,2]
    ) %>%
    arrange(Zone, Year) %>%
    group_by(Zone) %>%
    mutate(delta = euclid_step(NMDS1, NMDS2, lag(NMDS1), lag(NMDS2))) %>%
    filter(!is.na(delta)) %>%
    ungroup() %>%
    mutate(draw_id = d) %>%
    select(Zone, Year, draw_id, delta)
  
  turnover_draws[[d]] <- tmp_turn
  
  if (d %% 25 == 0) cat("Completed draw", d, "of", n_draws_nmds, "\n")
}

turnover_draws_long <- bind_rows(turnover_draws)
site_draws_long     <- bind_rows(site_draws)

# Save long draws
write.csv(turnover_draws_long, out_turn_draws_csv, row.names = FALSE)
cat("\nSaved turnover draws:", out_turn_draws_csv, "\n")

write.csv(site_draws_long, out_sites_draws_csv, row.names = FALSE)
cat("Saved aligned site draws:", out_sites_draws_csv, "\n")

# ----------------------------
# 6) SUMMARISE TURNOVER (median + 95% CrI) by Zone-Year
# ----------------------------
turnover_summary <- turnover_draws_long %>%
  group_by(Zone, Year) %>%
  summarise(
    delta_median = median(delta, na.rm = TRUE),
    delta_low    = as.numeric(quantile(delta, 0.025, na.rm = TRUE)),
    delta_high   = as.numeric(quantile(delta, 0.975, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(Zone, Year)

write.csv(turnover_summary, out_turn_sum_csv, row.names = FALSE)
cat("Saved turnover summary:", out_turn_sum_csv, "\n")

turnover_zone_summary <- turnover_draws_long %>%
  group_by(Zone, draw_id) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    sd_delta   = sd(delta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Zone) %>%
  summarise(
    mean_delta_median = median(mean_delta, na.rm = TRUE),
    mean_delta_low    = as.numeric(quantile(mean_delta, 0.025, na.rm = TRUE)),
    mean_delta_high   = as.numeric(quantile(mean_delta, 0.975, na.rm = TRUE)),
    sd_delta_median   = median(sd_delta, na.rm = TRUE),
    sd_delta_low      = as.numeric(quantile(sd_delta, 0.025, na.rm = TRUE)),
    sd_delta_high     = as.numeric(quantile(sd_delta, 0.975, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(Zone)

# ----------------------------
# 7) SUMMARISE COMMUNITY SPACE (median + 95% CrI) per Zone-Year
# ----------------------------
site_summary <- site_draws_long %>%
  group_by(Zone, Year) %>%
  summarise(
    NMDS1_med  = median(NMDS1, na.rm = TRUE),
    NMDS1_low  = as.numeric(quantile(NMDS1, 0.025, na.rm = TRUE)),
    NMDS1_high = as.numeric(quantile(NMDS1, 0.975, na.rm = TRUE)),
    NMDS2_med  = median(NMDS2, na.rm = TRUE),
    NMDS2_low  = as.numeric(quantile(NMDS2, 0.025, na.rm = TRUE)),
    NMDS2_high = as.numeric(quantile(NMDS2, 0.975, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(Zone, Year)

write.csv(site_summary, out_sites_sum_csv, row.names = FALSE)
cat("Saved NMDS site summary:", out_sites_sum_csv, "\n")

# ----------------------------
# 8) EXPORT XLSX (short sheet names to avoid truncation)
# ----------------------------
writexl::write_xlsx(
  list(
    "Turnover_ZoneYear" = turnover_summary,
    "Turnover_ZoneSum"  = turnover_zone_summary,
    "NMDS_Sites_Sum"    = site_summary
  ),
  path = out_xlsx
)
cat("Saved XLSX:", out_xlsx, "\n")

# ----------------------------
# 9) PLOT A: INTERANNUAL TURNOVER WITH UNCERTAINTY
# ----------------------------
ymax_plot <- max(turnover_summary$delta_high, na.rm = TRUE) * 1.10

# ----------------------------
# ADD-ON: Overall turnover label per zone (median + 95% CrI across draws)
# ----------------------------

overall_turnover <- turnover_draws_long %>%
  group_by(Zone, draw_id) %>%
  summarise(mean_delta = mean(delta, na.rm = TRUE), .groups = "drop") %>%
  group_by(Zone) %>%
  summarise(
    mean_delta_med  = median(mean_delta, na.rm = TRUE),
    mean_delta_low  = as.numeric(quantile(mean_delta, 0.025, na.rm = TRUE)),
    mean_delta_high = as.numeric(quantile(mean_delta, 0.975, na.rm = TRUE)),
    .groups = "drop"
  )

# Pick a consistent in-panel position for the label (top-left of each facet)
label_pos <- turnover_summary %>%
  group_by(Zone) %>%
  summarise(
    x = start_year + 1.2,
    y = max(delta_high, na.rm = TRUE) * 0.95,
    .groups = "drop"
  )

overall_turnover_labels <- overall_turnover %>%
  left_join(label_pos, by = "Zone") %>%
  mutate(
    label = paste0(
      "Mean annual turnover:\n",
      sprintf("%.3f", mean_delta_med),
      " (", sprintf("%.3f", mean_delta_low), "–", sprintf("%.3f", mean_delta_high), ")"
    )
  )


g_turn <- ggplot(turnover_summary, aes(x = Year, y = delta_median, color = Zone, fill = Zone)) +
  geom_ribbon(aes(ymin = delta_low, ymax = delta_high),
              alpha = 0.25, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_fill_colors) +
  facet_wrap(~Zone, ncol = 1) +
  coord_cartesian(xlim = c(start_year + 1, end_year), ylim = c(0, ymax_plot)) +
  labs(
    x = "Year",
    y = "Δ NMDS distance (aligned)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "gray30"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

print(g_turn)
ggsave(out_turn_png, g_turn, width = 10, height = 8, dpi = 900)
ggsave(out_turn_jpg, g_turn, width = 10, height = 8, dpi = 900)

# ----------------------------
# 10) PLOT B: COMMUNITY SPACE WITH UNCERTAINTY (median + 95% CrI crosshairs)
#     + ellipses by Zone on median positions (same style as reference)
# ----------------------------
g_unc <- ggplot(site_summary, aes(x = NMDS1_med, y = NMDS2_med, color = Zone)) +
  # 95% CrI “crosshairs”
  geom_segment(aes(x = NMDS1_low, xend = NMDS1_high, y = NMDS2_med, yend = NMDS2_med),
               alpha = 0.35, linewidth = 0.7) +
  geom_segment(aes(x = NMDS1_med, xend = NMDS1_med, y = NMDS2_low, yend = NMDS2_high),
               alpha = 0.35, linewidth = 0.7) +
  geom_point(size = 2.6, alpha = 0.85) +
  stat_ellipse(type = "t", linetype = 2, linewidth = 1, alpha = 0.9) +
  geom_path(aes(group = Zone), alpha = 0.35, linewidth = 0.8) +
  scale_color_manual(values = zone_colors) +
  labs(
    x = "NMDS1 (median, aligned)",
    y = "NMDS2 (median, aligned)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(color = "gray30"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    legend.position = "none"
  )

print(g_unc)
ggsave(out_unc_ord_png, g_unc, width = 8, height = 6, dpi = 900)
ggsave(out_unc_ord_jpg, g_unc, width = 8, height = 6, dpi = 900)

# ----------------------------
# 11) PLOT C: DRAW CLOUD (diagnostic)
# ----------------------------
set.seed(1)
cloud_plot_df <- site_draws_long %>%
  sample_frac(cloud_frac)

g_cloud <- ggplot(cloud_plot_df, aes(x = NMDS1, y = NMDS2, color = Zone)) +
  geom_point(alpha = 0.05, size = 0.9) +
  geom_point(
    data = site_summary,
    aes(x = NMDS1_med, y = NMDS2_med, color = Zone),
    inherit.aes = FALSE,
    size = 2.1, alpha = 0.9
  ) +
  stat_ellipse(
    data = site_summary,
    aes(x = NMDS1_med, y = NMDS2_med, color = Zone),
    type = "t", linetype = 2, linewidth = 1,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = zone_colors) +
  labs(
    x = "NMDS1 (aligned)",
    y = "NMDS2 (aligned)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

print(g_cloud)
ggsave(out_cloud_png, g_cloud, width = 8, height = 6, dpi = 900)
ggsave(out_cloud_jpg, g_cloud, width = 8, height = 6, dpi = 900)

# ----------------------------
# 12) OPTIONAL: PERMANOVA on reference (posterior means) matrix
# ----------------------------
set.seed(123)
perm_zone <- adonis2(comm_ref ~ Zone, data = metadata_ref, method = "bray", permutations = 999)
perm_zone_year <- adonis2(comm_ref ~ Zone + Year, data = metadata_ref, method = "bray", permutations = 999)
perm_margin <- adonis2(comm_ref ~ Zone + Year, data = metadata_ref, method = "bray", permutations = 999, by = "margin")

cat("\n=== PERMANOVA on posterior MEAN community matrix ===\n")
cat("\nCommunity composition ~ Zone\n"); print(perm_zone)
cat("\nCommunity composition ~ Zone + Year\n"); print(perm_zone_year)
cat("\nCommunity composition ~ Zone + Year (marginal)\n"); print(perm_margin)

cat("\nDONE.\n")






# year gradient #
g_unc_year <- ggplot(site_summary, aes(x = NMDS1_med, y = NMDS2_med)) +
  geom_point(aes(color = Year), size = 2.6, alpha = 0.9) +
  geom_path(aes(group = Zone), alpha = 0.35, linewidth = 0.8) +
  facet_wrap(~Zone, ncol = 3) +
  scale_color_viridis_c() +
  theme_minimal(base_size = 14) +
  labs(title = "NMDS trajectories by year (median, aligned)", x = "NMDS1", y = "NMDS2")
print(g_unc_year)


# end labels #
label_df <- site_summary %>%
  filter(Year %in% c(start_year, end_year)) %>%
  mutate(lbl = paste0(Zone, " ", Year))

g_unc_labels <- ggplot(site_summary, aes(x = NMDS1_med, y = NMDS2_med, color = Zone)) +
  geom_point(size = 2.3, alpha = 0.85) +
  geom_path(aes(group = Zone), alpha = 0.35, linewidth = 0.8) +
  geom_text(data = label_df, aes(label = lbl), vjust = -0.8, show.legend = FALSE) +
  scale_color_manual(values = zone_colors) +
  theme_minimal(base_size = 14) +
  labs(title = "NMDS trajectories (median, aligned) with start/end labels", x = "NMDS1", y = "NMDS2")
print(g_unc_labels)





















############################################
# 4_FINAL_NMDS_Turnover_BetaDiversity_Outhwaite
# Streamlined standalone "final figures" script (run top-to-bottom)
#
# INPUTS (assumes these exist from earlier steps in your workflow):
#   - site_draws_long  : draw_id, Zone, Year, NMDS1, NMDS2   (aligned draw coordinates)
#   - site_summary     : Zone, Year, NMDS1_med, NMDS2_med    (median trajectory coords)
#   - turnover_summary : Zone, Year, delta_median, delta_low, delta_high
#   - turnover_draws_long : Zone, Year, draw_id, delta
#   - bc_summary       : Zone, bray_median, bray_low, bray_high  (start→end BC)
#
# OUTPUTS (only final patchwork saved; no intermediate plots saved):
#   - 4_FINAL_Patchwork_NMDS_Turnover_BrayCurtis_2000_2023.png/.jpg
#   - 4_NMDS_Turnover_ZoneDiagnostics_2000_2023.csv
#   - 4_NMDS_Turnover_PairwiseDifferences_2000_2023.csv
#
# FIGURE LAYOUT:
#   Top row:  NMDS cloud+trajectory+labels  |  Turnover (with uncertainty + text labels)
#   Bottom:   Bray–Curtis start→end % change by zone
#
# IMPORTANT:
# - This script does NOT print intermediate plots.
# - It only prints console diagnostics and saves the final patchwork figure.
############################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)     # alpha()
  library(patchwork)
})

# ----------------------------
# 0) USER SETTINGS
# ----------------------------
start_year <- 2000
end_year   <- 2023

zones_keep <- c("Core","Buffer","Outside")

zone_colors <- c("Core"="#b2182b","Buffer"="black","Outside"="grey70")
zone_fill_colors <- c("Core"=alpha("#b2182b",0.20),
                      "Buffer"=alpha("black",0.20),
                      "Outside"=alpha("grey40",0.30))

# ----------------------------
# 1) SANITY CHECKS (objects must exist)
# ----------------------------
req_objs <- c("site_draws_long","site_summary","turnover_summary","turnover_draws_long")
missing <- req_objs[!sapply(req_objs, exists)]
if (length(missing) > 0) {
  stop("Missing required objects: ", paste(missing, collapse = ", "),
       "\nRun your NMDS/turnover sections first to create these.")
}

# ----------------------------
# 1b) If bc_summary is missing, compute Bray–Curtis start→end from BUGS summary
# ----------------------------
if (!exists("bc_summary")) {
  
  # --- REQUIRED for BC computation:
  # You must point to your BUGS summary csv
  if (!exists("bugs_file")) {
    stop("bc_summary is missing AND bugs_file is not set.\n",
         "Set bugs_file to your Combined BUGS summary CSV.")
  }
  
  # if not already defined in your script:
  if (!exists("year_index_offset")) year_index_offset <- 1999
  if (!exists("psi_zone_regex")) psi_zone_regex <- "^psi\\.fs\\.r_(Core|Buffer|Outside)\\[(\\d+)\\]$"
  if (!exists("eps")) eps <- 1e-6
  
  # beta helpers (only define if not already present)
  if (!exists("beta_ab_from_mean_sd")) {
    beta_ab_from_mean_sd <- function(mu, sd) {
      if (!is.finite(mu) || !is.finite(sd)) return(list(a=NA_real_, b=NA_real_))
      if (mu <= 0 || mu >= 1) return(list(a=NA_real_, b=NA_real_))
      if (sd <= 0) return(list(a=NA_real_, b=NA_real_))
      v <- sd^2
      vmax <- mu * (1 - mu)
      if (!is.finite(v) || !is.finite(vmax) || vmax <= 0) return(list(a=NA_real_, b=NA_real_))
      if (v >= vmax) v <- 0.95 * vmax
      k <- mu * (1 - mu) / v - 1
      if (!is.finite(k) || k <= 0) return(list(a=NA_real_, b=NA_real_))
      a <- mu * k
      b <- (1 - mu) * k
      if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(list(a=NA_real_, b=NA_real_))
      list(a=a, b=b)
    }
  }
  
  if (!exists("draw_row")) {
    draw_row <- function(mu, sd, a, b, n = 1) {
      if (is.finite(a) && is.finite(b) && a > 0 && b > 0) {
        d <- rbeta(n, a, b)
      } else {
        sduse <- ifelse(is.finite(sd) && sd > 0, sd, 0.05)
        d <- rnorm(n, mean = mu, sd = sduse)
      }
      pmin(1 - eps, pmax(eps, d))
    }
  }
  
  # n_draws for BC: match NMDS draws if present, else default
  if (!exists("n_draws_nmds")) n_draws_nmds <- 300
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(vegan)
  })
  
  cat("\n[bc_summary missing] Computing Bray–Curtis start→end from BUGS summary...\n")
  
  bugs <- read.csv(bugs_file, stringsAsFactors = FALSE)
  stopifnot(all(c("Parameter","Species","mean","sd") %in% names(bugs)))
  
  psi_long <- bugs %>%
    filter(str_detect(Parameter, "^psi\\.fs\\.r_")) %>%
    mutate(
      Zone = str_match(Parameter, psi_zone_regex)[,2],
      YearIndex = suppressWarnings(as.integer(str_match(Parameter, psi_zone_regex)[,3]))
    ) %>%
    filter(!is.na(Zone), !is.na(YearIndex)) %>%
    mutate(
      Year = YearIndex + year_index_offset,
      Zone = factor(Zone, levels = zones_keep),
      Species = as.character(Species),
      psi_mean = suppressWarnings(as.numeric(mean)),
      psi_sd   = suppressWarnings(as.numeric(sd))
    ) %>%
    filter(Year %in% c(start_year, end_year)) %>%
    filter(!is.na(Species), Species != "", is.finite(psi_mean), is.finite(psi_sd)) %>%
    mutate(
      psi_mean = pmin(pmax(psi_mean, eps), 1 - eps),
      psi_sd   = pmax(psi_sd, 1e-12)
    ) %>%
    select(Species, Zone, Year, psi_mean, psi_sd)
  
  psi_par_se <- psi_long %>%
    rowwise() %>%
    mutate(
      ab = list(beta_ab_from_mean_sd(psi_mean, psi_sd)),
      a = ab$a, b = ab$b
    ) %>%
    ungroup() %>%
    select(-ab)
  
  # species universe: only those present in BOTH years within each zone
  species_levels <- psi_par_se %>%
    group_by(Zone) %>%
    summarise(sp = list(intersect(Species[Year == start_year], Species[Year == end_year])), .groups = "drop") %>%
    pull(sp) %>% unlist() %>% unique()
  
  psi_indexed_se <- psi_par_se %>%
    mutate(species_id = match(Species, species_levels)) %>%
    filter(!is.na(species_id)) %>%
    select(Zone, Year, species_id, psi_mean, psi_sd, a, b)
  
  build_comm_vec <- function(df_zone_year) {
    v <- rep(0, length(species_levels))
    if (nrow(df_zone_year) == 0) return(v)
    vals <- numeric(nrow(df_zone_year))
    for (i in seq_len(nrow(df_zone_year))) {
      vals[i] <- draw_row(df_zone_year$psi_mean[i], df_zone_year$psi_sd[i],
                          df_zone_year$a[i], df_zone_year$b[i], n = 1)
    }
    v[df_zone_year$species_id] <- vals
    v
  }
  
  set.seed(202)
  bc_list <- vector("list", n_draws_nmds)
  
  for (d in seq_len(n_draws_nmds)) {
    per_zone <- lapply(zones_keep, function(zn) {
      sdat <- psi_indexed_se %>% filter(Zone == zn, Year == start_year)
      fdat <- psi_indexed_se %>% filter(Zone == zn, Year == end_year)
      
      v_start <- build_comm_vec(sdat)
      v_end   <- build_comm_vec(fdat)
      
      bc <- as.numeric(vegdist(rbind(v_start, v_end), method = "bray")[1])
      
      tibble(draw_id = d,
             Zone = factor(zn, levels = zones_keep),
             bray_curtis = bc)
    })
    bc_list[[d]] <- bind_rows(per_zone)
  }
  
  bc_draws <- bind_rows(bc_list)
  
  bc_summary <- bc_draws %>%
    group_by(Zone) %>%
    summarise(
      bray_median = median(bray_curtis, na.rm = TRUE),
      bray_low    = as.numeric(quantile(bray_curtis, 0.025, na.rm = TRUE)),
      bray_high   = as.numeric(quantile(bray_curtis, 0.975, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(Zone)
  
  cat("[bc_summary created] Done.\n")
}


# Cloud thinning fraction (keep small for speed + nice density)
cloud_frac <- 0.15

# Ordination plot limits (set what you want; clip=off so labels can sit outside panel)
ord_xlim <- c(-0.25, 0.30)
ord_ylim <- c(-0.175, 0.2)

# Turnover label position in each facet
turnover_label_y <- 0.075

# Output
out_final_png <- paste0("4_FINAL_Patchwork_NMDS_Turnover_BrayCurtis_", start_year, "_", end_year, ".png")
out_final_jpg <- paste0("4_FINAL_Patchwork_NMDS_Turnover_BrayCurtis_", start_year, "_", end_year, ".jpg")

out_diag_csv <- paste0("4_NMDS_Turnover_ZoneDiagnostics_", start_year, "_", end_year, ".csv")
out_pair_csv <- paste0("4_NMDS_Turnover_PairwiseDifferences_", start_year, "_", end_year, ".csv")

# ----------------------------
# 1) SANITY CHECKS (objects must exist)
# ----------------------------
req_objs <- c("site_draws_long","site_summary","turnover_summary","turnover_draws_long","bc_summary")
missing <- req_objs[!sapply(req_objs, exists)]
if (length(missing) > 0) {
  stop("Missing required objects: ", paste(missing, collapse = ", "),
       "\nRun your NMDS/turnover + Bray–Curtis sections first to create these.")
}

# Lock factor order
site_summary <- site_summary %>%
  mutate(Zone = factor(Zone, levels = zones_keep))
site_draws_long <- site_draws_long %>%
  mutate(Zone = factor(Zone, levels = zones_keep))
turnover_summary <- turnover_summary %>%
  mutate(Zone = factor(Zone, levels = zones_keep))
turnover_draws_long <- turnover_draws_long %>%
  mutate(Zone = factor(Zone, levels = zones_keep))
bc_summary <- bc_summary %>%
  mutate(Zone = factor(Zone, levels = zones_keep))

# ----------------------------
# 2) BUILD ORDINATION CLOUD PLOT WITH START+END LABELS (same format)
# ----------------------------

set.seed(1)
cloud_plot_df2 <- site_draws_long %>%
  sample_frac(cloud_frac)

traj_med <- site_summary %>%
  arrange(Zone, Year)

# helper for Euclidean distance (used later for diagnostics)
euclid_step <- function(x1, y1, x0, y0) sqrt((x1 - x0)^2 + (y1 - y0)^2)

# Start labels (unified format: "Core (2000)")
# ---- Push START labels outward (radially from the zone centroid)
start_push <- 0.05  # <- increase for more distance (try 0.08 if needed)
label_gap <- 0.02   # distance to leave between line end and label

zone_centroids <- cloud_plot_df2 %>%
  group_by(Zone) %>%
  summarise(cx = median(NMDS1, na.rm = TRUE),
            cy = median(NMDS2, na.rm = TRUE),
            .groups = "drop")

start_labels <- traj_med %>%
  filter(Year == start_year) %>%
  left_join(zone_centroids, by = "Zone") %>%
  mutate(
    label = paste0(Zone, " (", start_year, ")"),
    # direction from centroid -> start point
    vx = NMDS1_med - cx,
    vy = NMDS2_med - cy,
    vlen = sqrt(vx^2 + vy^2),
    # unit vector (guard against vlen==0)
    ux = ifelse(is.finite(vlen) & vlen > 0, vx / vlen, 0),
    uy = ifelse(is.finite(vlen) & vlen > 0, vy / vlen, 1),
    # push label further along that direction
    x_lab = NMDS1_med + start_push * ux,
    y_lab = NMDS2_med + start_push * uy
  ) %>%
  select(Zone, Year, NMDS1_med, NMDS2_med, label, x_lab, y_lab)


# End points + end labels (unified format: "Core (2023)")
end_pts <- traj_med %>%
  filter(Year == end_year) %>%
  transmute(Zone, NMDS1_end = NMDS1_med, NMDS2_end = NMDS2_med)

prev_end_pts <- traj_med %>%
  group_by(Zone) %>%
  filter(Year == max(Year[Year < end_year], na.rm = TRUE)) %>%
  ungroup() %>%
  transmute(Zone, NMDS1_prev = NMDS1_med, NMDS2_prev = NMDS2_med)

end_connect <- end_pts %>%
  left_join(prev_end_pts, by = "Zone") %>%
  filter(is.finite(NMDS1_prev), is.finite(NMDS2_prev))

# label nudges for end labels (tweak if needed)
label_offsets <- tibble(
  Zone = factor(zones_keep, levels = zones_keep),
  dx = c(-0.015, 0.015, 0.020),
  dy = c( 0.015, 0.000, 0.015)
)

# ---- Push END labels outward (radially from the zone centroid)
end_push <- 0.06   # <- increase for more distance (Outside often needs more)

end_labels_pos <- end_pts %>%
  left_join(zone_centroids, by = "Zone") %>%
  mutate(
    label = paste0(Zone, " (", end_year, ")"),
    vx = NMDS1_end - cx,
    vy = NMDS2_end - cy,
    vlen = sqrt(vx^2 + vy^2),
    ux = ifelse(is.finite(vlen) & vlen > 0, vx / vlen, 1),
    uy = ifelse(is.finite(vlen) & vlen > 0, vy / vlen, 0),
    x_lab = NMDS1_end + end_push * ux,
    y_lab = NMDS2_end + end_push * uy
  ) %>%
  select(Zone, NMDS1_end, NMDS2_end, label, x_lab, y_lab)


g_cloud_labeled <- ggplot() +
  #geom_point(
  #  data = cloud_plot_df2,
  #  aes(x = NMDS1, y = NMDS2, color = Zone),
  #  alpha = 0.05, size = 0.9
  #) +
  # --- Add zone ellipses (draw cloud)
  stat_ellipse(
    data = cloud_plot_df2,
    aes(x = NMDS1, y = NMDS2, color = Zone),
    type = "t",
    linetype = 2,
    linewidth = 0.9,
    level = 0.95,
    show.legend = FALSE
  ) +
  geom_path(
    data = traj_med,
    aes(x = NMDS1_med, y = NMDS2_med, color = Zone, group = Zone),
    alpha = 0.55, linewidth = 1.0
  ) +
  geom_point(
    data = traj_med,
    aes(x = NMDS1_med, y = NMDS2_med, color = Zone),
    size = 1.3, alpha = 0.9
  ) +
  # Start markers + leaders + labels
  geom_point(
    data = start_labels,
    aes(x = NMDS1_med, y = NMDS2_med, color = Zone),
    shape = 21, fill = "white", stroke = 1.2, size = 3,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = start_labels,
    aes(
      x = NMDS1_med,
      y = NMDS2_med,
      xend = x_lab - label_gap * (x_lab - NMDS1_med) /
        sqrt((x_lab - NMDS1_med)^2 + (y_lab - NMDS2_med)^2),
      yend = y_lab - label_gap * (y_lab - NMDS2_med) /
        sqrt((x_lab - NMDS1_med)^2 + (y_lab - NMDS2_med)^2),
      color = Zone
    ),
    linewidth = 0.6,
    inherit.aes = FALSE, 
    lineend = "round"
  )+
  geom_text(
    data = start_labels,
    aes(x = x_lab, y = y_lab, label = label, color = Zone),
    size = 4, fontface = "bold",
    hjust = 0.5, vjust = 0.5,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  # Connector into end + end marker + leader + end label
  geom_segment(
    data = end_connect,
    aes(x = NMDS1_prev, y = NMDS2_prev, xend = NMDS1_end, yend = NMDS2_end, color = Zone),
    linewidth = 1.3, alpha = 0.9
  ) +
  geom_point(
    data = end_pts,
    aes(x = NMDS1_end, y = NMDS2_end, color = Zone),
    shape = 21, fill = "white", stroke = 1.4, size = 3.0
  ) +
  geom_segment(
    data = end_labels_pos,
    aes(
      x = NMDS1_end,
      y = NMDS2_end,
      xend = x_lab - label_gap * (x_lab - NMDS1_end) /
        sqrt((x_lab - NMDS1_end)^2 + (y_lab - NMDS2_end)^2),
      yend = y_lab - label_gap * (y_lab - NMDS2_end) /
        sqrt((x_lab - NMDS1_end)^2 + (y_lab - NMDS2_end)^2),
      color = Zone
    ),
    linewidth = 0.9,
    lineend = "round",
    alpha = 0.9
  )+
  geom_text(
    data = end_labels_pos,
    aes(x = x_lab, y = y_lab, label = label, color = Zone),
    fontface = "bold", size = 4.2, hjust = 0, vjust = 0.5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = zone_colors) +
  coord_cartesian(xlim = ord_xlim, ylim = ord_ylim, clip = "off") +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 45, 10, 10)  # extra right margin for end labels
  )

g_cloud_labeled

# ----------------------------
# 3) TURNOVER PLOT (with coloured text diagnostics per zone)
# ----------------------------

# Overall mean annual turnover per zone per draw -> summaries
mean_turn_by_draw <- turnover_draws_long %>%
  group_by(Zone, draw_id) %>%
  summarise(mean_delta = mean(delta, na.rm = TRUE), .groups = "drop")

overall_turnover <- mean_turn_by_draw %>%
  group_by(Zone) %>%
  summarise(
    mean_delta_med  = median(mean_delta, na.rm = TRUE),
    mean_delta_low  = as.numeric(quantile(mean_delta, 0.025, na.rm = TRUE)),
    mean_delta_high = as.numeric(quantile(mean_delta, 0.975, na.rm = TRUE)),
    .groups = "drop"
  )

overall_turnover_labels <- overall_turnover %>%
  mutate(
    x = start_year + 1.2,
    y = turnover_label_y,
    label = paste0(
      "Mean annual turnover = ",
      sprintf("%.3f", mean_delta_med),
      " (", sprintf("%.3f", mean_delta_low), "–", sprintf("%.3f", mean_delta_high), ")"
    )
  )

ymax_plot <- max(turnover_summary$delta_high, na.rm = TRUE) * 1.10

g_turn_labeled <- ggplot(turnover_summary, aes(x = Year, y = delta_median, color = Zone, fill = Zone)) +
  geom_ribbon(aes(ymin = delta_low, ymax = delta_high),
              alpha = 0.25, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_text(
    data = overall_turnover_labels,
    aes(x = x, y = 0.11, label = label, color = Zone),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    fontface = "bold",
    size = 4.3,
    show.legend = FALSE
  ) +
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_fill_colors) +
  facet_wrap(~Zone, ncol = 1) +
  coord_cartesian(xlim = c(start_year + 1, end_year), ylim = c(0, 0.125)) +
  labs(x = "Year", y = "Δ NMDS distance") +
  theme_minimal(base_size = 15) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

# ----------------------------
# 4) BRAY–CURTIS START→END PLOT (with end caps)
# ----------------------------

bc_plot_df <- bc_summary %>%
  mutate(
    Zone = factor(Zone, levels = zones_keep),
    pct_median = 100 * bray_median,
    pct_low    = 100 * bray_low,
    pct_high   = 100 * bray_high
  )

g_bc <- ggplot(bc_plot_df, aes(x = Zone, y = pct_median, color = Zone)) +
  geom_errorbar(aes(ymin = pct_low, ymax = pct_high),
                width = 0.25, linewidth = 1.4, show.legend = FALSE) +
  geom_point(size = 5.2, show.legend = FALSE) +
  scale_color_manual(values = zone_colors) +
  labs(x = NULL, y = "% community change") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(35, 45))

bc_pct_summary <- bc_summary %>%
  mutate(
    pct_median = 100 * bray_median,
    pct_low    = 100 * bray_low,
    pct_high   = 100 * bray_high
  ) %>%
  arrange(Zone)

out_bc_csv <- paste0("4_BrayCurtis_StartEnd_PercentChange_", start_year, "_", end_year, ".csv")
write.csv(bc_pct_summary, out_bc_csv, row.names = FALSE)
cat("\nSaved Bray–Curtis % community change summary:\n  ", out_bc_csv, "\n")
print(bc_pct_summary)


# ============================================================
# 7) ADD-ON: % LOSS and % GAIN (per-species start→end, uncertainty-aware)
#    - Two bars per zone: loss (negative) and gain (positive)
#    - Uses per-species Beta-draw occupancy for start/end years
#    - Adds plot to patchwork to the RIGHT of Bray–Curtis panel
# ============================================================

# ---- REQUIREMENTS:
# This section needs per-species BUGS-derived parameters already in memory:
#   - psi_par (Species, Zone, Year, psi_mean, psi_sd, a, b)  OR equivalent
#   - species_levels (character vector of species columns used in NMDS reference)
#   - draw_row() helper (same as earlier; returns 1 draw with eps truncation)
# If these are missing, stop with a clear message.

need_objs <- c("psi_par", "species_levels", "draw_row")
missing2 <- need_objs[!sapply(need_objs, exists)]
if (length(missing2) > 0) {
  stop(
    "To compute % loss/% gain you need these objects in memory: ",
    paste(missing2, collapse = ", "),
    "\nRun the earlier BUGS extraction + indexing sections that create psi_par/species_levels/draw_row."
  )
}

# Build per-species index table for start/end years using the same species set as NMDS
psi_indexed_se <- psi_par %>%
  filter(Year %in% c(start_year, end_year)) %>%
  mutate(
    Zone = factor(Zone, levels = zones_keep),
    species_id = match(Species, species_levels)
  ) %>%
  filter(!is.na(species_id)) %>%
  select(Zone, Year, species_id, psi_mean, psi_sd, a, b)

# Helper: build a numeric vector (length = n species) for ONE zone-year for ONE draw
build_comm_vec_draw <- function(df_zone_year) {
  v <- numeric(length(species_levels))
  if (nrow(df_zone_year) == 0) return(v)
  
  vals <- numeric(nrow(df_zone_year))
  for (i in seq_len(nrow(df_zone_year))) {
    vals[i] <- draw_row(
      mu = df_zone_year$psi_mean[i],
      sd = df_zone_year$psi_sd[i],
      a  = df_zone_year$a[i],
      b  = df_zone_year$b[i],
      n  = 1
    )
  }
  v[df_zone_year$species_id] <- vals
  v
}

# Compute loss/gain per draw per zone
set.seed(303)
lg_list <- vector("list", n_draws_nmds)

for (d in seq_len(n_draws_nmds)) {
  tmp_zone <- vector("list", length(zones_keep))
  
  for (z in seq_along(zones_keep)) {
    zn <- zones_keep[z]
    
    sdat <- psi_indexed_se %>% filter(Zone == zn, Year == start_year)
    fdat <- psi_indexed_se %>% filter(Zone == zn, Year == end_year)
    
    s <- build_comm_vec_draw(sdat)
    f <- build_comm_vec_draw(fdat)
    
    denom <- sum(s, na.rm = TRUE)
    denom <- ifelse(is.finite(denom) && denom > 0, denom, NA_real_)
    
    loss_raw <- sum(pmax(0, s - f), na.rm = TRUE)
    gain_raw <- sum(pmax(0, f - s), na.rm = TRUE)
    
    pct_loss <- 100 * (loss_raw / denom)
    pct_gain <- 100 * (gain_raw / denom)
    
    tmp_zone[[z]] <- tibble(
      draw_id = d,
      Zone = factor(zn, levels = zones_keep),
      pct_loss = pct_loss,
      pct_gain = pct_gain
    )
  }
  
  lg_list[[d]] <- bind_rows(tmp_zone)
}

lg_draws <- bind_rows(lg_list)

# Summarise (median + 95% CrI) by zone
lg_summary <- lg_draws %>%
  group_by(Zone) %>%
  summarise(
    loss_median = median(pct_loss, na.rm = TRUE),
    loss_low    = as.numeric(quantile(pct_loss, 0.025, na.rm = TRUE)),
    loss_high   = as.numeric(quantile(pct_loss, 0.975, na.rm = TRUE)),
    gain_median = median(pct_gain, na.rm = TRUE),
    gain_low    = as.numeric(quantile(pct_gain, 0.025, na.rm = TRUE)),
    gain_high   = as.numeric(quantile(pct_gain, 0.975, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(Zone)

# Long format for plotting (loss negative, gain positive)
lg_plot_df <- lg_summary %>%
  transmute(
    Zone,
    type = "Loss",
    median = -loss_median, low = -loss_high, high = -loss_low
  ) %>%
  bind_rows(
    lg_summary %>%
      transmute(
        Zone,
        type = "Gain",
        median = gain_median, low = gain_low, high = gain_high
      )
  ) %>%
  mutate(
    type = factor(type, levels = c("Loss","Gain")),
    Zone = factor(Zone, levels = zones_keep)
  )

# ---- Make zone-tinted colours for Loss/Gain
zone_loss_cols <- setNames(scales::alpha(zone_colors, 0.80), names(zone_colors))  # darker
zone_gain_cols <- setNames(scales::alpha(zone_colors, 0.35), names(zone_colors))  # lighter

lg_plot_df2 <- lg_plot_df %>%
  mutate(
    Zone = factor(Zone, levels = zones_keep),
    fill_key = paste0(as.character(Zone), "_", as.character(type))
  )

fill_map <- c(
  setNames(zone_loss_cols, paste0(names(zone_loss_cols), "_Loss")),
  setNames(zone_gain_cols, paste0(names(zone_gain_cols), "_Gain"))
)


# Plot: two bars per zone, centered on 0
g_loss_gain <- ggplot(lg_plot_df2, aes(x = Zone, y = median, fill = fill_key)) +
  geom_hline(yintercept = 0, linewidth = 0.8, color = "grey40") +
  geom_col(position = position_dodge(width = 0.68), width = 0.62) +
  geom_errorbar(
    aes(ymin = low, ymax = high),
    width = 0.18,
    linewidth = 1.0,
    position = position_dodge(width = 0.68),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = fill_map) +
  labs(
    x = NULL,
    y = "% loss / % gain ",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 13)
  )



# (Optional) export these summaries
out_lg_csv <- paste0("4_LossGain_Summary_", start_year, "_", end_year, ".csv")
write.csv(lg_summary, out_lg_csv, row.names = FALSE)
cat("\nSaved loss/gain summary:", out_lg_csv, "\n")
print(lg_summary)


# ----------------------------
# 5) CONSOLE DIAGNOSTICS + CSV EXPORTS
# ----------------------------

# A) Median-trajectory diagnostics: total path length, net displacement, directionality
path_len <- traj_med %>%
  group_by(Zone) %>%
  arrange(Year) %>%
  summarise(
    total_path_length = sum(sqrt((NMDS1_med - lag(NMDS1_med))^2 + (NMDS2_med - lag(NMDS2_med))^2), na.rm = TRUE),
    .groups = "drop"
  )

disp <- traj_med %>%
  group_by(Zone) %>%
  summarise(
    NMDS1_start = NMDS1_med[Year == start_year][1],
    NMDS2_start = NMDS2_med[Year == start_year][1],
    NMDS1_end   = NMDS1_med[Year == end_year][1],
    NMDS2_end   = NMDS2_med[Year == end_year][1],
    net_displacement = euclid_step(NMDS1_end, NMDS2_end, NMDS1_start, NMDS2_start),
    .groups = "drop"
  )

directionality <- path_len %>%
  inner_join(disp %>% select(Zone, net_displacement), by = "Zone") %>%
  mutate(directionality_ratio = ifelse(total_path_length > 0, net_displacement / total_path_length, NA_real_))

zone_diagnostics <- overall_turnover %>%
  left_join(directionality, by = "Zone") %>%
  arrange(desc(mean_delta_med))

# B) Pairwise differences in mean annual turnover (draw-based)
mean_turn_wide <- mean_turn_by_draw %>%
  select(Zone, draw_id, mean_delta) %>%
  pivot_wider(names_from = Zone, values_from = mean_delta)

zone_pairs <- t(combn(zones_keep, 2))
pairwise_list <- vector("list", nrow(zone_pairs))

for (i in seq_len(nrow(zone_pairs))) {
  A <- zone_pairs[i, 1]
  B <- zone_pairs[i, 2]
  dvec <- mean_turn_wide[[A]] - mean_turn_wide[[B]]
  dvec <- dvec[is.finite(dvec)]
  low  <- as.numeric(quantile(dvec, 0.025, na.rm = TRUE))
  high <- as.numeric(quantile(dvec, 0.975, na.rm = TRUE))
  
  pairwise_list[[i]] <- tibble(
    Zone_A = A,
    Zone_B = B,
    diff_median = median(dvec, na.rm = TRUE),
    diff_low = low,
    diff_high = high,
    Pr_A_greater_B = mean(dvec > 0, na.rm = TRUE),
    Credible_Different_95CI = (low > 0) | (high < 0)
  )
}

pairwise_turnover_diffs <- bind_rows(pairwise_list) %>%
  arrange(desc(abs(diff_median)))

# Console prints (only diagnostics)
cat("\n========================================================\n")
cat("OVERALL TURNOVER SPEED (mean annual ΔNMDS distance; med + 95% CrI)\n")
cat("Higher = faster average year-to-year compositional turnover\n")
cat("========================================================\n")
print(overall_turnover)

cat("\n========================================================\n")
cat("MAIN CHANGE DIAGNOSTICS (median trajectories)\n")
cat(" - total_path_length: cumulative movement through NMDS space\n")
cat(" - net_displacement: straight-line start→end shift\n")
cat(" - directionality_ratio: net / total (closer to 1 = more directional)\n")
cat("========================================================\n")
print(directionality)

cat("\n========================================================\n")
cat("PAIRWISE DIFFERENCES in mean annual turnover (draw-based)\n")
cat(" - Credible_Different_95CI = TRUE means 95% CrI excludes 0\n")
cat("========================================================\n")
print(pairwise_turnover_diffs)

# Export CSVs
write.csv(zone_diagnostics, out_diag_csv, row.names = FALSE)
write.csv(pairwise_turnover_diffs, out_pair_csv, row.names = FALSE)
cat("\nSaved:\n  ", out_diag_csv, "\n  ", out_pair_csv, "\n")

# ----------------------------
# 8) FINAL PATCHWORK (UPDATED): bottom row has BC (left) and Loss/Gain (right)
# ----------------------------

top_row <- g_cloud_labeled + g_turn_labeled +
  plot_layout(widths = c(1, 1))

bottom_row <- g_bc + g_loss_gain +
  plot_layout(widths = c(1, 1))

final_fig <- top_row / bottom_row +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    tag_levels = "A"   # <-- adds A, B, C, D automatically
  ) &
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0.01, 0.99)   # top-left of each panel
  )

print(final_fig)

ggsave(out_final_png, final_fig, width = 10, height = 8, dpi = 1000)
ggsave(out_final_jpg, final_fig, width = 10, height = 8, dpi = 1000)

cat("\nSaved FINAL patchwork figure:\n  ", out_final_png, "\n  ", out_final_jpg, "\n")
























############################################
# ADD-ON: annotate NMDS with stress + PERMANOVA R2, then re-plot main figure
############################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# --- safety checks ---
needed_objs <- c("g_cloud_labeled", "ord_xlim", "ord_ylim")
missing <- needed_objs[!sapply(needed_objs, exists)]
if (length(missing) > 0) {
  stop("Missing objects needed for annotation: ",
       paste(missing, collapse = ", "),
       "\nRun your FINAL NMDS/turnover patchwork script first.")
}

if (!exists("nmds_ref")) {
  warning("Object 'nmds_ref' (reference NMDS) not found; stress will be reported as NA.")
}

if (!exists("perm_margin") && !exists("perm_zone_year")) {
  warning("PERMANOVA objects 'perm_margin'/'perm_zone_year' not found; PERMANOVA R² text will be NA.")
}

# --- 1) Extract stress ---
stress_val <- if (exists("nmds_ref")) nmds_ref$stress else NA_real_
stress_lab <- paste0("Stress = ", sprintf("%.3f", stress_val))

# --- 2) Extract PERMANOVA R² for Zone and Year ---
zone_r2 <- NA_real_
year_r2 <- NA_real_

if (exists("perm_margin")) {
  pm <- as.data.frame(perm_margin)
  if (all(c("Zone","Year") %in% rownames(pm)) && "R2" %in% colnames(pm)) {
    zone_r2 <- pm["Zone", "R2"]
    year_r2 <- pm["Year", "R2"]
  }
} else if (exists("perm_zone_year")) {
  # fall back: use overall model R² split by terms (less ideal, but better than nothing)
  pz <- as.data.frame(perm_zone_year)
  if (nrow(pz) >= 2 && "R2" %in% colnames(pz)) {
    zone_r2 <- pz[1, "R2"]
    year_r2 <- pz[2, "R2"]
  }
}

perm_lab <- paste0(
  "PERMANOVA (Bray–Curtis)\n",
  "R² Zone = ", sprintf("%.3f", zone_r2),
  "; R² Year = ", sprintf("%.3f", year_r2)
)

# --- 3) Add annotations to the ordination plot ---

g_cloud_labeled_annot <- g_cloud_labeled +
  # PERMANOVA label (top-left)
  annotate(
    "text",
    x = ord_xlim[2], y = ord_ylim[1],
    hjust = 0.6, vjust = -1,
    label = perm_lab,
    size = 3.5,
    fontface = "bold"
  ) +
  # stress label (bottom-right)
  annotate(
    "text",
    x = ord_xlim[2], y = ord_ylim[1],
    hjust = 1, vjust = 0,
    label = stress_lab,
    size = 3.5,
    fontface = "bold"
  )

# --- 4) Rebuild the patchwork with the updated ordination panel ---

# These objects already exist from your FINAL script:
#   g_turn_labeled, g_bc, g_loss_gain, start_year, end_year,
#   out_final_png, out_final_jpg
needed_fig_parts <- c("g_turn_labeled", "g_bc", "g_loss_gain",
                      "start_year", "end_year",
                      "out_final_png", "out_final_jpg")

missing2 <- needed_fig_parts[!sapply(needed_fig_parts, exists)]
if (length(missing2) > 0) {
  stop("Missing objects needed to rebuild patchwork: ",
       paste(missing2, collapse = ", "),
       "\nRe-run your FINAL patchwork script, then this add-on.")
}

top_row_new <- g_cloud_labeled_annot + g_turn_labeled +
  plot_layout(widths = c(1, 1))

bottom_row_new <- g_bc + g_loss_gain +
  plot_layout(widths = c(1, 1))

final_fig_annot <- top_row_new / bottom_row_new +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0.01, 0.99)
  )

# Print to device
print(final_fig_annot)

# Overwrite existing files with annotated version
ggsave(out_final_png, final_fig_annot, width = 10, height = 8, dpi = 1000)
ggsave(out_final_jpg, final_fig_annot, width = 10, height = 8, dpi = 1000)

cat("\nRe-saved FINAL patchwork with stress + PERMANOVA annotations:\n  ",
    out_final_png, "\n  ", out_final_jpg, "\n")




























############################################
# ADD-ON: Species-level increases & decreases
#         (counts + mean magnitude per zone)
############################################

suppressPackageStartupMessages({
  library(dplyr)
})

# ---- Requirements: same objects as loss/gain section ----
need_objs <- c("psi_par", "species_levels", "draw_row",
               "zones_keep", "start_year", "end_year", "n_draws_nmds")
missing3 <- need_objs[!sapply(need_objs, exists)]
if (length(missing3) > 0) {
  stop(
    "To compute species-level increases/decreases you need: ",
    paste(missing3, collapse = ", "),
    "\nRun your NMDS / BUGS extraction sections first."
  )
}

# If psi_indexed_se / build_comm_vec_draw not already defined, create them
if (!exists("psi_indexed_se")) {
  psi_indexed_se <- psi_par %>%
    filter(Year %in% c(start_year, end_year)) %>%
    mutate(
      Zone = factor(Zone, levels = zones_keep),
      species_id = match(Species, species_levels)
    ) %>%
    filter(!is.na(species_id)) %>%
    select(Zone, Year, species_id, psi_mean, psi_sd, a, b)
}

if (!exists("build_comm_vec_draw")) {
  build_comm_vec_draw <- function(df_zone_year) {
    v <- numeric(length(species_levels))
    if (nrow(df_zone_year) == 0) return(v)
    
    vals <- numeric(nrow(df_zone_year))
    for (i in seq_len(nrow(df_zone_year))) {
      vals[i] <- draw_row(
        mu = df_zone_year$psi_mean[i],
        sd = df_zone_year$psi_sd[i],
        a  = df_zone_year$a[i],
        b  = df_zone_year$b[i],
        n  = 1
      )
    }
    v[df_zone_year$species_id] <- vals
    v
  }
}

# ----------------------------
# 1) Per-draw species change stats
# ----------------------------
set.seed(404)  # separate seed for reproducibility
change_list <- vector("list", n_draws_nmds)

for (d in seq_len(n_draws_nmds)) {
  zone_list <- vector("list", length(zones_keep))
  
  for (z in seq_along(zones_keep)) {
    zn <- zones_keep[z]
    
    # Start & end per-species occupancy draws for this zone
    sdat <- psi_indexed_se %>% filter(Zone == zn, Year == start_year)
    fdat <- psi_indexed_se %>% filter(Zone == zn, Year == end_year)
    
    s <- build_comm_vec_draw(sdat)  # start occupancies
    f <- build_comm_vec_draw(fdat)  # end occupancies
    
    diff <- f - s  # positive = increase, negative = decrease
    
    inc_idx <- which(diff > 0)
    dec_idx <- which(diff < 0)
    
    n_inc <- length(inc_idx)
    n_dec <- length(dec_idx)
    
    mean_inc_mag <- if (n_inc > 0) mean(diff[inc_idx]) else NA_real_
    mean_dec_mag <- if (n_dec > 0) mean(abs(diff[dec_idx])) else NA_real_
    
    zone_list[[z]] <- tibble(
      draw_id      = d,
      Zone         = factor(zn, levels = zones_keep),
      n_increase   = n_inc,
      n_decrease   = n_dec,
      mean_inc_mag = mean_inc_mag,     # average positive change
      mean_dec_mag = mean_dec_mag      # average |negative change|
    )
  }
  
  change_list[[d]] <- bind_rows(zone_list)
}

species_change_draws <- bind_rows(change_list)

# ----------------------------
# 2) Summarise across draws (median + 95% CrI) per zone
# ----------------------------
species_change_summary <- species_change_draws %>%
  group_by(Zone) %>%
  summarise(
    # counts
    n_increase_median = median(n_increase, na.rm = TRUE),
    n_increase_low    = as.numeric(quantile(n_increase, 0.025, na.rm = TRUE)),
    n_increase_high   = as.numeric(quantile(n_increase, 0.975, na.rm = TRUE)),
    
    n_decrease_median = median(n_decrease, na.rm = TRUE),
    n_decrease_low    = as.numeric(quantile(n_decrease, 0.025, na.rm = TRUE)),
    n_decrease_high   = as.numeric(quantile(n_decrease, 0.975, na.rm = TRUE)),
    
    # magnitudes (on occupancy scale)
    mean_inc_mag_med  = median(mean_inc_mag, na.rm = TRUE),
    mean_inc_mag_low  = as.numeric(quantile(mean_inc_mag, 0.025, na.rm = TRUE)),
    mean_inc_mag_high = as.numeric(quantile(mean_inc_mag, 0.975, na.rm = TRUE)),
    
    mean_dec_mag_med  = median(mean_dec_mag, na.rm = TRUE),
    mean_dec_mag_low  = as.numeric(quantile(mean_dec_mag, 0.025, na.rm = TRUE)),
    mean_dec_mag_high = as.numeric(quantile(mean_dec_mag, 0.975, na.rm = TRUE)),
    
    .groups = "drop"
  ) %>%
  arrange(Zone)

# Optional: save to CSV
out_spchange_draws_csv <- paste0("4_SpeciesChange_Draws_", start_year, "_", end_year, ".csv")
out_spchange_sum_csv   <- paste0("4_SpeciesChange_Summary_", start_year, "_", end_year, ".csv")

write.csv(species_change_draws,  out_spchange_draws_csv, row.names = FALSE)
write.csv(species_change_summary, out_spchange_sum_csv,  row.names = FALSE)

cat("\nSaved species change draws:\n  ", out_spchange_draws_csv, "\n")
cat("Saved species change summary:\n  ", out_spchange_sum_csv, "\n")
print(species_change_summary)
