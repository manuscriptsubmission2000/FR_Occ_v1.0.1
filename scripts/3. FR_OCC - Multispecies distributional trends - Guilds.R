############################################
# GUILD composite trends + growth rates (2000–2023)
############################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(forcats)
  library(grid)     # unit()
  library(scales)   # alpha()
})

# ----------------------------
# 0) USER SETTINGS
# ----------------------------
setwd("anon")

bugs_file  <- "Combined_BUGS_Data_Outhwaite_finalclustv5_32000iterations_FINAL_WITH_TAXCORR.csv"
guild_path <- "combined_all_parks_with_guild_TAXCORR_ADDED.csv"

start_year <- 2000
end_year   <- 2023
years_n    <- end_year - start_year  # denominator for annual growth
year_index_offset <- 1999            # Year = YearIndex + 1999
zones_keep <- c("Core","Buffer","Outside")

# Optional: drop any guilds
#remove_guilds <- c("Ants")  # set to character(0) to keep all
remove_guilds <- character(0)  # set to character(0) to keep all

# Monte Carlo settings
set.seed(123)
n_draws <- 999
eps     <- 1e-6

# Regex to match exactly: psi.fs.r_Core[10], etc.
psi_zone_regex <- "^psi\\.fs\\.r_(Core|Buffer|Outside)\\[(\\d+)\\]$"

# Outputs
out_csv_clean  <- "BUGS_psi_fs_r_byZone_byYear_2000_2023.csv"
out_csv_trend  <- "guild_composite_trends_2000_2023.csv"
out_csv_growth <- "guild_composite_growth_rates_2000_2023.csv"

dir_out <- "GuildCompositeTrendPlots"
out_png_plot <- file.path(dir_out, paste0("guild_composite_trends_", start_year, "_", end_year, "_endLabels_%perYear.png"))
out_jpg_plot <- file.path(dir_out, paste0("guild_composite_trends_", start_year, "_", end_year, "_endLabels_%perYear.jpg"))

# Plot styling (Zones)
zone_colors <- c("Core"="#b2182b","Buffer"="black","Outside"="grey65")
zone_fill_colors <- c("Core"=alpha("#b2182b",0.20),
                      "Buffer"=alpha("black",0.20),
                      "Outside"=alpha("grey40",0.30))

# ----------------------------
# 1) HELPERS
# ----------------------------
std_key <- function(x) str_squish(str_to_lower(as.character(x)))

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

geo_mean <- function(x) {
  x <- x[is.finite(x) & !is.na(x) & x > 0]
  if (length(x) == 0) return(NA_real_)
  exp(mean(log(x)))
}

growth_rate_pct <- function(s, f, y) {
  if (!is.finite(s) || !is.finite(f) || s <= 0 || f <= 0) return(NA_real_)
  (((f / s)^(1 / y)) - 1) * 100
}

draw_row <- function(mu, sd, a, b, n = n_draws) {
  # Prefer Beta if feasible; fallback to truncated normal if Beta params invalid
  if (is.finite(a) && is.finite(b) && a > 0 && b > 0) {
    d <- rbeta(n, a, b)
  } else {
    sduse <- ifelse(is.finite(sd) && sd > 0, sd, 0.05)
    d <- rnorm(n, mean = mu, sd = sduse)
  }
  pmin(1 - eps, pmax(eps, d))
}

# Non-overlap label helper (per facet)
disperse_labels <- function(df, gap = 0.06) {
  # expects: label_base, ymin, ymax, pad
  n <- nrow(df)
  low  <- df$ymin[1] + df$pad[1]
  high <- df$ymax[1] - df$pad[1]
  avail <- high - low
  
  if (n <= 1 || !is.finite(avail) || avail <= 0) {
    df$label_y_adj <- pmin(pmax(df$label_base, low), high)
    return(df)
  }
  
  gap <- min(gap, avail / max(1, n - 1))
  
  ord <- order(df$label_base)
  y0  <- df$label_base[ord]
  
  y <- numeric(n)
  y[1] <- max(y0[1], low)
  for (i in 2:n) y[i] <- max(y0[i], y[i - 1] + gap)
  
  overflow <- y[n] - high
  if (overflow > 0) {
    y <- y - overflow
    if (y[1] < low) {
      y <- seq(low, high, length.out = n)
    }
  }
  
  df$label_y_adj <- NA_real_
  df$label_y_adj[ord] <- y
  df
}

# ----------------------------
# 2) LOAD BUGS + EXTRACT psi.fs.r_{Zone}[YearIndex]
# ----------------------------
bugs <- read.csv(bugs_file, stringsAsFactors = FALSE)

psi_zone <- bugs %>%
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
    sp_key = std_key(Species),
    psi_mean = suppressWarnings(as.numeric(mean)),
    psi_sd   = suppressWarnings(as.numeric(sd))
  ) %>%
  filter(Year >= start_year, Year <= end_year) %>%
  select(Species, sp_key, Zone, Year, psi_mean, psi_sd)

write.csv(psi_zone, out_csv_clean, row.names = FALSE)

cat("\n========================================================\n")
cat("Loaded psi.fs.r_{Zone}[YearIndex]\n")
cat("========================================================\n")
cat("Rows:", nrow(psi_zone), "\n")
cat("Species:", n_distinct(psi_zone$sp_key), "\n")
cat("Years:", paste(range(psi_zone$Year), collapse="–"), "\n")
cat("Zones:\n"); print(table(psi_zone$Zone, useNA="ifany"))
cat("Saved:", out_csv_clean, "\n")

# ----------------------------
# 3) LOAD + JOIN GUILD
# ----------------------------
guild_df <- read.csv(guild_path, stringsAsFactors = FALSE)

species_candidates <- c("species","Species","species_id","speciesID","species_name",
                        "Species_name","taxon","Taxon","scientificName","ScientificName")
guild_candidates   <- c("Guild","guild","GUILD")

species_col <- intersect(species_candidates, names(guild_df))[1]
guild_col   <- intersect(guild_candidates,   names(guild_df))[1]
if (is.na(species_col) || is.na(guild_col)) {
  stop("Could not detect Species and/or Guild columns in guild file. ",
       "Make sure it has a species column and a Guild column.")
}



guild_lookup <- guild_df %>%
  transmute(
    sp_key = std_key(.data[[species_col]]),
    Guild  = as.character(.data[[guild_col]])
  ) %>%
  filter(!is.na(sp_key), sp_key != "", !is.na(Guild), Guild != "") %>%
  distinct(sp_key, Guild)

psi_guild <- psi_zone %>%
  left_join(guild_lookup, by = "sp_key") %>%
  filter(!is.na(Guild)) %>%
  filter(!(Guild %in% remove_guilds)) %>%
  mutate(
    Guild = as.character(Guild)
  )

cat("\n========================================================\n")
cat("Guild join summary\n")
cat("========================================================\n")
cat("Species in psi_zone:", n_distinct(psi_zone$sp_key), "\n")
cat("Species with Guild  :", n_distinct(psi_guild$sp_key), "\n")
cat("Guilds:", n_distinct(psi_guild$Guild), "\n")

# Add beta parameters
psi_par <- psi_guild %>%
  rowwise() %>%
  mutate(
    ab = list(beta_ab_from_mean_sd(
      mu = pmin(pmax(psi_mean, eps), 1 - eps),
      sd = psi_sd
    )),
    a = ab$a,
    b = ab$b
  ) %>%
  ungroup() %>%
  select(-ab)

# ----------------------------
# ADD-ON A) Extract Rhat for psi.fs.r_{Zone}[YearIndex] and compute convergence %
# ----------------------------

# 1) Pull Rhat for the same psi.fs.r_ parameters you already extracted
psi_rhat <- bugs %>%
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
    sp_key = std_key(Species),
    Rhat = suppressWarnings(as.numeric(Rhat)),
    converged = !is.na(Rhat) & is.finite(Rhat) & (Rhat <= 1.1)
  ) %>%
  filter(Year >= start_year, Year <= end_year) %>%
  select(sp_key, Zone, Year, Rhat, converged)

# 2) Attach convergence flag onto psi_par (so it stays aligned with your guild filtering)
psi_par <- psi_par %>%
  left_join(psi_rhat %>% select(sp_key, Zone, Year, converged),
            by = c("sp_key","Zone","Year"))

# 3) Compute "% convergence per species in each year" (across zones)
#    - For each species-year: mean(converged across zones with data) * 100
species_year_conv <- psi_par %>%
  filter(Zone %in% zones_keep) %>%
  group_by(Guild, sp_key, Year) %>%
  summarise(
    conv_pct_species_year = 100 * mean(converged, na.rm = TRUE),
    .groups = "drop"
  )

# 4) Apply to guilds: average across species-year values (unweighted)
guild_conv <- species_year_conv %>%
  group_by(Guild) %>%
  summarise(
    Convergence_Pct = mean(conv_pct_species_year, na.rm = TRUE),
    .groups = "drop"
  )

# Optional quick check
cat("\n========================================================\n")
cat("Convergence summary (Rhat <= 1.1)\n")
cat("========================================================\n")
print(guild_conv %>% arrange(desc(Convergence_Pct)))


# ----------------------------
# 4) COMPOSITE TRENDS by Guild × Zone × Year (median + 95% CrI)
# ----------------------------
guilds_vec <- sort(unique(as.character(psi_par$Guild)))
years_vec  <- start_year:end_year

composite_list <- list()

set.seed(111)
for (gd in guilds_vec) {
  for (zn in zones_keep) {
    sub_gz <- psi_par %>% filter(Guild == gd, Zone == zn)
    if (nrow(sub_gz) == 0) next
    
    for (yr in years_vec) {
      sub <- sub_gz %>% filter(Year == yr)
      if (nrow(sub) == 0) next
      
      n_sp <- nrow(sub)
      draws_mat <- matrix(NA_real_, nrow = n_sp, ncol = n_draws)
      
      for (i in seq_len(n_sp)) {
        draws_mat[i, ] <- draw_row(sub$psi_mean[i], sub$psi_sd[i], sub$a[i], sub$b[i], n_draws)
      }
      
      comp_draws <- apply(draws_mat, 2, geo_mean)
      
      composite_list[[length(composite_list) + 1]] <- tibble(
        Guild = gd,
        Zone  = zn,
        Year  = yr,
        n_species = n_sp,
        comp_median = median(comp_draws, na.rm = TRUE),
        comp_low    = as.numeric(quantile(comp_draws, 0.025, na.rm = TRUE)),
        comp_high   = as.numeric(quantile(comp_draws, 0.975, na.rm = TRUE))
      )
    }
  }
}

guild_composite_trends <- bind_rows(composite_list) %>%
  mutate(
    Zone = factor(Zone, levels = zones_keep)
  ) %>%
  arrange(Guild, Zone, Year)

write.csv(guild_composite_trends, out_csv_trend, row.names = FALSE)
cat("\nSaved:", out_csv_trend, "\n")

# ----------------------------
# 5) GUILD × ZONE growth rate (% per year), uncertainty on draws
#     Growth labels will use ONLY the median (%/yr)
# ----------------------------
set.seed(222)

growth_list <- list()

for (gd in guilds_vec) {
  for (zn in zones_keep) {
    
    sub_gz <- psi_par %>% filter(Guild == gd, Zone == zn)
    if (nrow(sub_gz) == 0) next
    
    sdat <- sub_gz %>% filter(Year == start_year)
    fdat <- sub_gz %>% filter(Year == end_year)
    
    # Species common to start & end (Cooke/Outhwaite style)
    common_sp <- intersect(sdat$sp_key, fdat$sp_key)
    sdat <- sdat %>% filter(sp_key %in% common_sp)
    fdat <- fdat %>% filter(sp_key %in% common_sp)
    if (nrow(sdat) == 0 || nrow(fdat) == 0) next
    
    sim_comp <- function(df) {
      n_sp <- nrow(df)
      M <- matrix(NA_real_, nrow = n_sp, ncol = n_draws)
      for (i in seq_len(n_sp)) {
        M[i, ] <- draw_row(df$psi_mean[i], df$psi_sd[i], df$a[i], df$b[i], n_draws)
      }
      apply(M, 2, geo_mean)
    }
    
    comp_s <- sim_comp(sdat)
    comp_f <- sim_comp(fdat)
    
    gr_draws <- 100 * ((comp_f / comp_s)^(1 / years_n) - 1)
    
    growth_list[[length(growth_list) + 1]] <- tibble(
      Guild = gd,
      Zone  = zn,
      n_species_common = length(common_sp),
      growth_median = median(gr_draws, na.rm = TRUE),
      growth_low    = as.numeric(quantile(gr_draws, 0.025, na.rm = TRUE)),
      growth_high   = as.numeric(quantile(gr_draws, 0.975, na.rm = TRUE))
    )
  }
}

guild_composite_growth <- bind_rows(growth_list) %>%
  mutate(Zone = factor(Zone, levels = zones_keep)) %>%
  arrange(Guild, Zone)

write.csv(guild_composite_growth, "guild_composite_growth_rates_2000_2023.csv", row.names = FALSE)
cat("Saved: guild_composite_growth_rates_2000_2023.csv\n")


# ----------------------------
# 6) PLOT PREP: facets ordered by n species
#    Endpoint labels show ONLY annual growth (%/yr), no CI
# ----------------------------
dir.create(dir_out, showWarnings = FALSE)

# Facet label: "Guild (n = ...)" based on unique species in that guild
guild_counts <- psi_par %>%
  distinct(Guild, sp_key) %>%
  count(Guild, name = "Unique_Species") %>%
  mutate(Guild_Label = paste0(Guild, " (", Unique_Species, " sp.)")) %>%
  arrange(desc(Unique_Species)) %>%
  mutate(Guild_Label = factor(Guild_Label, levels = Guild_Label))  # lock order

plot_df <- guild_composite_trends %>%
  left_join(guild_counts %>% select(Guild, Guild_Label), by = "Guild") %>%
  mutate(
    Guild_Label = factor(Guild_Label, levels = levels(guild_counts$Guild_Label)),
    Zone = factor(Zone, levels = zones_keep)
  )

# ----------------------------
# ADD-ON B) Build new facet labels: "Guild (N species - X% convergence)"
# ----------------------------

guild_counts2 <- guild_counts %>%
  mutate(Guild_Label_Base = as.character(Guild_Label)) %>%
  select(-Guild_Label)

guild_counts2 <- guild_counts2 %>%
  left_join(guild_conv, by = "Guild") %>%
  mutate(
    Convergence_Pct = ifelse(is.finite(Convergence_Pct), Convergence_Pct, NA_real_),
    Guild_Label2 = paste0(
      Guild, " (", Unique_Species, " species - ",
      ifelse(is.na(Convergence_Pct), "NA", sprintf("%.0f", Convergence_Pct)),
      "% convergence)"
    )
  ) %>%
  arrange(desc(Unique_Species)) %>%
  mutate(Guild_Label2 = factor(Guild_Label2, levels = Guild_Label2))

plot_df2 <- guild_composite_trends %>%
  left_join(guild_counts2 %>% select(Guild, Guild_Label2), by = "Guild") %>%
  mutate(
    Guild_Label2 = factor(Guild_Label2, levels = levels(guild_counts2$Guild_Label2)),
    Zone = factor(Zone, levels = zones_keep)
  )


# Endpoint year per Guild × Zone (prefer end_year else latest)
guild_latest <- plot_df %>%
  filter(!is.na(comp_median)) %>%
  group_by(Guild, Guild_Label, Zone) %>%
  summarise(
    latest_year = if (any(Year == end_year)) end_year else max(Year),
    comp_latest = comp_median[Year == (if (any(Year == end_year)) end_year else max(Year))][1],
    .groups = "drop"
  )

# Per-facet y-limits from composite + CI (prevents leaders expanding panels)
facet_limits <- plot_df %>%
  group_by(Guild_Label) %>%
  summarise(
    ymin = min(c(comp_median, comp_low),  na.rm = TRUE),
    ymax = max(c(comp_median, comp_high), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    yrange = ymax - ymin,
    pad = pmax(0.005, 0.05 * yrange)
  )

# Endpoint labels = annual growth rate median (%/yr) for THAT Guild × Zone
guild_change <- guild_latest %>%
  left_join(
    guild_composite_growth %>% select(Guild, Zone, growth_median),
    by = c("Guild","Zone")
  ) %>%
  left_join(facet_limits, by = "Guild_Label") %>%
  mutate(
    label_text = paste0(sprintf("%.2f", growth_median), "%/yr"),
    label_base = comp_latest
  ) %>%
  group_by(Guild_Label) %>%
  group_modify(~disperse_labels(.x, gap = 0.08 * unique(.x$yrange))) %>%
  ungroup() %>%
  mutate(
    x_label = latest_year + 2.3,
    x_curve = latest_year + 2.5
  )

# Sanity check: should be EXACTLY 1 row per Guild_Label × Zone
dup_check <- guild_change %>%
  count(Guild_Label, Zone) %>%
  summarise(max_n = max(n), .groups = "drop") %>%
  pull(max_n)

cat("\nMax duplicates per Guild_Label×Zone in guild_change =", dup_check, "\n")
stopifnot(dup_check == 1)

# ----------------------------
# 7) PLOT: Guild composite trends by zone
# ----------------------------
p <- ggplot(plot_df, aes(x = Year, y = comp_median, color = Zone, fill = Zone)) +
  geom_ribbon(aes(ymin = comp_low, ymax = comp_high),
              alpha = 0.40, color = NA, show.legend = FALSE) +
  #geom_line(aes(y = comp_low),
  #          linewidth = 0.6, linetype = "dashed", alpha = 0.85, show.legend = FALSE) +
  #geom_line(aes(y = comp_high),
  #          linewidth = 0.6, linetype = "dashed", alpha = 0.85, show.legend = FALSE) +
  geom_line(linewidth = 1.05, na.rm = TRUE) +
  geom_point(size = 1.2, alpha = 0.55, na.rm = TRUE) +
  
  geom_vline(xintercept = end_year, linetype = "dashed",
             color = "grey40", linewidth = 0.8) +
  
  geom_point(
    data = guild_change,
    aes(x = latest_year, y = comp_latest, color = Zone),
    shape = 21, fill = "white", stroke = 1.4, size = 2.8,
    inherit.aes = FALSE
  ) +
  geom_curve(
    data = guild_change,
    aes(x = latest_year, y = comp_latest,
        xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.75, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = guild_change,
    aes(x = x_label, y = label_y_adj, label = label_text, color = Zone),
    hjust = -0.05, size = 4, fontface = "bold", show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_fill_colors) +
  guides(
    fill = "none",
    color = guide_legend(
      title = NULL,
      nrow = 1,
      override.aes = list(size = 5, linewidth = 1.4, alpha = 1)
    )
  ) +
  facet_wrap(~ Guild_Label, ncol = 3, scales = "free_y") +
  coord_cartesian(xlim = c(start_year, end_year + 13)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    legend.key.size = unit(30, "pt"),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Year",
    y = "Geometric Mean Occupancy"
  )

print(p)

ggsave(p, filename="guild_composite_trends_2000_2023.png", width = 12, height = 10, dpi = 1000)
ggsave(p, filename="guild_composite_trends_2000_2023.jpg", width = 12, height = 10, dpi = 1000)





# DONE





























# Guild plotting with converged dataset ####

############################################
# GUILD COMPOSITE TRENDS – CONVERGED SPECIES ONLY
#   (same plotting style as full guild analysis)
############################################

# ---- Convergence settings (as in Rare/Common converged script) ----
rhat_max_filtered  <- 1.10
neff_min_filtered  <- 0          # raise if you want n.eff threshold
require_all_zones_in_both_years <- TRUE

# 1) Rebuild full psi.fs.r_ table with Rhat + n.eff
psi_zone_full_conv <- bugs %>%
  dplyr::filter(stringr::str_detect(Parameter, "^psi\\.fs\\.r_")) %>%
  dplyr::mutate(
    Zone      = stringr::str_match(Parameter, psi_zone_regex)[, 2],
    YearIndex = suppressWarnings(
      as.integer(stringr::str_match(Parameter, psi_zone_regex)[, 3])
    )
  ) %>%
  dplyr::filter(!is.na(Zone), !is.na(YearIndex)) %>%
  dplyr::mutate(
    Year   = YearIndex + year_index_offset,
    Zone   = factor(Zone, levels = zones_keep),
    Species = as.character(Species),
    sp_key  = std_key(Species),
    psi     = suppressWarnings(as.numeric(mean)),
    sd      = suppressWarnings(as.numeric(sd)),
    Rhat_num = suppressWarnings(as.numeric(Rhat)),
    neff_num = suppressWarnings(as.numeric(n.eff))
  ) %>%
  dplyr::filter(Year >= start_year, Year <= end_year)

# 2) Decide which species are "converged" at BOTH endpoints in ALL zones
psi_zone_for_filter <- psi_zone_full_conv %>%
  dplyr::filter(Year %in% c(start_year, end_year),
                Zone %in% zones_keep) %>%
  dplyr::mutate(
    ok_conv = is.finite(Rhat_num) & (Rhat_num <= rhat_max_filtered) &
      is.finite(neff_num) & (neff_num >= neff_min_filtered) &
      is.finite(psi) & is.finite(sd) &
      !is.na(psi) & !is.na(sd)
  )

if (require_all_zones_in_both_years) {
  species_keep_filtered <- psi_zone_for_filter %>%
    dplyr::group_by(Species, Year) %>%
    dplyr::summarise(
      n_zones_present = dplyr::n_distinct(Zone),
      n_zones_ok      = dplyr::n_distinct(Zone[ok_conv]),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      ok_year = (n_zones_present == length(zones_keep)) &
        (n_zones_ok      == length(zones_keep))
    ) %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(ok_both_years = all(ok_year), .groups = "drop") %>%
    dplyr::filter(ok_both_years) %>%
    dplyr::pull(Species)
} else {
  species_keep_filtered <- psi_zone_for_filter %>%
    dplyr::group_by(Species, Year) %>%
    dplyr::summarise(any_ok = any(ok_conv), .groups = "drop") %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(ok_both_years = all(any_ok), .groups = "drop") %>%
    dplyr::filter(ok_both_years) %>%
    dplyr::pull(Species)
}

cat("\n==================== GUILD CONVERGED FILTER (endpoints) ====================\n")
cat("Total species (all guilds):", dplyr::n_distinct(psi_zone_full_conv$Species), "\n")
cat("Retained converged species :", length(species_keep_filtered), "\n")

# 3) Restrict to converged species only
psi_zone_full_conv <- psi_zone_full_conv %>%
  dplyr::filter(Species %in% species_keep_filtered)

# Optional: save what was actually used
write.csv(
  psi_zone_full_conv,
  "BUGS_psi_fs_r_byZone_byYear_2000_2023_GUILDS_CONVERGEDONLY_FIRSTLAST.csv",
  row.names = FALSE
)

# 4) Join Guilds and build Beta parameters (converged-only)
psi_guild_conv <- psi_zone_full_conv %>%
  dplyr::left_join(guild_lookup, by = "sp_key") %>%
  dplyr::filter(!is.na(Guild)) %>%
  dplyr::filter(!(Guild %in% remove_guilds)) %>%
  dplyr::mutate(Guild = as.character(Guild))

psi_par_conv <- psi_guild_conv %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    psi_mean = pmin(pmax(psi, eps), 1 - eps),
    psi_sd   = pmax(sd, 1e-12),
    ab       = list(beta_ab_from_mean_sd(psi_mean, psi_sd)),
    a        = ab$a,
    b        = ab$b
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-ab)

# 5) Convergence % per Guild (for facet labels, as before)
#    Use psi_guild_conv, which already has Guild attached.
psi_rhat_conv <- psi_guild_conv %>%
  dplyr::mutate(
    converged = is.finite(Rhat_num) & (Rhat_num <= rhat_max_filtered)
  )

# species-year convergence (% of zones converged for that species in that year)
species_year_conv_conv <- psi_rhat_conv %>%
  dplyr::filter(Zone %in% zones_keep) %>%
  dplyr::group_by(Guild, sp_key, Year) %>%
  dplyr::summarise(
    conv_pct_species_year = 100 * mean(converged, na.rm = TRUE),
    .groups = "drop"
  )

# average across species × years within each guild
guild_conv_conv <- species_year_conv_conv %>%
  dplyr::group_by(Guild) %>%
  dplyr::summarise(
    Convergence_Pct = mean(conv_pct_species_year, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nConvergence summary (converged species subset; Rhat <=", rhat_max_filtered, ")\n")
print(guild_conv_conv %>% dplyr::arrange(dplyr::desc(Convergence_Pct)))


# 6) COMPOSITE TRENDS (converged-only)
guilds_vec_conv <- sort(unique(as.character(psi_par_conv$Guild)))
years_vec <- start_year:end_year

composite_list_conv <- list()

set.seed(111)
for (gd in guilds_vec_conv) {
  for (zn in zones_keep) {
    sub_gz <- psi_par_conv %>% dplyr::filter(Guild == gd, Zone == zn)
    if (nrow(sub_gz) == 0) next
    
    for (yr in years_vec) {
      sub <- sub_gz %>% dplyr::filter(Year == yr)
      if (nrow(sub) == 0) next
      
      n_sp <- nrow(sub)
      draws_mat <- matrix(NA_real_, nrow = n_sp, ncol = n_draws)
      
      for (i in seq_len(n_sp)) {
        draws_mat[i, ] <- draw_row(sub$psi_mean[i], sub$psi_sd[i],
                                   sub$a[i], sub$b[i], n_draws)
      }
      
      comp_draws <- apply(draws_mat, 2, geo_mean)
      
      composite_list_conv[[length(composite_list_conv) + 1]] <- tibble(
        Guild      = gd,
        Zone       = zn,
        Year       = yr,
        n_species  = n_sp,
        comp_median = median(comp_draws, na.rm = TRUE),
        comp_low    = as.numeric(quantile(comp_draws, 0.025, na.rm = TRUE)),
        comp_high   = as.numeric(quantile(comp_draws, 0.975, na.rm = TRUE))
      )
    }
  }
}

guild_composite_trends_conv <- dplyr::bind_rows(composite_list_conv) %>%
  dplyr::mutate(Zone = factor(Zone, levels = zones_keep)) %>%
  dplyr::arrange(Guild, Zone, Year)

out_csv_trend_conv <- "guild_composite_trends_2000_2023_CONVERGEDONLY_FIRSTLAST.csv"
write.csv(guild_composite_trends_conv, out_csv_trend_conv, row.names = FALSE)
cat("\nSaved:", out_csv_trend_conv, "\n")

# 7) GROWTH RATES (%/yr, converged-only, same logic as before)
set.seed(222)
growth_list_conv <- list()

for (gd in guilds_vec_conv) {
  for (zn in zones_keep) {
    sub_gz <- psi_par_conv %>% dplyr::filter(Guild == gd, Zone == zn)
    if (nrow(sub_gz) == 0) next
    
    sdat <- sub_gz %>% dplyr::filter(Year == start_year)
    fdat <- sub_gz %>% dplyr::filter(Year == end_year)
    
    common_sp <- intersect(sdat$sp_key, fdat$sp_key)
    sdat <- sdat %>% dplyr::filter(sp_key %in% common_sp)
    fdat <- fdat %>% dplyr::filter(sp_key %in% common_sp)
    if (nrow(sdat) == 0 || nrow(fdat) == 0) next
    
    sim_comp <- function(df) {
      n_sp <- nrow(df)
      M <- matrix(NA_real_, nrow = n_sp, ncol = n_draws)
      for (i in seq_len(n_sp)) {
        M[i, ] <- draw_row(df$psi_mean[i], df$psi_sd[i],
                           df$a[i], df$b[i], n_draws)
      }
      apply(M, 2, geo_mean)
    }
    
    comp_s <- sim_comp(sdat)
    comp_f <- sim_comp(fdat)
    
    gr_draws <- 100 * ((comp_f / comp_s)^(1 / years_n) - 1)
    
    growth_list_conv[[length(growth_list_conv) + 1]] <- tibble(
      Guild = gd,
      Zone  = zn,
      n_species_common = length(common_sp),
      growth_median = median(gr_draws, na.rm = TRUE),
      growth_low    = as.numeric(quantile(gr_draws, 0.025, na.rm = TRUE)),
      growth_high   = as.numeric(quantile(gr_draws, 0.975, na.rm = TRUE))
    )
  }
}

guild_composite_growth_conv <- dplyr::bind_rows(growth_list_conv) %>%
  dplyr::mutate(Zone = factor(Zone, levels = zones_keep)) %>%
  dplyr::arrange(Guild, Zone)

out_csv_growth_conv <- "guild_composite_growth_rates_2000_2023_CONVERGEDONLY_FIRSTLAST.csv"
write.csv(guild_composite_growth_conv, out_csv_growth_conv, row.names = FALSE)
cat("Saved:", out_csv_growth_conv, "\n")

# 8) Facet labels with Unique species + convergence %
guild_counts_conv <- psi_par_conv %>%
  dplyr::distinct(Guild, sp_key) %>%
  dplyr::count(Guild, name = "Unique_Species") %>%
  dplyr::left_join(guild_conv_conv, by = "Guild") %>%
  dplyr::mutate(
    Guild_Label_Conv = paste0(
      Guild, " (", Unique_Species, " sp. - ",
      ifelse(is.na(Convergence_Pct), "NA",
             sprintf("%.0f", Convergence_Pct)),
      "% conv.)"
    )
  ) %>%
  dplyr::arrange(dplyr::desc(Unique_Species)) %>%
  dplyr::mutate(Guild_Label_Conv = factor(Guild_Label_Conv,
                                          levels = Guild_Label_Conv))

plot_df_conv <- guild_composite_trends_conv %>%
  dplyr::left_join(
    guild_counts_conv %>% dplyr::select(Guild, Guild_Label_Conv),
    by = "Guild"
  ) %>%
  dplyr::mutate(
    Guild_Label_Conv = factor(Guild_Label_Conv,
                              levels = levels(guild_counts_conv$Guild_Label_Conv)),
    Zone = factor(Zone, levels = zones_keep)
  )

# 9) Endpoint labels (using converged growth medians)
guild_latest_conv <- plot_df_conv %>%
  dplyr::filter(!is.na(comp_median)) %>%
  dplyr::group_by(Guild, Guild_Label_Conv, Zone) %>%
  dplyr::summarise(
    latest_year = if (any(Year == end_year)) end_year else max(Year),
    comp_latest = comp_median[Year == (if (any(Year == end_year)) end_year else max(Year))][1],
    .groups = "drop"
  )

facet_limits_conv <- plot_df_conv %>%
  dplyr::group_by(Guild_Label_Conv) %>%
  dplyr::summarise(
    ymin = min(c(comp_median, comp_low),  na.rm = TRUE),
    ymax = max(c(comp_median, comp_high), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    yrange = ymax - ymin,
    pad = pmax(0.005, 0.05 * yrange)
  )

guild_change_conv <- guild_latest_conv %>%
  dplyr::left_join(
    guild_composite_growth_conv %>% dplyr::select(Guild, Zone, growth_median),
    by = c("Guild","Zone")
  ) %>%
  dplyr::left_join(facet_limits_conv, by = "Guild_Label_Conv") %>%
  dplyr::mutate(
    label_text = paste0(sprintf("%.2f", growth_median), "%/yr"),
    label_base = comp_latest
  ) %>%
  dplyr::group_by(Guild_Label_Conv) %>%
  dplyr::group_modify(~disperse_labels(.x, gap = 0.08 * unique(.x$yrange))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    x_label = latest_year + 2.3,
    x_curve = latest_year + 2.5
  )

dup_check_conv <- guild_change_conv %>%
  dplyr::count(Guild_Label_Conv, Zone) %>%
  dplyr::summarise(max_n = max(n), .groups = "drop") %>%
  dplyr::pull(max_n)

cat("\nMax duplicates per Guild_Label_Conv×Zone =", dup_check_conv, "\n")
stopifnot(dup_check_conv == 1)

# 10) Plot – same style as full guild plot, but converged-only
p_conv <- ggplot(plot_df_conv,
                 aes(x = Year, y = comp_median, color = Zone, fill = Zone)) +
  geom_ribbon(aes(ymin = comp_low, ymax = comp_high),
              alpha = 0.40, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.05, na.rm = TRUE) +
  geom_point(size = 1.2, alpha = 0.55, na.rm = TRUE) +
  geom_vline(xintercept = end_year, linetype = "dashed",
             color = "grey40", linewidth = 0.8) +
  geom_point(
    data = guild_change_conv,
    aes(x = latest_year, y = comp_latest, color = Zone),
    shape = 21, fill = "white", stroke = 1.4, size = 2.8,
    inherit.aes = FALSE
  ) +
  geom_curve(
    data = guild_change_conv,
    aes(x = latest_year, y = comp_latest,
        xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.75, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = guild_change_conv,
    aes(x = x_label, y = label_y_adj, label = label_text, color = Zone),
    hjust = -0.05, size = 4, fontface = "bold", show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_fill_colors) +
  guides(
    fill = "none",
    color = guide_legend(
      title = NULL,
      nrow = 1,
      override.aes = list(size = 5, linewidth = 1.4, alpha = 1)
    )
  ) +
  facet_wrap(~ Guild_Label_Conv, ncol = 3, scales = "free_y") +
  coord_cartesian(xlim = c(start_year, end_year + 13)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Year",
    y = "Geometric Mean Occupancy",
    title = "Guild composite trends (converged species only)"
  )

print(p_conv)

safe_suffix <- "perYear"   # instead of "%perYear"

out_png_conv <- file.path(
  dir_out,
  paste0("guild_composite_trends_",
         start_year, "_", end_year,
         "_CONVERGEDONLY_endLabels_", safe_suffix, "_FIRSTLAST.png")
)

out_jpg_conv <- file.path(
  dir_out,
  paste0("guild_composite_trends_",
         start_year, "_", end_year,
         "_CONVERGEDONLY_endLabels_", safe_suffix, "_FIRSTLAST.jpg")
)


ggsave(p_conv, filename = out_png_conv, width = 12, height = 10, dpi = 1000)
ggsave(p_conv, filename = out_jpg_conv, width = 12, height = 10, dpi = 1000)

cat("\nSaved converged-only guild plot:\n  ", out_png_conv, "\n  ", out_jpg_conv, "\n", sep = "")






