############################################
# MASTER SCRIPT (unified): Rare vs Common plots + final patchwork
# Creates: g1, g2, g5, g4  -> patchworked into one column (one row each)
############################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(grid)      # unit()
  library(lubridate) # parse_date_time
  library(patchwork)
})

# ----------------------------
# 0) USER SETTINGS (shared)
# ----------------------------
setwd("C:/Users/georg/OneDrive - University of Reading/George Allen - PhD Master Folder/Year Three/Chapter 2 - Occupancy Modelling/Analysis Nov 2025/Nov_Outhwaite_Outputs")

bugs_file <- "Combined_BUGS_Data_Outhwaite_finalclustv5_32000iterations_FINAL_WITH_TAXCORR.csv"
allparks_path <- "C:/Users/georg/OneDrive - University of Reading/George Allen - PhD Master Folder/Year One/Statistics/Mapping/GBIF Sink Source AllParks/R analysis - 25.11.2024/R final/AllParks_records_FinalFiltered_Better2.csv"

start_year <- 2000
end_year   <- 2023
years_vec  <- start_year:end_year
year_index_offset <- 1999
zones_keep <- c("Core", "Buffer", "Outside")
psi_zone_regex <- "^psi\\.fs\\.r_(Core|Buffer|Outside)\\[(\\d+)\\]$"

# Colors
zone_colors <- c("Core"="#b2182b","Buffer"="black","Outside"="grey65")
zone_fill_colors <- c("Core"=alpha("#b2182b", 0.20),
                      "Buffer"=alpha("black", 0.20),
                      "Outside" = alpha("grey40",  0.30))   # darker than grey65 @ 0.20
                      
# Monte Carlo
n_draws_occ <- 999     # for g1/g2 (occupancy composites + growth)
n_draws_rich <- 999    # for g4 (AllParks richness MC)
eps <- 1e-6

# All output files
out_csv_clean   <- "BUGS_psi_fs_r_byZone_byYear_2000_2023.csv"
out_csv_trend   <- "Cooke_style_RareCommon_composite_geo_occ_median_CrI_byZoneYear.csv"
out_csv_growth  <- "Cooke_style_RareCommon_growth_rate_byZone_median_CrI.csv"

out_png_trend   <- "Cooke_style_RareCommon_composite_trends_2000_2023.png"
out_jpg_trend   <- "Cooke_style_RareCommon_composite_trends_2000_2023.jpg"
out_png_growth  <- "Cooke_style_RareCommon_growth_rate_2000_2023.png"
out_jpg_growth  <- "Cooke_style_RareCommon_growth_rate_2000_2023.jpg"

out_png_richness_trend <- "Cooke_style_RareCommon_richness_trends_2000_2023_NO_MC.png"
out_jpg_richness_trend <- "Cooke_style_RareCommon_richness_trends_2000_2023_NO_MC.jpg"

out_png_rich_mc <- "richness_by_zone_RareCommon.png"

out_png_patch <- "rarevscommon_patchwork.png"
out_jpg_patch <- "rarevscommon_patchwork.jpg"

# ----------------------------
# 1) LOAD BUGS ONCE
# ----------------------------
bugs <- read.csv(bugs_file, stringsAsFactors = FALSE)

# ----------------------------
# 2) Shared helpers
# ----------------------------
disperse_labels <- function(df, gap = 0.06) {
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
    if (y[1] < low) y <- seq(low, high, length.out = n)
  }
  
  df$label_y_adj <- NA_real_
  df$label_y_adj[ord] <- y
  df
}

beta_draws_from_mean_sd <- function(m, s, n = 1000) {
  m <- pmin(pmax(m, 1e-6), 1 - 1e-6)
  s <- pmax(s, 1e-12)
  v <- s^2
  vmax <- m * (1 - m)
  
  if (!is.finite(v) || v <= 0 || !is.finite(vmax) || vmax <= 0) v <- 0.01 * vmax
  if (v >= vmax) v <- 0.95 * vmax
  
  k <- (m * (1 - m) / v) - 1
  a <- m * k
  b <- (1 - m) * k
  
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
    v <- 0.5 * vmax
    k <- (m * (1 - m) / v) - 1
    a <- max(m * k, 1e-3)
    b <- max((1 - m) * k, 1e-3)
  }
  
  rbeta(n, shape1 = a, shape2 = b)
}

std_key <- function(x) str_squish(str_to_lower(as.character(x)))

hdi_vec <- function(x, credMass = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) return(c(low = NA_real_, high = NA_real_))
  x <- sort(x)
  m <- floor(credMass * n)
  if (m < 1) return(c(low = NA_real_, high = NA_real_))
  lows  <- x[1:(n - m)]
  highs <- x[(1 + m):n]
  widths <- highs - lows
  j <- which.min(widths)
  c(low = lows[j], high = highs[j])
}

fit_beta_from_ci <- function(l, u) {
  l <- pmin(1 - eps, pmax(eps, l))
  u <- pmin(1 - eps, pmax(eps, u))
  if (!is.finite(l) || !is.finite(u) || l >= u) return(c(NA_real_, NA_real_))
  
  obj <- function(par) {
    a <- exp(par[1]); b <- exp(par[2])
    ql <- qbeta(0.025, a, b)
    qu <- qbeta(0.975, a, b)
    (ql - l)^2 + (qu - u)^2
  }
  
  opt <- try(optim(log(c(2,2)), obj, method = "Nelder-Mead"), silent = TRUE)
  if (inherits(opt, "try-error")) return(c(NA_real_, NA_real_))
  
  a <- exp(opt$par[1]); b <- exp(opt$par[2])
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(c(NA_real_, NA_real_))
  c(a, b)
}

beta_from_mean_sd <- function(mu, sd) {
  mu <- pmin(1 - eps, pmax(eps, mu))
  if (!is.finite(mu) || !is.finite(sd) || sd <= 0) return(c(NA_real_, NA_real_))
  v <- sd^2
  k <- mu*(1-mu)/v - 1
  if (!is.finite(k) || k <= 0) return(c(NA_real_, NA_real_))
  c(mu*k, (1-mu)*k)
}


# ---- Helper: median + 95% HDI for a vector ----
median_hdi <- function(x, credMass = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(tibble(
      median  = NA_real_,
      hdi_low = NA_real_,
      hdi_high = NA_real_,
      n = length(x)
    ))
  }
  h <- HDInterval::hdi(x, credMass = credMass)
  tibble(
    median  = median(x),
    hdi_low = unname(h[1]),
    hdi_high = unname(h[2]),
    n = length(x)
  )
}

# ---- Helper: group summaries (median + 95% HDI) ----
group_summaries_hdi <- function(df, group_cols, value_col, credMass = 0.95) {
  df %>%
    dplyr::group_by(dplyr::across(all_of(group_cols))) %>%
    dplyr::summarise(
      tmp = list(median_hdi(.data[[value_col]], credMass = credMass)),
      .groups = "drop"
    ) %>%
    tidyr::unnest(tmp)
}


# ----------------------------
# 3) Extract psi.fs.r_* once (used by g1/g2 + g5)
# ----------------------------
psi_zone_full <- bugs %>%
  filter(str_detect(Parameter, "^psi\\.fs\\.r_")) %>%
  mutate(
    Zone = str_match(Parameter, psi_zone_regex)[, 2],
    YearIndex = suppressWarnings(as.integer(str_match(Parameter, psi_zone_regex)[, 3]))
  ) %>%
  filter(!is.na(Zone), !is.na(YearIndex)) %>%
  mutate(
    Year = YearIndex + year_index_offset,
    Zone = factor(Zone, levels = zones_keep),
    psi  = suppressWarnings(as.numeric(mean)),
    sd   = suppressWarnings(as.numeric(sd))
  ) %>%
  filter(Year >= start_year, Year <= end_year) %>%
  select(Species, Zone, Year, YearIndex, Parameter,
         psi, sd, X2.5., X25., X50., X75., X97.5., Rhat, n.eff)

write.csv(psi_zone_full, out_csv_clean, row.names = FALSE)

# ----------------------------
# 4) g1 + g2: Cooke-style multispecies trends (Rare vs Common; MC; geom mean; growth)
# ----------------------------
set.seed(1)

# Rare/Common definition (BUGS; Cooke-style):
# - For each species: compute annual occupancy (averaged across zones), then median across years
# - Rare  = lower quartile (Q1) of species medians
# - Common = upper quartile (Q3) of species medians

species_year_occ <- psi_zone_full %>%
  filter(is.finite(psi)) %>%
  group_by(Species, Year) %>%
  summarise(psi_year = mean(psi, na.rm = TRUE), .groups = "drop")  # avg across zones within year

species_median_occ <- species_year_occ %>%
  group_by(Species) %>%
  summarise(median_psi = median(psi_year, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(median_psi))

q1_psi <- quantile(species_median_occ$median_psi, 0.25, na.rm = TRUE)
q3_psi <- quantile(species_median_occ$median_psi, 0.75, na.rm = TRUE)

rarity_lookup <- species_median_occ %>%
  mutate(
    Rarity = case_when(
      median_psi <= q1_psi ~ "Rare",
      median_psi >= q3_psi ~ "Common",
      TRUE ~ NA_character_
    ),
    Rarity = factor(Rarity, levels = c("Rare", "Common"))
  ) %>%
  filter(!is.na(Rarity)) %>%
  select(Species, Rarity)

psi_zone_rc <- psi_zone_full %>%
  left_join(rarity_lookup, by = "Species") %>%
  filter(!is.na(Rarity)) %>%
  filter(!is.na(psi), !is.na(sd), is.finite(psi), is.finite(sd)) %>%
  mutate(
    psi_mean = pmin(pmax(psi, 1e-6), 1 - 1e-6),
    psi_sd   = pmax(sd, 1e-12)
  )



# Monte Carlo draws per species×zone×year
draws_long <- psi_zone_rc %>%
  rowwise() %>%
  mutate(draw = list(beta_draws_from_mean_sd(psi_mean, psi_sd, n = n_draws_occ))) %>%
  ungroup() %>%
  select(Rarity, Species, Zone, Year, draw) %>%
  unnest_longer(draw, values_to = "psi_draw") %>%
  group_by(Rarity, Species, Zone, Year) %>%
  mutate(draw_id = row_number()) %>%
  ungroup() %>%
  mutate(psi_draw = pmin(pmax(psi_draw, 1e-12), 1 - 1e-12))

# Composite geometric mean per draw -> median + 95% CrI
zone_year_draw_gm <- draws_long %>%
  group_by(Rarity, Zone, Year, draw_id) %>%
  summarise(gm = exp(mean(log(psi_draw))), .groups = "drop")

zone_trend <- zone_year_draw_gm %>%
  group_by(Rarity, Zone, Year) %>%
  summarise(
    gm_median = median(gm),
    gm_low    = as.numeric(quantile(gm, 0.025)),
    gm_high   = as.numeric(quantile(gm, 0.975)),
    .groups = "drop"
  ) %>%
  arrange(Rarity, Zone, Year)

write.csv(zone_trend, out_csv_trend, row.names = FALSE)

tapply(rarity_lookup$Species, list(rarity_lookup$Rarity), length) # number of rare and common

# 1) Species-level mean, min, max occupancy (across zones × years)
species_occ_summary <- psi_zone_full %>%
  filter(is.finite(psi)) %>%
  group_by(Species) %>%
  summarise(
    mean_psi = mean(psi, na.rm = TRUE),
    min_psi  = min(psi, na.rm = TRUE),
    max_psi  = max(psi, na.rm = TRUE),
    .groups = "drop"
  )

# 2) Add rarity class (Rare / Common)
species_occ_rc <- species_occ_summary %>%
  left_join(rarity_lookup, by = "Species") %>%
  filter(!is.na(Rarity))

# 3) Summarise within each rarity group
rarity_occ_stats <- species_occ_rc %>%
  group_by(Rarity) %>%
  summarise(
    mean_of_means = mean(mean_psi),
    min_of_mins   = min(min_psi),
    max_of_maxs   = max(max_psi),
    .groups = "drop"
  )

rarity_occ_stats

q1_psi <- quantile(species_median_occ$median_psi, 0.25, na.rm = TRUE)
q3_psi <- quantile(species_median_occ$median_psi, 0.75, na.rm = TRUE)

q1_psi
q3_psi

range_Q1_km2 <- q1_psi * 20000
range_Q3_km2 <- q3_psi * 20000

cat("Range-size threshold for RARE species (Q1):  ", round(range_Q1_km2, 1), " km²\n")
cat("Range-size threshold for COMMON species (Q3):", round(range_Q3_km2, 1), " km²\n")





# Converting to correct spatial coverage (km^2)

realised_terrestrial_km2 <- 41642.908159  # from your QGIS zone sums

range_Q1_km2 <- q1_psi * realised_terrestrial_km2
range_Q1_km2
range_Q3_km2 <- q3_psi * realised_terrestrial_km2
range_Q3_km2


cell_area_km2 <- 4
realised_terrestrial_km2 <- 41642.908159

n_cells_terrestrial_equiv <- realised_terrestrial_km2 / cell_area_km2
n_cells_terrestrial_equiv

range_Q1_cells <- q1_psi * n_cells_terrestrial_equiv
range_Q3_cells <- q3_psi * n_cells_terrestrial_equiv
range_Q1_cells
range_Q3_cells








# Growth rate (%/yr) from start->end on draws
y_years <- end_year - start_year

zone_start_end_draws <- zone_year_draw_gm %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Rarity, Zone, Year, draw_id, gm) %>%
  pivot_wider(names_from = Year, values_from = gm, names_prefix = "Year_") %>%
  filter(!is.na(.data[[paste0("Year_", start_year)]]),
         !is.na(.data[[paste0("Year_", end_year)]])) %>%
  mutate(
    s = pmax(.data[[paste0("Year_", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Year_", end_year)]], 1e-12),
    growth_rate_pct = 100 * ((f / s)^(1 / y_years) - 1)
  )

growth_summary <- zone_start_end_draws %>%
  group_by(Rarity, Zone) %>%
  summarise(
    growth_median = median(growth_rate_pct),
    growth_low    = as.numeric(quantile(growth_rate_pct, 0.025)),
    growth_high   = as.numeric(quantile(growth_rate_pct, 0.975)),
    .groups = "drop"
  ) %>%
  arrange(Rarity, Zone)

write.csv(growth_summary, out_csv_growth, row.names = FALSE)

# Endpoint label placement (show ONLY growth_median)
facet_limits_occ <- zone_trend %>%
  group_by(Rarity) %>%
  summarise(
    ymin = min(gm_low,  na.rm = TRUE),
    ymax = max(gm_high, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    yrange = ymax - ymin,
    pad = pmax(0.005, 0.05 * yrange)
  )

zone_end_occ <- zone_trend %>%
  filter(Year == end_year) %>%
  left_join(growth_summary %>% 
              select(Rarity, Zone, growth_median),
            by = c("Rarity","Zone")) %>%
  left_join(facet_limits_occ, by = "Rarity") %>%
  mutate(
    label_text = paste0(sprintf("%.2f", growth_median), "%/yr"),
    off = pmax(0.01, 0.18 * yrange),
    label_base = case_when(
      Zone == "Core"    ~ gm_median - off,
      Zone == "Buffer"  ~ gm_median,
      Zone == "Outside" ~ gm_median + off,
      TRUE              ~ gm_median
    )
  ) %>%
  group_by(Rarity) %>%
  group_modify(~disperse_labels(.x, gap = 0.10 * unique(.x$yrange))) %>%
  ungroup() %>%
  mutate(
    x_label = end_year + 3,
    x_curve = end_year + 2.5
  )



# g1
g1 <- ggplot(zone_trend, aes(x = Year, y = gm_median, color = Zone)) +
  geom_ribbon(aes(ymin = gm_low, ymax = gm_high, fill = Zone),
              alpha = 0.40, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = end_year, linetype = "dashed",
             color = "grey30", linewidth = 0.8) +
  geom_point(
    data = zone_end_occ,
    aes(x = end_year, y = gm_median, color = Zone),
    shape = 21, fill = "white", stroke = 1.5, size = 3.2,
    inherit.aes = FALSE
  ) +
  geom_curve(
    data = zone_end_occ,
    aes(x = end_year, y = gm_median,
        xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.7, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = zone_end_occ,
    aes(x = x_label, y = label_y_adj,
        label = label_text, color = Zone),
    hjust = 0, size = 4.2, fontface = "bold",
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  facet_wrap(~Rarity, nrow = 1, ncol = 2, scales = "free_y") +
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_fill_colors) +
  theme_minimal(base_size = 12) +
  labs(x = "Year", y = "Geometric mean occupancy", color = "Zone") +
  theme(
    legend.position = "right",
    legend.key.size = unit(2.2, "lines"),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  ) +
  coord_cartesian(xlim = c(start_year, end_year + 9))

g1

ggsave(g1, filename = out_png_trend, width = 10, height = 3, dpi = 1000)
ggsave(g1, filename = out_jpg_trend, width = 10, height = 3, dpi = 1000)

# g2
g2 <- ggplot(growth_summary, aes(x = Zone, y = growth_median, color = Zone)) +
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = "dashed", color = "grey30") +
  geom_errorbar(aes(ymin = growth_low, ymax = growth_high),
                width = 0.25, linewidth = 0.9, show.legend = FALSE) +
  geom_point(size = 4, stroke = 1.2, show.legend = FALSE) +
  facet_wrap(~Rarity, nrow = 1) +
  scale_color_manual(values = zone_colors) +
  theme_minimal(base_size = 12) +
  labs(x = "Zone", y = "Annual growth rate (% per year)") +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text  = element_text(size = 12)
  )

ggsave(g2, filename = out_png_growth, width = 10, height = 3, dpi = 1000)
ggsave(g2, filename = out_jpg_growth, width = 10, height = 3, dpi = 1000)

# ----------------------------
# 5) g5: Cooke-style richness trends (NO MC) (Rare/Common via median mean psi)
# ----------------------------
psi_zone_rich <- psi_zone_full %>%
  mutate(psi = suppressWarnings(as.numeric(psi))) %>%
  filter(is.finite(psi)) %>%
  select(Species, Zone, Year, psi) %>%
  left_join(rarity_lookup, by = "Species") %>%
  filter(!is.na(Rarity))

richness_trend <- psi_zone_rich %>%
  group_by(Rarity, Zone, Year) %>%
  summarise(
    richness = sum(psi, na.rm = TRUE),
    n_species = n_distinct(Species),
    .groups = "drop"
  ) %>%
  arrange(Rarity, Zone, Year)

# ----------------------------
# 5) g5: Cooke-style richness trends (NO MC) (Rare/Common via median mean psi)
# ----------------------------
psi_zone_rich <- psi_zone_full %>%
  mutate(psi = suppressWarnings(as.numeric(psi))) %>%
  filter(is.finite(psi)) %>%
  select(Species, Zone, Year, psi) %>%
  left_join(rarity_lookup, by = "Species") %>%
  filter(!is.na(Rarity))

richness_trend <- psi_zone_rich %>%
  group_by(Rarity, Zone, Year) %>%
  summarise(
    richness  = sum(psi, na.rm = TRUE),
    n_species = n_distinct(Species),
    .groups   = "drop"
  ) %>%
  arrange(Rarity, Zone, Year)

# ---- BOOTSTRAP % CHANGE 2000–2023 (this will feed BOTH g5 labels and HDI Excel) ----
set.seed(777)
B <- 999

make_mat <- function(df_group) {
  mat <- df_group %>%
    select(Species, Year, psi) %>%
    pivot_wider(names_from = Year, values_from = psi) %>%
    arrange(Species)
  
  yrs <- sort(unique(df_group$Year))
  M   <- as.matrix(mat[, as.character(yrs), drop = FALSE])
  rownames(M) <- mat$Species
  colnames(M) <- as.character(yrs)
  list(M = M, yrs = yrs)
}

groups_g5 <- psi_zone_rich %>%
  distinct(Rarity, Zone) %>%
  arrange(Rarity, Zone)

g5_boot_list <- vector("list", nrow(groups_g5))

for (k in seq_len(nrow(groups_g5))) {
  rr <- groups_g5$Rarity[k]
  zn <- groups_g5$Zone[k]
  
  df_g <- psi_zone_rich %>% filter(Rarity == rr, Zone == zn)
  obj  <- make_mat(df_g)
  M    <- obj$M
  
  nsp <- nrow(M)
  idx <- replicate(B, sample.int(nsp, size = nsp, replace = TRUE))
  
  sums <- apply(idx, 2, function(ii) colSums(M[ii, , drop = FALSE], na.rm = TRUE))
  
  s <- sums[as.character(start_year), ]
  f <- sums[as.character(end_year),   ]
  pct_change <- 100 * (f - s) / pmax(s, 1e-12)
  
  g5_boot_list[[k]] <- tibble(
    Rarity  = rr,
    Zone    = zn,
    boot_id = seq_len(B),
    value   = pct_change
  )
}

g5_change_boot <- bind_rows(g5_boot_list)

g5_groups <- group_summaries_hdi(g5_change_boot, c("Rarity", "Zone"), "value")
# g5_groups now has: Rarity, Zone, median (pct change), hdi_low, hdi_high, n


facet_limits_rich <- richness_trend %>%
  group_by(Rarity) %>%
  summarise(
    ymin = min(richness, na.rm = TRUE),
    ymax = max(richness, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    yrange = ymax - ymin,
    pad = pmax(0.5, 0.05 * yrange)
  )

# g5_groups: Rarity, Zone, median (= pct_change), hdi_low, hdi_high

zone_end_rich <- richness_trend %>%
  filter(Year == end_year) %>%
  left_join(
    g5_groups %>%
      rename(pct_change_median = median),
    by = c("Rarity", "Zone")
  ) %>%
  left_join(facet_limits_rich, by = "Rarity") %>%
  mutate(
    label_text = paste0(sprintf("%.1f", pct_change_median), "%"),
    off = pmax(0.5, 0.18 * yrange),
    label_base = case_when(
      Zone == "Core"    ~ richness - off,
      Zone == "Buffer"  ~ richness,
      Zone == "Outside" ~ richness + off,
      TRUE              ~ richness
    )
  ) %>%
  group_by(Rarity) %>%
  group_modify(~disperse_labels(.x, gap = 0.10 * unique(.x$yrange))) %>%
  ungroup() %>%
  mutate(
    x_label = end_year + 3,
    x_curve = end_year + 2.5
  )


g5 <- ggplot(richness_trend, aes(x = Year, y = richness, color = Zone)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  geom_vline(xintercept = end_year, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  geom_point(
    data = zone_end_rich,
    aes(x = end_year, y = richness, color = Zone),
    shape = 21, fill = "white", stroke = 1.4, size = 3.0,
    inherit.aes = FALSE
  ) +
  geom_curve(
    data = zone_end_rich,
    aes(x = end_year, y = richness, xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.7, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = zone_end_rich,
    aes(x = x_label, y = label_y_adj, label = label_text, color = Zone),
    hjust = 0, size = 4.2, fontface = "bold",
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  facet_wrap(~Rarity, scales = "free_y", nrow = 1) +
  scale_color_manual(values = zone_colors) +
  coord_cartesian(xlim = c(start_year, end_year + 9)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.key.size = unit(2.2, "lines"),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Year", y = "Species richness", color = "Zone")

ggsave(g5, filename = out_png_richness_trend, width = 10, height = 3, dpi = 1000)
ggsave(g5, filename = out_jpg_richness_trend, width = 10, height = 3, dpi = 1000)

# ----------------------------
# 6) g4: Cooke et al. 2023 - Rare/Common richness quartiles (BUGS definition) (MC cloud)
#     Rare/Common = lower/upper quartile of species' median occupancy across years (from BUGS)
#     (Optional filter retained): only species with >= min_records in AllParks
# ----------------------------
set.seed(123)

min_records <- 50

# A) Build AllParks species filter (>= min_records) -- no rarity from AllParks
rec_raw <- read.csv(allparks_path, stringsAsFactors = FALSE)

records0 <- rec_raw %>%
  select(species, eventDate, layer_grid) %>%
  distinct() %>%
  filter(!is.na(species), species != "", !is.na(layer_grid), layer_grid != "") %>%
  mutate(
    eventDate = suppressWarnings(parse_date_time(
      eventDate,
      orders = c("dmy", "ymd", "Ymd HMS", "Ymd", "mdy", "dmy HMS", "ymd HMS")
    )),
    year = as.integer(format(eventDate, "%Y"))
  ) %>%
  filter(!is.na(year), year >= start_year, year <= end_year)

sp_counts <- records0 %>%
  count(species, name = "n_records") %>%
  filter(n_records >= min_records) %>%
  mutate(sp_key = std_key(species))

# B) Extract BUGS psi summaries (mean + CI) and join BUGS rarity (quartiles)
psi_zone_ci <- bugs %>%
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
    l95 = suppressWarnings(as.numeric(X2.5.)),
    u95 = suppressWarnings(as.numeric(X97.5.))
  ) %>%
  filter(Year >= start_year, Year <= end_year) %>%
  filter(is.finite(psi_mean), is.finite(l95), is.finite(u95)) %>%
  select(Species, sp_key, Zone, Year, psi_mean, l95, u95)

psi_labeled <- psi_zone_ci %>%
  inner_join(sp_counts %>% select(sp_key), by = "sp_key") %>%   # >= min_records filter retained
  inner_join(rarity_lookup %>%
               mutate(sp_key = std_key(Species)) %>%
               select(sp_key, rarity = Rarity),
             by = "sp_key") %>%                                  # BUGS quartile rarity
  mutate(
    Zone = factor(Zone, levels = zones_keep),
    rarity = factor(rarity, levels = c("Rare","Common"))
  )

# C) Monte Carlo richness per rarity × zone × year
rich_draws_list <- list()

for (rr in levels(psi_labeled$rarity)) {
  for (zn in zones_keep) {
    for (yr in years_vec) {
      
      sub <- psi_labeled %>% filter(rarity == rr, Zone == zn, Year == yr)
      if (nrow(sub) == 0) next
      
      richness_draws <- numeric(n_draws_rich)
      
      for (i in seq_len(nrow(sub))) {
        l <- sub$l95[i]; u <- sub$u95[i]
        ab <- fit_beta_from_ci(l, u)
        
        if (anyNA(ab)) {
          sd_approx <- (u - l) / (2 * 1.96)
          ab <- beta_from_mean_sd(sub$psi_mean[i], sd_approx)
        }
        
        if (anyNA(ab)) {
          draws <- rnorm(n_draws_rich, mean = sub$psi_mean[i], sd = 0.02)
          draws <- pmin(1 - eps, pmax(eps, draws))
        } else {
          draws <- rbeta(n_draws_rich, ab[1], ab[2])
        }
        
        richness_draws <- richness_draws + draws
      }
      
      rich_draws_list[[length(rich_draws_list) + 1]] <- tibble(
        rarity = rr,
        Zone = zn,
        Year = yr,
        draw = seq_len(n_draws_rich),
        richness = richness_draws
      )
    }
  }
}

richness_draws_year <- bind_rows(rich_draws_list) %>%
  mutate(
    rarity = factor(rarity, levels = c("Rare","Common")),
    Zone = factor(Zone, levels = zones_keep)
  )

# Period mean richness across years per rarity × zone × draw
richness_period_draws <- richness_draws_year %>%
  group_by(rarity, Zone, draw) %>%
  summarise(richness_mean = mean(richness, na.rm = TRUE), .groups = "drop")

richness_zone_summary <- richness_period_draws %>%
  group_by(rarity, Zone) %>%
  summarise(
    median = median(richness_mean, na.rm = TRUE),
    hdi_low  = hdi_vec(richness_mean, 0.95)["low"],
    hdi_high = hdi_vec(richness_mean, 0.95)["high"],
    .groups = "drop"
  )

# g4 (same plotting style)
g4 <- ggplot(richness_period_draws, aes(x = Zone, y = richness_mean, color = Zone)) +
  geom_jitter(width = 0.4, alpha = 0.20, size = 1.5, show.legend = FALSE) +
  geom_hline(
    data = richness_zone_summary,
    aes(yintercept = median, color = Zone),
    linewidth = 0.7,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  facet_wrap(~rarity, nrow = 1, scales = "free_y") +
  scale_color_manual(values = zone_colors) +
  theme_minimal(base_size = 14) +
  labs(x = NULL, y = "Species richness (∑ years)") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  ) +
  coord_cartesian(clip = "off")

ggsave(out_png_rich_mc, g4, width = 10, height = 3, dpi = 1000)




# ----------------------------
# 7) FINAL PATCHWORK (one column, one plot per row; numbered; no legends)
# ----------------------------
g_all <- (g1 / g5 / g4 / g2) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "none")

ggsave(g_all, filename = out_png_patch, width = 10, height = 12, dpi = 1000)
ggsave(g_all, filename = out_jpg_patch, width = 10, height = 12, dpi = 1000)

print(g_all)











# Absolute effect sizes ####

############################################################
# ADD-ON: Absolute effect sizes (median Δ + 95% credible interval)
# Method: HDI via HDInterval::hdi (Meredith & Kruschke, 2018)
# Outputs Excel: one sheet per plot (g1, g2, g5, g4)
############################################################

# ---- Packages ----
if (!requireNamespace("HDInterval", quietly = TRUE)) install.packages("HDInterval")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

library(HDInterval)
library(openxlsx)
library(dplyr)
library(tidyr)
library(purrr)

# ---- Helper: median + 95% HDI for a vector ----
median_hdi <- function(x, credMass = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(tibble(median = NA_real_, hdi_low = NA_real_, hdi_high = NA_real_, n = length(x)))
  }
  h <- HDInterval::hdi(x, credMass = credMass)
  tibble(
    median = median(x),
    hdi_low = unname(h[1]),
    hdi_high = unname(h[2]),
    n = length(x)
  )
}

# ---- Helper: group summaries (median + 95% HDI) ----
group_summaries_hdi <- function(df, group_cols, value_col, credMass = 0.95) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      tmp = list(median_hdi(.data[[value_col]], credMass = credMass)),
      .groups = "drop"
    ) %>%
    tidyr::unnest(tmp)
}

# ---- Helper: absolute pairwise differences (A - B) with median + 95% HDI ----
# Expects df has ONE row per draw/replicate per group (e.g. draw_id, boot, draw)
pairwise_abs_effects <- function(df, strata_cols, treat_col, draw_col, value_col,
                                 treat_levels = NULL, credMass = 0.95) {
  
  d <- df %>%
    filter(is.finite(.data[[value_col]])) %>%
    mutate(.treat = as.character(.data[[treat_col]]))
  
  if (!is.null(treat_levels)) d <- d %>% mutate(.treat = factor(.treat, levels = treat_levels))
  lvl <- levels(factor(d$.treat))
  comps <- combn(lvl, 2, simplify = FALSE)
  
  # Wide per stratum + draw so we can subtract draw-wise
  wide <- d %>%
    select(all_of(strata_cols), !!sym(draw_col), .treat, !!sym(value_col)) %>%
    pivot_wider(names_from = .treat, values_from = !!sym(value_col))
  
  strata_keys <- wide %>% distinct(across(all_of(strata_cols)))
  
  out <- list()
  
  for (i in seq_len(nrow(strata_keys))) {
    key <- strata_keys[i, , drop = FALSE]
    subw <- wide %>% semi_join(key, by = strata_cols)
    
    for (cc in comps) {
      a <- as.character(cc[1]); b <- as.character(cc[2])
      
      if (!(a %in% names(subw)) || !(b %in% names(subw))) next
      
      delta <- subw[[a]] - subw[[b]]
      summ <- median_hdi(delta, credMass = credMass) %>%
        mutate(comparison = paste0(a, " - ", b))
      
      out[[length(out) + 1]] <- bind_cols(key, summ)
    }
  }
  
  bind_rows(out) %>%
    relocate(all_of(strata_cols), comparison)
}

# ==========================================================
# 1) g1 absolute effects
# Data: zone_year_draw_gm -> gm draws per Rarity×Zone×Year×draw_id
# For absolute effects in your paper:
#   - either use period-mean richness style (if g1 is "richness by zone"),
#   - OR if g1 is a time trend, use change 2000->2023 per draw as your effect variable.
#
# Below I do BOTH, so you can choose:
#   (A) Period-mean across years (draw-wise)
#   (B) Change 2000->2023 (draw-wise)
# ==========================================================

# (A) Period-mean (mean across years) per draw
g1_period_draw <- zone_year_draw_gm %>%
  group_by(Rarity, Zone, draw_id) %>%
  summarise(value = mean(gm, na.rm = TRUE), .groups = "drop")

g1_period_groups <- group_summaries_hdi(g1_period_draw, c("Rarity", "Zone"), "value")

g1_period_zone_abs <- pairwise_abs_effects(
  df = g1_period_draw,
  strata_cols = c("Rarity"),
  treat_col = "Zone",
  draw_col  = "draw_id",
  value_col = "value",
  treat_levels = zones_keep
)

g1_period_rarity_abs <- pairwise_abs_effects(
  df = g1_period_draw,
  strata_cols = c("Zone"),
  treat_col = "Rarity",
  draw_col  = "draw_id",
  value_col = "value",
  treat_levels = c("Rare", "Common")
)

# (B) Change 2000 -> 2023 per draw (percent change)
g1_change_draws <- zone_year_draw_gm %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Rarity, Zone, Year, draw_id, gm) %>%
  pivot_wider(names_from = Year, values_from = gm, names_prefix = "Y") %>%
  mutate(
    s = pmax(.data[[paste0("Y", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Y", end_year)]],   1e-12),
    pct_change = 100 * (f - s) / s
  ) %>%
  select(Rarity, Zone, draw_id, pct_change)

g1_change_groups <- group_summaries_hdi(g1_change_draws, c("Rarity", "Zone"), "pct_change")

g1_change_zone_abs <- pairwise_abs_effects(
  df = g1_change_draws,
  strata_cols = c("Rarity"),
  treat_col = "Zone",
  draw_col  = "draw_id",
  value_col = "pct_change",
  treat_levels = zones_keep
)

g1_change_rarity_abs <- pairwise_abs_effects(
  df = g1_change_draws,
  strata_cols = c("Zone"),
  treat_col = "Rarity",
  draw_col  = "draw_id",
  value_col = "pct_change",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# 2) g2 absolute effects (growth rate %/yr)
# Data: zone_start_end_draws has growth_rate_pct per draw_id
# ==========================================================
g2_draws <- zone_start_end_draws %>%
  select(Rarity, Zone, draw_id, growth_rate_pct) %>%
  rename(value = growth_rate_pct)

g2_groups <- group_summaries_hdi(g2_draws, c("Rarity", "Zone"), "value")

g2_zone_abs <- pairwise_abs_effects(
  df = g2_draws,
  strata_cols = c("Rarity"),
  treat_col = "Zone",
  draw_col  = "draw_id",
  value_col = "value",
  treat_levels = zones_keep
)

g2_rarity_abs <- pairwise_abs_effects(
  df = g2_draws,
  strata_cols = c("Zone"),
  treat_col = "Rarity",
  draw_col  = "draw_id",
  value_col = "value",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# 3) g5 absolute effects (reuse g5_change_boot + g5_groups from section 5)
# ==========================================================
g5_zone_abs <- pairwise_abs_effects(
  df          = g5_change_boot,
  strata_cols = c("Rarity"),
  treat_col   = "Zone",
  draw_col    = "boot_id",
  value_col   = "value",
  treat_levels = zones_keep
)

g5_rarity_abs <- pairwise_abs_effects(
  df          = g5_change_boot,
  strata_cols = c("Zone"),
  treat_col   = "Rarity",
  draw_col    = "boot_id",
  value_col   = "value",
  treat_levels = c("Rare", "Common")
)


# ==========================================================
# 4) g4 absolute effects (AllParks richness, period-mean across years)
# Data: richness_period_draws already has richness_mean per draw
# ==========================================================
g4_draws <- richness_period_draws %>%
  rename(Rarity = rarity, draw_id = draw) %>%
  select(Rarity, Zone, draw_id, richness_mean) %>%
  rename(value = richness_mean)

g4_groups <- group_summaries_hdi(g4_draws, c("Rarity", "Zone"), "value")

g4_zone_abs <- pairwise_abs_effects(
  df = g4_draws,
  strata_cols = c("Rarity"),
  treat_col = "Zone",
  draw_col  = "draw_id",
  value_col = "value",
  treat_levels = zones_keep
)

g4_rarity_abs <- pairwise_abs_effects(
  df = g4_draws,
  strata_cols = c("Zone"),
  treat_col = "Rarity",
  draw_col  = "draw_id",
  value_col = "value",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# 5) EXPORT: Excel workbook
# One sheet per plot/metric, placed after your Hedge's g workbook
# ==========================================================
wb_abs <- createWorkbook()

write_abs_sheet <- function(wb, sheet, groups_df, zone_abs_df, rarity_abs_df) {
  addWorksheet(wb, sheet)
  
  writeData(wb, sheet, "Group summaries (median + 95% HDI)", startRow = 1, startCol = 1)
  writeData(wb, sheet, groups_df, startRow = 2, startCol = 1)
  
  r2 <- 4 + nrow(groups_df)
  writeData(wb, sheet, "Absolute effect sizes: Zone contrasts within each Rarity (A - B; median + 95% HDI)",
            startRow = r2, startCol = 1)
  writeData(wb, sheet, zone_abs_df, startRow = r2 + 1, startCol = 1)
  
  r3 <- r2 + 3 + nrow(zone_abs_df)
  writeData(wb, sheet, "Absolute effect sizes: Rare vs Common within each Zone (Rare - Common; median + 95% HDI)",
            startRow = r3, startCol = 1)
  writeData(wb, sheet, rarity_abs_df, startRow = r3 + 1, startCol = 1)
  
  setColWidths(wb, sheet, cols = 1:40, widths = "auto")
}

# g1: two sheets so you can choose what you report
write_abs_sheet(wb_abs, "g1_period_mean",   g1_period_groups, g1_period_zone_abs, g1_period_rarity_abs)
write_abs_sheet(wb_abs, "g1_change_00_23",  g1_change_groups, g1_change_zone_abs,  g1_change_rarity_abs)

# g2, g5, g4
write_abs_sheet(wb_abs, "g2_growth_rate",   g2_groups, g2_zone_abs, g2_rarity_abs)
write_abs_sheet(wb_abs, "g5_change_00_23",  g5_groups, g5_zone_abs, g5_rarity_abs)
write_abs_sheet(wb_abs, "g4_period_rich",   g4_groups, g4_zone_abs, g4_rarity_abs)

excel_out_abs <- "2_AbsoluteEffectSizes_MedianDiff_HDI_byPlot.xlsx"
saveWorkbook(wb_abs, excel_out_abs, overwrite = TRUE)

cat("\nSaved Excel workbook (absolute effects):\n  ", excel_out_abs, "\n", sep = "")


















# Hedge's G CIs and relative effect sizes ####

############################################################
# ADD-ON: Hedge's g (95% CI) treatment differences per plot
# Exports to Excel: one sheet per plot (g1, g2, g5, g4)
############################################################

# ---- Packages ----
if (!requireNamespace("effsize", quietly = TRUE)) install.packages("effsize")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")

library(effsize)
library(openxlsx)
library(purrr)
library(dplyr)
library(tidyr)

# ---- Helper: Hedge's g + 95% CI via effsize::cohen.d ----
hedges_g_ci <- function(x, y, conf.level = 0.95) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  if (length(x) < 2 || length(y) < 2) {
    return(tibble(
      n1 = length(x), n2 = length(y),
      g = NA_real_, ci_low = NA_real_, ci_high = NA_real_
    ))
  }
  
  d <- c(x, y)
  f <- factor(c(rep("A", length(x)), rep("B", length(y))), levels = c("A", "B"))
  
  cd <- effsize::cohen.d(
    d = d,
    f = f,
    hedges.correction = TRUE,
    conf.level = conf.level
  )
  
  message(sprintf(
    "Hedge's g (95%% CI): g = %.3f [%.3f, %.3f] (n1 = %d, n2 = %d)",
    unname(cd$estimate),
    unname(cd$conf.int[1]),
    unname(cd$conf.int[2]),
    length(x),
    length(y)
  ))
  
  tibble(
    n1 = length(x),
    n2 = length(y),
    g = unname(cd$estimate),
    ci_low = unname(cd$conf.int[1]),
    ci_high = unname(cd$conf.int[2])
  )
}




# ---- Helper: group descriptives (mean/sd/n) ----
group_stats <- function(df, group_cols, value_col) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n = sum(is.finite(.data[[value_col]])),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd   = sd(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

# ---- Helper: Hedge's g for all pairwise comparisons within a stratum ----
hedges_pairwise <- function(df, strata_cols, treat_col, value_col, treat_levels = NULL) {
  d <- df %>%
    filter(is.finite(.data[[value_col]])) %>%
    mutate(.treat = as.character(.data[[treat_col]]))
  
  if (!is.null(treat_levels)) {
    d <- d %>% mutate(.treat = factor(.treat, levels = treat_levels))
  } else {
    d <- d %>% mutate(.treat = factor(.treat))
  }
  
  lvl <- levels(d$.treat)
  comps <- combn(lvl, 2, simplify = FALSE)
  
  map_dfr(comps, function(cc) {
    a <- as.character(cc[1]); b <- as.character(cc[2])
    
    d_strata <- d %>% select(all_of(strata_cols)) %>% distinct()
    
    map_dfr(seq_len(nrow(d_strata)), function(i) {
      key <- d_strata[i, , drop = FALSE]
      
      sub <- d %>% semi_join(key, by = strata_cols)
      x <- sub %>% filter(.treat == a) %>% pull(.data[[value_col]])
      y <- sub %>% filter(.treat == b) %>% pull(.data[[value_col]])
      
      out <- hedges_g_ci(x, y) %>%
        mutate(comparison = paste0(a, " - ", b))
      
      bind_cols(key, out)
    })
  }) %>%
    relocate(all_of(strata_cols), comparison)
}

# ==========================================================
# g1: Composite occupancy trend % change 2000 -> 2023 (MC draws)
# Requires: zone_year_draw_gm from your script
# ==========================================================
g1_change_draws <- zone_year_draw_gm %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Rarity, Zone, Year, draw_id, gm) %>%
  pivot_wider(names_from = Year, values_from = gm, names_prefix = "Y") %>%
  mutate(
    s = pmax(.data[[paste0("Y", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Y", end_year)]],   1e-12),
    pct_change = 100 * (f - s) / s
  ) %>%
  select(Rarity, Zone, draw_id, pct_change)

g1_stats <- group_stats(g1_change_draws, c("Rarity", "Zone"), "pct_change")

# Zone comparisons within each Rarity
g1_g_zone_within_rarity <- hedges_pairwise(
  g1_change_draws,
  strata_cols = c("Rarity"),
  treat_col   = "Zone",
  value_col   = "pct_change",
  treat_levels = zones_keep
)

# Rare vs Common within each Zone (Rare - Common is direction in comparison label)
g1_g_rarity_within_zone <- hedges_pairwise(
  g1_change_draws,
  strata_cols = c("Zone"),
  treat_col   = "Rarity",
  value_col   = "pct_change",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# g2: Growth rate (%/yr) (MC draws)
# Requires: zone_start_end_draws from your script
# ==========================================================
g2_draws <- zone_start_end_draws %>%
  select(Rarity, Zone, draw_id, growth_rate_pct)

g2_stats <- group_stats(g2_draws, c("Rarity", "Zone"), "growth_rate_pct")

g2_g_zone_within_rarity <- hedges_pairwise(
  g2_draws,
  strata_cols = c("Rarity"),
  treat_col   = "Zone",
  value_col   = "growth_rate_pct",
  treat_levels = zones_keep
)

g2_g_rarity_within_zone <- hedges_pairwise(
  g2_draws,
  strata_cols = c("Zone"),
  treat_col   = "Rarity",
  value_col   = "growth_rate_pct",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# g5: Richness trend % change 2000 -> 2023 (NO-MC -> bootstrap)
# Requires: psi_zone_rich from your script (Species, Zone, Year, psi, Rarity)
# ==========================================================
set.seed(777)
B <- 2000  # bootstrap replicates

make_mat <- function(df_group) {
  mat <- df_group %>%
    select(Species, Year, psi) %>%
    pivot_wider(names_from = Year, values_from = psi) %>%
    arrange(Species)
  
  yrs <- sort(unique(df_group$Year))
  M <- as.matrix(mat[, as.character(yrs), drop = FALSE])
  rownames(M) <- mat$Species
  colnames(M) <- as.character(yrs)
  list(M = M, yrs = yrs)
}

groups_g5 <- psi_zone_rich %>%
  distinct(Rarity, Zone) %>%
  arrange(Rarity, Zone)

g5_boot_list <- vector("list", nrow(groups_g5))

for (k in seq_len(nrow(groups_g5))) {
  rr <- groups_g5$Rarity[k]
  zn <- groups_g5$Zone[k]
  
  df_g <- psi_zone_rich %>% filter(Rarity == rr, Zone == zn)
  obj <- make_mat(df_g)
  M <- obj$M
  
  nsp <- nrow(M)
  idx <- replicate(B, sample.int(nsp, size = nsp, replace = TRUE))
  
  sums <- apply(idx, 2, function(ii) colSums(M[ii, , drop = FALSE], na.rm = TRUE))
  
  s <- sums[as.character(start_year), ]
  f <- sums[as.character(end_year),   ]
  pct_change <- 100 * (f - s) / pmax(s, 1e-12)
  
  g5_boot_list[[k]] <- tibble(
    Rarity = rr,
    Zone   = zn,
    boot   = seq_len(B),
    pct_change = pct_change
  )
}

g5_change_boot <- bind_rows(g5_boot_list)

g5_stats <- group_stats(g5_change_boot, c("Rarity", "Zone"), "pct_change")

g5_g_zone_within_rarity <- hedges_pairwise(
  g5_change_boot,
  strata_cols = c("Rarity"),
  treat_col   = "Zone",
  value_col   = "pct_change",
  treat_levels = zones_keep
)

g5_g_rarity_within_zone <- hedges_pairwise(
  g5_change_boot,
  strata_cols = c("Zone"),
  treat_col   = "Rarity",
  value_col   = "pct_change",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# g4: AllParks-based richness (period mean richness draws)
# Requires: richness_period_draws from your script (rarity, Zone, draw, richness_mean)
# ==========================================================
g4_draws <- richness_period_draws %>%
  rename(Rarity = rarity, draw_id = draw, value = richness_mean) %>%
  select(Rarity, Zone, draw_id, value)

g4_stats <- group_stats(g4_draws, c("Rarity", "Zone"), "value")

g4_g_zone_within_rarity <- hedges_pairwise(
  g4_draws,
  strata_cols = c("Rarity"),
  treat_col   = "Zone",
  value_col   = "value",
  treat_levels = zones_keep
)

g4_g_rarity_within_zone <- hedges_pairwise(
  g4_draws,
  strata_cols = c("Zone"),
  treat_col   = "Rarity",
  value_col   = "value",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# EXPORT: Excel workbook with a sheet per plot
# Each sheet contains:
# 1) group mean/sd/n
# 2) Hedge's g (zones within rarity)
# 3) Hedge's g (rarity within zone)
# ==========================================================
wb <- createWorkbook()

write_plot_sheet <- function(sheet, stats_df, g_zone_df, g_rarity_df) {
  addWorksheet(wb, sheet)
  
  writeData(wb, sheet, "Group summaries (mean, sd, n)", startRow = 1, startCol = 1)
  writeData(wb, sheet, stats_df, startRow = 2, startCol = 1)
  
  r2 <- 4 + nrow(stats_df)
  writeData(wb, sheet, "Hedge's g (95% CI): Zone comparisons within each Rarity", startRow = r2, startCol = 1)
  writeData(wb, sheet, g_zone_df, startRow = r2 + 1, startCol = 1)
  
  r3 <- r2 + 3 + nrow(g_zone_df)
  writeData(wb, sheet, "Hedge's g (95% CI): Rare vs Common within each Zone", startRow = r3, startCol = 1)
  writeData(wb, sheet, g_rarity_df, startRow = r3 + 1, startCol = 1)
  
  setColWidths(wb, sheet, cols = 1:30, widths = "auto")
}

write_plot_sheet("g1_change_00_23", g1_stats, g1_g_zone_within_rarity, g1_g_rarity_within_zone)
write_plot_sheet("g2_growth_rate",  g2_stats, g2_g_zone_within_rarity, g2_g_rarity_within_zone)
write_plot_sheet("g5_change_00_23", g5_stats, g5_g_zone_within_rarity, g5_g_rarity_within_zone)
write_plot_sheet("g4_period_rich",  g4_stats, g4_g_zone_within_rarity, g4_g_rarity_within_zone)

excel_out <- "2_HedgesG_Treatment_Differences_byPlot.xlsx"
saveWorkbook(wb, excel_out, overwrite = TRUE)

cat("\nSaved Excel workbook:\n  ", excel_out, "\n", sep = "")

















































# Repeating this with the converged dataset only ####


############################################
# MASTER SCRIPT (unified): Rare vs Common plots + final patchwork
# Creates: g1, g2, g5, g4  -> patchworked into one column (one row each)
############################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(grid)      # unit()
  library(lubridate) # parse_date_time
  library(patchwork)
})

# ----------------------------
# 0) USER SETTINGS (shared)
# ----------------------------
setwd("C:/Users/georg/OneDrive - University of Reading/George Allen - PhD Master Folder/Year Three/Chapter 2 - Occupancy Modelling/Analysis Nov 2025/Nov_Outhwaite_Outputs")

bugs_file <- "Combined_BUGS_Data_Outhwaite_finalclustv5_32000iterations_FINAL_WITH_TAXCORR.csv"
allparks_path <- "C:/Users/georg/OneDrive - University of Reading/George Allen - PhD Master Folder/Year One/Statistics/Mapping/GBIF Sink Source AllParks/R analysis - 25.11.2024/R final/AllParks_records_FinalFiltered_Better2.csv"

start_year <- 2000
end_year   <- 2023
years_vec  <- start_year:end_year
year_index_offset <- 1999
zones_keep <- c("Core", "Buffer", "Outside")
psi_zone_regex <- "^psi\\.fs\\.r_(Core|Buffer|Outside)\\[(\\d+)\\]$"

# Colors
zone_colors <- c("Core"="#b2182b","Buffer"="black","Outside"="grey65")
zone_fill_colors <- c("Core"=alpha("#b2182b", 0.20),
                      "Buffer"=alpha("black", 0.20),
                      "Outside" = alpha("grey40",  0.30))   # darker than grey65 @ 0.20

# Monte Carlo
n_draws_occ <- 999     # for g1/g2 (occupancy composites + growth)
n_draws_rich <- 999    # for g4 (AllParks richness MC)
eps <- 1e-6

# All output files
out_csv_clean   <- "BUGS_psi_fs_r_byZone_byYear_2000_2023.csv"
out_csv_trend   <- "Cooke_style_RareCommon_composite_geo_occ_median_CrI_byZoneYear.csv"
out_csv_growth  <- "Cooke_style_RareCommon_growth_rate_byZone_median_CrI.csv"

out_png_trend   <- "Cooke_style_RareCommon_composite_trends_2000_2023.png"
out_jpg_trend   <- "Cooke_style_RareCommon_composite_trends_2000_2023.jpg"
out_png_growth  <- "Cooke_style_RareCommon_growth_rate_2000_2023.png"
out_jpg_growth  <- "Cooke_style_RareCommon_growth_rate_2000_2023.jpg"

out_png_richness_trend <- "Cooke_style_RareCommon_richness_trends_2000_2023_NO_MC.png"
out_jpg_richness_trend <- "Cooke_style_RareCommon_richness_trends_2000_2023_NO_MC.jpg"

out_png_rich_mc <- "richness_by_zone_RareCommon.png"

out_png_patch <- "rarevscommon_patchwork.png"
out_jpg_patch <- "rarevscommon_patchwork.jpg"


# Append _FIRSTLAST to all filenames
append_FL <- function(x) sub("(\\.\\w+)$", "_FIRSTLAST\\1", x)

out_csv_clean  <- append_FL(out_csv_clean)
out_csv_trend  <- append_FL(out_csv_trend)
out_csv_growth <- append_FL(out_csv_growth)

out_png_trend <- append_FL(out_png_trend)
out_jpg_trend <- append_FL(out_jpg_trend)
out_png_growth <- append_FL(out_png_growth)
out_jpg_growth <- append_FL(out_jpg_growth)

out_png_richness_trend <- append_FL(out_png_richness_trend)
out_jpg_richness_trend <- append_FL(out_jpg_richness_trend)
out_png_rich_mc <- append_FL(out_png_rich_mc)

out_png_patch <- append_FL(out_png_patch)
out_jpg_patch <- append_FL(out_jpg_patch)

excel_out_abs <- append_FL("2_AbsoluteEffectSizes_MedianDiff_HDI_byPlot.xlsx")
excel_out <- append_FL("2_HedgesG_Treatment_Differences_byPlot.xlsx")

# Also rename the converged-only export
converged_csv <- append_FL("BUGS_psi_fs_r_byZone_byYear_2000_2023_CONVERGEDONLY.csv")


# ----------------------------
# 1) LOAD BUGS ONCE
# ----------------------------
bugs <- read.csv(bugs_file, stringsAsFactors = FALSE)

# ----------------------------
# 2) Shared helpers
# ----------------------------
disperse_labels <- function(df, gap = 0.06) {
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
    if (y[1] < low) y <- seq(low, high, length.out = n)
  }
  
  df$label_y_adj <- NA_real_
  df$label_y_adj[ord] <- y
  df
}

beta_draws_from_mean_sd <- function(m, s, n = 1000) {
  m <- pmin(pmax(m, 1e-6), 1 - 1e-6)
  s <- pmax(s, 1e-12)
  v <- s^2
  vmax <- m * (1 - m)
  
  if (!is.finite(v) || v <= 0 || !is.finite(vmax) || vmax <= 0) v <- 0.01 * vmax
  if (v >= vmax) v <- 0.95 * vmax
  
  k <- (m * (1 - m) / v) - 1
  a <- m * k
  b <- (1 - m) * k
  
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
    v <- 0.5 * vmax
    k <- (m * (1 - m) / v) - 1
    a <- max(m * k, 1e-3)
    b <- max((1 - m) * k, 1e-3)
  }
  
  rbeta(n, shape1 = a, shape2 = b)
}

std_key <- function(x) str_squish(str_to_lower(as.character(x)))

hdi_vec <- function(x, credMass = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) return(c(low = NA_real_, high = NA_real_))
  x <- sort(x)
  m <- floor(credMass * n)
  if (m < 1) return(c(low = NA_real_, high = NA_real_))
  lows  <- x[1:(n - m)]
  highs <- x[(1 + m):n]
  widths <- highs - lows
  j <- which.min(widths)
  c(low = lows[j], high = highs[j])
}

fit_beta_from_ci <- function(l, u) {
  l <- pmin(1 - eps, pmax(eps, l))
  u <- pmin(1 - eps, pmax(eps, u))
  if (!is.finite(l) || !is.finite(u) || l >= u) return(c(NA_real_, NA_real_))
  
  obj <- function(par) {
    a <- exp(par[1]); b <- exp(par[2])
    ql <- qbeta(0.025, a, b)
    qu <- qbeta(0.975, a, b)
    (ql - l)^2 + (qu - u)^2
  }
  
  opt <- try(optim(log(c(2,2)), obj, method = "Nelder-Mead"), silent = TRUE)
  if (inherits(opt, "try-error")) return(c(NA_real_, NA_real_))
  
  a <- exp(opt$par[1]); b <- exp(opt$par[2])
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(c(NA_real_, NA_real_))
  c(a, b)
}

beta_from_mean_sd <- function(mu, sd) {
  mu <- pmin(1 - eps, pmax(eps, mu))
  if (!is.finite(mu) || !is.finite(sd) || sd <= 0) return(c(NA_real_, NA_real_))
  v <- sd^2
  k <- mu*(1-mu)/v - 1
  if (!is.finite(k) || k <= 0) return(c(NA_real_, NA_real_))
  c(mu*k, (1-mu)*k)
}

# ----------------------------
# 3) Extract psi.fs.r_* and then FILTER TO CONVERGED SPECIES ONLY
# ----------------------------

# First: full psi.fs.r_ table
psi_zone_full <- bugs %>%
  filter(str_detect(Parameter, "^psi\\.fs\\.r_")) %>%
  mutate(
    Zone = str_match(Parameter, psi_zone_regex)[, 2],
    YearIndex = suppressWarnings(as.integer(str_match(Parameter, psi_zone_regex)[, 3]))
  ) %>%
  filter(!is.na(Zone), !is.na(YearIndex)) %>%
  mutate(
    Year = YearIndex + year_index_offset,
    Zone = factor(Zone, levels = zones_keep),
    psi  = suppressWarnings(as.numeric(mean)),
    sd   = suppressWarnings(as.numeric(sd))
  ) %>%
  filter(Year >= start_year, Year <= end_year) %>%
  select(Species, Zone, Year, YearIndex, Parameter,
         psi, sd, X2.5., X25., X50., X75., X97.5., Rhat, n.eff)

# ---- CONVERGENCE FILTER (endpoint years, all zones must converge) ----
rhat_max_filtered  <- 1.10
neff_min_filtered  <- 0        # you can set >0 if you want n.eff threshold
require_all_zones_in_both_years <- TRUE

psi_zone_for_filter <- psi_zone_full %>%
  filter(Year %in% c(start_year, end_year),
         Zone %in% zones_keep) %>%
  mutate(
    Rhat_num = suppressWarnings(as.numeric(Rhat)),
    neff_num = suppressWarnings(as.numeric(n.eff)),
    ok_conv  = is.finite(Rhat_num) & (Rhat_num <= rhat_max_filtered) &
      is.finite(neff_num) & (neff_num >= neff_min_filtered) &
      is.finite(psi) & is.finite(sd) & !is.na(psi) & !is.na(sd)
  )

if (require_all_zones_in_both_years) {
  species_keep_filtered <- psi_zone_for_filter %>%
    group_by(Species, Year) %>%
    summarise(
      n_zones_present = n_distinct(Zone),
      n_zones_ok      = n_distinct(Zone[ok_conv]),
      .groups = "drop"
    ) %>%
    mutate(
      ok_year = (n_zones_present == length(zones_keep)) &
        (n_zones_ok      == length(zones_keep))
    ) %>%
    group_by(Species) %>%
    summarise(ok_both_years = all(ok_year), .groups = "drop") %>%
    filter(ok_both_years) %>%
    pull(Species)
} else {
  # looser: at least one converged zone in each endpoint year
  species_keep_filtered <- psi_zone_for_filter %>%
    group_by(Species, Year) %>%
    summarise(any_ok = any(ok_conv), .groups = "drop") %>%
    group_by(Species) %>%
    summarise(ok_both_years = all(any_ok), .groups = "drop") %>%
    filter(ok_both_years) %>%
    pull(Species)
}

cat("\n==================== CONVERGED FILTER (endpoints) ====================\n")
cat("Total species:", n_distinct(psi_zone_full$Species), "\n")
cat("Retained converged species:", length(species_keep_filtered), "\n")

# NOW restrict psi_zone_full to converged species ONLY
psi_zone_full <- psi_zone_full %>%
  filter(Species %in% species_keep_filtered)

# Save the converged-only table (so you know what this script actually used)
write.csv(psi_zone_full, converged_csv, row.names = FALSE)
cat("Saved converged-only table:", converged_csv, "\n")


# ----------------------------
# 4) g1 + g2: Cooke-style multispecies trends (Rare vs Common; MC; geom mean; growth)
# ----------------------------
set.seed(1)

# Rare/Common definition (BUGS; Cooke-style):
# - For each species: compute annual occupancy (averaged across zones), then median across years
# - Rare  = lower quartile (Q1) of species medians
# - Common = upper quartile (Q3) of species medians

species_year_occ <- psi_zone_full %>%
  filter(is.finite(psi)) %>%
  group_by(Species, Year) %>%
  summarise(psi_year = mean(psi, na.rm = TRUE), .groups = "drop")  # avg across zones within year

species_median_occ <- species_year_occ %>%
  group_by(Species) %>%
  summarise(median_psi = median(psi_year, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(median_psi))

q1_psi <- quantile(species_median_occ$median_psi, 0.25, na.rm = TRUE)
q3_psi <- quantile(species_median_occ$median_psi, 0.75, na.rm = TRUE)

rarity_lookup <- species_median_occ %>%
  mutate(
    Rarity = case_when(
      median_psi <= q1_psi ~ "Rare",
      median_psi >= q3_psi ~ "Common",
      TRUE ~ NA_character_
    ),
    Rarity = factor(Rarity, levels = c("Rare", "Common"))
  ) %>%
  filter(!is.na(Rarity)) %>%
  select(Species, Rarity)

psi_zone_rc <- psi_zone_full %>%
  left_join(rarity_lookup, by = "Species") %>%
  filter(!is.na(Rarity)) %>%
  filter(!is.na(psi), !is.na(sd), is.finite(psi), is.finite(sd)) %>%
  mutate(
    psi_mean = pmin(pmax(psi, 1e-6), 1 - 1e-6),
    psi_sd   = pmax(sd, 1e-12)
  )



# Monte Carlo draws per species×zone×year
draws_long <- psi_zone_rc %>%
  rowwise() %>%
  mutate(draw = list(beta_draws_from_mean_sd(psi_mean, psi_sd, n = n_draws_occ))) %>%
  ungroup() %>%
  select(Rarity, Species, Zone, Year, draw) %>%
  unnest_longer(draw, values_to = "psi_draw") %>%
  group_by(Rarity, Species, Zone, Year) %>%
  mutate(draw_id = row_number()) %>%
  ungroup() %>%
  mutate(psi_draw = pmin(pmax(psi_draw, 1e-12), 1 - 1e-12))

# Composite geometric mean per draw -> median + 95% CrI
zone_year_draw_gm <- draws_long %>%
  group_by(Rarity, Zone, Year, draw_id) %>%
  summarise(gm = exp(mean(log(psi_draw))), .groups = "drop")

zone_trend <- zone_year_draw_gm %>%
  group_by(Rarity, Zone, Year) %>%
  summarise(
    gm_median = median(gm),
    gm_low    = as.numeric(quantile(gm, 0.025)),
    gm_high   = as.numeric(quantile(gm, 0.975)),
    .groups = "drop"
  ) %>%
  arrange(Rarity, Zone, Year)

write.csv(zone_trend, out_csv_trend, row.names = FALSE)

# Growth rate (%/yr) from start->end on draws
y_years <- end_year - start_year

zone_start_end_draws <- zone_year_draw_gm %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Rarity, Zone, Year, draw_id, gm) %>%
  pivot_wider(names_from = Year, values_from = gm, names_prefix = "Year_") %>%
  filter(!is.na(.data[[paste0("Year_", start_year)]]),
         !is.na(.data[[paste0("Year_", end_year)]])) %>%
  mutate(
    s = pmax(.data[[paste0("Year_", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Year_", end_year)]], 1e-12),
    growth_rate_pct = 100 * ((f / s)^(1 / y_years) - 1)
  )

growth_summary <- zone_start_end_draws %>%
  group_by(Rarity, Zone) %>%
  summarise(
    growth_median = median(growth_rate_pct),
    growth_low    = as.numeric(quantile(growth_rate_pct, 0.025)),
    growth_high   = as.numeric(quantile(growth_rate_pct, 0.975)),
    .groups = "drop"
  ) %>%
  arrange(Rarity, Zone)

write.csv(growth_summary, out_csv_growth, row.names = FALSE)

# Endpoint label placement (show ONLY growth_median)
facet_limits_occ <- zone_trend %>%
  group_by(Rarity) %>%
  summarise(
    ymin = min(gm_low,  na.rm = TRUE),
    ymax = max(gm_high, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    yrange = ymax - ymin,
    pad = pmax(0.005, 0.05 * yrange)
  )

zone_end_occ <- zone_trend %>%
  filter(Year == end_year) %>%
  left_join(growth_summary %>% select(Rarity, Zone, growth_median), by = c("Rarity","Zone")) %>%
  left_join(facet_limits_occ, by = "Rarity") %>%
  mutate(
    label_text = paste0(sprintf("%.2f", growth_median), "%/yr"),
    off = pmax(0.01, 0.18 * yrange),
    label_base = case_when(
      Zone == "Core"    ~ gm_median - off,
      Zone == "Buffer"  ~ gm_median,
      Zone == "Outside" ~ gm_median + off,
      TRUE              ~ gm_median
    )
  ) %>%
  group_by(Rarity) %>%
  group_modify(~disperse_labels(.x, gap = 0.10 * unique(.x$yrange))) %>%
  ungroup() %>%
  mutate(
    x_label = end_year + 3,
    x_curve = end_year + 2.5
  )

# g1
g1 <- ggplot(zone_trend, aes(x = Year, y = gm_median, color = Zone)) +
  geom_ribbon(aes(ymin = gm_low, ymax = gm_high, fill = Zone),
              alpha = 0.40, color = NA, show.legend = FALSE) +
  #geom_line(aes(y = gm_high, group = Zone), linetype = "dashed", linewidth = 0.7, show.legend = FALSE) +
  #geom_line(aes(y = gm_low,  group = Zone), linetype = "dashed", linewidth = 0.7, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = end_year, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  geom_point(
    data = zone_end_occ,
    aes(x = end_year, y = gm_median, color = Zone),
    shape = 21, fill = "white", stroke = 1.5, size = 3.2,
    inherit.aes = FALSE
  ) +
  geom_curve(
    data = zone_end_occ,
    aes(x = end_year, y = gm_median, xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.7, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = zone_end_occ,
    aes(x = x_label, y = label_y_adj, label = label_text, color = Zone),
    hjust = 0, size = 4.2, fontface = "bold",
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  facet_wrap(~Rarity, nrow = 1, ncol = 2, scales = "free_y") +
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_fill_colors) +
  theme_minimal(base_size = 12) +
  labs(x = "Year", y = "Geometric mean occupancy", color = "Zone") +
  theme(
    legend.position = "right",
    legend.key.size = unit(2.2, "lines"),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  ) +
  coord_cartesian(xlim = c(start_year, end_year + 9))

g1

ggsave(g1, filename = out_png_trend, width = 10, height = 3, dpi = 1000)
ggsave(g1, filename = out_jpg_trend, width = 10, height = 3, dpi = 1000)

# g2
g2 <- ggplot(growth_summary, aes(x = Zone, y = growth_median, color = Zone)) +
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = "dashed", color = "grey30") +
  geom_errorbar(aes(ymin = growth_low, ymax = growth_high),
                width = 0.25, linewidth = 0.9, show.legend = FALSE) +
  geom_point(size = 4, stroke = 1.2, show.legend = FALSE) +
  facet_wrap(~Rarity, nrow = 1) +
  scale_color_manual(values = zone_colors) +
  theme_minimal(base_size = 12) +
  labs(x = "Zone", y = "Annual growth rate (% per year)") +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text  = element_text(size = 12)
  )

ggsave(g2, filename = out_png_growth, width = 10, height = 3, dpi = 1000)
ggsave(g2, filename = out_jpg_growth, width = 10, height = 3, dpi = 1000)

# ----------------------------
# 5) g5: Cooke-style richness trends (NO MC) (Rare/Common via median mean psi)
# ----------------------------
psi_zone_rich <- psi_zone_full %>%
  mutate(psi = suppressWarnings(as.numeric(psi))) %>%
  filter(is.finite(psi)) %>%
  select(Species, Zone, Year, psi) %>%
  left_join(rarity_lookup, by = "Species") %>%
  filter(!is.na(Rarity))

richness_trend <- psi_zone_rich %>%
  group_by(Rarity, Zone, Year) %>%
  summarise(
    richness = sum(psi, na.rm = TRUE),
    n_species = n_distinct(Species),
    .groups = "drop"
  ) %>%
  arrange(Rarity, Zone, Year)

change_df <- richness_trend %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Rarity, Zone, Year, richness) %>%
  pivot_wider(names_from = Year, values_from = richness, names_prefix = "Y") %>%
  mutate(
    start_val = .data[[paste0("Y", start_year)]],
    end_val   = .data[[paste0("Y", end_year)]],
    pct_change = 100 * (end_val - start_val) / start_val
  ) %>%
  select(Rarity, Zone, start_val, end_val, pct_change)

facet_limits_rich <- richness_trend %>%
  group_by(Rarity) %>%
  summarise(
    ymin = min(richness, na.rm = TRUE),
    ymax = max(richness, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    yrange = ymax - ymin,
    pad = pmax(0.5, 0.05 * yrange)
  )

zone_end_rich <- richness_trend %>%
  filter(Year == end_year) %>%
  left_join(change_df, by = c("Rarity","Zone")) %>%
  left_join(facet_limits_rich, by = "Rarity") %>%
  mutate(
    label_text = paste0(sprintf("%.1f", pct_change), "%"),
    off = pmax(0.5, 0.18 * yrange),
    label_base = case_when(
      Zone == "Core"    ~ richness - off,
      Zone == "Buffer"  ~ richness,
      Zone == "Outside" ~ richness + off,
      TRUE              ~ richness
    )
  ) %>%
  group_by(Rarity) %>%
  group_modify(~disperse_labels(.x, gap = 0.10 * unique(.x$yrange))) %>%
  ungroup() %>%
  mutate(
    x_label = end_year + 3,
    x_curve = end_year + 2.5
  )

g5 <- ggplot(richness_trend, aes(x = Year, y = richness, color = Zone)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  geom_vline(xintercept = end_year, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  geom_point(
    data = zone_end_rich,
    aes(x = end_year, y = richness, color = Zone),
    shape = 21, fill = "white", stroke = 1.4, size = 3.0,
    inherit.aes = FALSE
  ) +
  geom_curve(
    data = zone_end_rich,
    aes(x = end_year, y = richness, xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.7, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = zone_end_rich,
    aes(x = x_label, y = label_y_adj, label = label_text, color = Zone),
    hjust = 0, size = 4.2, fontface = "bold",
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  facet_wrap(~Rarity, scales = "free_y", nrow = 1) +
  scale_color_manual(values = zone_colors) +
  coord_cartesian(xlim = c(start_year, end_year + 9)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.key.size = unit(2.2, "lines"),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Year", y = "Species richness", color = "Zone")

ggsave(g5, filename = out_png_richness_trend, width = 10, height = 3, dpi = 1000)
ggsave(g5, filename = out_jpg_richness_trend, width = 10, height = 3, dpi = 1000)

# ----------------------------
# 6) g4: Cooke et al. 2023 - Rare/Common richness quartiles (BUGS definition) (MC cloud)
#     Rare/Common = lower/upper quartile of species' median occupancy across years (from BUGS)
#     (Optional filter retained): only species with >= min_records in AllParks
# ----------------------------
set.seed(123)

min_records <- 50

# A) Build AllParks species filter (>= min_records) -- no rarity from AllParks
rec_raw <- read.csv(allparks_path, stringsAsFactors = FALSE)

records0 <- rec_raw %>%
  select(species, eventDate, layer_grid) %>%
  distinct() %>%
  filter(!is.na(species), species != "", !is.na(layer_grid), layer_grid != "") %>%
  mutate(
    eventDate = suppressWarnings(parse_date_time(
      eventDate,
      orders = c("dmy", "ymd", "Ymd HMS", "Ymd", "mdy", "dmy HMS", "ymd HMS")
    )),
    year = as.integer(format(eventDate, "%Y"))
  ) %>%
  filter(!is.na(year), year >= start_year, year <= end_year)

sp_counts <- records0 %>%
  count(species, name = "n_records") %>%
  filter(n_records >= min_records) %>%
  mutate(sp_key = std_key(species))

# B) Extract BUGS psi summaries (mean + CI) and join BUGS rarity (quartiles)
psi_zone_ci <- bugs %>%
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
    l95 = suppressWarnings(as.numeric(X2.5.)),
    u95 = suppressWarnings(as.numeric(X97.5.))
  ) %>%
  filter(Year >= start_year, Year <= end_year) %>%
  filter(is.finite(psi_mean), is.finite(l95), is.finite(u95)) %>%
  select(Species, sp_key, Zone, Year, psi_mean, l95, u95)

psi_labeled <- psi_zone_ci %>%
  inner_join(sp_counts %>% select(sp_key), by = "sp_key") %>%   # >= min_records filter retained
  inner_join(rarity_lookup %>%
               mutate(sp_key = std_key(Species)) %>%
               select(sp_key, rarity = Rarity),
             by = "sp_key") %>%                                  # BUGS quartile rarity
  mutate(
    Zone = factor(Zone, levels = zones_keep),
    rarity = factor(rarity, levels = c("Rare","Common"))
  )

# C) Monte Carlo richness per rarity × zone × year
rich_draws_list <- list()

for (rr in levels(psi_labeled$rarity)) {
  for (zn in zones_keep) {
    for (yr in years_vec) {
      
      sub <- psi_labeled %>% filter(rarity == rr, Zone == zn, Year == yr)
      if (nrow(sub) == 0) next
      
      richness_draws <- numeric(n_draws_rich)
      
      for (i in seq_len(nrow(sub))) {
        l <- sub$l95[i]; u <- sub$u95[i]
        ab <- fit_beta_from_ci(l, u)
        
        if (anyNA(ab)) {
          sd_approx <- (u - l) / (2 * 1.96)
          ab <- beta_from_mean_sd(sub$psi_mean[i], sd_approx)
        }
        
        if (anyNA(ab)) {
          draws <- rnorm(n_draws_rich, mean = sub$psi_mean[i], sd = 0.02)
          draws <- pmin(1 - eps, pmax(eps, draws))
        } else {
          draws <- rbeta(n_draws_rich, ab[1], ab[2])
        }
        
        richness_draws <- richness_draws + draws
      }
      
      rich_draws_list[[length(rich_draws_list) + 1]] <- tibble(
        rarity = rr,
        Zone = zn,
        Year = yr,
        draw = seq_len(n_draws_rich),
        richness = richness_draws
      )
    }
  }
}

richness_draws_year <- bind_rows(rich_draws_list) %>%
  mutate(
    rarity = factor(rarity, levels = c("Rare","Common")),
    Zone = factor(Zone, levels = zones_keep)
  )

# Period mean richness across years per rarity × zone × draw
richness_period_draws <- richness_draws_year %>%
  group_by(rarity, Zone, draw) %>%
  summarise(richness_mean = mean(richness, na.rm = TRUE), .groups = "drop")

richness_zone_summary <- richness_period_draws %>%
  group_by(rarity, Zone) %>%
  summarise(
    median = median(richness_mean, na.rm = TRUE),
    hdi_low  = hdi_vec(richness_mean, 0.95)["low"],
    hdi_high = hdi_vec(richness_mean, 0.95)["high"],
    .groups = "drop"
  )

# g4 (same plotting style)
g4 <- ggplot(richness_period_draws, aes(x = Zone, y = richness_mean, color = Zone)) +
  geom_jitter(width = 0.4, alpha = 0.20, size = 1.5, show.legend = FALSE) +
  geom_hline(
    data = richness_zone_summary,
    aes(yintercept = median, color = Zone),
    linewidth = 0.7,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  facet_wrap(~rarity, nrow = 1, scales = "free_y") +
  scale_color_manual(values = zone_colors) +
  theme_minimal(base_size = 14) +
  labs(x = NULL, y = "Species richness (∑ years)") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  ) +
  coord_cartesian(clip = "off")

ggsave(out_png_rich_mc, g4, width = 10, height = 3, dpi = 1000)




# ----------------------------
# 7) FINAL PATCHWORK (one column, one plot per row; numbered; no legends)
# ----------------------------
g_all <- (g1 / g5 / g4 / g2) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "none")

ggsave(g_all, filename = out_png_patch, width = 10, height = 12, dpi = 1000)
ggsave(g_all, filename = out_jpg_patch, width = 10, height = 12, dpi = 1000)

print(g_all)














# Absolute effect sizes ####

############################################################
# ADD-ON: Absolute effect sizes (median Δ + 95% credible interval)
# Method: HDI via HDInterval::hdi (Meredith & Kruschke, 2018)
# Outputs Excel: one sheet per plot (g1, g2, g5, g4)
############################################################

# ---- Packages ----
if (!requireNamespace("HDInterval", quietly = TRUE)) install.packages("HDInterval")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

library(HDInterval)
library(openxlsx)
library(dplyr)
library(tidyr)
library(purrr)

# ---- Helper: median + 95% HDI for a vector ----
median_hdi <- function(x, credMass = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(tibble(median = NA_real_, hdi_low = NA_real_, hdi_high = NA_real_, n = length(x)))
  }
  h <- HDInterval::hdi(x, credMass = credMass)
  tibble(
    median = median(x),
    hdi_low = unname(h[1]),
    hdi_high = unname(h[2]),
    n = length(x)
  )
}

# ---- Helper: group summaries (median + 95% HDI) ----
group_summaries_hdi <- function(df, group_cols, value_col, credMass = 0.95) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      tmp = list(median_hdi(.data[[value_col]], credMass = credMass)),
      .groups = "drop"
    ) %>%
    tidyr::unnest(tmp)
}

# ---- Helper: absolute pairwise differences (A - B) with median + 95% HDI ----
# Expects df has ONE row per draw/replicate per group (e.g. draw_id, boot, draw)
pairwise_abs_effects <- function(df, strata_cols, treat_col, draw_col, value_col,
                                 treat_levels = NULL, credMass = 0.95) {
  
  d <- df %>%
    filter(is.finite(.data[[value_col]])) %>%
    mutate(.treat = as.character(.data[[treat_col]]))
  
  if (!is.null(treat_levels)) d <- d %>% mutate(.treat = factor(.treat, levels = treat_levels))
  lvl <- levels(factor(d$.treat))
  comps <- combn(lvl, 2, simplify = FALSE)
  
  # Wide per stratum + draw so we can subtract draw-wise
  wide <- d %>%
    select(all_of(strata_cols), !!sym(draw_col), .treat, !!sym(value_col)) %>%
    pivot_wider(names_from = .treat, values_from = !!sym(value_col))
  
  strata_keys <- wide %>% distinct(across(all_of(strata_cols)))
  
  out <- list()
  
  for (i in seq_len(nrow(strata_keys))) {
    key <- strata_keys[i, , drop = FALSE]
    subw <- wide %>% semi_join(key, by = strata_cols)
    
    for (cc in comps) {
      a <- as.character(cc[1]); b <- as.character(cc[2])
      
      if (!(a %in% names(subw)) || !(b %in% names(subw))) next
      
      delta <- subw[[a]] - subw[[b]]
      summ <- median_hdi(delta, credMass = credMass) %>%
        mutate(comparison = paste0(a, " - ", b))
      
      out[[length(out) + 1]] <- bind_cols(key, summ)
    }
  }
  
  bind_rows(out) %>%
    relocate(all_of(strata_cols), comparison)
}

# ==========================================================
# 1) g1 absolute effects
# Data: zone_year_draw_gm -> gm draws per Rarity×Zone×Year×draw_id
# For absolute effects in your paper:
#   - either use period-mean richness style (if g1 is "richness by zone"),
#   - OR if g1 is a time trend, use change 2000->2023 per draw as your effect variable.
#
# Below I do BOTH, so you can choose:
#   (A) Period-mean across years (draw-wise)
#   (B) Change 2000->2023 (draw-wise)
# ==========================================================

# (A) Period-mean (mean across years) per draw
g1_period_draw <- zone_year_draw_gm %>%
  group_by(Rarity, Zone, draw_id) %>%
  summarise(value = mean(gm, na.rm = TRUE), .groups = "drop")

g1_period_groups <- group_summaries_hdi(g1_period_draw, c("Rarity", "Zone"), "value")

g1_period_zone_abs <- pairwise_abs_effects(
  df = g1_period_draw,
  strata_cols = c("Rarity"),
  treat_col = "Zone",
  draw_col  = "draw_id",
  value_col = "value",
  treat_levels = zones_keep
)

g1_period_rarity_abs <- pairwise_abs_effects(
  df = g1_period_draw,
  strata_cols = c("Zone"),
  treat_col = "Rarity",
  draw_col  = "draw_id",
  value_col = "value",
  treat_levels = c("Rare", "Common")
)

# (B) Change 2000 -> 2023 per draw (percent change)
g1_change_draws <- zone_year_draw_gm %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Rarity, Zone, Year, draw_id, gm) %>%
  pivot_wider(names_from = Year, values_from = gm, names_prefix = "Y") %>%
  mutate(
    s = pmax(.data[[paste0("Y", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Y", end_year)]],   1e-12),
    pct_change = 100 * (f - s) / s
  ) %>%
  select(Rarity, Zone, draw_id, pct_change)

g1_change_groups <- group_summaries_hdi(g1_change_draws, c("Rarity", "Zone"), "pct_change")

g1_change_zone_abs <- pairwise_abs_effects(
  df = g1_change_draws,
  strata_cols = c("Rarity"),
  treat_col = "Zone",
  draw_col  = "draw_id",
  value_col = "pct_change",
  treat_levels = zones_keep
)

g1_change_rarity_abs <- pairwise_abs_effects(
  df = g1_change_draws,
  strata_cols = c("Zone"),
  treat_col = "Rarity",
  draw_col  = "draw_id",
  value_col = "pct_change",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# 2) g2 absolute effects (growth rate %/yr)
# Data: zone_start_end_draws has growth_rate_pct per draw_id
# ==========================================================
g2_draws <- zone_start_end_draws %>%
  select(Rarity, Zone, draw_id, growth_rate_pct) %>%
  rename(value = growth_rate_pct)

g2_groups <- group_summaries_hdi(g2_draws, c("Rarity", "Zone"), "value")

g2_zone_abs <- pairwise_abs_effects(
  df = g2_draws,
  strata_cols = c("Rarity"),
  treat_col = "Zone",
  draw_col  = "draw_id",
  value_col = "value",
  treat_levels = zones_keep
)

g2_rarity_abs <- pairwise_abs_effects(
  df = g2_draws,
  strata_cols = c("Zone"),
  treat_col = "Rarity",
  draw_col  = "draw_id",
  value_col = "value",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# 3) g5 absolute effects (richness trend, NO-MC -> bootstrap)
# Data: psi_zone_rich used to build bootstrap change 2000->2023
# ==========================================================
set.seed(777)
B <- 2000

make_mat <- function(df_group) {
  mat <- df_group %>%
    select(Species, Year, psi) %>%
    pivot_wider(names_from = Year, values_from = psi) %>%
    arrange(Species)
  
  yrs <- sort(unique(df_group$Year))
  M <- as.matrix(mat[, as.character(yrs), drop = FALSE])
  rownames(M) <- mat$Species
  colnames(M) <- as.character(yrs)
  list(M = M, yrs = yrs)
}

groups_g5 <- psi_zone_rich %>%
  distinct(Rarity, Zone) %>%
  arrange(Rarity, Zone)

g5_boot_list <- vector("list", nrow(groups_g5))

for (k in seq_len(nrow(groups_g5))) {
  rr <- groups_g5$Rarity[k]
  zn <- groups_g5$Zone[k]
  
  df_g <- psi_zone_rich %>% filter(Rarity == rr, Zone == zn)
  obj <- make_mat(df_g)
  M <- obj$M
  
  nsp <- nrow(M)
  idx <- replicate(B, sample.int(nsp, size = nsp, replace = TRUE))
  
  sums <- apply(idx, 2, function(ii) colSums(M[ii, , drop = FALSE], na.rm = TRUE))
  
  s <- sums[as.character(start_year), ]
  f <- sums[as.character(end_year),   ]
  pct_change <- 100 * (f - s) / pmax(s, 1e-12)
  
  g5_boot_list[[k]] <- tibble(
    Rarity = rr,
    Zone   = zn,
    boot_id = seq_len(B),
    value = pct_change
  )
}

g5_change_boot <- bind_rows(g5_boot_list)

g5_groups <- group_summaries_hdi(g5_change_boot, c("Rarity", "Zone"), "value")

g5_zone_abs <- pairwise_abs_effects(
  df = g5_change_boot,
  strata_cols = c("Rarity"),
  treat_col = "Zone",
  draw_col  = "boot_id",
  value_col = "value",
  treat_levels = zones_keep
)

g5_rarity_abs <- pairwise_abs_effects(
  df = g5_change_boot,
  strata_cols = c("Zone"),
  treat_col = "Rarity",
  draw_col  = "boot_id",
  value_col = "value",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# 4) g4 absolute effects (AllParks richness, period-mean across years)
# Data: richness_period_draws already has richness_mean per draw
# ==========================================================
g4_draws <- richness_period_draws %>%
  rename(Rarity = rarity, draw_id = draw) %>%
  select(Rarity, Zone, draw_id, richness_mean) %>%
  rename(value = richness_mean)

g4_groups <- group_summaries_hdi(g4_draws, c("Rarity", "Zone"), "value")

g4_zone_abs <- pairwise_abs_effects(
  df = g4_draws,
  strata_cols = c("Rarity"),
  treat_col = "Zone",
  draw_col  = "draw_id",
  value_col = "value",
  treat_levels = zones_keep
)

g4_rarity_abs <- pairwise_abs_effects(
  df = g4_draws,
  strata_cols = c("Zone"),
  treat_col = "Rarity",
  draw_col  = "draw_id",
  value_col = "value",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# 5) EXPORT: Excel workbook
# One sheet per plot/metric, placed after your Hedge's g workbook
# ==========================================================
wb_abs <- createWorkbook()

write_abs_sheet <- function(wb, sheet, groups_df, zone_abs_df, rarity_abs_df) {
  addWorksheet(wb, sheet)
  
  writeData(wb, sheet, "Group summaries (median + 95% HDI)", startRow = 1, startCol = 1)
  writeData(wb, sheet, groups_df, startRow = 2, startCol = 1)
  
  r2 <- 4 + nrow(groups_df)
  writeData(wb, sheet, "Absolute effect sizes: Zone contrasts within each Rarity (A - B; median + 95% HDI)",
            startRow = r2, startCol = 1)
  writeData(wb, sheet, zone_abs_df, startRow = r2 + 1, startCol = 1)
  
  r3 <- r2 + 3 + nrow(zone_abs_df)
  writeData(wb, sheet, "Absolute effect sizes: Rare vs Common within each Zone (Rare - Common; median + 95% HDI)",
            startRow = r3, startCol = 1)
  writeData(wb, sheet, rarity_abs_df, startRow = r3 + 1, startCol = 1)
  
  setColWidths(wb, sheet, cols = 1:40, widths = "auto")
}

# g1: two sheets so you can choose what you report
write_abs_sheet(wb_abs, "g1_period_mean",   g1_period_groups, g1_period_zone_abs, g1_period_rarity_abs)
write_abs_sheet(wb_abs, "g1_change_00_23",  g1_change_groups, g1_change_zone_abs,  g1_change_rarity_abs)

# g2, g5, g4
write_abs_sheet(wb_abs, "g2_growth_rate",   g2_groups, g2_zone_abs, g2_rarity_abs)
write_abs_sheet(wb_abs, "g5_change_00_23",  g5_groups, g5_zone_abs, g5_rarity_abs)
write_abs_sheet(wb_abs, "g4_period_rich",   g4_groups, g4_zone_abs, g4_rarity_abs)

excel_out_abs <- "2_AbsoluteEffectSizes_MedianDiff_HDI_byPlot.xlsx"
saveWorkbook(wb_abs, excel_out_abs, overwrite = TRUE)

cat("\nSaved Excel workbook (absolute effects):\n  ", excel_out_abs, "\n", sep = "")


















# Hedge's G CIs and relative effect sizes ####

############################################################
# ADD-ON: Hedge's g (95% CI) treatment differences per plot
# Exports to Excel: one sheet per plot (g1, g2, g5, g4)
############################################################

# ---- Packages ----
if (!requireNamespace("effsize", quietly = TRUE)) install.packages("effsize")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")

library(effsize)
library(openxlsx)
library(purrr)
library(dplyr)
library(tidyr)

# ---- Helper: Hedge's g + 95% CI via effsize::cohen.d ----
hedges_g_ci <- function(x, y, conf.level = 0.95) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  if (length(x) < 2 || length(y) < 2) {
    return(tibble(
      n1 = length(x), n2 = length(y),
      g = NA_real_, ci_low = NA_real_, ci_high = NA_real_
    ))
  }
  
  d <- c(x, y)
  f <- factor(c(rep("A", length(x)), rep("B", length(y))), levels = c("A", "B"))
  
  cd <- effsize::cohen.d(
    d = d,
    f = f,
    hedges.correction = TRUE,
    conf.level = conf.level
  )
  
  message(sprintf(
    "Hedge's g (95%% CI): g = %.3f [%.3f, %.3f] (n1 = %d, n2 = %d)",
    unname(cd$estimate),
    unname(cd$conf.int[1]),
    unname(cd$conf.int[2]),
    length(x),
    length(y)
  ))
  
  tibble(
    n1 = length(x),
    n2 = length(y),
    g = unname(cd$estimate),
    ci_low = unname(cd$conf.int[1]),
    ci_high = unname(cd$conf.int[2])
  )
}




# ---- Helper: group descriptives (mean/sd/n) ----
group_stats <- function(df, group_cols, value_col) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n = sum(is.finite(.data[[value_col]])),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd   = sd(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

# ---- Helper: Hedge's g for all pairwise comparisons within a stratum ----
hedges_pairwise <- function(df, strata_cols, treat_col, value_col, treat_levels = NULL) {
  d <- df %>%
    filter(is.finite(.data[[value_col]])) %>%
    mutate(.treat = as.character(.data[[treat_col]]))
  
  if (!is.null(treat_levels)) {
    d <- d %>% mutate(.treat = factor(.treat, levels = treat_levels))
  } else {
    d <- d %>% mutate(.treat = factor(.treat))
  }
  
  lvl <- levels(d$.treat)
  comps <- combn(lvl, 2, simplify = FALSE)
  
  map_dfr(comps, function(cc) {
    a <- as.character(cc[1]); b <- as.character(cc[2])
    
    d_strata <- d %>% select(all_of(strata_cols)) %>% distinct()
    
    map_dfr(seq_len(nrow(d_strata)), function(i) {
      key <- d_strata[i, , drop = FALSE]
      
      sub <- d %>% semi_join(key, by = strata_cols)
      x <- sub %>% filter(.treat == a) %>% pull(.data[[value_col]])
      y <- sub %>% filter(.treat == b) %>% pull(.data[[value_col]])
      
      out <- hedges_g_ci(x, y) %>%
        mutate(comparison = paste0(a, " - ", b))
      
      bind_cols(key, out)
    })
  }) %>%
    relocate(all_of(strata_cols), comparison)
}

# ==========================================================
# g1: Composite occupancy trend % change 2000 -> 2023 (MC draws)
# Requires: zone_year_draw_gm from your script
# ==========================================================
g1_change_draws <- zone_year_draw_gm %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Rarity, Zone, Year, draw_id, gm) %>%
  pivot_wider(names_from = Year, values_from = gm, names_prefix = "Y") %>%
  mutate(
    s = pmax(.data[[paste0("Y", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Y", end_year)]],   1e-12),
    pct_change = 100 * (f - s) / s
  ) %>%
  select(Rarity, Zone, draw_id, pct_change)

g1_stats <- group_stats(g1_change_draws, c("Rarity", "Zone"), "pct_change")

# Zone comparisons within each Rarity
g1_g_zone_within_rarity <- hedges_pairwise(
  g1_change_draws,
  strata_cols = c("Rarity"),
  treat_col   = "Zone",
  value_col   = "pct_change",
  treat_levels = zones_keep
)

# Rare vs Common within each Zone (Rare - Common is direction in comparison label)
g1_g_rarity_within_zone <- hedges_pairwise(
  g1_change_draws,
  strata_cols = c("Zone"),
  treat_col   = "Rarity",
  value_col   = "pct_change",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# g2: Growth rate (%/yr) (MC draws)
# Requires: zone_start_end_draws from your script
# ==========================================================
g2_draws <- zone_start_end_draws %>%
  select(Rarity, Zone, draw_id, growth_rate_pct)

g2_stats <- group_stats(g2_draws, c("Rarity", "Zone"), "growth_rate_pct")

g2_g_zone_within_rarity <- hedges_pairwise(
  g2_draws,
  strata_cols = c("Rarity"),
  treat_col   = "Zone",
  value_col   = "growth_rate_pct",
  treat_levels = zones_keep
)

g2_g_rarity_within_zone <- hedges_pairwise(
  g2_draws,
  strata_cols = c("Zone"),
  treat_col   = "Rarity",
  value_col   = "growth_rate_pct",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# g5: Richness trend % change 2000 -> 2023 (NO-MC -> bootstrap)
# Requires: psi_zone_rich from your script (Species, Zone, Year, psi, Rarity)
# ==========================================================
set.seed(777)
B <- 2000  # bootstrap replicates

make_mat <- function(df_group) {
  mat <- df_group %>%
    select(Species, Year, psi) %>%
    pivot_wider(names_from = Year, values_from = psi) %>%
    arrange(Species)
  
  yrs <- sort(unique(df_group$Year))
  M <- as.matrix(mat[, as.character(yrs), drop = FALSE])
  rownames(M) <- mat$Species
  colnames(M) <- as.character(yrs)
  list(M = M, yrs = yrs)
}

groups_g5 <- psi_zone_rich %>%
  distinct(Rarity, Zone) %>%
  arrange(Rarity, Zone)

g5_boot_list <- vector("list", nrow(groups_g5))

for (k in seq_len(nrow(groups_g5))) {
  rr <- groups_g5$Rarity[k]
  zn <- groups_g5$Zone[k]
  
  df_g <- psi_zone_rich %>% filter(Rarity == rr, Zone == zn)
  obj <- make_mat(df_g)
  M <- obj$M
  
  nsp <- nrow(M)
  idx <- replicate(B, sample.int(nsp, size = nsp, replace = TRUE))
  
  sums <- apply(idx, 2, function(ii) colSums(M[ii, , drop = FALSE], na.rm = TRUE))
  
  s <- sums[as.character(start_year), ]
  f <- sums[as.character(end_year),   ]
  pct_change <- 100 * (f - s) / pmax(s, 1e-12)
  
  g5_boot_list[[k]] <- tibble(
    Rarity = rr,
    Zone   = zn,
    boot   = seq_len(B),
    pct_change = pct_change
  )
}

g5_change_boot <- bind_rows(g5_boot_list)

g5_stats <- group_stats(g5_change_boot, c("Rarity", "Zone"), "pct_change")

g5_g_zone_within_rarity <- hedges_pairwise(
  g5_change_boot,
  strata_cols = c("Rarity"),
  treat_col   = "Zone",
  value_col   = "pct_change",
  treat_levels = zones_keep
)

g5_g_rarity_within_zone <- hedges_pairwise(
  g5_change_boot,
  strata_cols = c("Zone"),
  treat_col   = "Rarity",
  value_col   = "pct_change",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# g4: AllParks-based richness (period mean richness draws)
# Requires: richness_period_draws from your script (rarity, Zone, draw, richness_mean)
# ==========================================================
g4_draws <- richness_period_draws %>%
  rename(Rarity = rarity, draw_id = draw, value = richness_mean) %>%
  select(Rarity, Zone, draw_id, value)

g4_stats <- group_stats(g4_draws, c("Rarity", "Zone"), "value")

g4_g_zone_within_rarity <- hedges_pairwise(
  g4_draws,
  strata_cols = c("Rarity"),
  treat_col   = "Zone",
  value_col   = "value",
  treat_levels = zones_keep
)

g4_g_rarity_within_zone <- hedges_pairwise(
  g4_draws,
  strata_cols = c("Zone"),
  treat_col   = "Rarity",
  value_col   = "value",
  treat_levels = c("Rare", "Common")
)

# ==========================================================
# EXPORT: Excel workbook with a sheet per plot
# Each sheet contains:
# 1) group mean/sd/n
# 2) Hedge's g (zones within rarity)
# 3) Hedge's g (rarity within zone)
# ==========================================================
wb <- createWorkbook()

write_plot_sheet <- function(sheet, stats_df, g_zone_df, g_rarity_df) {
  addWorksheet(wb, sheet)
  
  writeData(wb, sheet, "Group summaries (mean, sd, n)", startRow = 1, startCol = 1)
  writeData(wb, sheet, stats_df, startRow = 2, startCol = 1)
  
  r2 <- 4 + nrow(stats_df)
  writeData(wb, sheet, "Hedge's g (95% CI): Zone comparisons within each Rarity", startRow = r2, startCol = 1)
  writeData(wb, sheet, g_zone_df, startRow = r2 + 1, startCol = 1)
  
  r3 <- r2 + 3 + nrow(g_zone_df)
  writeData(wb, sheet, "Hedge's g (95% CI): Rare vs Common within each Zone", startRow = r3, startCol = 1)
  writeData(wb, sheet, g_rarity_df, startRow = r3 + 1, startCol = 1)
  
  setColWidths(wb, sheet, cols = 1:30, widths = "auto")
}

write_plot_sheet("g1_change_00_23", g1_stats, g1_g_zone_within_rarity, g1_g_rarity_within_zone)
write_plot_sheet("g2_growth_rate",  g2_stats, g2_g_zone_within_rarity, g2_g_rarity_within_zone)
write_plot_sheet("g5_change_00_23", g5_stats, g5_g_zone_within_rarity, g5_g_rarity_within_zone)
write_plot_sheet("g4_period_rich",  g4_stats, g4_g_zone_within_rarity, g4_g_rarity_within_zone)

excel_out <- "2_HedgesG_Treatment_Differences_byPlot.xlsx"
saveWorkbook(wb, excel_out, overwrite = TRUE)

cat("\nSaved Excel workbook:\n  ", excel_out, "\n", sep = "")


























# Richness change overall - No rare/common split ####


# ----------------------------
# ADD-ON: Overall richness trends (all species, no Rare/Common)
# ----------------------------

# Optional: output filenames
out_png_richness_all <- "Cooke_style_ALLSPECIES_richness_trends_2000_2023_NO_MC.png"
out_jpg_richness_all <- "Cooke_style_ALLSPECIES_richness_trends_2000_2023_NO_MC.jpg"

# 1) Total richness per Zone × Year (sum of psi across all species)
richness_trend_all <- psi_zone_full %>%
  mutate(psi = suppressWarnings(as.numeric(psi))) %>%
  filter(is.finite(psi)) %>%
  group_by(Zone, Year) %>%
  summarise(
    richness  = sum(psi, na.rm = TRUE),
    n_species = n_distinct(Species),
    .groups   = "drop"
  ) %>%
  arrange(Zone, Year)

# 2) % change 2000 -> 2023 per Zone
change_all_df <- richness_trend_all %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Zone, Year, richness) %>%
  tidyr::pivot_wider(
    names_from  = Year,
    values_from = richness,
    names_prefix = "Y"
  ) %>%
  mutate(
    start_val  = .data[[paste0("Y", start_year)]],
    end_val    = .data[[paste0("Y", end_year)]],
    pct_change = 100 * (end_val - start_val) / pmax(start_val, 1e-12)
  ) %>%
  select(Zone, start_val, end_val, pct_change)

# 3) Panel-wise y limits for label placement
facet_limits_all <- richness_trend_all %>%
  summarise(
    ymin = min(richness, na.rm = TRUE),
    ymax = max(richness, na.rm = TRUE)
  ) %>%
  mutate(
    yrange = ymax - ymin,
    pad    = pmax(0.5, 0.05 * yrange)
  )

# 4) Endpoint labels (2023 richness + % change)
zone_end_rich_all <- richness_trend_all %>%
  filter(Year == end_year) %>%
  left_join(change_all_df, by = "Zone") %>%
  crossing(facet_limits_all) %>%  # to grab ymin/ymax/yrange/pad
  mutate(
    label_text = paste0(sprintf("%.1f", pct_change), "%"),
    off = pmax(0.5, 0.18 * yrange),
    label_base = case_when(
      Zone == "Core"    ~ richness - off,
      Zone == "Buffer"  ~ richness,
      Zone == "Outside" ~ richness + off,
      TRUE              ~ richness
    )
  ) %>%
  disperse_labels(gap = 0.10 *10) %>%
  mutate(
    x_label = end_year + 3,
    x_curve = end_year + 2.5
  )

# 5) Plot: overall richness trend (no Rare/Common facets)
g5_all <- ggplot(richness_trend_all, aes(x = Year, y = richness, color = Zone)) +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  geom_vline(xintercept = end_year, linetype = "dashed",
             color = "grey30", linewidth = 0.8) +
  geom_point(
    data = zone_end_rich_all,
    aes(x = end_year, y = richness, color = Zone),
    shape = 21, fill = "white", stroke = 1.4, size = 3.0,
    inherit.aes = FALSE
  ) +
  geom_curve(
    data = zone_end_rich_all,
    aes(x = end_year, y = richness,
        xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.7, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = zone_end_rich_all,
    aes(x = x_label, y = label_y_adj, label = label_text, color = Zone),
    hjust = 0, size = 4.2, fontface = "bold",
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  scale_color_manual(values = zone_colors) +
  coord_cartesian(xlim = c(start_year, end_year + 9)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position      = "right",
    legend.key.size      = unit(2.2, "lines"),
    legend.text          = element_text(size = 16),
    strip.text           = element_text(size = 14, face = "bold"),
    axis.title           = element_text(size = 13),
    axis.text            = element_text(size = 12),
    panel.grid.major.y   = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x   = element_blank()
  ) +
  labs(
    x     = "Year",
    y     = "Species richness (all species)",
    color = "Zone",
    title = "Change in overall richness (2000–2023)"
  )

g5_all

ggsave(g5_all, filename = out_png_richness_all, width = 10, height = 3, dpi = 1000)
ggsave(g5_all, filename = out_jpg_richness_all, width = 10, height = 3, dpi = 1000)

































# RICHNESS Total richness per zone (all years, all species) ####

# ==========================================================
# ADD-ON: Jitter plot of overall richness by Zone (ALL species)
# ==========================================================

# 1) Prepare psi summaries with CIs (all species, all zones)
psi_zone_ci_all <- psi_zone_full %>%
  mutate(
    psi_mean = suppressWarnings(as.numeric(psi)),
    l95      = suppressWarnings(as.numeric(X2.5.)),
    u95      = suppressWarnings(as.numeric(X97.5.))
  ) %>%
  filter(
    is.finite(psi_mean), is.finite(l95), is.finite(u95),
    Zone %in% zones_keep,
    Year >= start_year, Year <= end_year
  ) %>%
  select(Species, Zone, Year, psi_mean, l95, u95)

# 2) Monte Carlo richness draws per Zone × Year (ALL species)
rich_draws_list_all <- list()

for (zn in zones_keep) {
  for (yr in years_vec) {
    
    sub <- psi_zone_ci_all %>% filter(Zone == zn, Year == yr)
    if (nrow(sub) == 0) next
    
    richness_draws <- numeric(n_draws_rich)
    
    for (i in seq_len(nrow(sub))) {
      l  <- sub$l95[i]
      u  <- sub$u95[i]
      ab <- fit_beta_from_ci(l, u)
      
      if (anyNA(ab)) {
        sd_approx <- (u - l) / (2 * 1.96)
        ab <- beta_from_mean_sd(sub$psi_mean[i], sd_approx)
      }
      
      if (anyNA(ab)) {
        draws <- rnorm(n_draws_rich, mean = sub$psi_mean[i], sd = 0.02)
        draws <- pmin(1 - eps, pmax(eps, draws))
      } else {
        draws <- rbeta(n_draws_rich, ab[1], ab[2])
      }
      
      richness_draws <- richness_draws + draws
    }
    
    rich_draws_list_all[[length(rich_draws_list_all) + 1]] <- tibble(
      Zone     = zn,
      Year     = yr,
      draw     = seq_len(n_draws_rich),
      richness = richness_draws
    )
  }
}

richness_draws_year_all <- bind_rows(rich_draws_list_all) %>%
  mutate(Zone = factor(Zone, levels = zones_keep))

# 3) Period-mean richness across years per Zone × draw
richness_period_draws_all <- richness_draws_year_all %>%
  group_by(Zone, draw) %>%
  summarise(richness_mean = mean(richness, na.rm = TRUE), .groups = "drop")

# 4) Zone-wise medians + 95% HDIs (for dashed lines, if wanted)
richness_zone_summary_all <- richness_period_draws_all %>%
  group_by(Zone) %>%
  summarise(
    median   = median(richness_mean, na.rm = TRUE),
    hdi_low  = hdi_vec(richness_mean, 0.95)["low"],
    hdi_high = hdi_vec(richness_mean, 0.95)["high"],
    .groups = "drop"
  )

# 5) Jitter plot (like g4, but no rarity facet)
g_rich_all_jitter <- ggplot(richness_period_draws_all,
                            aes(x = Zone, y = richness_mean, color = Zone)) +
  geom_jitter(width = 0.4, alpha = 0.20, size = 1.5, show.legend = FALSE) +
  geom_hline(
    data = richness_zone_summary_all,
    aes(yintercept = median, color = Zone),
    linewidth = 0.7,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = zone_colors) +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Species richness (mean over 2000–2023,\nall species)",
    title = "Overall richness by zone (ALL species, MC)"
  ) +
  theme(
    legend.position    = "none",
    strip.text         = element_text(size = 14, face = "bold"),
    axis.title         = element_text(size = 13),
    axis.text          = element_text(size = 12),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  ) +
  coord_cartesian(clip = "off")

g_rich_all_jitter

# Optional save:
ggsave("richness_by_zone_ALLSPECIES_jitter.png",
       g_rich_all_jitter, width = 8, height = 3, dpi = 1000)
ggsave("richness_by_zone_ALLSPECIES_jitter.jpg",
       g_rich_all_jitter, width = 8, height = 3, dpi = 1000)






############################################################
# HDI add-on for overall richness jitter (ALL species)
# Uses: richness_period_draws_all (Zone, draw, richness_mean)
############################################################

# 1) Group summaries: median + 95% HDI per Zone
g_all_groups <- group_summaries_hdi(
  df        = richness_period_draws_all,
  group_cols = c("Zone"),
  value_col = "richness_mean"
)

# 2) Pairwise absolute differences between Zones (A - B)
#    Trick: add a dummy stratum "All" so pairwise_abs_effects works
g_all_zone_abs <- pairwise_abs_effects(
  df = richness_period_draws_all %>% dplyr::mutate(All = "All"),
  strata_cols   = c("All"),
  treat_col     = "Zone",
  draw_col      = "draw",
  value_col     = "richness_mean",
  treat_levels  = zones_keep
) %>%
  dplyr::select(-All)   # drop dummy column

# 3) Optional: export to a standalone Excel workbook
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

wb_all <- createWorkbook()

addWorksheet(wb_all, "ALLSPECIES_richness")

writeData(wb_all, "ALLSPECIES_richness",
          "Group summaries (median + 95% HDI) by Zone",
          startRow = 1, startCol = 1)
writeData(wb_all, "ALLSPECIES_richness",
          g_all_groups,
          startRow = 2, startCol = 1)

r2 <- 4 + nrow(g_all_groups)
writeData(wb_all, "ALLSPECIES_richness",
          "Absolute effect sizes: Zone contrasts (A - B; median Δ + 95% HDI)",
          startRow = r2, startCol = 1)
writeData(wb_all, "ALLSPECIES_richness",
          g_all_zone_abs,
          startRow = r2 + 1, startCol = 1)

setColWidths(wb_all, "ALLSPECIES_richness", cols = 1:30, widths = "auto")

excel_out_all <- "AbsoluteEffectSizes_ALLSPECIES_richness_jitter.xlsx"
saveWorkbook(wb_all, excel_out_all, overwrite = TRUE)

cat("\nSaved Excel workbook for ALL-species richness jitter:\n  ", excel_out_all, "\n", sep = "")








###############################################################################
# PROPORTIONAL REPRESENTATION OF RARE vs COMMON SPECIES IN EACH ZONE (MC)
###############################################################################

# 1) Prepare psi summaries with rarity labels
psi_zone_ci_rc <- psi_zone_full %>%
  left_join(rarity_lookup, by = "Species") %>%  # adds Rarity = Rare / Common
  filter(!is.na(Rarity)) %>%
  mutate(
    psi_mean = suppressWarnings(as.numeric(psi)),
    l95      = suppressWarnings(as.numeric(X2.5.)),
    u95      = suppressWarnings(as.numeric(X97.5.))
  ) %>%
  filter(
    is.finite(psi_mean), is.finite(l95), is.finite(u95),
    Zone %in% zones_keep,
    Year >= start_year, Year <= end_year
  ) %>%
  select(Species, Rarity, Zone, Year, psi_mean, l95, u95)


# 2) Monte Carlo richness draws per Zone × Year × Rarity
rich_rc_list <- list()

for (rr in c("Rare", "Common")) {
  for (zn in zones_keep) {
    for (yr in years_vec) {
      
      sub <- psi_zone_ci_rc %>% filter(Rarity == rr, Zone == zn, Year == yr)
      if (nrow(sub) == 0) next
      
      richness_draws <- numeric(n_draws_rich)
      
      for (i in seq_len(nrow(sub))) {
        l  <- sub$l95[i]
        u  <- sub$u95[i]
        ab <- fit_beta_from_ci(l, u)
        
        if (anyNA(ab)) {
          sd_approx <- (u - l) / (2 * 1.96)
          ab <- beta_from_mean_sd(sub$psi_mean[i], sd_approx)
        }
        
        if (anyNA(ab)) {
          draws <- rnorm(n_draws_rich, mean = sub$psi_mean[i], sd = 0.02)
          draws <- pmin(1 - eps, pmax(eps, draws))
        } else {
          draws <- rbeta(n_draws_rich, ab[1], ab[2])
        }
        
        richness_draws <- richness_draws + draws
      }
      
      rich_rc_list[[length(rich_rc_list) + 1]] <- tibble(
        Rarity   = rr,
        Zone     = zn,
        Year     = yr,
        draw     = seq_len(n_draws_rich),
        richness = richness_draws
      )
    }
  }
}

richness_rc_year <- bind_rows(rich_rc_list) %>%
  mutate(
    Zone   = factor(Zone,   levels = zones_keep),
    Rarity = factor(Rarity, levels = c("Rare","Common"))
  )


# 3) Period-mean richness across years per Zone × Rarity × draw
richness_rc_period <- richness_rc_year %>%
  group_by(Rarity, Zone, draw) %>%
  summarise(richness_mean = mean(richness, na.rm = TRUE), .groups = "drop")


# 4) Join with ALL-species richness (computed earlier)
#    richness_period_draws_all: Zone, draw, richness_mean (ALL species)

prop_df <- richness_rc_period %>%
  rename(richness_rc = richness_mean) %>%
  left_join(
    richness_period_draws_all %>% rename(richness_all = richness_mean),
    by = c("Zone","draw")
  ) %>%
  mutate(
    prop = richness_rc / pmax(richness_all, 1e-12)
  )


# 5) Summaries: median + 95% HDI (per Zone × Rarity)
prop_summary <- prop_df %>%
  group_by(Zone, Rarity) %>%
  summarise(
    median   = median(prop),
    hdi_low  = hdi_vec(prop)["low"],
    hdi_high = hdi_vec(prop)["high"],
    .groups = "drop"
  )


# 6) Print results
prop_summary
cat("\nProportional representation of Rare/Common per zone (median + 95% HDI):\n")
print(prop_summary)


# 7) Optional Excel export
wb_prop <- createWorkbook()
addWorksheet(wb_prop, "Proportion_RareCommon")

writeData(wb_prop, "Proportion_RareCommon",
          "Proportional representation (Rare/Common) by Zone", 
          startRow = 1, startCol = 1)

writeData(wb_prop, "Proportion_RareCommon",
          prop_summary,
          startRow = 2, startCol = 1)

setColWidths(wb_prop, "Proportion_RareCommon", cols = 1:20, widths = "auto")

excel_prop_out <- "ProportionalRepresentation_RareCommon_byZone.xlsx"
saveWorkbook(wb_prop, excel_prop_out, overwrite = TRUE)

cat("\nSaved:", excel_prop_out, "\n")





###############################################################################
# STACKED BAR: Proportional representation of Rare vs Common per Zone
###############################################################################

# Starting from prop_df (Zone, Rarity, draw, richness_rc, richness_all, prop)
# and prop_summary (Zone, Rarity, median, hdi_low, hdi_high)

# 1) Use median proportion per Zone × Rarity
prop_med <- prop_df %>%
  group_by(Zone, Rarity) %>%
  summarise(
    prop_med_all = median(prop, na.rm = TRUE),
    .groups = "drop"
  )

# 2) Rescale within each Zone so Rare + Common = 1
#    (i.e. proportion within the Rare+Common subset)
prop_med_stack <- prop_med %>%
  group_by(Zone) %>%
  mutate(
    prop_within_rc = prop_med_all / sum(prop_med_all, na.rm = TRUE)
  ) %>%
  ungroup()

# 3) A simple color palette for Rarity
rarity_colors <- c("Rare" = "#2166ac",  # blue
                   "Common" = "#b2182b") # red (matches your Core-ish red)

# 4) Stacked bar plot
g_prop_stack <- ggplot(prop_med_stack,
                       aes(x = Zone,
                           y = prop_within_rc,
                           fill = Rarity)) +
  geom_col(color = "white", width = 0.7) +
  scale_fill_manual(values = rarity_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Proportion of species (Rare vs Common)\nwithin Rare+Common subset",
    fill = "Rarity",
    title = "Proportional representation of Rare vs Common species by zone"
  ) +
  coord_cartesian(ylim=c(0.9,1))+
  theme(
    legend.position    = "right",
    legend.key.size    = unit(1.8, "lines"),
    legend.text        = element_text(size = 12),
    axis.title         = element_text(size = 13),
    axis.text          = element_text(size = 12),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  )

g_prop_stack

# Optional: save to file
ggsave("Proportional_RareCommon_byZone_STACKED.png",
       g_prop_stack, width = 6, height = 4, dpi = 1000)
ggsave("Proportional_RareCommon_byZone_STACKED.jpg",
       g_prop_stack, width = 6, height = 4, dpi = 1000)



############################################################
# Proportional representation (Rare/Common) by Zone + HDIs
# Uses richness_period_draws (Rare/Common) and
# richness_period_draws_all (ALL species)
############################################################

# 1) Build draw-wise proportions of richness for Rare/Common per Zone
prop_draws_zone_rarity <- richness_period_draws %>%
  dplyr::rename(Rarity = rarity) %>%
  dplyr::left_join(
    richness_period_draws_all %>%
      dplyr::rename(rich_all = richness_mean),
    by = c("Zone", "draw")
  ) %>%
  dplyr::mutate(
    prop = richness_mean / pmax(rich_all, eps),
    Rarity = factor(Rarity, levels = c("Rare", "Common")),
    Zone   = factor(Zone,   levels = zones_keep)
  )

# 2) Group summaries: median + 95% HDI of proportional representation
prop_groups_zone_rarity <- group_summaries_hdi(
  df         = prop_draws_zone_rarity,
  group_cols = c("Zone", "Rarity"),
  value_col  = "prop"
)

# 3) Zonewise absolute differences within each Rarity (Core–Buffer, etc.)
prop_zone_abs_within_rarity <- pairwise_abs_effects(
  df           = prop_draws_zone_rarity,
  strata_cols  = c("Rarity"),
  treat_col    = "Zone",
  draw_col     = "draw",
  value_col    = "prop",
  treat_levels = zones_keep
)

# 4) Rare vs Common within each Zone (Rare – Common)
prop_rarity_abs_within_zone <- pairwise_abs_effects(
  df           = prop_draws_zone_rarity,
  strata_cols  = c("Zone"),
  treat_col    = "Rarity",
  draw_col     = "draw",
  value_col    = "prop",
  treat_levels = c("Rare", "Common")
)

# 5) Export to Excel
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

wb_prop <- createWorkbook()

addWorksheet(wb_prop, "Prop_RareCommon_byZone")

## (i) Group summaries
writeData(
  wb_prop, "Prop_RareCommon_byZone",
  "Proportional representation (Rare/Common) by Zone",
  startRow = 1, startCol = 1
)
writeData(
  wb_prop, "Prop_RareCommon_byZone",
  prop_groups_zone_rarity,
  startRow = 2, startCol = 1
)

## (ii) Zone contrasts within each Rarity
r2 <- 4 + nrow(prop_groups_zone_rarity)
writeData(
  wb_prop, "Prop_RareCommon_byZone",
  "Absolute effect sizes: Zone contrasts within each Rarity (A - B; median + 95% HDI)",
  startRow = r2, startCol = 1
)
writeData(
  wb_prop, "Prop_RareCommon_byZone",
  prop_zone_abs_within_rarity,
  startRow = r2 + 1, startCol = 1
)

## (iii) Rare vs Common within each Zone
r3 <- r2 + 3 + nrow(prop_zone_abs_within_rarity)
writeData(
  wb_prop, "Prop_RareCommon_byZone",
  "Absolute effect sizes: Rare vs Common within each Zone (Rare - Common; median + 95% HDI)",
  startRow = r3, startCol = 1
)
writeData(
  wb_prop, "Prop_RareCommon_byZone",
  prop_rarity_abs_within_zone,
  startRow = r3 + 1, startCol = 1
)

setColWidths(wb_prop, "Prop_RareCommon_byZone", cols = 1:40, widths = "auto")

excel_out_prop <- "2_ProportionalRepresentation_RareCommon_byZone.xlsx"
saveWorkbook(wb_prop, excel_out_prop, overwrite = TRUE)

cat("\nSaved Excel workbook (proportional representation):\n  ",
    excel_out_prop, "\n", sep = "")

