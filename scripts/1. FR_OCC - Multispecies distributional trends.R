############################################
# multispecies trends
############################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(grid) # for unit()
  library(ggrepel)
  library(openxlsx)
})

# ----------------------------
# 0) USER SETTINGS
# ----------------------------
setwd("anon")

bugs_file <- "Combined_BUGS_Data_Outhwaite_finalclustv5_32000iterations_FINAL_WITH_TAXCORR.csv"

start_year <- 2000
end_year   <- 2023
year_index_offset <- 1999

zones_keep <- c("Core", "Buffer", "Outside")

# Monte Carlo settings
set.seed(1)
n_draws <- 999   

# Outputs
out_csv_clean   <- "BUGS_psi_fs_r_byZone_byYear_2000_2023.csv"
out_csv_comp    <- "composite_geo_occ_byZone_byYear_median_CrI.csv"
out_csv_growth  <- "growth_rate_byZone_draws_median_CrI.csv"

out_png_occ  <- "Cooke_style_composite_geo_occ_2000_2023.png"
out_jpg_occ  <- "Cooke_style_composite_geo_occ_2000_2023.jpg"
out_png_gr   <- "Cooke_style_growth_rate_2000_2023.png"
out_jpg_gr   <- "Cooke_style_growth_rate_2000_2023.jpg"

psi_zone_regex <- "^psi\\.fs\\.r_(Core|Buffer|Outside)\\[(\\d+)\\]$"

# Plot colours
zone_colors <- c("Core"="#b2182b","Buffer"="black","Outside"="grey65")
zone_fill_colors <- c("Core"=alpha("#b2182b", 0.20),
                      "Buffer"=alpha("black", 0.20),
                      "Outside"=alpha("grey40", 0.30))

# ----------------------------
# 1) LOAD
# ----------------------------
bugs <- read.csv(bugs_file, stringsAsFactors = FALSE)

# ----------------------------
# 2) EXTRACT psi.fs.r_{Zone}[YearIndex]
# ----------------------------
psi_zone <- bugs %>%
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

# ----------------------------
# 3) SANITY CHECKS + SAVE CLEAN TABLE
# ----------------------------
cat("\n==================== psi_zone sanity checks ====================\n")
cat("Rows:", nrow(psi_zone), "\n")
cat("Zones:\n"); print(table(psi_zone$Zone, useNA = "ifany"))
cat("Year range:\n"); print(range(psi_zone$Year, na.rm = TRUE))
cat("Unique species:", n_distinct(psi_zone$Species), "\n")

check_counts <- psi_zone %>%
  distinct(Species, Year, Zone) %>%
  count(Species, Year, name = "n_zones") %>%
  summarise(min_zones = min(n_zones), max_zones = max(n_zones), mean_zones = mean(n_zones))
print(check_counts)

write.csv(psi_zone, out_csv_clean, row.names = FALSE)
cat("\nSaved cleaned table:", out_csv_clean, "\n")

# ----------------------------
# 4) BETA MONTE CARLO DRAWS PER SPECIES-ZONE-YEAR
# ----------------------------
# We approximate each species-year-zone posterior for occupancy with a Beta distribution
# parameterised from mean and sd.

beta_draws_from_mean_sd <- function(m, s, n = 1000) {
  # m in (0,1), s >= 0
  m <- pmin(pmax(m, 1e-6), 1 - 1e-6)
  s <- pmax(s, 1e-12)
  v <- s^2
  
  vmax <- m * (1 - m)
  # Guardrails: make variance feasible
  if (!is.finite(v) || v <= 0 || !is.finite(vmax) || vmax <= 0) v <- 0.01 * vmax
  if (v >= vmax) v <- 0.95 * vmax
  
  k <- (m * (1 - m) / v) - 1
  a <- m * k
  b <- (1 - m) * k
  
  # Last-resort shrink if numerically unstable
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
    v <- 0.5 * vmax
    k <- (m * (1 - m) / v) - 1
    a <- max(m * k, 1e-3)
    b <- max((1 - m) * k, 1e-3)
  }
  
  rbeta(n, shape1 = a, shape2 = b)
}

plot_df <- psi_zone %>%
  filter(!is.na(psi), !is.na(sd), is.finite(psi), is.finite(sd)) %>%
  mutate(
    psi_mean = pmin(pmax(psi, 1e-6), 1 - 1e-6),
    psi_sd   = pmax(sd, 1e-12)
  )

# Generate draws
draws_long <- plot_df %>%
  rowwise() %>%
  mutate(draw = list(beta_draws_from_mean_sd(psi_mean, psi_sd, n = n_draws))) %>%
  ungroup() %>%
  select(Species, Zone, Year, draw) %>%
  unnest_longer(draw, values_to = "psi_draw") %>%
  group_by(Species, Zone, Year) %>%
  mutate(draw_id = row_number()) %>%
  ungroup() %>%
  mutate(psi_draw = pmin(pmax(psi_draw, 1e-12), 1 - 1e-12))  # safe for log

# ----------------------------
# 5) COMPOSITE MULTISPECIES GEOMETRIC MEAN BY ZONE-YEAR (PER DRAW)
# propagate uncertainty then summarise as median + 95% CrI
# ----------------------------
zone_year_draw_gm <- draws_long %>%
  group_by(Zone, Year, draw_id) %>%
  summarise(gm = exp(mean(log(psi_draw))), .groups = "drop")

zone_trend <- zone_year_draw_gm %>%
  group_by(Zone, Year) %>%
  summarise(
    gm_median = median(gm),
    gm_low    = as.numeric(quantile(gm, 0.025)),
    gm_high   = as.numeric(quantile(gm, 0.975)),
    .groups = "drop"
  ) %>%
  arrange(Zone, Year)

write.csv(zone_trend, out_csv_comp, row.names = FALSE)
cat("\nSaved composite occupancy summaries:", out_csv_comp, "\n")

# ----------------------------
# 6) CHANGE (START -> END) + LABEL DATA (endpoint annotation)
# ----------------------------
zone_change <- zone_trend %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Zone, Year, gm_median, gm_low, gm_high) %>%
  pivot_wider(
    names_from  = Year,
    values_from = c(gm_median, gm_low, gm_high),
    names_glue  = "{.value}_{Year}"
  ) %>%
  mutate(
    abs_change = .data[[paste0("gm_median_", end_year)]] - .data[[paste0("gm_median_", start_year)]],
    pct_change = 100 * abs_change / .data[[paste0("gm_median_", start_year)]]
  ) %>%
  arrange(Zone)

cat("\n==================== Composite change (start -> end) ====================\n")
print(zone_change)
# ----------------------------
# ENDPOINT LABELS (no boxes) + curved leaders
# - uses your disperse_labels() approach to prevent overlap
# ----------------------------

# Per-plot y-limits (based on ribbon + median line)
plot_limits <- zone_trend %>%
  summarise(
    ymin = min(gm_low,  na.rm = TRUE),
    ymax = max(gm_high, na.rm = TRUE)
  ) %>%
  mutate(
    yrange = ymax - ymin,
    pad = pmax(0.005, 0.05 * yrange)
  )

# Your helper, adapted to a single panel (3 zones)
disperse_labels <- function(df, gap = 0.06) {
  # expects df has: label_base, ymin, ymax, pad
  n <- nrow(df)
  low  <- df$ymin[1] + df$pad[1]
  high <- df$ymax[1] - df$pad[1]
  avail <- high - low
  
  if (n <= 1 || !is.finite(avail) || avail <= 0) {
    df$label_y_adj <- pmin(pmax(df$label_base, low), high)
    return(df)
  }
  
  # interpret gap as absolute y units if provided as absolute; we’ll pass absolute gap
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

# zone_end must exist; if not, recreate it cleanly:
zone_end <- zone_trend %>%
  filter(Year == end_year) %>%
  left_join(
    zone_change %>%
      transmute(
        Zone,
        start_val = .data[[paste0("gm_median_", start_year)]],
        end_val   = .data[[paste0("gm_median_", end_year)]],
        pct_change
      ),
    by = "Zone"
  ) %>%
  mutate(
    # label text (no boxes)
    label_text = paste0(round(end_val, 3), " (", round(pct_change, 1), "%)")
  ) %>%
  tidyr::crossing(plot_limits) %>%  # adds ymin/ymax/yrange/pad to each row
  mutate(
    # anchor labels near endpoint but biased by zone (same aesthetic as your guild plot)
    off = pmax(0.01, 0.18 * yrange),
    label_base = case_when(
      Zone == "Core"    ~ gm_median - off,
      Zone == "Buffer"  ~ gm_median,
      Zone == "Outside" ~ gm_median + off,
      TRUE              ~ gm_median
    )
  ) %>%
  # disperse within this single panel
  disperse_labels(gap = 0.10 * unique(plot_limits$yrange)) %>%
  mutate(
    # keep curve/label positions consistent
    x_label = end_year + 2.3,
    x_curve = end_year + 2.5
  )


# ----------------------------
# 7) PLOT 1: COMPOSITE OCCUPANCY (median + 95% CrI)
#   - ribbon + dashed CI edges
#   - endpoint markers
#   - curved leaders to endpoint labels
# ----------------------------
# Curved leader lines: draw a curve from endpoint to label position
# (geom_curve supports curvature; adjust curvature/sign if you want different bend)

g_occ <- ggplot(zone_trend, aes(x = Year, y = gm_median, color = Zone)) +
  geom_ribbon(aes(ymin = gm_low, ymax = gm_high, fill = Zone),
              alpha = 0.4, color = NA, show.legend = FALSE) +
  #geom_line(aes(y = gm_high, group = Zone), linetype = "dashed",
  #          linewidth = 0.8, show.legend = FALSE) +
  #geom_line(aes(y = gm_low, group = Zone), linetype = "dashed",
  #          linewidth = 0.8, show.legend = FALSE) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_fill_colors) +
  geom_vline(xintercept = end_year, linetype = "dashed",
             color = "grey30", linewidth = 0.8) +
  geom_point(
    aes(fill = Zone),          # fill mapped to Zone as well
    shape = 16,                # filled circle
    stroke = 1.5,
    size  = 2
  ) +
  # endpoint markers
  geom_point(
    data = zone_end,
    aes(x = end_year, y = gm_median, color = Zone),
    shape = 21, fill = "white", stroke = 1.5, size = 3.5,
    inherit.aes = FALSE
  ) +
  
  # curved leaders (like your guild plot)
  geom_curve(
    data = zone_end,
    aes(x = end_year, y = gm_median,
        xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.8, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  
  # endpoint labels (NO boxes)
  geom_text(
    data = zone_end,
    aes(x = x_label, y = label_y_adj, label = label_text, color = Zone),
    hjust = -0.05, size = 4, fontface = "bold",
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  theme_minimal(base_size = 10) +
  labs(
    x = "Year", y = "Geometric mean occupancy", color = "Zone"
  ) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(2.6, "lines"),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13, color = "gray30"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  ) +
  coord_cartesian(xlim = c(start_year, end_year + 10), ylim=c(0,0.2))

print(g_occ)

ggsave(g_occ, filename = out_png_occ, width = 8, height = 6, dpi = 1000)
ggsave(g_occ, filename = out_jpg_occ, width = 8, height = 6, dpi = 1000)

cat("\nSaved occupancy plots:\n  ", out_png_occ, "\n  ", out_jpg_occ, "\n")




# ----------------------------
# 8) COOKE/OUTHWAITE-STYLE GROWTH RATE (% PER YEAR), WITH UNCERTAINTY
# growth_rate = 100 * ( (f/s)^(1/y) - 1 )
# computed on composite draws (per zone) so uncertainty is propagated.
# ----------------------------
y_years <- end_year - start_year

# Get start/end composite per draw, per zone
zone_start_end_draws <- zone_year_draw_gm %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Zone, Year, draw_id, gm) %>%
  pivot_wider(names_from = Year, values_from = gm, names_prefix = "Year_") %>%
  filter(!is.na(.data[[paste0("Year_", start_year)]]),
         !is.na(.data[[paste0("Year_", end_year)]])) %>%
  mutate(
    s = pmax(.data[[paste0("Year_", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Year_", end_year)]], 1e-12),
    growth_rate_pct = 100 * ((f / s)^(1 / y_years) - 1)
  )

growth_summary <- zone_start_end_draws %>%
  group_by(Zone) %>%
  summarise(
    growth_median = median(growth_rate_pct),
    growth_low    = as.numeric(quantile(growth_rate_pct, 0.025)),
    growth_high   = as.numeric(quantile(growth_rate_pct, 0.975)),
    .groups = "drop"
  ) %>%
  arrange(Zone)

write.csv(growth_summary, out_csv_growth, row.names = FALSE)
cat("\nSaved growth rate summary:", out_csv_growth, "\n")
cat("\n==================== Growth rate (% per year) ====================\n")
print(growth_summary)

# ----------------------------
# 9) PLOT 2: GROWTH RATE BY ZONE
# ----------------------------
g_gr <- ggplot(growth_summary, aes(x = Zone, y = growth_median, color = Zone)) +
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(ymin = growth_low, ymax = growth_high),
                  linewidth = 1.0, fatten = 2.5, show.legend = FALSE) +
  scale_color_manual(values = zone_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Annual growth rate of composite occupancy (start → end)",
    subtitle = paste0("Growth rate = 100 * ((f/s)^(1/y) - 1), where y = ", y_years,
                      " years; points are medians, bars are 95% credible intervals"),
    x = "Zone", y = "Annual growth rate of occupancy (% per year)"
  ) +
  theme(
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13, color = "gray30")
  )

print(g_gr)

ggsave(g_gr, filename = out_png_gr, width = 10, height = 7, dpi = 1000)
ggsave(g_gr, filename = out_jpg_gr, width = 10, height = 7, dpi = 1000)

cat("\nSaved growth rate plots:\n  ", out_png_gr, "\n  ", out_jpg_gr, "\n")
cat("\nDONE.\n")














# trends with annual growth rates as endpoint labels ####

# ==========================================================
# ADD-ON: Endpoint labels show composite annual growth rate
# (from Section 8: growth_summary)
# ==========================================================

# Make sure growth_summary exists before running this block!
# (i.e., run Section 8 first)

zone_end <- zone_end %>%
  left_join(growth_summary, by = "Zone") %>%
  mutate(
    label_text = paste0(
      sprintf("%.2f", growth_median), "%/yr [",
      sprintf("%.2f", growth_low), ", ",
      sprintf("%.2f", growth_high), "]"
    )
  )

# (Optional) If you want a slightly shorter label, use this instead:
# zone_end <- zone_end %>%
#   mutate(label_text = paste0(sprintf("%.2f", gr_med), "%/yr"))




g_occ <- ggplot(zone_trend, aes(x = Year, y = gm_median, color = Zone)) +
  geom_ribbon(aes(ymin = gm_low, ymax = gm_high, fill = Zone),
              alpha = 0.30, color = NA, show.legend = FALSE) +
  #geom_line(aes(y = gm_high, group = Zone), linetype = "dashed",
  #          linewidth = 0.8, show.legend = FALSE) +
  #geom_line(aes(y = gm_low, group = Zone), linetype = "dashed",
  #          linewidth = 0.8, show.legend = FALSE) +
  geom_point(
    aes(fill = Zone),          # fill mapped to Zone as well
    shape = 16,                # filled circle
    stroke = 1.5,
    size  = 2
  ) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_fill_colors) +
  geom_vline(xintercept = end_year, linetype = "dashed",
             color = "grey30", linewidth = 0.8) +
  
  # endpoint markers
  geom_point(
    data = zone_end,
    aes(x = end_year, y = gm_median, color = Zone),
    shape = 21, fill = "white", stroke = 1.5, size = 3.5,
    inherit.aes = FALSE
  ) +
  
  # curved leaders (like your guild plot)
  geom_curve(
    data = zone_end,
    aes(x = end_year, y = gm_median,
        xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.8, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  
  # endpoint labels (NO boxes)
  geom_text(
    data = zone_end,
    aes(x = x_label, y = label_y_adj, label = label_text, color = Zone),
    hjust = 0, size = 4, fontface = "bold",
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  theme_minimal(base_size = 10) +
  labs(
    x = "Year", y = "Geometric mean occupancy (median)", color = "Zone"
  ) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(2.6, "lines"),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13, color = "gray30"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  ) +
  coord_cartesian(xlim = c(start_year, end_year + 9), ylim=c(0,0.2))

print(g_occ)

# filenames for the ANNUAL-GROWTH-LABELLED trend plot (so you don't overwrite)
out_png_occ_gr <- "Cooke_style_composite_geo_occ_2000_2023_AnnualGrowthLabels.png"
out_jpg_occ_gr <- "Cooke_style_composite_geo_occ_2000_2023_AnnualGrowthLabels.jpg"

ggsave(g_occ, filename = out_png_occ_gr, width = 8, height = 6, dpi = 1000)
ggsave(g_occ, filename = out_jpg_occ_gr, width = 8, height = 6, dpi = 1000)

cat("\nSaved ANNUAL-GROWTH-labelled occupancy plots:\n  ",
    out_png_occ_gr, "\n  ", out_jpg_occ_gr, "\n")




















############################################
# Absolute and Relative effect sizes and CIs
############################################


# ---- Packages ----
if (!requireNamespace("HDInterval", quietly = TRUE)) install.packages("HDInterval")
if (!requireNamespace("effsize", quietly = TRUE)) install.packages("effsize")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")

library(HDInterval)
library(effsize)
library(openxlsx)
library(purrr)
library(dplyr)
library(tidyr)

# ---- Settings ----
credMass <- 0.95

out_csv_ts_hdi     <- "1_Trends_composite_geo_occ_byZone_byYear_median_HDI.csv"
out_csv_change_hdi <- "endpoint_change_byZone_median_HDI.csv"
out_csv_growth_hdi <- "growth_rate_byZone_median_HDI.csv"
out_csv_hedges_g   <- "hedges_g_zone_contrasts.csv"
out_xlsx_hdi_g     <- "Composite_HDI_and_HedgesG_outputs.xlsx"

# ==========================================================
# Helpers
# ==========================================================
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

group_summaries_hdi <- function(df, group_cols, value_col, credMass = 0.95) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(tmp = list(median_hdi(.data[[value_col]], credMass = credMass)), .groups = "drop") %>%
    tidyr::unnest(tmp)
}

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
    d = d, f = f,
    hedges.correction = TRUE,
    conf.level = conf.level
  )
  
  tibble(
    n1 = length(x),
    n2 = length(y),
    g = unname(cd$estimate),
    ci_low = unname(cd$conf.int[1]),
    ci_high = unname(cd$conf.int[2])
  )
}

# ==========================================================
# 1) TIME-SERIES HDIs: composite gm by Zone-Year
# Uses zone_year_draw_gm (Zone, Year, draw_id, gm)
# ==========================================================
ts_hdi <- group_summaries_hdi(zone_year_draw_gm, c("Zone", "Year"), "gm", credMass = credMass) %>%
  rename(
    gm_median = median,
    gm_hdi_low = hdi_low,
    gm_hdi_high = hdi_high
  ) %>%
  arrange(Zone, Year)

write.csv(ts_hdi, out_csv_ts_hdi, row.names = FALSE)
cat("\nSaved time-series HDIs:", out_csv_ts_hdi, "\n")

# ==========================================================
# 2) ENDPOINT CHANGE HDIs (abs + %): draw-wise then summarise
# ==========================================================
endpoint_draws <- zone_year_draw_gm %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Zone, Year, draw_id, gm) %>%
  pivot_wider(names_from = Year, values_from = gm, names_prefix = "Y") %>%
  filter(
    is.finite(.data[[paste0("Y", start_year)]]),
    is.finite(.data[[paste0("Y", end_year)]])
  ) %>%
  mutate(
    s = pmax(.data[[paste0("Y", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Y", end_year)]],   1e-12),
    abs_change = f - s,
    pct_change = 100 * (f - s) / s
  ) %>%
  select(Zone, draw_id, abs_change, pct_change)

change_abs_hdi <- group_summaries_hdi(endpoint_draws, "Zone", "abs_change", credMass = credMass) %>%
  rename(abs_median = median, abs_hdi_low = hdi_low, abs_hdi_high = hdi_high)

change_pct_hdi <- group_summaries_hdi(endpoint_draws, "Zone", "pct_change", credMass = credMass) %>%
  rename(pct_median = median, pct_hdi_low = hdi_low, pct_hdi_high = hdi_high)

change_hdi <- left_join(change_abs_hdi, change_pct_hdi, by = "Zone") %>%
  arrange(Zone)

write.csv(change_hdi, out_csv_change_hdi, row.names = FALSE)
cat("Saved endpoint change HDIs:", out_csv_change_hdi, "\n")

# ==========================================================
# 3) GROWTH RATE HDIs (%/yr): you already compute zone_start_end_draws
# Uses zone_start_end_draws (Zone, draw_id, growth_rate_pct)
# ==========================================================
growth_hdi <- group_summaries_hdi(zone_start_end_draws, "Zone", "growth_rate_pct", credMass = credMass) %>%
  rename(growth_median = median, growth_hdi_low = hdi_low, growth_hdi_high = hdi_high) %>%
  arrange(Zone)

write.csv(growth_hdi, out_csv_growth_hdi, row.names = FALSE)
cat("Saved growth rate HDIs:", out_csv_growth_hdi, "\n")

# ==========================================================
# 4) HEDGE'S g (95% CI) ZONE CONTRASTS (using posterior draws)
#
# Here we treat the Monte Carlo draws of composite quantities
# as posterior samples of the underlying parameters, following
# Outhwaite et al. (2019) and Cooke et al. (2023).
# ==========================================================

zone_levels <- zones_keep
zone_pairs  <- combn(zone_levels, 2, simplify = FALSE)

# ---- (A) Composite % change (start → end): Hedge's g on pct_change draws ----
g_change_draws <- map_dfr(zone_pairs, function(cc) {
  a <- cc[1]; b <- cc[2]
  
  xa <- endpoint_draws %>%
    filter(Zone == a) %>%
    pull(pct_change)
  
  xb <- endpoint_draws %>%
    filter(Zone == b) %>%
    pull(pct_change)
  
  hedges_g_ci(xa, xb) %>%
    mutate(
      metric         = paste0("Composite % change ", start_year, "→", end_year),
      comparison     = paste0(a, " - ", b),
      replicate_unit = "Posterior draw",
      note           = "Uses Monte Carlo draws of % change between zones."
    )
})

# ---- (B) Annual growth rate (%/yr): Hedge's g on growth_rate_pct draws ----
g_growth_draws <- map_dfr(zone_pairs, function(cc) {
  a <- cc[1]; b <- cc[2]
  
  xa <- zone_start_end_draws %>%
    filter(Zone == a) %>%
    pull(growth_rate_pct)
  
  xb <- zone_start_end_draws %>%
    filter(Zone == b) %>%
    pull(growth_rate_pct)
  
  hedges_g_ci(xa, xb) %>%
    mutate(
      metric         = "Annual growth rate (%/yr)",
      comparison     = paste0(a, " - ", b),
      replicate_unit = "Posterior draw",
      note           = "Uses Monte Carlo draws of composite growth rate between zones."
    )
})

hedges_g_out <- bind_rows(g_change_draws, g_growth_draws) %>%
  arrange(metric, comparison)

write.csv(hedges_g_out, out_csv_hedges_g, row.names = FALSE)
cat("Saved Hedge's g contrasts:", out_csv_hedges_g, "\n")
























# Converged only dataset - Multispecies Trends ####

# this an add-on, run the above scripts first 

############################################################
# ADD-ON (FILTERED): Repeat Cooke-style trends + effect sizes
# for ONLY species that converged in BOTH start & end years
# (2000 AND 2023), then redo:
#   - composite geometric mean occupancy by Zone-Year (median + 95% CrI)
#   - endpoint change (abs + %)
#   - annual growth rate (%/yr) (median + 95% CrI)
#   - Cooke-style plots (occupancy + growth rate)
#   - HDIs + Hedge’s g exports (same as your “Absolute and Relative…” section)
#
# Output files are the same formats but suffixed with "_filtered".
############################################################

# ----------------------------
# FILTER SETTINGS
# ----------------------------
# Convergence rule (edit if your project uses something else)
rhat_max_filtered  <- 1.10
neff_min_filtered  <- 0  # conservative; raise if you want stricter filtering (e.g. 1000)

# Filter requires convergence in BOTH 2000 and 2023.
# By default, require it for ALL zones_keep in BOTH years.
require_all_zones_in_both_years <- TRUE

# Filtered outputs (avoid overwriting)
out_csv_clean_f   <- "1_BUGS_psi_fs_r_filteredFirstLast.csv"
out_csv_comp_f    <- "1_composite_geo_occ_filteredFirstLast.csv"
out_csv_growth_f  <- "1_growth_rate_filteredFirstLast.csv"

out_png_occ_f     <- "1_Cooke_style_composite_geo_occ_2000_2023_filtered.png"
out_jpg_occ_f     <- "1_Cooke_style_composite_geo_occ_2000_2023_filtered.jpg"
out_png_gr_f      <- "1_Cooke_style_growth_rate_2000_2023_filtered.png"
out_jpg_gr_f      <- "1_Cooke_style_growth_rate_2000_2023_filtered.jpg"

# HDI/effect-size outputs
out_csv_ts_hdi_f     <- "1_Trends_composite_geo_occ_byZone_byYear_median_HDI_filtered.csv"
out_csv_change_hdi_f <- "1_endpoint_change_byZone_median_HDI_filtered.csv"
out_csv_growth_hdi_f <- "1_growth_rate_byZone_median_HDI_filtered.csv"
out_csv_hedges_g_f   <- "1_hedges_g_zone_contrasts_filtered.csv"
out_xlsx_hdi_g_f     <- "1_Composite_HDI_and_HedgesG_outputs_filtered.xlsx"

cat("\n\n==================== FILTERED ADD-ON ====================\n")
cat("Convergence rule:\n")
cat("  Rhat <=", rhat_max_filtered, "\n")
cat("  n.eff >=", neff_min_filtered, "\n")
cat("  Years:", start_year, "and", end_year, "\n")
cat("  require_all_zones_in_both_years:", require_all_zones_in_both_years, "\n")

# ----------------------------
# 1) Identify species that converged in BOTH endpoint years
# ----------------------------
psi_zone_for_filter <- psi_zone %>%
  filter(
    Year %in% c(start_year, end_year),
    Zone %in% zones_keep
  ) %>%
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
        (n_zones_ok == length(zones_keep))
    ) %>%
    group_by(Species) %>%
    summarise(ok_both_years = all(ok_year), .groups = "drop") %>%
    filter(ok_both_years) %>%
    pull(Species)
} else {
  # Less strict: species must have at least one converged zone in both years
  species_keep_filtered <- psi_zone_for_filter %>%
    group_by(Species, Year) %>%
    summarise(any_ok = any(ok_conv), .groups = "drop") %>%
    group_by(Species) %>%
    summarise(ok_both_years = all(any_ok), .groups = "drop") %>%
    filter(ok_both_years) %>%
    pull(Species)
}

cat("\nFiltered species count:", length(species_keep_filtered), "\n")

# ----------------------------
# 2) Filter full psi_zone to those species, then save cleaned filtered table
# ----------------------------
psi_zone_f <- psi_zone %>%
  filter(Species %in% species_keep_filtered)

write.csv(psi_zone_f, out_csv_clean_f, row.names = FALSE)
cat("Saved FILTERED cleaned table:", out_csv_clean_f, "\n")

# ----------------------------
# 3) FILTERED Beta MC draws per species-zone-year (same process)
# ----------------------------
plot_df_f <- psi_zone_f %>%
  filter(!is.na(psi), !is.na(sd), is.finite(psi), is.finite(sd)) %>%
  mutate(
    psi_mean = pmin(pmax(psi, 1e-6), 1 - 1e-6),
    psi_sd   = pmax(sd, 1e-12)
  )

draws_long_f <- plot_df_f %>%
  rowwise() %>%
  mutate(draw = list(beta_draws_from_mean_sd(psi_mean, psi_sd, n = n_draws))) %>%
  ungroup() %>%
  select(Species, Zone, Year, draw) %>%
  unnest_longer(draw, values_to = "psi_draw") %>%
  group_by(Species, Zone, Year) %>%
  mutate(draw_id = row_number()) %>%
  ungroup() %>%
  mutate(psi_draw = pmin(pmax(psi_draw, 1e-12), 1 - 1e-12))

# ----------------------------
# 4) FILTERED Composite multispecies geometric mean by Zone-Year per draw
# ----------------------------
zone_year_draw_gm_f <- draws_long_f %>%
  group_by(Zone, Year, draw_id) %>%
  summarise(gm = exp(mean(log(psi_draw))), .groups = "drop")

zone_trend_f <- zone_year_draw_gm_f %>%
  group_by(Zone, Year) %>%
  summarise(
    gm_median = median(gm),
    gm_low    = as.numeric(quantile(gm, 0.025)),
    gm_high   = as.numeric(quantile(gm, 0.975)),
    .groups = "drop"
  ) %>%
  arrange(Zone, Year)

write.csv(zone_trend_f, out_csv_comp_f, row.names = FALSE)
cat("Saved FILTERED composite occupancy summaries:", out_csv_comp_f, "\n")

# ----------------------------
# 5) FILTERED Endpoint change (for annotation / reporting)
# ----------------------------
zone_change_f <- zone_trend_f %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Zone, Year, gm_median, gm_low, gm_high) %>%
  pivot_wider(
    names_from  = Year,
    values_from = c(gm_median, gm_low, gm_high),
    names_glue  = "{.value}_{Year}"
  ) %>%
  mutate(
    abs_change = .data[[paste0("gm_median_", end_year)]] - .data[[paste0("gm_median_", start_year)]],
    pct_change = 100 * abs_change / .data[[paste0("gm_median_", start_year)]]
  ) %>%
  arrange(Zone)

cat("\n==================== FILTERED Composite change (start -> end) ====================\n")
print(zone_change_f)

# ----------------------------
# 6) FILTERED labels (same approach)
# ----------------------------
plot_limits_f <- zone_trend_f %>%
  summarise(
    ymin = min(gm_low,  na.rm = TRUE),
    ymax = max(gm_high, na.rm = TRUE)
  ) %>%
  mutate(
    yrange = ymax - ymin,
    pad = pmax(0.005, 0.05 * yrange)
  )

zone_end_f <- zone_trend_f %>%
  filter(Year == end_year) %>%
  left_join(
    zone_change_f %>%
      transmute(
        Zone,
        start_val = .data[[paste0("gm_median_", start_year)]],
        end_val   = .data[[paste0("gm_median_", end_year)]],
        pct_change
      ),
    by = "Zone"
  ) %>%
  mutate(
    label_text = paste0(round(end_val, 3), " (", round(pct_change, 1), "%)")
  ) %>%
  tidyr::crossing(plot_limits_f) %>%
  mutate(
    off = pmax(0.01, 0.18 * yrange),
    label_base = case_when(
      Zone == "Core"    ~ gm_median - off,
      Zone == "Buffer"  ~ gm_median,
      Zone == "Outside" ~ gm_median + off,
      TRUE              ~ gm_median
    )
  ) %>%
  disperse_labels(gap = 0.10 * unique(plot_limits_f$yrange)) %>%
  mutate(
    x_label = end_year + 2.3,
    x_curve = end_year + 2.5
  )

# ----------------------------
# 7) FILTERED Plot 1: composite occupancy
# ----------------------------
g_occ_f <- ggplot(zone_trend_f, aes(x = Year, y = gm_median, color = Zone)) +
  geom_ribbon(aes(ymin = gm_low, ymax = gm_high, fill = Zone),
              alpha = 0.30, color = NA, show.legend = FALSE) +
  geom_line(aes(y = gm_high, group = Zone), linetype = "dashed",
            linewidth = 0.8, show.legend = FALSE) +
  geom_line(aes(y = gm_low, group = Zone), linetype = "dashed",
            linewidth = 0.8, show.legend = FALSE) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_fill_colors) +
  geom_vline(xintercept = end_year, linetype = "dashed",
             color = "grey30", linewidth = 0.8) +
  geom_point(
    data = zone_end_f,
    aes(x = end_year, y = gm_median, color = Zone),
    shape = 21, fill = "white", stroke = 1.5, size = 3.5,
    inherit.aes = FALSE
  ) +
  geom_curve(
    data = zone_end_f,
    aes(x = end_year, y = gm_median,
        xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.8, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = zone_end_f,
    aes(x = x_label, y = label_y_adj, label = label_text, color = Zone),
    hjust = 0, size = 6, fontface = "bold",
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  theme_minimal(base_size = 10) +
  labs(
    x = "Year",
    y = "Geometric mean occupancy (median)",
    color = "Zone",
    title = "Composite geometric mean occupancy (FILTERED species)"
  ) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(2.6, "lines"),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13, color = "gray30"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  ) +
  coord_cartesian(xlim = c(start_year, end_year + 9), ylim = c(0, 0.2))

print(g_occ_f)

ggsave(g_occ_f, filename = out_png_occ_f, width = 8, height = 6, dpi = 1000)
ggsave(g_occ_f, filename = out_jpg_occ_f, width = 8, height = 6, dpi = 1000)

cat("\nSaved FILTERED occupancy plots:\n  ", out_png_occ_f, "\n  ", out_jpg_occ_f, "\n")

# ----------------------------
# 8) FILTERED Growth rate (% per year) on composite draws
# ----------------------------
y_years <- end_year - start_year

zone_start_end_draws_f <- zone_year_draw_gm_f %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Zone, Year, draw_id, gm) %>%
  pivot_wider(names_from = Year, values_from = gm, names_prefix = "Year_") %>%
  filter(!is.na(.data[[paste0("Year_", start_year)]]),
         !is.na(.data[[paste0("Year_", end_year)]])) %>%
  mutate(
    s = pmax(.data[[paste0("Year_", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Year_", end_year)]], 1e-12),
    growth_rate_pct = 100 * ((f / s)^(1 / y_years) - 1)
  )

growth_summary_f <- zone_start_end_draws_f %>%
  group_by(Zone) %>%
  summarise(
    growth_median = median(growth_rate_pct),
    growth_low    = as.numeric(quantile(growth_rate_pct, 0.025)),
    growth_high   = as.numeric(quantile(growth_rate_pct, 0.975)),
    .groups = "drop"
  ) %>%
  arrange(Zone)

write.csv(growth_summary_f, out_csv_growth_f, row.names = FALSE)
cat("\nSaved FILTERED growth rate summary:", out_csv_growth_f, "\n")
cat("\n==================== FILTERED Growth rate (% per year) ====================\n")
print(growth_summary_f)

# ----------------------------
# 9) FILTERED Plot 2: growth rate by zone
# ----------------------------
g_gr_f <- ggplot(growth_summary_f, aes(x = Zone, y = growth_median, color = Zone)) +
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(ymin = growth_low, ymax = growth_high),
                  linewidth = 1.0, fatten = 2.5, show.legend = FALSE) +
  scale_color_manual(values = zone_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Annual growth rate of composite occupancy (FILTERED species; start → end)",
    subtitle = paste0("Growth rate = 100 * ((f/s)^(1/y) - 1), where y = ", y_years,
                      " years; points are medians, bars are 95% credible intervals"),
    x = "Zone", y = "Annual growth rate of occupancy (% per year)"
  ) +
  theme(
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13, color = "gray30")
  )

print(g_gr_f)

ggsave(g_gr_f, filename = out_png_gr_f, width = 10, height = 7, dpi = 1000)
ggsave(g_gr_f, filename = out_jpg_gr_f, width = 10, height = 7, dpi = 1000)

cat("\nSaved FILTERED growth rate plots:\n  ", out_png_gr_f, "\n  ", out_jpg_gr_f, "\n")

# ==========================================================
# FILTERED “Absolute and Relative effect sizes and CIs”
# - HDIs: time series, endpoint change, growth rate
# - Hedge’s g: (A) time-series uses yearly medians to avoid pseudo-replication
#              (B) growth rate uses posterior draws
# - Export Excel workbook
# ==========================================================

# ---- HDIs: composite gm by Zone-Year ----
ts_hdi_f <- group_summaries_hdi(zone_year_draw_gm_f, c("Zone", "Year"), "gm", credMass = credMass) %>%
  rename(
    gm_median = median,
    gm_hdi_low = hdi_low,
    gm_hdi_high = hdi_high
  ) %>%
  arrange(Zone, Year)

write.csv(ts_hdi_f, out_csv_ts_hdi_f, row.names = FALSE)
cat("\nSaved FILTERED time-series HDIs:", out_csv_ts_hdi_f, "\n")

# ---- Endpoint change HDIs (abs + %) ----
endpoint_draws_f <- zone_year_draw_gm_f %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Zone, Year, draw_id, gm) %>%
  pivot_wider(names_from = Year, values_from = gm, names_prefix = "Y") %>%
  filter(
    is.finite(.data[[paste0("Y", start_year)]]),
    is.finite(.data[[paste0("Y", end_year)]])
  ) %>%
  mutate(
    s = pmax(.data[[paste0("Y", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Y", end_year)]],   1e-12),
    abs_change = f - s,
    pct_change = 100 * (f - s) / s
  ) %>%
  select(Zone, draw_id, abs_change, pct_change)

change_abs_hdi_f <- group_summaries_hdi(endpoint_draws_f, "Zone", "abs_change", credMass = credMass) %>%
  rename(abs_median = median, abs_hdi_low = hdi_low, abs_hdi_high = hdi_high)

change_pct_hdi_f <- group_summaries_hdi(endpoint_draws_f, "Zone", "pct_change", credMass = credMass) %>%
  rename(pct_median = median, pct_hdi_low = hdi_low, pct_hdi_high = hdi_high)

change_hdi_f <- left_join(change_abs_hdi_f, change_pct_hdi_f, by = "Zone") %>%
  arrange(Zone)

write.csv(change_hdi_f, out_csv_change_hdi_f, row.names = FALSE)
cat("Saved FILTERED endpoint change HDIs:", out_csv_change_hdi_f, "\n")

# ---- Growth rate HDIs ----
growth_hdi_f <- group_summaries_hdi(zone_start_end_draws_f, "Zone", "growth_rate_pct", credMass = credMass) %>%
  rename(growth_median = median, growth_hdi_low = hdi_low, growth_hdi_high = hdi_high) %>%
  arrange(Zone)

write.csv(growth_hdi_f, out_csv_growth_hdi_f, row.names = FALSE)
cat("Saved FILTERED growth rate HDIs:", out_csv_growth_hdi_f, "\n")

# ---- Hedge’s g contrasts (same logic as your unfiltered section) ----
yearly_medians_f <- zone_trend_f %>%
  select(Zone, Year, gm_median) %>%
  filter(is.finite(gm_median))

zone_levels <- zones_keep
zone_pairs <- combn(zone_levels, 2, simplify = FALSE)

g_occ_years_f <- purrr::map_dfr(zone_pairs, function(cc) {
  a <- cc[1]; b <- cc[2]
  xa <- yearly_medians_f %>% filter(Zone == a) %>% pull(gm_median)
  xb <- yearly_medians_f %>% filter(Zone == b) %>% pull(gm_median)
  
  hedges_g_ci(xa, xb) %>%
    mutate(
      metric = "Composite occupancy (gm_median by Year) [FILTERED]",
      comparison = paste0(a, " - ", b),
      replicate_unit = "Year",
      note = "Uses yearly medians (one value per Zone-Year) to avoid pseudo-replication."
    )
})

g_growth_draws_f <- purrr::map_dfr(zone_pairs, function(cc) {
  a <- cc[1]; b <- cc[2]
  xa <- zone_start_end_draws_f %>% filter(Zone == a) %>% pull(growth_rate_pct)
  xb <- zone_start_end_draws_f %>% filter(Zone == b) %>% pull(growth_rate_pct)
  
  hedges_g_ci(xa, xb) %>%
    mutate(
      metric = "Annual growth rate (%/yr) [FILTERED]",
      comparison = paste0(a, " - ", b),
      replicate_unit = "Posterior draw",
      note = "Uses MC composite draws as posterior samples of growth rate."
    )
})

hedges_g_out_f <- bind_rows(g_occ_years_f, g_growth_draws_f) %>%
  arrange(metric, comparison)

write.csv(hedges_g_out_f, out_csv_hedges_g_f, row.names = FALSE)
cat("Saved FILTERED Hedge's g contrasts:", out_csv_hedges_g_f, "\n")

# ---- Export: Excel workbook (filtered) ----
wb_hdi_f <- createWorkbook()

write_sheet(
  wb_hdi_f, "TimeSeries_HDI_filtered",
  ts_hdi_f,
  header_text = paste0("FILTERED: Composite occupancy by Zone-Year: median + ", credMass*100, "% HDI (from MC draws).")
)

write_sheet(
  wb_hdi_f, "EndpointChange_HDI_filtered",
  change_hdi_f,
  header_text = paste0("FILTERED: Endpoint change ", start_year, "→", end_year, ": abs + % change; median + ", credMass*100, "% HDI.")
)

write_sheet(
  wb_hdi_f, "GrowthRate_HDI_filtered",
  growth_hdi_f,
  header_text = paste0("FILTERED: Annual growth rate (%/yr): median + ", credMass*100, "% HDI (computed on composite draws).")
)

write_sheet(
  wb_hdi_f, "HedgesG_ZoneContrasts_filtered",
  hedges_g_out_f,
  header_text = "FILTERED: Hedge's g (95% CI) for Zone contrasts. See 'replicate_unit' + 'note' for interpretation."
)

saveWorkbook(wb_hdi_f, out_xlsx_hdi_g_f, overwrite = TRUE)
cat("\nSaved FILTERED Excel workbook:\n  ", out_xlsx_hdi_g_f, "\n", sep = "")

cat("\nFILTERED ADD-ON COMPLETE.\n")






# ==========================================================
# QUICK SUMMARY: FILTERED (First + Last year only)
# - counts retained species (converged in 2000 AND 2023)
# - convergence rates *within the retained subset* for:
#     (A) endpoint rows (2000 + 2023 only)
#     (B) endpoint Species-Year strict (all zones ok in each endpoint year)
# - also reports how many retained species have complete coverage across all years
#     (this is NOT required by the endpoint filter, but is useful context)
# ==========================================================

cat("\n==================== FILTERED SUBSET SUMMARY (FIRST+LAST YEARS) ====================\n")

# Basic counts
n_species_total_fl <- dplyr::n_distinct(psi_zone$Species)
n_species_kept_fl  <- length(species_keep_filtered)

cat("Total species in psi_zone:", n_species_total_fl, "\n")
cat("Retained species (converged in ", start_year, " AND ", end_year, "): ",
    n_species_kept_fl, "\n", sep = "")
cat("Retained proportion:", round(100 * n_species_kept_fl / n_species_total_fl, 2), "%\n")

# Rebuild the same ok_conv flag used in the endpoint filter (includes n.eff here)
psi_zone_filter_scope_fl <- psi_zone %>%
  dplyr::filter(
    Year %in% c(start_year, end_year),
    Zone %in% zones_keep
  ) %>%
  dplyr::mutate(
    Rhat_num = suppressWarnings(as.numeric(Rhat)),
    neff_num = suppressWarnings(as.numeric(n.eff)),
    ok_conv  = is.finite(Rhat_num) & (Rhat_num <= rhat_max_filtered) &
      is.finite(neff_num) & (neff_num >= neff_min_filtered) &
      is.finite(psi) & is.finite(sd) & !is.na(psi) & !is.na(sd)
  )

# Retained subset only (endpoints)
psi_zone_kept_scope_fl <- psi_zone_filter_scope_fl %>%
  dplyr::filter(Species %in% species_keep_filtered)

# (A) Row-level convergence rate across endpoint Zone-Year rows
row_conv_rate_fl <- mean(psi_zone_kept_scope_fl$ok_conv, na.rm = TRUE)

cat("\n--- Convergence rates within RETAINED subset (endpoints only) ---\n")
cat("Row-level convergence rate (endpoint Zone-Year rows ok):",
    round(100 * row_conv_rate_fl, 2), "%\n")

# (B) Species-Year strict at endpoints: all zones present + ok in that year
species_year_strict_fl <- psi_zone_kept_scope_fl %>%
  dplyr::group_by(Species, Year) %>%
  dplyr::summarise(
    n_zones_present = dplyr::n_distinct(Zone),
    n_zones_ok      = dplyr::n_distinct(Zone[ok_conv]),
    ok_year_strict  = (n_zones_present == length(zones_keep)) &
      (n_zones_ok == length(zones_keep)),
    .groups = "drop"
  )

species_year_strict_rate_fl <- mean(species_year_strict_fl$ok_year_strict, na.rm = TRUE)
cat("Species-Year strict convergence rate at endpoints (all zones ok):",
    round(100 * species_year_strict_rate_fl, 2), "%\n")

# Optional: breakdown of any failures within endpoint Species-Year strictness
fail_breakdown_fl <- species_year_strict_fl %>%
  dplyr::mutate(
    reason = dplyr::case_when(
      n_zones_present < length(zones_keep) ~ "Missing zone(s)",
      n_zones_ok      < length(zones_keep) ~ "Rhat/n.eff/NA fail in >=1 zone",
      TRUE ~ "OK"
    )
  ) %>%
  dplyr::count(reason) %>%
  dplyr::mutate(pct = 100 * n / sum(n))

cat("\nEndpoint Species-Year strict status breakdown (retained subset):\n")
print(fail_breakdown_fl)

# Extra context (NOT a requirement of endpoint filtering):
# how many retained species have complete Year coverage across the full range?
species_year_coverage_fl <- psi_zone %>%
  dplyr::filter(Zone %in% zones_keep, Species %in% species_keep_filtered,
                Year >= start_year, Year <= end_year) %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    n_years_present = dplyr::n_distinct(Year),
    complete_years  = (n_years_present == length(seq(start_year, end_year))),
    .groups = "drop"
  )

cat("\nRetained species with complete Year coverage (", start_year, "–", end_year, "): ",
    sum(species_year_coverage_fl$complete_years),
    "/", n_species_kept_fl, "\n", sep = "")

cat("\n==============================================================\n")




















# Filtered to ONLY FULLY CONVERGED SPECIES in all years ####

############################################################
# ADD-ON (FULLY CONVERGED): Repeat Cooke-style trends + effect sizes
# for ONLY species that converged in EVERY year (2000–2023).
#
# NOTE: This section is self-contained in terms of outputs:
#   - All plot + table outputs start with "1_" and end with "FULLYCONVERGED"
#   - Filtering uses ONLY Rhat threshold (n.eff ignored; set to 0)
############################################################

library(openxlsx)

# helper: add a sheet with an optional header text and a data frame
write_sheet <- function(wb, sheet_name, data, header_text = NULL) {
  # create sheet
  openxlsx::addWorksheet(wb, sheetName = sheet_name)
  
  row <- 1
  
  # optional header text at the top
  if (!is.null(header_text)) {
    openxlsx::writeData(
      wb, sheet = sheet_name,
      x = header_text,
      startRow = row, startCol = 1
    )
    row <- row + 2  # blank line before the table
  }
  
  # write the data frame
  openxlsx::writeData(
    wb, sheet = sheet_name,
    x = data,
    startRow = row, startCol = 1
  )
}


# ----------------------------
# FULLY CONVERGED FILTER SETTINGS
# ----------------------------
rhat_max_filtered  <- 1.10
neff_min_filtered  <- 0     # ignored effectively; kept for compatibility with your prints
require_all_zones_each_year <- TRUE  # strict: all zones_keep must converge each year

# ----------------------------
# OUTPUTS: prefix "1_" + suffix "FULLYCONVERGED"
# ----------------------------
suffix_fc <- "FULLYCONVERGED"

out_csv_clean_fc   <- paste0("1_BUGS_psi_fs_r_", suffix_fc, ".csv")
out_csv_comp_fc    <- paste0("1_composite_geo_occ_", suffix_fc, ".csv")
out_csv_growth_fc  <- paste0("1_growth_rate_", suffix_fc, ".csv")

out_png_occ_fc     <- paste0("1_Cooke_style_composite_geo_occ_2000_2023_", suffix_fc, ".png")
out_jpg_occ_fc     <- paste0("1_Cooke_style_composite_geo_occ_2000_2023_", suffix_fc, ".jpg")
out_png_gr_fc      <- paste0("1_Cooke_style_growth_rate_2000_2023_", suffix_fc, ".png")
out_jpg_gr_fc      <- paste0("1_Cooke_style_growth_rate_2000_2023_", suffix_fc, ".jpg")

out_csv_ts_hdi_fc     <- paste0("1_Trends_composite_geo_occ_byZone_byYear_median_HDI_", suffix_fc, ".csv")
out_csv_change_hdi_fc <- paste0("1_endpoint_change_byZone_median_HDI_", suffix_fc, ".csv")
out_csv_growth_hdi_fc <- paste0("1_growth_rate_byZone_median_HDI_", suffix_fc, ".csv")
out_csv_hedges_g_fc   <- paste0("1_hedges_g_zone_contrasts_", suffix_fc, ".csv")
out_xlsx_hdi_g_fc     <- paste0("1_Composite_HDI_and_HedgesG_outputs_", suffix_fc, ".xlsx")

cat("\n\n==================== ADD-ON: FULLY CONVERGED (EVERY YEAR) ====================\n")
cat("Convergence rule:\n")
cat("  Rhat <=", rhat_max_filtered, "\n")
cat("  n.eff >=", neff_min_filtered, " (ignored)\n")
cat("  Years:", start_year, "to", end_year, "\n")
cat("  require_all_zones_each_year:", require_all_zones_each_year, "\n")

# ----------------------------
# 1) Identify species that converged in EVERY year
# ----------------------------
psi_zone_for_filter_fc <- psi_zone %>%
  dplyr::filter(
    Year >= start_year, Year <= end_year,
    Zone %in% zones_keep
  ) %>%
  dplyr::mutate(
    Rhat_num = suppressWarnings(as.numeric(Rhat)),
    neff_num = suppressWarnings(as.numeric(n.eff)),
    ok_conv  = is.finite(Rhat_num) & (Rhat_num <= rhat_max_filtered) &
      is.finite(psi) & is.finite(sd) & !is.na(psi) & !is.na(sd)
  )

years_expected_fc <- seq(start_year, end_year)

if (require_all_zones_each_year) {
  
  species_keep_fullyconverged <- psi_zone_for_filter_fc %>%
    dplyr::group_by(Species, Year) %>%
    dplyr::summarise(
      n_zones_present = dplyr::n_distinct(Zone),
      n_zones_ok      = dplyr::n_distinct(Zone[ok_conv]),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      ok_year = (n_zones_present == length(zones_keep)) &
        (n_zones_ok == length(zones_keep))
    ) %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(
      n_years_present = dplyr::n_distinct(Year),
      all_years_ok    = all(ok_year) & (n_years_present == length(years_expected_fc)),
      .groups = "drop"
    ) %>%
    dplyr::filter(all_years_ok) %>%
    dplyr::pull(Species)
  
} else {
  
  species_keep_fullyconverged <- psi_zone_for_filter_fc %>%
    dplyr::group_by(Species, Year) %>%
    dplyr::summarise(any_ok = any(ok_conv), .groups = "drop") %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(
      n_years_present = dplyr::n_distinct(Year),
      all_years_ok    = all(any_ok) & (n_years_present == length(years_expected_fc)),
      .groups = "drop"
    ) %>%
    dplyr::filter(all_years_ok) %>%
    dplyr::pull(Species)
  
}

cat("\nRetained species count (FULLYCONVERGED):", length(species_keep_fullyconverged), "\n")

# ----------------------------
# 2) Filter psi_zone to fully converged species + save
# ----------------------------
psi_zone_fc <- psi_zone %>%
  dplyr::filter(Species %in% species_keep_fullyconverged)

write.csv(psi_zone_fc, out_csv_clean_fc, row.names = FALSE)
cat("Saved CLEAN (FULLYCONVERGED):", out_csv_clean_fc, "\n")

# ----------------------------
# 3) Beta MC draws per species-zone-year (same process)
# ----------------------------
plot_df_fc <- psi_zone_fc %>%
  dplyr::filter(!is.na(psi), !is.na(sd), is.finite(psi), is.finite(sd)) %>%
  dplyr::mutate(
    psi_mean = pmin(pmax(psi, 1e-6), 1 - 1e-6),
    psi_sd   = pmax(sd, 1e-12)
  )

draws_long_fc <- plot_df_fc %>%
  dplyr::rowwise() %>%
  dplyr::mutate(draw = list(beta_draws_from_mean_sd(psi_mean, psi_sd, n = n_draws))) %>%
  dplyr::ungroup() %>%
  dplyr::select(Species, Zone, Year, draw) %>%
  tidyr::unnest_longer(draw, values_to = "psi_draw") %>%
  dplyr::group_by(Species, Zone, Year) %>%
  dplyr::mutate(draw_id = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(psi_draw = pmin(pmax(psi_draw, 1e-12), 1 - 1e-12))

# ----------------------------
# 4) Composite geometric mean by Zone-Year per draw
# ----------------------------
zone_year_draw_gm_fc <- draws_long_fc %>%
  dplyr::group_by(Zone, Year, draw_id) %>%
  dplyr::summarise(gm = exp(mean(log(psi_draw))), .groups = "drop")

zone_trend_fc <- zone_year_draw_gm_fc %>%
  dplyr::group_by(Zone, Year) %>%
  dplyr::summarise(
    gm_median = median(gm),
    gm_low    = as.numeric(quantile(gm, 0.025)),
    gm_high   = as.numeric(quantile(gm, 0.975)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Zone, Year)

write.csv(zone_trend_fc, out_csv_comp_fc, row.names = FALSE)
cat("Saved composite occupancy (FULLYCONVERGED):", out_csv_comp_fc, "\n")

# ----------------------------
# 5) Endpoint change (for annotation / reporting)
# ----------------------------
zone_change_fc <- zone_trend_fc %>%
  dplyr::filter(Year %in% c(start_year, end_year)) %>%
  dplyr::select(Zone, Year, gm_median, gm_low, gm_high) %>%
  tidyr::pivot_wider(
    names_from  = Year,
    values_from = c(gm_median, gm_low, gm_high),
    names_glue  = "{.value}_{Year}"
  ) %>%
  dplyr::mutate(
    abs_change = .data[[paste0("gm_median_", end_year)]] - .data[[paste0("gm_median_", start_year)]],
    pct_change = 100 * abs_change / .data[[paste0("gm_median_", start_year)]]
  ) %>%
  dplyr::arrange(Zone)

cat("\n==================== FULLYCONVERGED change (start -> end) ====================\n")
print(zone_change_fc)

# ----------------------------
# 6) Endpoint label placement (same approach)
# ----------------------------
plot_limits_fc <- zone_trend_fc %>%
  dplyr::summarise(
    ymin = min(gm_low,  na.rm = TRUE),
    ymax = max(gm_high, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    yrange = ymax - ymin,
    pad = pmax(0.005, 0.05 * yrange)
  )

zone_end_fc <- zone_trend_fc %>%
  dplyr::filter(Year == end_year) %>%
  dplyr::left_join(
    zone_change_fc %>%
      dplyr::transmute(
        Zone,
        start_val = .data[[paste0("gm_median_", start_year)]],
        end_val   = .data[[paste0("gm_median_", end_year)]],
        pct_change
      ),
    by = "Zone"
  ) %>%
  dplyr::mutate(
    label_text = paste0(round(end_val, 3), " (", round(pct_change, 1), "%)")
  ) %>%
  tidyr::crossing(plot_limits_fc) %>%
  dplyr::mutate(
    off = pmax(0.01, 0.18 * yrange),
    label_base = dplyr::case_when(
      Zone == "Core"    ~ gm_median - off,
      Zone == "Buffer"  ~ gm_median,
      Zone == "Outside" ~ gm_median + off,
      TRUE              ~ gm_median
    )
  ) %>%
  disperse_labels(gap = 0.10 * unique(plot_limits_fc$yrange)) %>%
  dplyr::mutate(
    x_label = end_year + 2.3,
    x_curve = end_year + 2.5
  )

# ----------------------------
# 7) Plot 1: composite occupancy (FULLYCONVERGED)
# ----------------------------
g_occ_fc <- ggplot(zone_trend_fc, aes(x = Year, y = gm_median, color = Zone)) +
  geom_ribbon(aes(ymin = gm_low, ymax = gm_high, fill = Zone),
              alpha = 0.30, color = NA, show.legend = FALSE) +
  geom_line(aes(y = gm_high, group = Zone), linetype = "dashed",
            linewidth = 0.8, show.legend = FALSE) +
  geom_line(aes(y = gm_low, group = Zone), linetype = "dashed",
            linewidth = 0.8, show.legend = FALSE) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_fill_colors) +
  geom_vline(xintercept = end_year, linetype = "dashed",
             color = "grey30", linewidth = 0.8) +
  geom_point(
    data = zone_end_fc,
    aes(x = end_year, y = gm_median, color = Zone),
    shape = 21, fill = "white", stroke = 1.5, size = 3.5,
    inherit.aes = FALSE
  ) +
  geom_curve(
    data = zone_end_fc,
    aes(x = end_year, y = gm_median,
        xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.8, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = zone_end_fc,
    aes(x = x_label, y = label_y_adj, label = label_text, color = Zone),
    hjust = 0, size = 6, fontface = "bold",
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  theme_minimal(base_size = 10) +
  labs(
    x = "Year",
    y = "Geometric mean occupancy",
    color = "Zone",
  ) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(2.6, "lines"),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13, color = "gray30"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  ) +
  coord_cartesian(xlim = c(start_year, end_year + 9), ylim = c(0.05, 0.3))

print(g_occ_fc)

ggsave(g_occ_fc, filename = out_png_occ_fc, width = 8, height = 6, dpi = 1000)
ggsave(g_occ_fc, filename = out_jpg_occ_fc, width = 8, height = 6, dpi = 1000)

cat("\nSaved occupancy plots (FULLYCONVERGED):\n  ",
    out_png_occ_fc, "\n  ", out_jpg_occ_fc, "\n")

# ----------------------------
# 8) Growth rate (%/yr) on composite draws (FULLYCONVERGED)
# ----------------------------
y_years <- end_year - start_year

zone_start_end_draws_fc <- zone_year_draw_gm_fc %>%
  dplyr::filter(Year %in% c(start_year, end_year)) %>%
  dplyr::select(Zone, Year, draw_id, gm) %>%
  tidyr::pivot_wider(names_from = Year, values_from = gm, names_prefix = "Year_") %>%
  dplyr::filter(!is.na(.data[[paste0("Year_", start_year)]]),
                !is.na(.data[[paste0("Year_", end_year)]])) %>%
  dplyr::mutate(
    s = pmax(.data[[paste0("Year_", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Year_", end_year)]], 1e-12),
    growth_rate_pct = 100 * ((f / s)^(1 / y_years) - 1)
  )

growth_summary_fc <- zone_start_end_draws_fc %>%
  dplyr::group_by(Zone) %>%
  dplyr::summarise(
    growth_median = median(growth_rate_pct),
    growth_low    = as.numeric(quantile(growth_rate_pct, 0.025)),
    growth_high   = as.numeric(quantile(growth_rate_pct, 0.975)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Zone)

write.csv(growth_summary_fc, out_csv_growth_fc, row.names = FALSE)
cat("\nSaved growth rate summary (FULLYCONVERGED):", out_csv_growth_fc, "\n")
cat("\n==================== FULLYCONVERGED Growth rate (%/yr) ====================\n")
print(growth_summary_fc)

# ----------------------------
# 9) Plot 2: growth rate by zone (FULLYCONVERGED)
# ----------------------------
g_gr_fc <- ggplot(growth_summary_fc, aes(x = Zone, y = growth_median, color = Zone)) +
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(ymin = growth_low, ymax = growth_high),
                  linewidth = 1.0, fatten = 2.5, show.legend = FALSE) +
  scale_color_manual(values = zone_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Annual growth rate of composite occupancy (FULLYCONVERGED species; start → end)",
    subtitle = paste0("Growth rate = 100 * ((f/s)^(1/y) - 1), where y = ", y_years,
                      " years; points are medians, bars are 95% credible intervals"),
    x = "Zone", y = "Annual growth rate of occupancy (% per year)"
  ) +
  theme(
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13, color = "gray30")
  )

print(g_gr_fc)

ggsave(g_gr_fc, filename = out_png_gr_fc, width = 10, height = 7, dpi = 1000)
ggsave(g_gr_fc, filename = out_jpg_gr_fc, width = 10, height = 7, dpi = 1000)

cat("\nSaved growth rate plots (FULLYCONVERGED):\n  ",
    out_png_gr_fc, "\n  ", out_jpg_gr_fc, "\n")

# ==========================================================
# Effect sizes + intervals (FULLYCONVERGED)
# ==========================================================

# ---- HDIs: composite gm by Zone-Year ----
ts_hdi_fc <- group_summaries_hdi(zone_year_draw_gm_fc, c("Zone", "Year"), "gm", credMass = credMass) %>%
  dplyr::rename(
    gm_median = median,
    gm_hdi_low = hdi_low,
    gm_hdi_high = hdi_high
  ) %>%
  dplyr::arrange(Zone, Year)

write.csv(ts_hdi_fc, out_csv_ts_hdi_fc, row.names = FALSE)
cat("\nSaved time-series HDIs (FULLYCONVERGED):", out_csv_ts_hdi_fc, "\n")

# ---- Endpoint change HDIs (abs + %) ----
endpoint_draws_fc <- zone_year_draw_gm_fc %>%
  dplyr::filter(Year %in% c(start_year, end_year)) %>%
  dplyr::select(Zone, Year, draw_id, gm) %>%
  tidyr::pivot_wider(names_from = Year, values_from = gm, names_prefix = "Y") %>%
  dplyr::filter(
    is.finite(.data[[paste0("Y", start_year)]]),
    is.finite(.data[[paste0("Y", end_year)]])
  ) %>%
  dplyr::mutate(
    s = pmax(.data[[paste0("Y", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Y", end_year)]],   1e-12),
    abs_change = f - s,
    pct_change = 100 * (f - s) / s
  ) %>%
  dplyr::select(Zone, draw_id, abs_change, pct_change)

change_abs_hdi_fc <- group_summaries_hdi(endpoint_draws_fc, "Zone", "abs_change", credMass = credMass) %>%
  dplyr::rename(abs_median = median, abs_hdi_low = hdi_low, abs_hdi_high = hdi_high)

change_pct_hdi_fc <- group_summaries_hdi(endpoint_draws_fc, "Zone", "pct_change", credMass = credMass) %>%
  dplyr::rename(pct_median = median, pct_hdi_low = hdi_low, pct_hdi_high = hdi_high)

change_hdi_fc <- dplyr::left_join(change_abs_hdi_fc, change_pct_hdi_fc, by = "Zone") %>%
  dplyr::arrange(Zone)

write.csv(change_hdi_fc, out_csv_change_hdi_fc, row.names = FALSE)
cat("Saved endpoint change HDIs (FULLYCONVERGED):", out_csv_change_hdi_fc, "\n")

# ---- Growth rate HDIs ----
growth_hdi_fc <- group_summaries_hdi(zone_start_end_draws_fc, "Zone", "growth_rate_pct", credMass = credMass) %>%
  dplyr::rename(growth_median = median, growth_hdi_low = hdi_low, growth_hdi_high = hdi_high) %>%
  dplyr::arrange(Zone)

write.csv(growth_hdi_fc, out_csv_growth_hdi_fc, row.names = FALSE)
cat("Saved growth rate HDIs (FULLYCONVERGED):", out_csv_growth_hdi_fc, "\n")

# ---- Hedge’s g contrasts ----
yearly_medians_fc <- zone_trend_fc %>%
  dplyr::select(Zone, Year, gm_median) %>%
  dplyr::filter(is.finite(gm_median))

zone_levels <- zones_keep
zone_pairs <- combn(zone_levels, 2, simplify = FALSE)

g_occ_years_fc <- purrr::map_dfr(zone_pairs, function(cc) {
  a <- cc[1]; b <- cc[2]
  xa <- yearly_medians_fc %>% dplyr::filter(Zone == a) %>% dplyr::pull(gm_median)
  xb <- yearly_medians_fc %>% dplyr::filter(Zone == b) %>% dplyr::pull(gm_median)
  
  hedges_g_ci(xa, xb) %>%
    dplyr::mutate(
      metric = "Composite occupancy (gm_median by Year) [FULLYCONVERGED]",
      comparison = paste0(a, " - ", b),
      replicate_unit = "Year",
      note = "Uses yearly medians (one value per Zone-Year) to avoid pseudo-replication."
    )
})

g_growth_draws_fc <- purrr::map_dfr(zone_pairs, function(cc) {
  a <- cc[1]; b <- cc[2]
  xa <- zone_start_end_draws_fc %>% dplyr::filter(Zone == a) %>% dplyr::pull(growth_rate_pct)
  xb <- zone_start_end_draws_fc %>% dplyr::filter(Zone == b) %>% dplyr::pull(growth_rate_pct)
  
  hedges_g_ci(xa, xb) %>%
    dplyr::mutate(
      metric = "Annual growth rate (%/yr) [FULLYCONVERGED]",
      comparison = paste0(a, " - ", b),
      replicate_unit = "Posterior draw",
      note = "Uses MC composite draws as posterior samples of growth rate."
    )
})

hedges_g_out_fc <- dplyr::bind_rows(g_occ_years_fc, g_growth_draws_fc) %>%
  dplyr::arrange(metric, comparison)

write.csv(hedges_g_out_fc, out_csv_hedges_g_fc, row.names = FALSE)
cat("Saved Hedge's g contrasts (FULLYCONVERGED):", out_csv_hedges_g_fc, "\n")

# ---- Export: Excel workbook (FULLYCONVERGED) ----
wb_hdi_fc <- createWorkbook()

# Excel worksheet names must be <= 31 chars
sheet_ts_fc     <- "TS_HDI_FULLYCONVERGED"        # 19
sheet_change_fc <- "EndptChg_HDI_FULLYCONVERGED"  # 25
sheet_growth_fc <- "Growth_HDI_FULLYCONVERGED"    # 22
sheet_g_fc      <- "HedgesG_FULLYCONVERGED"       # 20

write_sheet(
  wb_hdi_fc, sheet_ts_fc,
  ts_hdi_fc,
  header_text = paste0(
    "FULLYCONVERGED: Composite occupancy by Zone-Year: median + ",
    credMass*100, "% HDI (from MC draws)."
  )
)

write_sheet(
  wb_hdi_fc, sheet_change_fc,
  change_hdi_fc,
  header_text = paste0(
    "FULLYCONVERGED: Endpoint change ", start_year, "→", end_year,
    ": abs + % change; median + ", credMass*100, "% HDI."
  )
)

write_sheet(
  wb_hdi_fc, sheet_growth_fc,
  growth_hdi_fc,
  header_text = paste0(
    "FULLYCONVERGED: Annual growth rate (%/yr): median + ",
    credMass*100, "% HDI (computed on composite draws)."
  )
)

write_sheet(
  wb_hdi_fc, sheet_g_fc,
  hedges_g_out_fc,
  header_text = "FULLYCONVERGED: Hedge's g (95% CI) for Zone contrasts."
)

saveWorkbook(wb_hdi_fc, out_xlsx_hdi_g_fc, overwrite = TRUE)
cat("\nSaved Excel workbook (FULLYCONVERGED):\n  ", out_xlsx_hdi_g_fc, "\n", sep = "")











# ==========================================================
# QUICK SUMMARY (FULLYCONVERGED): Retained species + convergence rates
# Uses species_keep_fullyconverged (the ALL-years, ALL-zones filter)
# ==========================================================

cat("\n==================== FULLYCONVERGED SUBSET SUMMARY ====================\n")

n_species_total_fc <- dplyr::n_distinct(psi_zone$Species)
n_species_kept_fc  <- length(species_keep_fullyconverged)

cat("Total species in psi_zone:", n_species_total_fc, "\n")
cat("Retained species (FULLYCONVERGED):", n_species_kept_fc, "\n")
cat("Retained proportion:", round(100 * n_species_kept_fc / n_species_total_fc, 2), "%\n")

# Recreate the same ok_conv flag used in filtering
psi_zone_filter_scope_fc <- psi_zone %>%
  dplyr::filter(
    Year >= start_year, Year <= end_year,
    Zone %in% zones_keep
  ) %>%
  dplyr::mutate(
    Rhat_num = suppressWarnings(as.numeric(Rhat)),
    ok_conv  = is.finite(Rhat_num) & (Rhat_num <= rhat_max_filtered) &
      is.finite(psi) & is.finite(sd) & !is.na(psi) & !is.na(sd)
  )

psi_zone_kept_scope_fc <- psi_zone_filter_scope_fc %>%
  dplyr::filter(Species %in% species_keep_fullyconverged)

# Row-level convergence rate (Zone-Year rows ok) WITHIN retained subset
row_conv_rate_fc <- mean(psi_zone_kept_scope_fc$ok_conv, na.rm = TRUE)

cat("\n--- Convergence rates within FULLYCONVERGED subset ---\n")
cat("Row-level convergence rate (Zone-Year rows ok):",
    round(100 * row_conv_rate_fc, 2), "%\n")

# Species-Year strict: all zones ok in each year
species_year_strict_fc <- psi_zone_kept_scope_fc %>%
  dplyr::group_by(Species, Year) %>%
  dplyr::summarise(
    n_zones_present = dplyr::n_distinct(Zone),
    n_zones_ok      = dplyr::n_distinct(Zone[ok_conv]),
    ok_year_strict  = (n_zones_present == length(zones_keep)) &
      (n_zones_ok == length(zones_keep)),
    .groups = "drop"
  )

species_year_strict_rate_fc <- mean(species_year_strict_fc$ok_year_strict, na.rm = TRUE)
cat("Species-Year strict convergence rate (all zones ok):",
    round(100 * species_year_strict_rate_fc, 2), "%\n")

# Completeness check
years_expected_fc <- seq(start_year, end_year)
species_year_completeness_fc <- species_year_strict_fc %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    n_years_present = dplyr::n_distinct(Year),
    complete_years  = (n_years_present == length(years_expected_fc)),
    .groups = "drop"
  )

cat("Retained species with complete Year coverage:",
    sum(species_year_completeness_fc$complete_years),
    "/", n_species_kept_fc, "\n")

# Optional breakdown
fail_breakdown_fc <- species_year_strict_fc %>%
  dplyr::mutate(
    reason = dplyr::case_when(
      n_zones_present < length(zones_keep) ~ "Missing zone(s)",
      n_zones_ok      < length(zones_keep) ~ "Rhat/NA fail in >=1 zone",
      TRUE ~ "OK"
    )
  ) %>%
  dplyr::count(reason) %>%
  dplyr::mutate(pct = 100 * n / sum(n))

cat("\nSpecies-Year strict status breakdown (FULLYCONVERGED subset):\n")
print(fail_breakdown_fc)

cat("\n==============================================================\n")








# ==========================================================
# ADD-ON: Rarity groups + removal summary (FULLYCONVERGED filter)
# Uses the SAME rarity definition as your Rare/Common script:
#   Rare/Common based on median of species mean psi (from psi_zone_full / psi_zone)
#
# Requires objects already in your workflow:
#   - psi_zone (or psi_zone_full) containing Species + psi
#   - species_keep_fullyconverged (vector from FULLYCONVERGED filter)
#   - start_year, end_year, zones_keep, rhat_max_filtered, require_all_zones_each_year
# ==========================================================

cat("\n==================== FULLYCONVERGED SUBSET SUMMARY (WITH RARITY) ====================\n")

# ----------------------------
# 0) Build Rare/Common labels (same method as your master script)
# ----------------------------
# Choose psi source: if psi_zone_full exists (from master script), use it.
psi_for_rarity <- if (exists("psi_zone_full")) psi_zone_full else psi_zone

species_occ_all <- psi_for_rarity %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(mean_psi = mean(psi, na.rm = TRUE), .groups = "drop")

rarity_threshold_all <- median(species_occ_all$mean_psi, na.rm = TRUE)

rarity_lookup_all <- species_occ_all %>%
  dplyr::mutate(
    Rarity = ifelse(mean_psi <= rarity_threshold_all, "Rare", "Common"),
    Rarity = factor(Rarity, levels = c("Rare", "Common"))
  ) %>%
  dplyr::select(Species, Rarity)

# ----------------------------
# 1) Overall counts + kept/removed by rarity
# ----------------------------
all_species_vec <- species_occ_all$Species
kept_species_vec <- species_keep_fullyconverged
removed_species_vec <- setdiff(all_species_vec, kept_species_vec)

rarity_counts_all <- rarity_lookup_all %>%
  dplyr::count(Rarity, name = "n_total") %>%
  tidyr::complete(Rarity = factor(c("Rare","Common"), levels = c("Rare","Common")),
                  fill = list(n_total = 0))

rarity_counts_kept <- rarity_lookup_all %>%
  dplyr::filter(Species %in% kept_species_vec) %>%
  dplyr::count(Rarity, name = "n_kept") %>%
  tidyr::complete(Rarity = factor(c("Rare","Common"), levels = c("Rare","Common")),
                  fill = list(n_kept = 0))

rarity_counts_removed <- rarity_lookup_all %>%
  dplyr::filter(Species %in% removed_species_vec) %>%
  dplyr::count(Rarity, name = "n_removed") %>%
  tidyr::complete(Rarity = factor(c("Rare","Common"), levels = c("Rare","Common")),
                  fill = list(n_removed = 0))

rarity_summary <- rarity_counts_all %>%
  dplyr::left_join(rarity_counts_kept, by = "Rarity") %>%
  dplyr::left_join(rarity_counts_removed, by = "Rarity") %>%
  dplyr::mutate(
    pct_removed_within_rarity = 100 * n_removed / pmax(n_total, 1),
    pct_kept_within_rarity    = 100 * n_kept / pmax(n_total, 1)
  )

cat("Rarity threshold (median of species mean psi): ", signif(rarity_threshold_all, 4), "\n", sep = "")
cat("Total species (rarity-classified): ", length(all_species_vec), "\n", sep = "")
cat("Kept species (FULLYCONVERGED): ", length(kept_species_vec), "\n", sep = "")
cat("Removed species: ", length(removed_species_vec), "\n\n", sep = "")

cat("--- Rare/Common kept vs removed (counts + % within rarity) ---\n")
print(rarity_summary)

# ----------------------------
# 2) Existing FULLYCONVERGED quick summary (now with correct variables)
#    (counts + convergence rates within retained subset)
# ----------------------------

# Basic counts
n_species_total_fc <- dplyr::n_distinct(psi_zone$Species)
n_species_kept_fc  <- length(species_keep_fullyconverged)

cat("\n--- FULLYCONVERGED filter counts ---\n")
cat("Total species in psi_zone:", n_species_total_fc, "\n")
cat("Retained species (fully converged all years):", n_species_kept_fc, "\n")
cat("Retained proportion:", round(100 * n_species_kept_fc / n_species_total_fc, 2), "%\n")

# Build the ok_conv flag used in FULLYCONVERGED filtering scope (ALL years + zones_keep)
psi_zone_filter_scope_fc <- psi_zone %>%
  dplyr::filter(
    Year >= start_year, Year <= end_year,
    Zone %in% zones_keep
  ) %>%
  dplyr::mutate(
    Rhat_num = suppressWarnings(as.numeric(Rhat)),
    ok_conv  = is.finite(Rhat_num) & (Rhat_num <= rhat_max_filtered) &
      is.finite(psi) & is.finite(sd) & !is.na(psi) & !is.na(sd)
  )

# Kept subset only
psi_zone_kept_scope_fc <- psi_zone_filter_scope_fc %>%
  dplyr::filter(Species %in% species_keep_fullyconverged)

# Row-level convergence rate (Zone-Year rows)
row_conv_rate_fc <- mean(psi_zone_kept_scope_fc$ok_conv, na.rm = TRUE)

cat("\n--- Convergence rates within RETAINED subset ---\n")
cat("Row-level convergence rate (Zone-Year rows ok):",
    round(100 * row_conv_rate_fc, 2), "%\n")

# Species-Year strict convergence rate (all zones ok) (should be 100% if filter was strict)
species_year_strict_fc <- psi_zone_kept_scope_fc %>%
  dplyr::group_by(Species, Year) %>%
  dplyr::summarise(
    n_zones_present = dplyr::n_distinct(Zone),
    n_zones_ok      = dplyr::n_distinct(Zone[ok_conv]),
    ok_year_strict  = (n_zones_present == length(zones_keep)) &
      (n_zones_ok == length(zones_keep)),
    .groups = "drop"
  )

species_year_strict_rate_fc <- mean(species_year_strict_fc$ok_year_strict, na.rm = TRUE)

cat("Species-Year strict convergence rate (all zones ok):",
    round(100 * species_year_strict_rate_fc, 2), "%\n")

# Completeness check: do retained species have all years present?
years_expected_fc <- seq(start_year, end_year)
species_year_completeness_fc <- species_year_strict_fc %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    n_years_present = dplyr::n_distinct(Year),
    complete_years  = (n_years_present == length(years_expected_fc)),
    .groups = "drop"
  )

cat("Retained species with complete Year coverage:",
    sum(species_year_completeness_fc$complete_years),
    "/", n_species_kept_fc, "\n")

# Optional breakdown (kept subset) — if strict filter, should be all OK
fail_breakdown_fc <- species_year_strict_fc %>%
  dplyr::mutate(
    reason = dplyr::case_when(
      n_zones_present < length(zones_keep) ~ "Missing zone(s)",
      n_zones_ok      < length(zones_keep) ~ "Rhat/NA fail in >=1 zone",
      TRUE ~ "OK"
    )
  ) %>%
  dplyr::count(reason) %>%
  dplyr::mutate(pct = 100 * n / sum(n))

cat("\nSpecies-Year strict status breakdown (retained subset):\n")
print(fail_breakdown_fc)

cat("\n==============================================================\n")























# patchwork - reporting ALL and fully converged datasets together ####

## =========================================================
## ADD-ON: Full vs Fully Converged composite trends (patchwork)
##  - Endpoint labels = % change in composite occupancy (start → end)
##  - Two panels: Full Dataset and Converged Dataset
##  - Convergence reported at ROW level (species–zone–year)
## =========================================================

library(patchwork)

## ---------- 0) Helper for label text ----------
fmt_pct <- function(x, digits = 1) {
  paste0(ifelse(x >= 0, "+", ""), round(x, digits), "%")
}

## ---------- 1) FULL DATASET: % change labels + ROW-LEVEL convergence ----------

# % change between start & end for full dataset
zone_change_full <- zone_trend %>%
  dplyr::filter(Year %in% c(start_year, end_year)) %>%
  dplyr::select(Zone, Year, gm_median) %>%
  tidyr::pivot_wider(
    names_from  = Year,
    values_from = gm_median,
    names_glue  = "gm_{Year}"
  ) %>%
  dplyr::mutate(
    pct_change = 100 * (
      .data[[paste0("gm_", end_year)]] -
        .data[[paste0("gm_", start_year)]]
    ) /
      .data[[paste0("gm_", start_year)]]
  )

# y-limits for label positioning
plot_limits_full <- zone_trend %>%
  dplyr::summarise(
    ymin = min(gm_low,  na.rm = TRUE),
    ymax = max(gm_high, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    yrange = ymax - ymin,
    pad    = pmax(0.005, 0.05 * yrange)
  )

# endpoint rows + label positions
zone_end_full <- zone_trend %>%
  dplyr::filter(Year == end_year) %>%
  dplyr::left_join(
    zone_change_full %>%
      dplyr::transmute(
        Zone,
        end_val    = .data[[paste0("gm_", end_year)]],
        pct_change
      ),
    by = "Zone"
  ) %>%
  dplyr::mutate(
    label_text = paste0(
      sprintf("%.3f", end_val), " (", fmt_pct(pct_change), ")"
    )
  ) %>%
  tidyr::crossing(plot_limits_full) %>%
  dplyr::mutate(
    off = pmax(0.01, 0.18 * yrange),
    label_base = dplyr::case_when(
      Zone == "Core"    ~ gm_median - off,
      Zone == "Buffer"  ~ gm_median,
      Zone == "Outside" ~ gm_median + off,
      TRUE              ~ gm_median
    )
  ) %>%
  disperse_labels(gap = 0.10 * unique(plot_limits_full$yrange)) %>%
  dplyr::mutate(
    x_label = end_year + 2.3,
    x_curve = end_year + 2.5
  )

# ROW-LEVEL convergence for full dataset (species–zone–year)
rhat_max_full <- 1.10
neff_min_full <- 0

psi_zone_full_scope <- psi_zone %>%
  dplyr::filter(
    Year >= start_year, Year <= end_year,
    Zone %in% zones_keep
  ) %>%
  dplyr::mutate(
    Rhat_num = suppressWarnings(as.numeric(Rhat)),
    neff_num = suppressWarnings(as.numeric(n.eff)),
    ok_conv  = is.finite(Rhat_num) & (Rhat_num <= rhat_max_full) &
      is.finite(neff_num) & (neff_num >= neff_min_full) &
      is.finite(psi) & is.finite(sd) & !is.na(psi) & !is.na(sd)
  )

row_conv_rate_full <- mean(psi_zone_full_scope$ok_conv, na.rm = TRUE)
n_species_full     <- dplyr::n_distinct(psi_zone$Species)

title_full <- paste0(
  "Full Dataset (", n_species_full, " species; ",
  round(100 * row_conv_rate_full, 1), "% convergence)"
)

## ---------- 2) FULLY CONVERGED DATASET: % change labels + ROW-LEVEL convergence ----------

# % change between start & end for fully converged composite
zone_change_fc2 <- zone_trend_fc %>%
  dplyr::filter(Year %in% c(start_year, end_year)) %>%
  dplyr::select(Zone, Year, gm_median) %>%
  tidyr::pivot_wider(
    names_from  = Year,
    values_from = gm_median,
    names_glue  = "gm_{Year}"
  ) %>%
  dplyr::mutate(
    pct_change = 100 * (
      .data[[paste0("gm_", end_year)]] -
        .data[[paste0("gm_", start_year)]]
    ) /
      .data[[paste0("gm_", start_year)]]
  )

plot_limits_fc2 <- zone_trend_fc %>%
  dplyr::summarise(
    ymin = min(gm_low,  na.rm = TRUE),
    ymax = max(gm_high, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    yrange = ymax - ymin,
    pad    = pmax(0.005, 0.05 * yrange)
  )

zone_end_fc2 <- zone_trend_fc %>%
  dplyr::filter(Year == end_year) %>%
  dplyr::left_join(
    zone_change_fc2 %>%
      dplyr::transmute(
        Zone,
        end_val    = .data[[paste0("gm_", end_year)]],
        pct_change
      ),
    by = "Zone"
  ) %>%
  dplyr::mutate(
    label_text = paste0(
      sprintf("%.3f", end_val), " (", fmt_pct(pct_change), ")"
    )
  ) %>%
  tidyr::crossing(plot_limits_fc2) %>%
  dplyr::mutate(
    off = pmax(0.01, 0.18 * yrange),
    label_base = dplyr::case_when(
      Zone == "Core"    ~ gm_median - off,
      Zone == "Buffer"  ~ gm_median,
      Zone == "Outside" ~ gm_median + off,
      TRUE              ~ gm_median
    )
  ) %>%
  disperse_labels(gap = 0.10 * unique(plot_limits_fc2$yrange)) %>%
  dplyr::mutate(
    x_label = end_year + 2.3,
    x_curve = end_year + 2.5
  )

# ROW-LEVEL convergence for fully converged subset
psi_zone_fc_scope <- psi_zone_full_scope %>%
  dplyr::filter(Species %in% species_keep_fullyconverged)

row_conv_rate_fc <- mean(psi_zone_fc_scope$ok_conv, na.rm = TRUE)
n_species_fc     <- length(species_keep_fullyconverged)

title_fc <- paste0(
  "Converged Dataset (", n_species_fc, " species; ",
  round(100 * row_conv_rate_fc, 1), "% convergence)"
)

## ---------- 3) Rebuild g_occ and g_occ_fc with % change endpoint labels ----------

g_occ <- ggplot(zone_trend, aes(x = Year, y = gm_median, color = Zone)) +
  geom_ribbon(aes(ymin = gm_low, ymax = gm_high, fill = Zone),
              alpha = 0.4, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.4) +
  geom_point(aes(fill = Zone), shape = 16, stroke = 1.5, size = 2) +
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_fill_colors) +
  geom_vline(xintercept = end_year, linetype = "dashed",
             color = "grey30", linewidth = 0.8) +
  geom_point(
    data = zone_end_full,
    aes(x = end_year, y = gm_median, color = Zone),
    shape = 21, fill = "white", stroke = 1.5, size = 3.5,
    inherit.aes = FALSE
  ) +
  geom_curve(
    data = zone_end_full,
    aes(x = end_year, y = gm_median,
        xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.8, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = zone_end_full,
    aes(x = x_label, y = label_y_adj, label = label_text, color = Zone),
    hjust = -0.05, size = 6, fontface = "bold",
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(start_year, end_year + 10), ylim = c(0, 0.2)) +
  labs(
    title = title_full,
    x = "Year",
    y = "Geometric mean occupancy",
    color = "Zone"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",          # NO legend in top plot
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 16),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  )

g_occ_fc <- ggplot(zone_trend_fc, aes(x = Year, y = gm_median, color = Zone)) +
  geom_ribbon(aes(ymin = gm_low, ymax = gm_high, fill = Zone),
              alpha = 0.4, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.4) +
  geom_point(aes(fill = Zone), shape = 16, stroke = 1.5, size = 2) +
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_fill_colors) +
  geom_vline(xintercept = end_year, linetype = "dashed",
             color = "grey30", linewidth = 0.8) +
  geom_point(
    data = zone_end_fc2,
    aes(x = end_year, y = gm_median, color = Zone),
    shape = 21, fill = "white", stroke = 1.5, size = 3.5,
    inherit.aes = FALSE
  ) +
  geom_curve(
    data = zone_end_fc2,
    aes(x = end_year, y = gm_median,
        xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.8, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = zone_end_fc2,
    aes(x = x_label, y = label_y_adj, label = label_text, color = Zone),
    hjust = -0.05, size = 6, fontface = "bold",
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(start_year, end_year + 10), ylim = c(0.08, 0.28)) +
  labs(
    title = title_fc,
    x = "Year",
    y = "Geometric mean occupancy",
    color = "Zone"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position   = "bottom",           # ONE shared legend here
    legend.key.size   = unit(1.24, "lines"), # smaller legend
    legend.text       = element_text(size = 16),
    legend.title      = element_text(size = 16),
    plot.title        = element_text(size = 18, face = "bold"),
    axis.title        = element_text(size = 18),
    axis.text         = element_text(size = 16),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  )

## ---------- 4) Patchwork: 2 rows × 1 column, legend under second plot ----------
library(patchwork)

g_combined <- (g_occ / g_occ_fc +
                 plot_annotation(tag_levels = 'A')) &
  theme(plot.tag = element_text(size = 18, face = "bold"))

print(g_combined)

# Optional save
ggsave("Composite_geo_occ_Full_vs_Converged_patchwork.png",
       g_combined, width = 8, height = 10, dpi = 600)
ggsave("Composite_geo_occ_Full_vs_Converged_patchwork.jpg",
       g_combined, width = 8, height = 10, dpi = 600)










# absolute change ####

## =========================================================
## ADD-ON: change endpoint labels to show ABSOLUTE change + 95% CrI
##  e.g. -0.032 [-0.050, -0.015]
##  Uses zone_year_draw_gm and zone_year_draw_gm_fc
## =========================================================

# Helper for formatting numeric values
fmt_abs <- function(x, digits = 3) {
  round(x, digits)
}

## ---------- 1) FULL DATASET: absolute change HDIs ----------

endpoint_draws_full <- zone_year_draw_gm %>%
  dplyr::filter(Year %in% c(start_year, end_year)) %>%
  dplyr::select(Zone, Year, draw_id, gm) %>%
  tidyr::pivot_wider(
    names_from  = Year,
    values_from = gm,
    names_prefix = "Y"
  ) %>%
  dplyr::mutate(
    s = pmax(.data[[paste0("Y", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Y", end_year)]],   1e-12),
    abs_change = f - s
  )

change_ci_full <- endpoint_draws_full %>%
  dplyr::group_by(Zone) %>%
  dplyr::summarise(
    abs_median = median(abs_change),
    abs_low    = as.numeric(quantile(abs_change, 0.025)),
    abs_high   = as.numeric(quantile(abs_change, 0.975)),
    .groups = "drop"
  )

zone_end_full_ci <- zone_end_full %>%
  dplyr::select(-label_text) %>%  # drop old label
  dplyr::left_join(change_ci_full, by = "Zone") %>%
  dplyr::mutate(
    label_text = paste0(
      fmt_abs(abs_median, 3), " [",
      fmt_abs(abs_low, 3), ", ",
      fmt_abs(abs_high, 3), "]"
    )
  )

## ---------- 2) FULLY CONVERGED DATASET: absolute change HDIs ----------

endpoint_draws_fc <- zone_year_draw_gm_fc %>%
  dplyr::filter(Year %in% c(start_year, end_year)) %>%
  dplyr::select(Zone, Year, draw_id, gm) %>%
  tidyr::pivot_wider(
    names_from  = Year,
    values_from = gm,
    names_prefix = "Y"
  ) %>%
  dplyr::mutate(
    s = pmax(.data[[paste0("Y", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Y", end_year)]],   1e-12),
    abs_change = f - s
  )

change_ci_fc <- endpoint_draws_fc %>%
  dplyr::group_by(Zone) %>%
  dplyr::summarise(
    abs_median = median(abs_change),
    abs_low    = as.numeric(quantile(abs_change, 0.025)),
    abs_high   = as.numeric(quantile(abs_change, 0.975)),
    .groups = "drop"
  )

zone_end_fc2_ci <- zone_end_fc2 %>%
  dplyr::select(-label_text) %>%  # drop old label
  dplyr::left_join(change_ci_fc, by = "Zone") %>%
  dplyr::mutate(
    label_text = paste0(
      fmt_abs(abs_median, 3), " [",
      fmt_abs(abs_low, 3), ", ",
      fmt_abs(abs_high, 3), "]"
    )
  )

## ---------- 3) Rebuild g_occ and g_occ_fc with new labels ----------

g_occ <- ggplot(zone_trend, aes(x = Year, y = gm_median, color = Zone)) +
  geom_ribbon(aes(ymin = gm_low, ymax = gm_high, fill = Zone),
              alpha = 0.4, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.4) +
  geom_point(aes(fill = Zone), shape = 16, stroke = 1.5, size = 2) +
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_fill_colors) +
  geom_vline(xintercept = end_year, linetype = "dashed",
             color = "grey30", linewidth = 0.8) +
  geom_point(
    data = zone_end_full_ci,
    aes(x = end_year, y = gm_median, color = Zone),
    shape = 21, fill = "white", stroke = 1.5, size = 3.5,
    inherit.aes = FALSE
  ) +
  geom_curve(
    data = zone_end_full_ci,
    aes(x = end_year, y = gm_median,
        xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.8, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = zone_end_full_ci,
    aes(x = x_label, y = label_y_adj, label = label_text, color = Zone),
    hjust = -0.05, size = 4, fontface = "bold",
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(start_year, end_year + 10), ylim = c(0, 0.2)) +
  labs(
    title = title_full,
    x = "Year",
    y = "Geometric mean occupancy",
    color = "Zone"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 16),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  )

g_occ_fc <- ggplot(zone_trend_fc, aes(x = Year, y = gm_median, color = Zone)) +
  geom_ribbon(aes(ymin = gm_low, ymax = gm_high, fill = Zone),
              alpha = 0.4, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.4) +
  geom_point(aes(fill = Zone), shape = 16, stroke = 1.5, size = 2) +
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_fill_colors) +
  geom_vline(xintercept = end_year, linetype = "dashed",
             color = "grey30", linewidth = 0.8) +
  geom_point(
    data = zone_end_fc2_ci,
    aes(x = end_year, y = gm_median, color = Zone),
    shape = 21, fill = "white", stroke = 1.5, size = 3.5,
    inherit.aes = FALSE
  ) +
  geom_curve(
    data = zone_end_fc2_ci,
    aes(x = end_year, y = gm_median,
        xend = x_curve, yend = label_y_adj, color = Zone),
    curvature = -0.25, linewidth = 0.8, show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = zone_end_fc2_ci,
    aes(x = x_label, y = label_y_adj, label = label_text, color = Zone),
    hjust = -0.05, size = 4, fontface = "bold",
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(start_year, end_year + 10), ylim = c(0.08, 0.28)) +
  labs(
    title = title_fc,
    x = "Year",
    y = "Geometric mean occupancy",
    color = "Zone"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position   = "bottom",
    legend.key.size   = unit(0.9, "lines"),
    legend.text       = element_text(size = 10),
    legend.title      = element_text(size = 11),
    plot.title        = element_text(size = 18, face = "bold"),
    axis.title        = element_text(size = 18),
    axis.text         = element_text(size = 16),
    panel.grid.major.y = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.major.x = element_blank()
  )

## ---------- 4) Patchwork: 2 rows × 1 column, legend under second plot ----------

g_combined <- g_occ / g_occ_fc

print(g_combined)

ggsave("Composite_geo_occ_Full_vs_Converged_patchwork_ABSchangelabels.png",
       g_combined, width = 8, height = 10, dpi = 600)
ggsave("Composite_geo_occ_Full_vs_Converged_patchwork_ABSchangelabels.jpg",
       g_combined, width = 8, height = 10, dpi = 600)



















# Convergence summary ####

############################################################
# ADD-ON: Convergence summaries (psi_zone only)
# 1) Row-level convergence across ALL psi_zone rows
# 2) Species that converged in BOTH first and last years
# 3) Species that converged in AT LEAST ONE year
#
# Assumptions:
# - A psi_zone row "converged" if:
#       Rhat <= conv_rhat_max
#       n.eff >= conv_neff_min
#       psi, sd are finite and non-NA
# - A species-year "converged" if ALL its psi_zone rows
#   (Species–Zone combinations for that Year) satisfy the convergence rule.
#   Species-level summaries then aggregate over these species-year flags.
############################################################

# ----------------------------
# SETTINGS
# ----------------------------
conv_rhat_max <- 1.10
conv_neff_min <- 0   # change if you want a stricter n.eff threshold

cat("\n==================== PSI_ZONE CONVERGENCE SUMMARY ====================\n")
cat("Convergence rule:\n")
cat("  Rhat <=", conv_rhat_max, "\n")
cat("  n.eff >=", conv_neff_min, "\n\n")

# ==========================================================
# Helper: convergence flag for psi_zone rows
# (species–zone–year; used in all parts)
# ==========================================================
psi_zone_conv <- psi_zone %>%
  mutate(
    Rhat_num = suppressWarnings(as.numeric(Rhat)),
    neff_num = suppressWarnings(as.numeric(n.eff)),
    ok_conv  = is.finite(Rhat_num) & (Rhat_num <= conv_rhat_max) &
      is.finite(neff_num) & (neff_num >= conv_neff_min) &
      is.finite(psi) & is.finite(sd) &
      !is.na(psi) & !is.na(sd)
  )

n_rows_total    <- nrow(psi_zone_conv)
n_rows_conv     <- sum(psi_zone_conv$ok_conv, na.rm = TRUE)
prop_rows_conv  <- 100 * n_rows_conv / n_rows_total
n_species_total <- n_distinct(psi_zone_conv$Species)

# ==========================================================
# 1) ROW-LEVEL convergence in psi_zone DATASET
# ==========================================================
cat("1) Row-level convergence (ALL psi_zone rows; species–zone–year):\n")
cat("   Rows converged:", n_rows_conv, "of", n_rows_total,
    sprintf("(%.2f%%)\n\n", prop_rows_conv))


# ==========================================================
# 2) Species that CONVERGED in BOTH first & last years
#    (e.g. 2000 AND 2023)
#
# Definition:
# - For each Species-Year, we say that year "converged" if
#   ALL psi_zone rows (across zones) are ok_conv.
# - A species "converged in first & last years" if it has a
#   converged year at start_year AND at end_year.
# ==========================================================
psi_end_years <- psi_zone_conv %>%
  filter(Year %in% c(start_year, end_year),
         Zone %in% zones_keep)

species_year_conv_end <- psi_end_years %>%
  group_by(Species, Year) %>%
  summarise(
    year_converged = all(ok_conv),
    .groups = "drop"
  )

species_both_end <- species_year_conv_end %>%
  mutate(
    is_start = (Year == start_year) & year_converged,
    is_end   = (Year == end_year)   & year_converged
  ) %>%
  group_by(Species) %>%
  summarise(
    conv_start = any(is_start),
    conv_end   = any(is_end),
    conv_both  = conv_start & conv_end,
    .groups = "drop"
  )

n_species_both_conv <- sum(species_both_end$conv_both, na.rm = TRUE)
prop_species_both   <- 100 * n_species_both_conv / n_species_total

cat("2) Species that converged in BOTH first & last years (",
    start_year, " AND ", end_year, "):\n", sep = "")
cat("   Species converged in both:", n_species_both_conv, "of",
    n_species_total, sprintf("(%.2f%%)\n\n", prop_species_both))


# ==========================================================
# 3) Species that CONVERGED in AT LEAST ONE YEAR
#
# Definition:
# - For each Species-Year (within [start_year, end_year] & zones_keep),
#   year_converged = ALL its psi_zone rows are ok_conv.
# - A species is counted here if ANY of its years are converged.
# ==========================================================
psi_all_years <- psi_zone_conv %>%
  filter(Year >= start_year, Year <= end_year,
         Zone %in% zones_keep)

species_year_conv_all <- psi_all_years %>%
  group_by(Species, Year) %>%
  summarise(
    year_converged = all(ok_conv),
    .groups = "drop"
  )

species_any_year <- species_year_conv_all %>%
  group_by(Species) %>%
  summarise(
    conv_in_at_least_one_year = any(year_converged),
    .groups = "drop"
  )

n_species_any_conv <- sum(species_any_year$conv_in_at_least_one_year, na.rm = TRUE)
prop_species_any   <- 100 * n_species_any_conv / n_species_total

cat("3) Species that converged in AT LEAST ONE year (",
    start_year, "–", end_year, "):\n", sep = "")
cat("   Species converged in ≥1 year:", n_species_any_conv, "of",
    n_species_total, sprintf("(%.2f%%)\n", prop_species_any))

cat("\nPSI_ZONE CONVERGENCE SUMMARY COMPLETE.\n")




#
# ==========================================================
# 4c) Species that CONVERGED IN ALL YEARS
#
# Definition:
# - For each Species-Year (within [start_year, end_year] & zones_keep),
#   year_converged = ALL its psi_zone rows are ok_conv.
# - A species "converged in ALL years" if EVERY one of its years
#   between start_year and end_year (inclusive) are converged.
# ==========================================================

species_all_years <- species_year_conv_all %>%
  group_by(Species) %>%
  summarise(
    conv_in_all_years = all(year_converged),
    .groups = "drop"
  )

species_all_years_list <- species_all_years %>%
  filter(conv_in_all_years) %>%
  arrange(Species) %>%
  pull(Species)

n_species_all_conv <- length(species_all_years_list)
prop_species_all   <- 100 * n_species_all_conv / n_species_total

cat("4c) Species that converged in ALL years (",
    start_year, "–", end_year, "):\n", sep = "")
cat("   Species converged in ALL years:", n_species_all_conv, "of",
    n_species_total, sprintf("(%.2f%%)\n", prop_species_all))
cat("   Species list:\n")
print(species_all_years_list)














############################################
# ADD-ON: Overall decline across ALL species in ALL zones
# Uses draws_long (Species, Zone, Year, draw_id, psi_draw)
# to build a single global composite geometric mean
# across all species × zones.
#
# OUTPUTS:
#   - 1_Global_composite_geo_occ_byYear_median_CrI.csv
#   - 1_Global_change_startEnd_median_CrI.csv
#   - 1_Global_growth_rate_median_CrI.csv
#
# Console summary reports:
#   - Absolute change in global gm occupancy (start → end)
#   - % change (start → end)
#   - Annual growth rate (%/yr)
############################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# Safety check --------------------------------------------------------------
if (!exists("draws_long")) {
  stop("This add-on requires 'draws_long' (Species, Zone, Year, draw_id, psi_draw).\n",
       "Run the main Cooke-style script up to the Monte Carlo section first.")
}

if (!exists("start_year") || !exists("end_year")) {
  stop("Need 'start_year' and 'end_year' objects defined (as in your main script).")
}

# User output file names ----------------------------------------------------
out_csv_global_comp   <- "1_Global_composite_geo_occ_byYear_median_CrI.csv"
out_csv_global_change <- "1_Global_change_startEnd_median_CrI.csv"
out_csv_global_growth <- "1_Global_growth_rate_median_CrI.csv"

# 1) Global composite gm by Year & draw -------------------------------------
# Here we ignore zones and treat every Species × Zone × Year cell
# as one "species-unit" contributing to the global composite.
# For each Year & draw_id we compute:
#   gm_global = exp(mean(log(psi_draw)))

global_year_draw_gm <- draws_long %>%
  group_by(Year, draw_id) %>%
  summarise(
    gm_global = exp(mean(log(psi_draw))),
    .groups = "drop"
  )

# Summarise over draws: median + 95% CrI per Year
global_trend <- global_year_draw_gm %>%
  group_by(Year) %>%
  summarise(
    gm_median = median(gm_global, na.rm = TRUE),
    gm_low    = as.numeric(quantile(gm_global, 0.025, na.rm = TRUE)),
    gm_high   = as.numeric(quantile(gm_global, 0.975, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(Year)

write.csv(global_trend, out_csv_global_comp, row.names = FALSE)
cat("\n[GLOBAL] Saved overall composite occupancy by year:\n  ",
    out_csv_global_comp, "\n", sep = "")

# 2) Start → end change (absolute + %) --------------------------------------
y_years <- end_year - start_year

global_start_end <- global_year_draw_gm %>%
  filter(Year %in% c(start_year, end_year)) %>%
  select(Year, draw_id, gm_global) %>%
  pivot_wider(
    names_from  = Year,
    values_from = gm_global,
    names_glue  = "Year_{Year}"
  ) %>%
  filter(
    is.finite(.data[[paste0("Year_", start_year)]]),
    is.finite(.data[[paste0("Year_", end_year)]])
  ) %>%
  mutate(
    s = pmax(.data[[paste0("Year_", start_year)]], 1e-12),
    f = pmax(.data[[paste0("Year_", end_year)]],   1e-12),
    abs_change = f - s,
    pct_change = 100 * (f - s) / s,
    growth_rate_pct = 100 * ((f / s)^(1 / y_years) - 1)
  )

# Absolute & % change summaries
global_change_summary <- global_start_end %>%
  summarise(
    abs_median = median(abs_change, na.rm = TRUE),
    abs_low    = as.numeric(quantile(abs_change, 0.025, na.rm = TRUE)),
    abs_high   = as.numeric(quantile(abs_change, 0.975, na.rm = TRUE)),
    pct_median = median(pct_change, na.rm = TRUE),
    pct_low    = as.numeric(quantile(pct_change, 0.025, na.rm = TRUE)),
    pct_high   = as.numeric(quantile(pct_change, 0.975, na.rm = TRUE))
  )

write.csv(global_change_summary, out_csv_global_change, row.names = FALSE)
cat("[GLOBAL] Saved start→end change summary:\n  ",
    out_csv_global_change, "\n", sep = "")

# 3) Annual growth rate (%/yr) ----------------------------------------------
global_growth_summary <- global_start_end %>%
  summarise(
    growth_median = median(growth_rate_pct, na.rm = TRUE),
    growth_low    = as.numeric(quantile(growth_rate_pct, 0.025, na.rm = TRUE)),
    growth_high   = as.numeric(quantile(growth_rate_pct, 0.975, na.rm = TRUE))
  )

write.csv(global_growth_summary, out_csv_global_growth, row.names = FALSE)
cat("[GLOBAL] Saved annual growth rate summary:\n  ",
    out_csv_global_growth, "\n", sep = "")

# 4) Console report: "overall decline" text ---------------------------------
cat("\n==================== OVERALL DECLINE: ALL SPECIES × ALL ZONES ====================\n")
cat("Period: ", start_year, "–", end_year, " (", y_years, " years)\n\n", sep = "")

with(global_change_summary, {
  cat("Start → end absolute change in global composite occupancy:\n")
  cat("  Δ = ", sprintf("%.3f", abs_median),
      " [", sprintf("%.3f", abs_low), ", ", sprintf("%.3f", abs_high), "]\n\n", sep = "")
  
  cat("Start → end % change in global composite occupancy:\n")
  cat("  ", sprintf("%.1f", pct_median), "% [",
      sprintf("%.1f", pct_low), ", ", sprintf("%.1f", pct_high), "]\n\n", sep = "")
})

with(global_growth_summary, {
  cat("Annual growth rate of global composite occupancy:\n")
  cat("  ", sprintf("%.2f", growth_median), "%/yr [",
      sprintf("%.2f", growth_low), ", ", sprintf("%.2f", growth_high), "]\n\n", sep = "")
})

cat("Interpretation:\n")
cat("  • A NEGATIVE Δ and NEGATIVE % change indicate an overall decline\n")
cat("    in the global composite (geometric mean) occupancy across ALL\n")
cat("    species in ALL zones over the study period.\n")
cat("  • The annual growth rate is the average % change per year over\n")
cat("    ", y_years, " years, calculated on the composite draws.\n", sep = "")
cat("=============================================================================\n")

