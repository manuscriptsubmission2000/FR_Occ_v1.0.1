############################################################
# TAXONOMIC WORKFLOW: diagnosis → guild update → TAXCORR records → guild totals
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(lubridate)
})

# ----------------------------------------------------------
# 0) Paths, working directory, helpers
# ----------------------------------------------------------

setwd("anon")

records_path      <- "anon/AllParks_records_FinalFiltered_Better2.csv"
guild_path_in     <- "anon/combined_all_parks_with_guild.csv"
guild_path_out    <- "combined_all_parks_with_guild_TAXCORR_ADDED.csv"
bugs_final_path   <- "Combined_BUGS_Data_Outhwaite_finalclustv5_32000iterations_FINAL_WITH_TAXCORR.csv"

std_key <- function(x) str_squish(str_to_lower(as.character(x)))

# ----------------------------------------------------------
# 1) Define taxonomic maps: synonyms + cryptic lumps
# ----------------------------------------------------------

syn_map <- c(
  # Maculinea -> Phengaris
  "Maculinea arion" = "Phengaris arion",
  "Maculinea alcon" = "Phengaris alcon",
  
  # Agrodiaetus -> Polyommatus
  "Agrodiaetus damon" = "Polyommatus damon",
  "Agrodiaetus dolus" = "Polyommatus dolus",
  "Agrodiaetus humedasae" = "Polyommatus humedasae",
  
  # Plebicula -> Polyommatus
  "Plebicula dorylas" = "Polyommatus dorylas",
  "Plebicula escheri" = "Polyommatus escheri",
  
  # Cyaniris -> Polyommatus
  "Cyaniris semiargus" = "Polyommatus semiargus",
  
  # Vacciniina -> Agriades
  "Vacciniina optilete" = "Agriades optilete",
  
  # Palaeochrysophanus -> Lycaena
  "Palaeochrysophanus hippothoe" = "Lycaena hippothoe",
  
  # Clossiana -> Boloria
  "Clossiana dia" = "Boloria dia",
  "Clossiana euphrosyne" = "Boloria euphrosyne",
  
  # Mellicta -> Melitaea
  "Mellicta athalia" = "Melitaea athalia",
  "Mellicta parthenoides" = "Melitaea parthenoides",
  "Mellicta varia" = "Melitaea varia",
  
  # Damora -> Argynnis
  "Damora pandora" = "Argynnis pandora",
  
  # Ladoga -> Limenitis
  "Ladoga camilla" = "Limenitis camilla",
  
  # Various hairstreak genera -> Satyrium
  "Nordmannia ilicis" = "Satyrium ilicis",
  "Quercusia quercus" = "Satyrium quercus",
  "Fixsenia esculi"   = "Satyrium esculi"
)

lump_map <- c(
  # Cryptic Leptidea complex
  "Leptidea sinapis" = "Leptidea sinapis_complex",
  "Leptidea reali"   = "Leptidea sinapis_complex"
)

# Helper tables for tracing
syn_map_tbl <- tibble(
  from_syn   = names(syn_map),
  to_taxcorr = unname(syn_map)
)

lump_map_tbl <- tibble(
  from_lump  = names(lump_map),
  to_taxcorr = unname(lump_map)
)

# ----------------------------------------------------------
# 2) Load raw records and diagnose taxonomic issues
# ----------------------------------------------------------

dat_raw <- read_csv(records_path, show_col_types = FALSE)

stopifnot(all(c("order", "species", "eventDate", "layer_grid", "Code_Zone") %in% names(dat_raw)))

dat <- dat_raw %>%
  mutate(
    species_raw = species,
    species     = str_squish(as.character(species))
  )

# Apply synonym map to see what actually changes
dat_tax_diag <- dat %>%
  mutate(
    species_revised = recode(species, !!!syn_map, .default = species),
    species_model   = recode(species_revised, !!!lump_map, .default = species_revised),
    needs_revision  = species != species_revised,
    needs_lumping   = species_revised != species_model,
    needs_any_merge = needs_revision | needs_lumping
  )

revision_pairs_present <- dat_tax_diag %>%
  filter(needs_revision) %>%
  distinct(species, species_revised) %>%
  arrange(species_revised, species)

lump_pairs_present <- dat_tax_diag %>%
  filter(needs_lumping) %>%
  distinct(species_revised, species_model) %>%
  arrange(species_model, species_revised)

write_csv(revision_pairs_present, "taxonomic_revision_pairs_present_in_dataset.csv")
write_csv(lump_pairs_present,     "taxonomic_lump_pairs_present_in_dataset.csv")

cat("\n==================== TAXONOMIC DIAGNOSIS ====================\n")
cat("Total rows in records:", nrow(dat), "\n")
cat("Rows needing revision:", sum(dat_tax_diag$needs_revision, na.rm = TRUE), "\n")
cat("Rows needing lumping :", sum(dat_tax_diag$needs_lumping,  na.rm = TRUE), "\n")
cat("Distinct revision pairs present:", nrow(revision_pairs_present), "\n")
cat("Distinct lump pairs present    :", nrow(lump_pairs_present), "\n")

# ----------------------------------------------------------
# 3) Update guild lookup to include TAXCORR names
# ----------------------------------------------------------

guild_df <- read_csv(guild_path_in, show_col_types = FALSE)

species_candidates <- c("species","Species","species_id","speciesID",
                        "species_name","Species_name","taxon","Taxon",
                        "scientificName","ScientificName")
guild_candidates   <- c("Guild","guild","GUILD")

species_col <- intersect(species_candidates, names(guild_df))[1]
guild_col   <- intersect(guild_candidates,   names(guild_df))[1]

if (is.na(species_col) || is.na(guild_col)) {
  stop("Could not detect Species and/or Guild columns in the original guild file.")
}

guild_lookup_orig <- guild_df %>%
  transmute(
    species = str_squish(as.character(.data[[species_col]])),
    Guild   = str_squish(as.character(.data[[guild_col]]))
  ) %>%
  filter(!is.na(species), species != "", !is.na(Guild), Guild != "") %>%
  distinct()

# Add revised names (synonyms)
guild_added_revised <- guild_lookup_orig %>%
  inner_join(syn_map_tbl, by = c("species" = "from_syn")) %>%
  transmute(
    species = to_taxcorr,
    Guild   = Guild,
    source  = "added_from_syn_map"
  ) %>%
  distinct()

# Add lump labels for Leptidea complex
guild_added_lumps <- guild_lookup_orig %>%
  inner_join(lump_map_tbl, by = c("species" = "from_lump")) %>%
  transmute(
    species = to_taxcorr,
    Guild   = Guild,
    source  = "added_from_lump_map"
  ) %>%
  distinct()

guild_updated <- bind_rows(
  guild_lookup_orig %>% mutate(source = "original"),
  guild_added_revised,
  guild_added_lumps
) %>%
  distinct(species, Guild, .keep_all = TRUE)

guild_conflicts <- guild_updated %>%
  group_by(species) %>%
  summarise(n_guilds = n_distinct(Guild), .groups = "drop") %>%
  filter(n_guilds > 1)

cat("\n==================== GUILD UPDATE SUMMARY ====================\n")
cat("Original guild rows:", nrow(guild_lookup_orig), "\n")
cat("Added revised-name rows:", nrow(guild_added_revised), "\n")
cat("Added lump-label rows  :", nrow(guild_added_lumps), "\n")
cat("Final guild rows:", nrow(guild_updated), "\n")
cat("Conflicting species->guild mappings:", nrow(guild_conflicts), "\n")

if (nrow(guild_conflicts) > 0) {
  cat("\nWARNING: Some species map to multiple guilds. Inspect these:\n")
  print(
    guild_updated %>%
      semi_join(guild_conflicts, by = "species") %>%
      arrange(species, Guild)
  )
}

write_csv(guild_updated, guild_path_out)
write_csv(
  bind_rows(guild_added_revised, guild_added_lumps) %>% arrange(source, species),
  "guild_rows_added_from_taxcorr_maps.csv"
)

cat("\nSaved updated guild file:\n  ", guild_path_out, "\n", sep = "")

# Build key-based version for joins later
guild_lookup_taxcorr <- guild_updated %>%
  mutate(
    species_taxcorr = species,
    sp_key = std_key(species_taxcorr)
  ) %>%
  select(sp_key, species_taxcorr, Guild)

# ----------------------------------------------------------
# 4) Apply TAXCORR names to ALL records & filter to 2000–2023
# ----------------------------------------------------------

dat_taxcorr <- dat %>%
  mutate(
    species_revised = recode(species, !!!syn_map, .default = species),
    species_taxcorr = recode(species_revised, !!!lump_map, .default = species_revised),
    sp_key          = std_key(species_taxcorr)
  )

records_2000_2023 <- dat_taxcorr %>%
  mutate(
    eventDate_parsed = suppressWarnings(parse_date_time(
      eventDate,
      orders = c("dmy", "ymd", "Ymd HMS", "Ymd", "mdy")
    )),
    year = year(eventDate_parsed)
  ) %>%
  filter(!is.na(eventDate_parsed), !is.na(year),
         year >= 2000, year <= 2023) %>%
  filter(!is.na(species_taxcorr), species_taxcorr != "",
         !is.na(layer_grid), layer_grid != "")

cat("\n==================== RECORD FILTER SUMMARY ====================\n")
cat("Rows after 2000–2023 + non-missing species/grid:", nrow(records_2000_2023), "\n")
cat("Distinct TAXCORR species in 2000–2023:", n_distinct(records_2000_2023$species_taxcorr), "\n")

species_counts_2000_2023 <- records_2000_2023 %>%
  count(sp_key, species_taxcorr, name = "n_records_2000_2023")

# ----------------------------------------------------------
# 5) Load FINAL TAXCORR BUGS file and extract species list
# ----------------------------------------------------------

bugs_final <- read_csv(bugs_final_path, show_col_types = FALSE)

stopifnot("Species" %in% names(bugs_final))

bugs_species <- bugs_final %>%
  filter(!is.na(Species), Species != "") %>%
  mutate(
    Species = str_squish(as.character(Species)),
    sp_key  = std_key(Species)
  ) %>%
  distinct(sp_key, Species)

n_bugs_species <- n_distinct(bugs_species$sp_key)
cat("\nSpecies in FINAL TAXCORR BUGS file:", n_bugs_species, "\n")

# ----------------------------------------------------------
# 6) Trace origins of BUGS species wrt records & taxonomic changes
# ----------------------------------------------------------

syn_targets  <- unique(unname(syn_map))
lump_targets <- unique(unname(lump_map))

species_in_raw_taxcorr <- dat_taxcorr %>%
  distinct(sp_key, species_taxcorr)

species_trace <- bugs_species %>%
  # attach record counts (2000–2023)
  left_join(species_counts_2000_2023,
            by = c("sp_key", "Species" = "species_taxcorr")) %>%
  # attach guild info
  left_join(guild_lookup_taxcorr, by = "sp_key") %>%
  mutate(
    is_syn_target  = Species %in% syn_targets,
    is_lump_target = Species %in% lump_targets
  ) %>%
  # where did this name come from?
  left_join(
    syn_map_tbl %>%
      group_by(to_taxcorr) %>%
      summarise(from_synonyms = paste(unique(from_syn), collapse = "; "),
                .groups = "drop"),
    by = c("Species" = "to_taxcorr")
  ) %>%
  left_join(
    lump_map_tbl %>%
      group_by(to_taxcorr) %>%
      summarise(from_lump_sources = paste(unique(from_lump), collapse = "; "),
                .groups = "drop"),
    by = c("Species" = "to_taxcorr")
  ) %>%
  # presence in raw records at all
  left_join(
    species_in_raw_taxcorr %>%
      transmute(sp_key, present_in_raw_records = TRUE),
    by = "sp_key"
  ) %>%
  mutate(
    present_in_raw_records = ifelse(is.na(present_in_raw_records), FALSE, present_in_raw_records),
    n_records_2000_2023    = ifelse(is.na(n_records_2000_2023), 0L, n_records_2000_2023),
    meets_50_filter        = n_records_2000_2023 >= 50,
    has_guild              = !is.na(Guild)
  ) %>%
  mutate(
    origin_tag = case_when(
      !present_in_raw_records & is_syn_target  ~ "synonym_target_no_raw_records",
      !present_in_raw_records & is_lump_target ~ "lump_target_no_raw_records",
      !present_in_raw_records                  ~ "no_raw_records",
      present_in_raw_records & !meets_50_filter & is_syn_target  ~ "synonym_target_<50_records_2000_2023",
      present_in_raw_records & !meets_50_filter & is_lump_target ~ "lump_target_<50_records_2000_2023",
      present_in_raw_records & !meets_50_filter                  ~ "<50_records_2000_2023",
      present_in_raw_records & meets_50_filter & is_syn_target   ~ "model_species_synonym_target",
      present_in_raw_records & meets_50_filter & is_lump_target  ~ "model_species_lump_target",
      present_in_raw_records & meets_50_filter                   ~ "model_species_direct",
      TRUE                                                       ~ "unknown"
    )
  )

write_csv(species_trace, "species_trace_TAXCORR_vs_records.csv")

cat("\n==================== SPECIES TRACE SUMMARY ====================\n")
print(
  species_trace %>%
    count(origin_tag, has_guild) %>%
    arrange(desc(n))
)
cat("\nSaved species trace table: species_trace_TAXCORR_vs_records.csv\n")

# ----------------------------------------------------------
# 7) Final guild-wise totals for modelled species (what you used in results)
# ----------------------------------------------------------

# Restrict records to:
#  - 2000–2023
#  - species_taxcorr present in BUGS_FINAL (sp_key overlap)
records_for_guild <- records_2000_2023 %>%
  mutate(sp_key = std_key(species_taxcorr)) %>%
  semi_join(bugs_species %>% select(sp_key), by = "sp_key") %>%
  left_join(guild_lookup_taxcorr %>% select(sp_key, Guild), by = "sp_key") %>%
  filter(!is.na(Guild))

guild_table_final <- records_for_guild %>%
  group_by(Guild) %>%
  summarise(
    Number_of_Records   = n(),
    Number_of_Species   = n_distinct(sp_key),
    Grid_Cells_Occupied = n_distinct(layer_grid),
    .groups = "drop"
  ) %>%
  arrange(desc(Number_of_Records))

write_csv(guild_table_final, "Guild_summary_TAXCORR_2000_2023_from_610species.csv")

cat("\n==================== FINAL GUILD SUMMARY (MODELLED SPECIES) ====================\n")
print(guild_table_final)
cat("Saved: Guild_summary_TAXCORR_2000_2023_from_610species.csv\n")

total_species_sum <- sum(guild_table_final$Number_of_Species)
unique_species_guild <- n_distinct(records_for_guild$sp_key)

cat("\nTotal species (sum over guild rows): ", total_species_sum, "\n")
cat("Unique species across all guilds:    ", unique_species_guild, "\n")
cat("Species in BUGS_FINAL (total):       ", n_bugs_species, "\n")
cat("================================================================\n")















# # ----------------------------------------------------------
# 5b) Build psi_zone_full (Species × Zone × Year estimates) from BUGS_FINAL
# ----------------------------------------------------------

start_year <- 2000
end_year   <- 2023
year_index_offset <- 1999
zones_keep <- c("Core", "Buffer", "Outside")
psi_zone_regex <- "^psi\\.fs\\.r_(Core|Buffer|Outside)\\[(\\d+)\\]$"

psi_zone_full <- bugs_final %>%
  filter(str_detect(Parameter, "^psi\\.fs\\.r_")) %>%
  mutate(
    Zone = str_match(Parameter, psi_zone_regex)[, 2],
    YearIndex = suppressWarnings(as.integer(str_match(Parameter, psi_zone_regex)[, 3])),
    Year = YearIndex + year_index_offset,
    psi = suppressWarnings(as.numeric(mean)),
    sd  = suppressWarnings(as.numeric(sd))
  ) %>%
  filter(!is.na(Zone), !is.na(YearIndex), Year >= start_year, Year <= end_year) %>%
  mutate(
    Zone = factor(Zone, levels = zones_keep),
    Species = str_squish(as.character(Species))
  ) %>%
  select(Species, Zone, Year, YearIndex, Parameter, psi, sd, Rhat, n.eff)



expected_years <- length(start_year:end_year)   # 24
expected_zones <- length(zones_keep)            # 3
expected_n_est <- expected_years * expected_zones

species_est_counts <- psi_zone_full %>%
  filter(is.finite(psi)) %>%
  distinct(Species, Zone, Year) %>%
  group_by(Species) %>%
  summarise(
    n_species_year_zone_estimates = n(),
    n_years_with_any_estimate     = n_distinct(Year),
    n_zones_with_any_estimate     = n_distinct(Zone),
    expected_n_estimates          = expected_n_est,
    prop_expected                 = n_species_year_zone_estimates / expected_n_est,
    .groups = "drop"
  )

species_trace <- species_trace %>%
  left_join(species_est_counts, by = "Species") %>%
  mutate(
    missing_any_species_year_zone = n_species_year_zone_estimates < expected_n_est
  )


# ----------------------------------------------------------
# 8) Totals per guild: Species×Year×Zone estimate coverage
# ----------------------------------------------------------

guild_est_totals <- species_trace %>%
  filter(!is.na(Guild)) %>%
  group_by(Guild) %>%
  summarise(
    n_species = n_distinct(sp_key),
    
    # totals across species
    total_species_year_zone_estimates = sum(n_species_year_zone_estimates, na.rm = TRUE),
    expected_total_estimates          = sum(expected_n_estimates, na.rm = TRUE),
    
    # coverage diagnostics
    mean_prop_expected = mean(prop_expected, na.rm = TRUE),
    n_species_missing_any = sum(missing_any_species_year_zone, na.rm = TRUE),
    
    # optional: how many unique years/zones are represented across species
    median_years_with_any_estimate = median(n_years_with_any_estimate, na.rm = TRUE),
    median_zones_with_any_estimate = median(n_zones_with_any_estimate, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(desc(total_species_year_zone_estimates))

cat("\n==================== GUILD TOTALS (Species×Year×Zone estimates) ====================\n")
print(guild_est_totals)
cat("==========================================================================\n")

# (optional) save
write_csv(guild_est_totals, "Guild_totals_species_year_zone_estimates.csv")


