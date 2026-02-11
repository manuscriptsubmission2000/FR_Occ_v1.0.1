# Modelling Occupancy of insects in France's National Parks

options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Function to check and install CRAN packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing missing package:", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
}

# Uninstall and reinstall 'sparta' package from GitHub
if ("sparta" %in% rownames(installed.packages())) {
  message("Removing existing 'sparta' package...")
  remove.packages("sparta")
}

message("Installing 'sparta' from GitHub...")
remotes::install_github("biologicalrecordscentre/sparta", force = TRUE)

# Ensure 'R2Jags' is installed
install_if_missing("R2jags")

# List of required packages
packages <- c("R2jags", "sparta", "dplyr", "lubridate", "tidyr", "doParallel", "foreach")

# Function to check and install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Apply the function to each package
sapply(packages, install_if_missing)

# Load the packages
lapply(packages, library, character.only = TRUE)


# Read data
df <- read.csv("AllParks_records_FinalFiltered_Better2.csv")
myData <- df[, c("order", "species", "eventDate", "layer_grid", "Code_Zone")]
myData <- myData[!duplicated(myData), ]

myData$time_period <- myData$eventDate
myData$taxa <- myData$species
myData$site <- myData$layer_grid

# Remove rows with missing or empty species information
myData <- myData %>%
  filter(!is.na(species) & species != "")

# Parse dates
myData$eventDate <- parse_date_time(myData$eventDate, orders = c("dmy", "ymd", "Ymd HMS", "Ymd", "mdy"))
myData <- myData[!is.na(myData$eventDate), ]

# Extract year and filter data from 2000 to 2023
myData$year <- as.numeric(format(myData$eventDate, "%Y"))
myData <- myData %>% filter(year >= 2000 & year <= 2023)

# Count occurrences per species within each order
species_occurrence_counts <- myData %>%
  group_by(order, species) %>%
  summarise(record_count = n(), .groups = "drop")

# NEW STEP: Filter species with at least 50 records 
filtered_species_counts <- species_occurrence_counts %>%
  filter(record_count >= 50)  # Remove species with fewer than 50 records

# Rank species by occurrence count within each order and extract the top 5000
top_500_species <- filtered_species_counts %>%
  group_by(order) %>%
  arrange(desc(record_count)) %>%
  slice_head(n = 5000) %>%
  ungroup()

# Filter the original dataset to include only the top species per order
# >>> THIS IS THE DATASET WE WILL USE FOR MODELLING <<<
filtered_data <- myData %>%
  filter(species %in% top_500_species$species)

# Summarize the filtered dataset for inspection
top_species_summary <- filtered_data %>%
  group_by(order, species) %>%
  summarise(total_records = n(), .groups = "drop")

# (Optional) summaries by zone
core_zone_data <- filtered_data %>%
  filter(Code_Zone == "Core_Zone")

buffer_zone_data <- filtered_data %>%
  filter(Code_Zone == "Buffer_Zone")

outside_boundary_data <- filtered_data %>%
  filter(Code_Zone == "Outside_Park_Boundary")

cat("Core Zone Summary:\n")
print(core_zone_data %>% summarise(total_records = n(), unique_species = n_distinct(species)))

cat("\nBuffer Zone Summary:\n")
print(buffer_zone_data %>% summarise(total_records = n(), unique_species = n_distinct(species)))

cat("\nOutside Park Boundary Summary:\n")
print(outside_boundary_data %>% summarise(total_records = n(), unique_species = n_distinct(species)))



### Model 1 - All zones jointly, with Code_Zone as region ####

# >>> NEW <<< use the combined filtered_data (all zones) for modelling
myData_subset <- filtered_data

# Extract and parse dates from eventDate column and assign to time_period
myData_subset$time_period <- as.Date(parse_date_time(
  myData_subset$eventDate,
  orders = c("dmy HMS", "dmy", "ymd HMS", "ymd", "mdy")
))

# First format our data for sparta
formattedOccData <- formatOccData(
  taxa   = myData_subset$taxa,
  site   = myData_subset$site,
  survey = myData_subset$time_period
)

# This is a list of two elements
names(formattedOccData)

# Have a look at spp_vis
head(formattedOccData$spp_vis[, 1:5])

# Have a look at occDetdata
head(formattedOccData$occDetdata)

species_list <- unique(myData_subset$species)
species_list

formattedOccData$spp_vis

species_present <- species_list %in% formattedOccData$spp_vis
names(species_present) <- species_list
species_present

# Ensure the year column is extracted
myData_subset$year <- as.numeric(format(myData_subset$time_period, "%Y"))

# >>> NEW <<< keep filtered_data intact; use a new object for the summary
species_year_data <- myData_subset[myData_subset$taxa %in% species_list, ]

# Group by species and year, then count records
species_year_summary <- species_year_data %>%
  group_by(taxa, year) %>%
  summarise(record_count = n(), .groups = "drop")

# View the summary
print(species_year_summary)

# >>> NEW <<< Build regional_codes from Code_Zone for use as regional factor
# We create a site × region matrix: one row per site, columns for Core/Buffer/Outside
regional_codes <- myData_subset %>%
  select(site, Code_Zone) %>%
  distinct() %>%
  mutate(
    Core    = ifelse(Code_Zone == "Core_Zone", 1, NA),
    Buffer  = ifelse(Code_Zone == "Buffer_Zone", 1, NA),
    Outside = ifelse(Code_Zone == "Outside_Park_Boundary", 1, NA)
  ) %>%
  select(site, Core, Buffer, Outside)

# Optional: check
head(regional_codes)


# Detect available cores based on the cluster environment
num_cores <- ifelse(
  Sys.getenv("SLURM_CPUS_PER_TASK") != "",
  min(4, as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))),  # Limit to 4 cores max
  ifelse(Sys.getenv("PBS_NP") != "", 
         min(4, as.numeric(Sys.getenv("PBS_NP"))),  # Limit to 4 cores max
         min(4, detectCores() - 1))  # Default: Use up to 4 cores but leave 1 free
)

# Ensure at least 1 core is used
num_cores <- 80  # Manually set to 80 cores

# Register parallel backend with a safe number of workers
cl <- makeCluster(80, type = "PSOCK")
registerDoParallel(cl)

# Print detected cores
print(paste("Using", num_cores, "CPU cores for parallel processing."))


# Define the occupancy modeling function
occ_mod_function <- function(taxa_name) {
  library(sparta)  # Ensure package is loaded inside each worker
  
  occ_out <- occDetFunc(
    taxa_name    = taxa_name,
    occDetdata   = formattedOccData$occDetdata,
    spp_vis      = formattedOccData$spp_vis,
    n_chains     = 3,
    n_iterations = 40000,
    burnin       = 20000,
    thinning     = 3,
    nyr          = 2,
    modeltype    = c("ranwalk", "halfcauchy", "catlistlength"),
    write_results = TRUE,
    seed          = 123,
    # >>> NEW <<< use Code_Zone as regional factor
    regional_codes = regional_codes,
    region_aggs    = list(AllZones = c("Core", "Buffer", "Outside"))
  )
  
  return(occ_out)
}

# Run the model in parallel using foreach
system.time({
  para_out <- foreach(taxa_name = species_list, .packages = "sparta") %dopar% occ_mod_function(taxa_name)
})

# Name each element of the output by species
names(para_out) <- sapply(para_out, function(x) x$SPP_NAM)

# Stop the cluster after processing
stopCluster(cl)

# Print completion message
print("Occupancy modelling (all zones, regional by Code_Zone) completed successfully!")


# Combine all species data into a single data.frame
combined_data <- do.call(rbind, lapply(names(para_out), function(species) {
  if (!is.null(para_out[[species]]) && !is.null(para_out[[species]]$BUGSoutput$summary)) {
    # Extract the BUGS summary and add a species column
    species_data <- para_out[[species]]$BUGSoutput$summary
    species_data <- data.frame(Species = species,
                               Parameter = rownames(species_data),
                               species_data)
    return(species_data)
  } else {
    return(NULL)
  }
}))

# Save the combined data to a single CSV file
write.csv(combined_data,
          "Combined_BUGS_Data_Out_finalclustv5_32000iterations.csv",
          row.names = FALSE)

saveRDS(para_out, "Out_model_finalclustv5.rds")

