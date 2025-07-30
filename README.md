
<!-- README.md is generated from README.Rmd. Please edit that file -->

# biofetchR

<!-- badges: start -->
<!-- badges: end -->

`biofetchR` is an R package designed to streamline the batch import,
spatial processing, and standardized export of GBIF occurrence data
tailored for large-scale biodiversity, invasion biology, and
macroecological research.

## Key Features

Authenticated GBIF downloads across multiple species × country
combinations, with robust retry and error handling.

Pre-download presence verification to minimize redundant queries and
optimize data retrieval efficiency.

Flexible spatial joins with Global Administrative Areas (GADM) at
multiple administrative levels, and optional marine Exclusive Economic
Zone (EEZ) joins for marine taxa.

Spatial thinning algorithms to mitigate spatial sampling bias,
customizable to user-specified distances and diagnostics.

Time-series data handling, enabling occurrence filtering and mapping
across user-defined temporal blocks, including multi-decadal
aggregation.

Automated export of clean, standardized CSV files for each species ×
country (or marine zone) combination, facilitating downstream analyses.

Support for generating publication-quality temporal and spatial
visualizations of species richness and occurrence patterns using
ggplot2.

Cache management and modular workflow components for reproducible,
high-throughput processing.

Seamless integration with species trait and environmental data pipelines
to enable advanced ecological and biogeographic modeling.

biofetchR provides a robust, modular foundation for integrating spatial,
temporal, and taxonomic filters in ecological data workflows, promoting
reproducibility and scalability in biodiversity science.

## Installation

You can install the development version of biofetchR from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("Dstuart505/biofetchR")
```

## Example

This is a basic example

``` r

library(biofetchR)

# --- GBIF User Credentials ---
# Provide your GBIF username, password, and email for authenticated downloads
user <- "your_username"
pwd <- "your_password"
email <- "your_email_address"

# --- Define Output Directories ---
# Specify where processed GBIF data and outputs will be saved
output_dir_terrestrial <- "D:/datasets/output_gbif_terrestrial"
output_dir_marine <- "D:/datasets/output_gbif_marine"  # reserved for marine species pipeline

# --- Prepare Terrestrial Species × Country Data Frame ---
# Define a tibble with species names and ISO2 country codes for batch processing
# Includes diverse taxa: mammals, plants, insects, freshwater species, birds, herpetofauna, fungi, and pathogens
test_terrestrial <- tibble::tibble(
  species = c(
    "Rattus norvegicus", "Sus scrofa",               # mammals
    "Lantana camara", "Imperata cylindrica",         # plants
    "Harmonia axyridis", "Aedes aegypti",            # insects
    "Salmo trutta", "Cyprinus carpio", "Eichhornia crassipes", # freshwater species
    "Sturnus vulgaris", "Columba livia",             # birds
    "Rhinella marina", "Trachemys scripta elegans",  # amphibians and reptiles
    "Cryphonectria parasitica",                       # fungi
    "Batrachochytrium dendrobatidis"                  # amphibian pathogen
  ),
  iso2c = c(
    "US", "AU", "IN", "BR", "DE", "CO",
    "GB", "PL", "NG", "ZA", "FR", "AU",
    "US", "US", "CR"
  )
)

# --- Run Terrestrial GBIF Processing Pipeline ---
# Downloads occurrences, performs spatial joins with GADM polygons,
# applies spatial thinning at 5 km distance to reduce sampling bias,
# exports standardized CSVs, and returns results in memory.
all_data_gadm <- process_gbif_gadm_pipeline(
  df = test_terrestrial,
  output_dir = output_dir_terrestrial,
  user = user,
  pwd = pwd,
  email = email,
  apply_thinning = TRUE,        # Enable spatial thinning to reduce oversampling bias
  dist_km = 5,                  # Minimum distance between retained points (km)
  batch_size = 3,               # Number of species-country combos per GBIF download batch
  return_all_results = TRUE,    # Return combined processed data as an sf object
  export_summary = TRUE,        # Generate summary tables for each batch
  store_in_memory = TRUE,       # Keep results in memory for immediate use
  use_planar = FALSE            # Use spherical geometry (recommended for global data)
)

# --- Define Temporal Blocks for Richness Mapping ---
# Custom time blocks covering mid-20th century to present
time_blocks <- list(
  "1935_1964" = 1935:1964,
  "1965_1994" = 1965:1994,
  "1995_2024" = 1995:2024
)

# --- Generate Global Species Richness Maps by Time Block ---
# Produces PNG maps and optionally an animated GIF showing temporal richness changes
plot_global_gadm_richness_by_time(
  data = all_data_gadm,           # Combined species occurrence data
  start_year = 1930,             # Start year for filtering data
  end_year = 2020,               # End year for filtering data
  block_length = 10,             # Length of each temporal block (years)
  output_dir = output_dir_terrestrial,
  gadm_cache_dir = NULL,         # Default cache location for GADM shapefiles
  include_antarctica = TRUE,     # Include Antarctica polygons in maps
  annotate_blocks = TRUE,        # Overlay period labels on each map
  animated_gif = TRUE,           # Save an animated GIF of the temporal sequence (requires 'magick')
  quiet = FALSE                  # Display progress messages and alerts
)
```

## Temporal Species Richness Animation

The animation below visualizes invasive species richness changes across
global GADM Level 1 regions over time.

This GIF was generated using the `plot_global_gadm_richness_by_time()`
function.

<img src="inst/figures/richness_animation.gif" width="700px" style="display: block; margin: auto;" />

## License

`biofetchR` is licensed under the [MIT
License](https://opensource.org/licenses/MIT).  
See the [`LICENSE`](LICENSE) file for details.

------------------------------------------------------------------------

## Bug Reports

If you encounter any issues or have feature requests, please report them
on the GitHub issue tracker:

<https://github.com/Dstuart505/biofetchR/issues>

Before reporting a bug, please ensure you are using the latest version
of `biofetchR` and provide a minimal reproducible example if possible.

------------------------------------------------------------------------
