
# biofetchR

<!-- badges: start -->
<!-- badges: end -->

`biofetchR` is an R package designed to automate the batch import, spatial processing, and standardized export of GBIF occurrence data for ecological and biogeographic research. It is optimized for large-scale, cross-species applications in biodiversity science, invasion biology, and macroecology.

Key features include:

Authenticated GBIF downloads for multiple species × country combinations

Pre-download data presence checks to reduce redundant queries

Spatial joins with Global Administrative Area (GADM) polygons

Optional spatial thinning based on user-defined parameters to reduce spatial sampling bias

Export of standardized CSVs for each species × country combination

Support for time-series mapping of species occurrences across custom temporal windows

biofetchR facilitates reproducible, high-throughput processing of occurrence data and provides a modular foundation for integrating spatial, temporal, and taxonomic filters in ecological analyses.

## Installation

You can install the development version of biofetchR from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("Dstuart505/biofetchR")
```

## Features

Batch import of GBIF occurrence data for multiple species, taxa, or species × country combinations

Automated coordinate cleaning using the CoordinateCleaner package

Optional spatial thinning to reduce spatial clustering and autocorrelation

Rapid visualization of occurrence distributions with time and space filters

Designed for integration into reproducible trait–environment workflows

## Example

This is a basic example

```{r example}


library(biofetchR)

# Sample input: data frame with species and ISO2 country codes
input_df <- data.frame(
  species = c("Panthera leo", "Panthera tigris"),
  iso2c = c("ZA", "IN")
)

# Set your GBIF credentials (you must be a registered GBIF user)
user <- "your_gbif_username"
pwd <- "your_gbif_password"
email <- "your_email@example.com"

# Run the batch processing pipeline
process_gbif(
  df = input_df,
  output_dir = "data/processed_gbif",
  user = user,
  pwd = pwd,
  email = email,
  batch_size = 2,
  dist_km = 5,
  apply_thinning = TRUE
)

# OPTIONAL: Plot temporal trends in occurrences
# This assumes that plot_occurrences_over_time() reads CSVs from the output_dir
plot_occurrences_over_time(
  input_dir = "data/processed_gbif",
  species = "Panthera tigris",
  iso2c = "IN",
  start_year = 1950,
  end_year = 2022,
  interval = 10  # e.g., decadal bins

```{r}

## Citation
```{r}
citation("biofetchR")

```

## License

```{r}

MIT © Darren Stuart

```

## Bug Reports

Please open issues or feature requests at:
https://github.com/Dstuart505/biofetchR/issues
)

```

