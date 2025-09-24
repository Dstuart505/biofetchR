#' Export GBIF spatial data to CSV with optional thinning and coordinate restoration
#'
#' Applies spatial thinning (if enabled), restores coordinates, and writes to a standardized CSV.
#' Handles both terrestrial (GADM) and marine (EEZ/global) workflows.
#'
#' @param sf_df An `sf` object of GBIF occurrences joined to GADM or EEZ.
#' @param region_label Region identifier (e.g., `"ZA"` for terrestrial or `"global"` for marine).
#' @param species_name Scientific name of the species (e.g., `"Panthera leo"`).
#' @param output_dir Output directory for the resulting CSV file.
#' @param thinning Logical. Apply spatial thinning if `TRUE`. Default is `FALSE`.
#' @param dist_km Distance in km for thinning. Default is `5`.
#' @param store_in_memory Logical. If `TRUE`, returns a `data.frame` instead of writing to disk.
#' @param workflow Character. Either `"gadm"` (terrestrial) or `"eez"` (marine) for suffix handling.
#'
#' @return If `store_in_memory = TRUE`, returns a `data.frame`; otherwise returns output path (invisible).
#' @export
export_gbif_csv <- function(sf_df,
                            region_label,
                            species_name,
                            output_dir,
                            thinning = FALSE,
                            dist_km = 5,
                            store_in_memory = FALSE,
                            workflow = c("gadm", "eez")) {

  workflow <- match.arg(workflow)

  if (is.null(sf_df)) {
    cli::cli_alert_warning("⚠️ NULL input for {.val {species_name}} in {.val {region_label}} — skipping export.")
    return(if (store_in_memory) tibble::tibble() else NULL)
  }

  # Thinning
  processed_sf <- if (thinning) {
    thin_spatial_points(sf_df, dist_km = dist_km, return_all = FALSE)
  } else {
    sf_df
  }

  # Restore coordinates
  coords <- sf::st_coordinates(processed_sf)
  processed_sf$decimalLongitude <- coords[, 1]
  processed_sf$decimalLatitude  <- coords[, 2]
  output_df <- sf::st_drop_geometry(processed_sf)

  # Resolve region name (only for GADM)
  if (workflow == "gadm") {
    country_name <- countrycode::countrycode(region_label, origin = "iso2c", destination = "country.name", warn = TRUE)
    if (is.na(country_name)) {
      cli::cli_alert_danger("❌ Could not resolve country name for {.val {region_label}}")
      return(if (store_in_memory) tibble::tibble() else NULL)
    }
    region_clean <- gsub("[^a-zA-Z0-9]", "_", country_name)
  } else {
    region_clean <- gsub("[^a-zA-Z0-9]", "_", region_label)
  }

  species_clean <- gsub(" ", "_", species_name)
  suffix <- paste0("_gbif_", workflow, if (thinning) "_thinned", ".csv")
  out_path <- file.path(output_dir, paste0(species_clean, "_", region_label, "_", region_clean, suffix))

  if (store_in_memory) {
    cli::cli_alert_success("✅ Cleaned data for {.val {species_name}} in {.val {region_label}} [rows: {.val {nrow(output_df)}}]")
    return(output_df)
  } else {
    tryCatch({
      readr::write_csv(output_df, out_path)
      cli::cli_alert_success("✅ Saved: {.val {out_path}} [rows: {.val {nrow(output_df)}}]")
    }, error = function(e) {
      cli::cli_alert_danger("❌ Write failed for {.val {out_path}}: {e$message}")
    })
    return(invisible(out_path))
  }
}
