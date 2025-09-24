#' Check if GBIF has coordinate-based records for a species in a region
#'
#' Uses `rgbif::occ_search()` to test whether any GBIF records with coordinates
#' exist for a given species in a specific ISO2 country or globally (for marine).
#'
#' @param species Character. Scientific name.
#' @param region_id Character. ISO2 country code or NULL for global.
#' @param use_eez_for_marine Logical. If TRUE, skips country filter for marine species.
#' @param verbose Logical. If TRUE, prints a CLI message. Default is FALSE.
#'
#' @return Logical. TRUE if at least one record is found, FALSE otherwise.
#' @export
check_gbif_presence <- function(species, region_id = NULL, use_eez_for_marine = FALSE, verbose = FALSE) {
  out <- tryCatch({
    args <- list(scientificName = species, hasCoordinate = TRUE, limit = 1)
    if (!(use_eez_for_marine && is_marine_species(species))) {
      args$country <- region_id
    }

    res <- do.call(rgbif::occ_search, args)
    has_records <- is.data.frame(res$data) && nrow(res$data) > 0

    if (verbose) {
      msg <- if (has_records) "✅ Records found" else "❌ No records"
      cli::cli_alert_info("{.val {species}} in {.val {region_id %||% 'global'}}: {msg}")
    }

    has_records
  }, error = function(e) {
    if (verbose) cli::cli_alert_danger("{.val {species}}: Error - {e$message}")
    FALSE
  })

  return(out)
}

#' Retrieve the GBIF taxonKey for a species
#'
#' Uses `rgbif::name_backbone()` to fetch the GBIF taxonKey for a given species.
#'
#' @param species Character. Scientific name.
#' @param verbose Logical. If TRUE, prints status message.
#'
#' @return Integer taxonKey or NA if lookup fails.
#' @export
get_taxon_key <- function(species, verbose = FALSE) {
  key <- tryCatch({
    backbone <- rgbif::name_backbone(name = species)
    if (!is.null(backbone$usageKey)) {
      if (verbose) cli::cli_alert_success("{.val {species}}: taxonKey = {.val {backbone$usageKey}}")
      backbone$usageKey
    } else {
      if (verbose) cli::cli_alert_warning("{.val {species}}: ❌ No taxonKey found")
      NA_integer_
    }
  }, error = function(e) {
    if (verbose) cli::cli_alert_danger("{.val {species}}: Error - {e$message}")
    NA_integer_
  })

  return(key)
}

#' Wait for GBIF downloads to complete and import results
#'
#' Polls GBIF for download status and imports results when ready.
#'
#' @param keys A named list of GBIF download keys.
#' @param wait_time Seconds to wait between polling attempts (default = 30).
#' @param max_tries Maximum number of polling attempts per key (default = 60).
#'
#' @return A named list of `sf` objects (one per region ID), or empty list on failure.
#' @export
wait_and_import_gbif <- function(keys, wait_time = 30, max_tries = 60) {
  results <- list()

  for (region_id in names(keys)) {
    key <- keys[[region_id]]
    success <- FALSE

    for (i in seq_len(max_tries)) {
      status <- tryCatch({
        rgbif::occ_download_meta(key)$status
      }, error = function(e) NA_character_)

      if (!is.na(status)) {
        if (status == "SUCCEEDED") {
          success <- TRUE
          break
        } else if (status %in% c("KILLED", "CANCELLED")) {
          cli::cli_alert_danger("❌ Download {key} failed with status: {status}")
          break
        }
      }

      Sys.sleep(wait_time)
    }

    if (!success) {
      cli::cli_alert_warning("⏱️ Timeout waiting for download {.val {key}} — skipping.")
      next
    }

    occ_data <- tryCatch({
      zipfile <- rgbif::occ_download_get(key, overwrite = TRUE, path = tempdir())
      df <- rgbif::occ_download_import(zipfile)
      df <- dplyr::filter(df, !is.na(decimalLatitude), !is.na(decimalLongitude))
      sf::st_as_sf(df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
    }, error = function(e) {
      cli::cli_alert_danger("❌ Failed to import download {.val {key}}: {e$message}")
      NULL
    })

    if (!is.null(occ_data)) {
      results[[region_id]] <- occ_data
    }
  }

  return(results)
}

#' Harmonize Column Types Across GBIF sf Objects
#'
#' Ensures consistent column types across GBIF `sf` objects to prevent `bind_rows()` errors.
#'
#' @param sf_list A named list of `sf` objects.
#' @return A list of `sf` objects with harmonized column types.
#' @export
harmonize_column_types <- function(sf_list) {
  all_cols <- unique(unlist(lapply(sf_list, names)))

  col_classes <- lapply(all_cols, function(col) {
    classes <- sapply(sf_list, function(x) {
      if (!is.null(x) && col %in% names(x)) class(x[[col]])[1] else NA_character_
    }, USE.NAMES = FALSE)
    common_class <- na.omit(classes)
    if (length(common_class) == 0) {
      return("character")
    } else {
      return(names(sort(table(common_class), decreasing = TRUE))[1])
    }
  })
  names(col_classes) <- all_cols

  sf_list_fixed <- lapply(sf_list, function(x) {
    if (is.null(x)) return(NULL)
    for (col in all_cols) {
      if (!col %in% names(x)) {
        x[[col]] <- NA
      }

      target_class <- col_classes[[col]]

      if (length(target_class) != 1 || is.na(target_class)) {
        warning("Column ", col, " has ambiguous or missing class. Defaulting to character.")
        target_class <- "character"
      }

      x[[col]] <- switch(
        target_class,
        character = as.character(x[[col]]),
        integer   = as.integer(x[[col]]),
        numeric   = as.numeric(x[[col]]),
        factor    = as.factor(x[[col]]),
        Date      = as.Date(x[[col]]),
        POSIXct   = as.POSIXct(x[[col]]),
        x[[col]]
      )
    }
    return(x)
  })

  return(sf_list_fixed)
}
