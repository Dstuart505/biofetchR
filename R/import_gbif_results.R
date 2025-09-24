#' Import and wait for GBIF downloads to complete (GADM or EEZ)
#'
#' Waits for each GBIF download to complete, downloads the data,
#' and parses it as a spatial `sf` object. Handles both terrestrial and marine workflows.
#'
#' @param keys Named list of GBIF download keys (e.g., from `download_gbif_batch_gadm()` or `download_gbif_batch_eez()`).
#' @param crs Coordinate reference system for output (default: EPSG 4326).
#' @param start_year Optional: filter records from this year onwards.
#' @param end_year Optional: filter records up to and including this year.
#' @param max_wait_minutes Maximum wait time per download (default: 30 minutes).
#' @param poll_interval_sec Polling frequency for GBIF API (default: 30 seconds).
#'
#' @return Named list of `sf` objects per region. Adds a `"summary"` attribute (tibble).
#' @export
import_gbif_results <- function(keys,
                                crs = 4326,
                                start_year = NULL,
                                end_year = NULL,
                                max_wait_minutes = 30,
                                poll_interval_sec = 30) {
  results <- list()

  for (region in names(keys)) {
    key <- keys[[region]]
    cli::cli_alert_info("⏳ Waiting for GBIF download {.val {key}} for {.val {region}}")

    max_attempts <- ceiling((max_wait_minutes * 60) / poll_interval_sec)
    status <- NA_character_
    success <- FALSE

    for (i in seq_len(max_attempts)) {
      status <- tryCatch(rgbif::occ_download_meta(key)$status,
                         error = function(e) NA_character_)
      if (!is.na(status)) {
        if (status == "SUCCEEDED") {
          success <- TRUE
          break
        } else if (status %in% c("KILLED", "CANCELLED")) {
          cli::cli_alert_danger("❌ GBIF download {.val {key}} failed: {.val {status}}")
          break
        }
      }
      Sys.sleep(poll_interval_sec)
    }

    if (!success) {
      cli::cli_alert_danger("⏱️ Timeout or error for key {.val {key}} after {.val {max_wait_minutes}} minutes.")
      next
    }

    zipfile <- tryCatch(
      rgbif::occ_download_get(key, overwrite = TRUE, path = tempdir()),
      error = function(e) {
        cli::cli_alert_danger("⚠️ Failed to retrieve ZIP for {.val {region}}: {e$message}")
        NULL
      }
    )
    if (is.null(zipfile)) next

    occ <- tryCatch(
      rgbif::occ_download_import(zipfile),
      error = function(e) {
        cli::cli_alert_danger("⚠️ Failed to import data for {.val {region}}: {e$message}")
        NULL
      }
    )
    if (is.null(occ) || nrow(occ) == 0) {
      cli::cli_alert_warning("⚠️ No records retrieved for {.val {region}}.")
      next
    }

    # Filter bad coords
    occ <- dplyr::filter(occ, !is.na(decimalLatitude), !is.na(decimalLongitude))

    # Parse eventDate and apply temporal filtering
    occ$eventDate <- suppressWarnings(lubridate::ymd(occ$eventDate))
    occ$year <- lubridate::year(occ$eventDate)

    if (!is.null(start_year)) {
      occ <- dplyr::filter(occ, year >= start_year)
    }
    if (!is.null(end_year)) {
      occ <- dplyr::filter(occ, year <= end_year)
    }

    if (nrow(occ) == 0) {
      cli::cli_alert_warning("⚠️ No records after filtering for {.val {region}}.")
      next
    }

    # Drop all-NA columns
    occ <- occ[, colSums(!is.na(occ)) > 0]

    # Convert to sf
    occ_sf <- sf::st_as_sf(occ, coords = c("decimalLongitude", "decimalLatitude"),
                           crs = 4326, remove = FALSE)
    occ_sf <- sf::st_transform(occ_sf, crs = crs)

    results[[region]] <- occ_sf
    cli::cli_alert_success("✅ Imported {.val {nrow(occ_sf)}} records for {.val {region}}.")
  }

  attr(results, "summary") <- tibble::tibble(
    region = names(results),
    n_records = vapply(results, nrow, integer(1))
  )

  return(results)
}
