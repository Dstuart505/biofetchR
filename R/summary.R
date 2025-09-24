#' Initialize a GBIF Processing Summary Table
#'
#' Creates a typed empty tibble to track GBIF data processing results
#' for each species × region (GADM or EEZ).
#'
#' @return A tibble with columns: species, region_id, region_type, n_total, n_cleaned, n_thinned, output_file, status.
#' @export
initialize_summary <- function() {
  tibble::tibble(
    species     = character(),  # Scientific name (e.g., "Rattus norvegicus")
    region_id   = character(),  # GADM country code or EEZ name (e.g., "US", "Gulf of Guinea")
    region_type = character(),  # "GADM" (terrestrial) or "EEZ" (marine)
    n_total     = integer(),    # Number of raw GBIF records before filtering
    n_cleaned   = integer(),    # Number after coordinate cleaning
    n_thinned   = integer(),    # Number after spatial thinning
    output_file = character(),  # Path to exported CSV or NA if in-memory
    status      = character()   # Status of operation ("success", "no_data", "error", etc.)
  )
}

#' Append a Row to the GBIF Summary Table
#'
#' Records the outcome of processing a species × region combination.
#'
#' @param summary_tbl A tibble from `initialize_summary()`.
#' @param species Character. Scientific name of the species.
#' @param region_id Character. GADM ISO2 code or EEZ label.
#' @param region_type Character. Either "GADM" or "EEZ".
#' @param n_total Integer. Raw GBIF records before filtering.
#' @param n_cleaned Integer. Records after coordinate filtering.
#' @param n_thinned Integer. Records after spatial thinning.
#' @param output_file Character. File path to CSV or NA if stored in-memory.
#' @param status Character. Summary status message.
#' @param quiet Logical. Suppress console message if TRUE.
#'
#' @return An updated summary tibble with the new row appended.
#' @export
append_summary_row <- function(summary_tbl,
                               species,
                               region_id,
                               region_type,
                               n_total,
                               n_cleaned,
                               n_thinned,
                               output_file,
                               status,
                               quiet = FALSE) {
  stopifnot(is.data.frame(summary_tbl))

  new_row <- tibble::tibble(
    species     = species,
    region_id   = region_id,
    region_type = region_type,
    n_total     = as.integer(n_total),
    n_cleaned   = as.integer(n_cleaned),
    n_thinned   = as.integer(n_thinned),
    output_file = output_file,
    status      = status
  )

  if (!quiet) {
    cli::cli_alert_success(
      "\U0001F4CA Summary updated: {.val {species}} in {.val {region_id}} ({.val {region_type}}) — {.val {status}} ({.val {n_thinned}} / {.val {n_total}})"
    )
  }

  dplyr::bind_rows(summary_tbl, new_row)
}
