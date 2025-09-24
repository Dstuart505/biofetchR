#' Check if GBIF has records with coordinates for a species in a country
#'
#' @param species Scientific name (character)
#' @param cc ISO2 country code (character)
#' @param verbose Logical. If TRUE, prints a short message. Default is FALSE.
#'
#' @return Logical TRUE/FALSE
#' @export
check_gbif_presence <- function(species, cc, verbose = FALSE) {
  tryCatch({
    res <- rgbif::occ_search(
      scientificName = species,
      country = cc,
      hasCoordinate = TRUE,
      limit = 1
    )
    has_records <- !is.null(res$data) && is.data.frame(res$data) && nrow(res$data) > 0

    if (verbose) {
      msg <- if (has_records) "✅ Records found" else "❌ No records"
      cli::cli_alert_info("{species} in {cc}: {msg}")
    }

    has_records
  }, error = function(e) {
    if (verbose) cli::cli_alert_danger("{species} in {cc}: Error - {e$message}")
    FALSE
  })
}



#' Get GBIF taxonKey for a species name
#'
#' @param species Scientific name (character)
#' @param verbose Logical. If TRUE, prints status message. Default is FALSE.
#'
#' @return Integer GBIF taxonKey or NA if lookup fails
#' @export
get_taxon_key <- function(species, verbose = FALSE) {
  tryCatch({
    backbone <- rgbif::name_backbone(name = species)
    key <- backbone$usageKey
    if (verbose) {
      if (!is.null(key)) {
        cli::cli_alert_success("{species}: TaxonKey = {key}")
      } else {
        cli::cli_alert_warning("{species}: ❌ No taxonKey found")
      }
    }
    if (!is.null(key)) key else NA
  }, error = function(e) {
    if (verbose) cli::cli_alert_danger("{species}: Error - {e$message}")
    NA
  })
}

#' Check if GBIF has global records with coordinates for a marine species
#'
#' This version omits the country filter and checks globally.
#'
#' @param species Scientific name (character)
#' @param verbose Logical. If TRUE, prints a short message. Default is FALSE.
#'
#' @return Logical TRUE/FALSE
#' @export
check_gbif_presence_global <- function(species, verbose = FALSE) {
  tryCatch({
    res <- rgbif::occ_search(
      scientificName = species,
      hasCoordinate = TRUE,
      limit = 1
    )
    has_records <- !is.null(res$data) && is.data.frame(res$data) && nrow(res$data) > 0

    if (verbose) {
      msg <- if (has_records) "✅ Global records found" else "❌ No global records"
      cli::cli_alert_info("{species}: {msg}")
    }

    has_records
  }, error = function(e) {
    if (verbose) cli::cli_alert_danger("{species}: Error - {e$message}")
    FALSE
  })
}
