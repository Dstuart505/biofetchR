#' Submit GBIF downloads per country for a terrestrial species (GADM only)
#'
#' Submits one GBIF occurrence download per country using the 'country' predicate.
#'
#' @param species Character. Scientific name of the species to query.
#' @param iso2_codes Character vector of ISO 2-letter country codes.
#' @param user,pwd,email GBIF credentials.
#'
#' @return Named list of GBIF download keys (one per country).
#' @export
download_gbif_batch_gadm <- function(species, iso2_codes, user, pwd, email) {
  taxon_key <- get_taxon_key(species)
  if (is.na(taxon_key)) {
    cli::cli_alert_warning("❌ Could not get taxonKey for {.val {species}}")
    return(NULL)
  }

  keys <- list()
  for (cc in iso2_codes) {
    attempt <- 1
    success <- FALSE

    while (!success && attempt <= 3) {
      result <- tryCatch({
        rgbif::occ_download(
          rgbif::pred_and(
            rgbif::pred("taxonKey", taxon_key),
            rgbif::pred("hasCoordinate", TRUE),
            rgbif::pred("country", cc)
          ),
          format = "SIMPLE_CSV",
          user = user, pwd = pwd, email = email
        )
      }, error = function(e) NULL)

      key <- if (is.character(result)) result else result$key

      if (!is.null(key)) {
        keys[[cc]] <- key
        cli::cli_alert_success("✅ Submitted download for {.val {species}} in {.val {cc}}")
        success <- TRUE
      } else {
        attempt <- attempt + 1
        Sys.sleep(60)
      }
    }
  }

  if (length(keys) == 0) return(NULL)
  return(keys)
}

#' Submit global GBIF download for a marine species (EEZ only)
#'
#' Submits a single global download for marine species and assigns the result to "global".
#'
#' @param species Character. Scientific name of the species to query.
#' @param user,pwd,email GBIF credentials.
#'
#' @return Named list of GBIF download keys (one global key).
#' @export
download_gbif_batch_eez <- function(species, user, pwd, email) {
  taxon_key <- get_taxon_key(species)
  if (is.na(taxon_key)) {
    cli::cli_alert_warning("❌ Could not get taxonKey for {.val {species}}")
    return(NULL)
  }

  keys <- list()
  attempt <- 1
  success <- FALSE

  while (!success && attempt <= 3) {
    result <- tryCatch({
      rgbif::occ_download(
        rgbif::pred_and(
          rgbif::pred("taxonKey", taxon_key),
          rgbif::pred("hasCoordinate", TRUE)
        ),
        format = "SIMPLE_CSV",
        user = user, pwd = pwd, email = email
      )
    }, error = function(e) NULL)

    key <- if (is.character(result)) result else result$key

    if (!is.null(key)) {
      keys[["global"]] <- key
      cli::cli_alert_success("🌊 Submitted global marine download for {.val {species}}")
      success <- TRUE
    } else {
      attempt <- attempt + 1
      Sys.sleep(60)
    }
  }

  if (length(keys) == 0) return(NULL)
  return(keys)
}

#' Submit global GBIF downloads for marine species (EEZ only)
#'
#' Submits one global GBIF occurrence download per species.
#'
#' @param df Data frame with at least a column `species`.
#' @param user,pwd,email GBIF credentials.
#'
#' @return Named list of GBIF download keys (one per species).
#' @export
download_gbif_batch_eez <- function(df, user, pwd, email) {
  if (!"species" %in% names(df)) {
    stop("Marine dataframe must contain a 'species' column.")
  }

  keys <- list()
  for (i in seq_len(nrow(df))) {
    sp <- df$species[i]

    taxon_key <- get_taxon_key(sp)
    if (is.na(taxon_key)) {
      cli::cli_alert_warning("❌ Could not get taxonKey for {.val {sp}}")
      next
    }

    attempt <- 1
    success <- FALSE
    while (!success && attempt <= 3) {
      result <- tryCatch({
        rgbif::occ_download(
          rgbif::pred_and(
            rgbif::pred("taxonKey", taxon_key),
            rgbif::pred("hasCoordinate", TRUE)
          ),
          format = "SIMPLE_CSV",
          user = user, pwd = pwd, email = email
        )
      }, error = function(e) NULL)

      key <- if (is.character(result)) result else result$key
      if (!is.null(key)) {
        keys[[sp]] <- key
        cli::cli_alert_success("🌊 Submitted global marine download for {.val {sp}}")
        success <- TRUE
      } else {
        attempt <- attempt + 1
        Sys.sleep(60)
      }
    }
  }

  if (length(keys) == 0) return(NULL)
  return(keys)
}
