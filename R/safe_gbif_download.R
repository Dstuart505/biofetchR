#' Submit GBIF downloads per country for terrestrial species (GADM only)
#'
#' Submits one GBIF occurrence download per country using the 'country' filter.
#'
#' @param species Scientific name (character).
#' @param iso2_codes Character vector of ISO2 country codes.
#' @param user,pwd,email GBIF credentials.
#'
#' @return Named list of GBIF download keys (one per country).
#' @export
safe_gbif_download_gadm <- function(species, iso2_codes, user, pwd, email) {
  taxon_key <- get_taxon_key(species)
  if (is.na(taxon_key)) {
    cli::cli_alert_warning("❌ Could not resolve taxonKey for {.val {species}}")
    return(NULL)
  }

  keys <- list()
  count_active <- function(klist) {
    active_keys <- unlist(klist)
    if (length(active_keys) == 0) return(0)
    sum(purrr::map_lgl(active_keys, function(k) {
      status <- tryCatch(rgbif::occ_download_meta(k)$status, error = function(e) NA)
      status %in% c("PREPARING", "RUNNING")
    }))
  }

  for (cc in iso2_codes) {
    attempt <- 1
    success <- FALSE
    while (!success && attempt <= 3) {
      while (count_active(keys) >= 3) Sys.sleep(60)

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
        cli::cli_alert_success("📦 Submitted download for {.val {species}} in {.val {cc}}")
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

#' Submit global GBIF download for marine species (EEZ only)
#'
#' Submits a single global download (no country filter) and duplicates the key for bookkeeping.
#'
#' @param species Scientific name (character).
#' @param user,pwd,email GBIF credentials.
#'
#' @return Named list of GBIF download keys (single key with name "global").
#' @export
safe_gbif_download_eez <- function(species, user, pwd, email) {
  taxon_key <- get_taxon_key(species)
  if (is.na(taxon_key)) {
    cli::cli_alert_warning("❌ Could not resolve taxonKey for {.val {species}}")
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


