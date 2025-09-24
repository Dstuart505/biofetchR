#' List available native/non-native status presets
#'
#' Presets used by plotting/summarising functions to filter by origin/establishment.
#' @return data.frame with columns: preset, includes, notes
#' @export
list_status_presets <- function() {
  data.frame(
    preset = c(
      "all",
      "native",
      "non_native",         # alias for "alien"
      "alien",              # alias of non_native
      "introduced",
      "naturalized",
      "invasive",
      "native_plus_cryptogenic",
      "non_native_plus_cryptogenic"
    ),
    includes = c(
      "native, introduced, naturalized, invasive, cryptogenic, managed, unknown",
      "native",
      "introduced, naturalized, invasive",
      "introduced, naturalized, invasive",
      "introduced",
      "naturalized",
      "invasive",
      "native, cryptogenic",
      "introduced, naturalized, invasive, cryptogenic"
    ),
    notes = c(
      "No filtering",
      "Only taxa flagged as native/indigenous/endemic",
      "All alien/non-native categories (excl. cryptogenic)",
      "Alias for non_native",
      "Strictly introduced (often non-established)",
      "Established alien (naturalised)",
      "Alien flagged as invasive",
      "Treat cryptogenic with native (conservative)",
      "Treat cryptogenic with non-native (liberal)"
    ),
    stringsAsFactors = FALSE
  )
}

# ---- internal mappers (do not export) --------------------------------------
#' @keywords internal
.bf_normalize_status <- function(x) {
  z <- tolower(trimws(as.character(x)))

  # Treat empty/NA as unknown early
  z[is.na(z) | z %in% c("", "na", "n/a", "none", "unspecified")] <- "unknown"

  # Collapse punctuation/underscores and condense spaces
  z <- gsub("[^a-z]+", " ", z)
  z <- gsub("\\s+", " ", z)

  # canonical outputs: native | introduced | naturalized | invasive | cryptogenic | managed | unknown
  res <- ifelse(grepl("\\b(native|indigenous|endemic|resident)\\b", z), "native",
                ifelse(grepl("\\b(invasive)\\b", z), "invasive",
                       ifelse(grepl("\\b(naturaliz|naturalised|naturalized|establish(ed|ment)?)\\b", z), "naturalized",
                              ifelse(grepl("\\b(alien|non ?native|exotic|introduced|adventive|casual)\\b", z), "introduced",
                                     ifelse(grepl("\\b(cryptogenic)\\b", z), "cryptogenic",
                                            ifelse(grepl("\\b(unknown|uncertain|undetermined|unresolved|unverified|not specified)\\b", z), "unknown",
                                                   ifelse(grepl("\\b(managed|cultivat|captive|controlled)\\b", z), "managed", "unknown")))))))

res
}


#' @keywords internal
.bf_resolve_status_keep <- function(preset, custom_values = NULL) {
  if (!is.null(custom_values)) {
    canon <- c("native","introduced","naturalized","invasive","cryptogenic","managed","unknown")
    return(intersect(custom_values, canon))
  }
  preset <- tolower(preset)
  if (preset == "alien") preset <- "non_native"
  switch(preset,
         all                          = c("native","introduced","naturalized","invasive","cryptogenic","managed","unknown"),
         native                       = "native",
         non_native                   = c("introduced","naturalized","invasive"),
         introduced                   = "introduced",
         naturalized                  = "naturalized",
         invasive                     = "invasive",
         native_plus_cryptogenic      = c("native","cryptogenic"),
         non_native_plus_cryptogenic  = c("introduced","naturalized","invasive","cryptogenic"),
         c("native","introduced","naturalized","invasive","cryptogenic","managed","unknown")
  )
}

#' Filter/normalise records by native/alien status
#'
#' Normalises a status column to canonical classes and (optionally) filters rows.
#'
#' @param data A data.frame (or tibble) containing status information.
#' @param status_col Column name holding status labels. Defaults to `"alien_status"`.
#'   If not present, the function returns `data` unchanged (with a message if filtering requested).
#' @param preset One of [list_status_presets()]$preset. Default `"all"` (no filtering).
#' @param values Optional vector of canonical statuses to **keep**, e.g.
#'   `c("introduced","naturalized")`. Overrides `preset` if supplied.
#' @param normalize_only If `TRUE`, only add `status_norm` (no filtering).
#' @param verbose If `TRUE`, print a brief summary of counts kept/removed.
#' @return The input data with `status_norm` added/updated; filtered if requested.
#' @export
filter_by_status <- function(
    data,
    status_col = "alien_status",
    preset = "all",
    values = NULL,
    normalize_only = FALSE,
    verbose = TRUE
) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame-like object.")
  if (!status_col %in% names(data)) {
    if (!normalize_only && (!identical(tolower(preset), "all") || !is.null(values))) {
      warning("Status column '", status_col, "' not found; skipping status filter.")
    }
    return(data)
  }

  data$status_norm <- .bf_normalize_status(data[[status_col]])

  if (normalize_only) return(data)

  keep_vals <- .bf_resolve_status_keep(preset, values)
  n0 <- nrow(data)
  data <- data[data$status_norm %in% keep_vals, , drop = FALSE]
  n1 <- nrow(data)

  if (verbose) {
    kept_tab <- sort(table(data$status_norm), decreasing = TRUE)
    msg <- sprintf("Status filter: kept %d / %d rows (preset='%s').", n1, n0, preset)
    message(msg)
    if (length(kept_tab)) utils::capture.output(kept_tab) |> paste(collapse = "; ") |> message()
  }

  if (n1 == 0L) warning("After status filtering, no rows remain.")
  data
}

#' List available map projection keywords for biofetchR plots
#'
#' Returns a table of supported projection **keywords** and their corresponding
#' PROJ strings that can be passed via `crs_proj` in plotting functions.
#'
#' Use this to discover convenient, banding-safe global projections and a few
#' region-centered templates. You can always pass a full custom PROJ string to
#' `crs_proj` instead of a keyword.
#'
#' @return A `data.frame` with columns: `keyword`, `proj4`, and `description`.
#'
#' @examples
#' list_plot_projections()
#'
#' @export
list_plot_projections <- function() {
  data.frame(
    keyword = c(
      # Global (recommended)
      "wintri", "eqearth", "robin", "moll", "eck4",
      # Geographic / conformal
      "platecarree", "merc",
      # Lambert azimuthal equal-area presets
      "laea_world", "laea_north", "laea_south",
      # Orthographic hemispheric views (for demos)
      "ortho_atlantic", "ortho_pacific"
    ),
    proj4 = c(
      "+proj=wintri",
      "+proj=eqearth",
      "+proj=robin",
      "+proj=moll",
      "+proj=eck4",
      "+proj=longlat +datum=WGS84 +no_defs",
      "+proj=merc",
      "+proj=laea +lat_0=0  +lon_0=0",
      "+proj=laea +lat_0=90 +lon_0=0",
      "+proj=laea +lat_0=-90 +lon_0=0",
      "+proj=ortho +lat_0=30 +lon_0=-30",
      "+proj=ortho +lat_0=0  +lon_0=160"
    ),
    description = c(
      "Winkel Tripel (default; low seam issues)",
      "Equal Earth (equal-area; global)",
      "Robinson (compromise; banding-prone without wrapping)",
      "Mollweide (equal-area; banding-prone without wrapping)",
      "Eckert IV (equal-area; banding-prone without wrapping)",
      "Plate Carrée (geographic lon/lat)",
      "Mercator (conformal; strong area distortion)",
      "LAEA world (0/0 center)",
      "LAEA north polar",
      "LAEA south polar",
      "Orthographic centered on Atlantic (hemispheric view)",
      "Orthographic centered on Pacific (hemispheric view)"
    ),
    stringsAsFactors = FALSE
  )
}

#' Harmonize and bind a list of GBIF result data frames/sf objects
#' - No dependency on sf::st_geometry_name()
#' - Standardizes geometry column name to "geometry"
#' - Forces all non-geometry columns to character (avoids type conflicts)
#' - Ensures WGS84 (EPSG:4326)
harmonize_column_types <- function(x_list, coerce_all_to_character = TRUE) {
  # keep only non-empty elements
  x_list <- Filter(function(x) !is.null(x) && NROW(x) > 0, x_list)
  if (!length(x_list)) return(tibble::tibble())

  # helper: standardize to sf with geometry column "geometry"
  .to_sf_std <- function(x) {
    if (!inherits(x, "sf")) {
      if (all(c("decimalLongitude", "decimalLatitude") %in% names(x))) {
        x <- sf::st_as_sf(
          x, coords = c("decimalLongitude", "decimalLatitude"),
          crs = 4326, remove = FALSE
        )
      } else {
        # no coords, cannot make sf -> drop this chunk
        return(NULL)
      }
    }
    # ensure WGS84
    if (!is.na(sf::st_crs(x)) && sf::st_crs(x)$epsg != 4326) {
      x <- sf::st_transform(x, 4326)
    }
    # find the sfc column name from attributes/classes, not st_geometry_name()
    gcol <- attr(x, "sf_column", exact = TRUE)
    if (is.null(gcol) || !gcol %in% names(x)) {
      # locate first sfc column by class
      sfc_cols <- names(x)[vapply(x, function(col) inherits(col, "sfc"), logical(1))]
      if (length(sfc_cols)) {
        gcol <- sfc_cols[1]
      } else {
        return(NULL)  # no geometry found, drop
      }
    }
    # rename to "geometry" if needed and set sf attr
    if (!identical(gcol, "geometry")) {
      names(x)[names(x) == gcol] <- "geometry"
      attr(x, "sf_column") <- "geometry"
      sf::st_geometry(x) <- "geometry"
    }
    x
  }

  x_list <- lapply(x_list, .to_sf_std)
  x_list <- Filter(Negate(is.null), x_list)
  if (!length(x_list)) return(tibble::tibble())

  # union of attribute names (excluding geometry)
  name_union <- Reduce(union, lapply(x_list, function(x) setdiff(names(x), "geometry")))

  # align columns & coerce types
  x_list <- lapply(x_list, function(x) {
    missing_cols <- setdiff(name_union, names(x))
    if (length(missing_cols)) x[missing_cols] <- NA_character_
    x <- x[, c(name_union, "geometry"), drop = FALSE]
    if (isTRUE(coerce_all_to_character)) {
      for (nm in name_union) x[[nm]] <- as.character(x[[nm]])
    }
    x
  })

  # bind (sf rbind)
  out <- do.call(rbind, x_list)

  # if geometry types vary, try to cast to POINT (our GBIF outputs are points)
  try({
    gtyp <- unique(sf::st_geometry_type(out))
    if (length(gtyp) > 1) out <- sf::st_cast(out, "POINT", warn = FALSE)
  }, silent = TRUE)

  out
}

