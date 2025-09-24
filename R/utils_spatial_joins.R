#' Load GADM Geometries for a Vector of ISO Country Codes
#'
#' Downloads and returns GADM spatial data for a list of countries at the specified administrative level.
#' Useful for batch loading only the geometries needed for your analysis (e.g., GADM Level 0, 1, or 2).
#'
#' @param iso_codes Character vector of ISO 2-letter country codes (e.g., "US", "AU", "ZA").
#' @param level Integer. GADM administrative level (0 = country, 1 = first-level divisions, etc.).
#' @param path Character. Directory to cache downloaded GADM data. Default is `"gadm_cache"`.
#'
#' @return An `sf` object containing geometries for the specified countries and level.
#' @export
load_all_gadm <- function(iso_codes, level = 1, path = "gadm_cache") {
  cli::cli_alert_info("📦 Loading GADM Level {level} data for {length(iso_codes)} countries...")

  all_sf <- purrr::map_dfr(iso_codes, function(code) {
    tryCatch({
      gadm <- geodata::gadm(country = code, level = level, path = path)
      sf::st_as_sf(gadm)
    }, error = function(e) {
      cli::cli_alert_warning("⚠️ Failed to load GADM for {.val {code}}: {e$message}")
      NULL
    })
  })

  if (nrow(all_sf) == 0) stop("❌ No GADM data could be loaded.")

  cli::cli_alert_success("✅ Successfully loaded GADM Level {level} data for {length(unique(all_sf$GID_0))} countries.")
  return(all_sf)
}

# helper (internal): pick usable name+gid columns and ensure standardized fields exist
.pick_join_cols <- function(g, level) {
  nm <- names(g)
  # prefer standardized columns if already normalized
  if (all(c("GADM_name","GID") %in% nm)) return(list(name="GADM_name", gid="GID", g=g))

  name_col <- paste0("NAME_", level)
  gid_col  <- paste0("GID_",  level)

  # fallbacks for name column
  cand_name <- c(name_col, "NAME","name","Name","ADMIN","COUNTRY","COUNTRY_NAME","CNTRY_NAME","SOVEREIGNT")
  name_hit  <- cand_name[cand_name %in% nm][1]
  if (is.na(name_hit)) stop(sprintf("No suitable name column found for level %d. Available: %s", level, paste(nm, collapse=", ")))

  # fallbacks for gid column
  cand_gid <- c(gid_col, "GID")
  gid_hit  <- cand_gid[cand_gid %in% nm][1]
  if (is.na(gid_hit)) stop(sprintf("No suitable GID column found for level %d. Available: %s", level, paste(nm, collapse=", ")))

  # create standardized columns for downstream
  g$GADM_name <- as.character(g[[name_hit]])
  g$GID       <- g[[gid_hit]]

  list(name="GADM_name", gid="GID", g=g)
}

#' Spatially Join GBIF Points to GADM Polygons (robust to level 0/1/2 schemas)
#'
#' @param results_list Named list of GBIF `sf` point data by ISO2 code.
#' @param gadm_data    Named list of GADM `sf` polygons (from `load_gadm()`).
#' @param level        GADM administrative level (0 = country, 1, 2, …).
#' @param use_planar   Logical; use planar CRS for joins (default: TRUE).
#' @return Named list of joined `sf` objects with `GADM_name` column.
#' @export
gadm_join <- function(results_list,
                      gadm_data,
                      level = 1,
                      use_planar = TRUE) {

  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")
  if (use_planar) sf::sf_use_s2(FALSE) else sf::sf_use_s2(TRUE)

  out <- vector("list", length(results_list))
  names(out) <- names(results_list)

  for (cc in names(results_list)) {
    gbif_sf <- results_list[[cc]]
    if (is.null(gbif_sf) || nrow(gbif_sf) == 0) {
      # keep structure consistent
      if (!"GADM_name" %in% names(gbif_sf)) gbif_sf$GADM_name <- NA_character_
      out[[cc]] <- gbif_sf
      next
    }

    gadm_sf <- gadm_data[[cc]]
    if (is.null(gadm_sf)) {
      message(sprintf("❌ No GADM data available for %s", cc))
      if (!"GADM_name" %in% names(gbif_sf)) gbif_sf$GADM_name <- NA_character_
      out[[cc]] <- gbif_sf
      next
    }

    # normalize / pick columns
    picked <- .pick_join_cols(gadm_sf, level)
    gadm_sf <- picked$g
    join_cols <- c(picked$name, picked$gid)

    # transform to a common CRS
    crs_target <- if (use_planar) 3857 else sf::st_crs(gadm_sf)
    gbif_sf <- sf::st_transform(gbif_sf, crs_target)
    gadm_sf <- sf::st_transform(gadm_sf, crs_target)

    # do the spatial join using the standardized columns
    gsel <- gadm_sf[, join_cols, drop = FALSE]

    joined <- tryCatch(
      sf::st_join(gbif_sf, gsel, join = sf::st_within),
      error = function(e) {
        message(sprintf("❌ st_within join failed for %s: %s", cc, e$message))
        gbif_sf
      }
    )

    # nearest fallback if within returns no matches
    if (!picked$name %in% names(joined) || all(is.na(joined[[picked$name]]))) {
      message(sprintf("⚠️ No matches for %s — using nearest feature fallback", cc))
      joined <- tryCatch(
        sf::st_join(gbif_sf, gsel, join = sf::st_nearest_feature),
        error = function(e) {
          message(sprintf("❌ st_nearest_feature failed for %s: %s", cc, e$message))
          if (!picked$name %in% names(gbif_sf)) gbif_sf[[picked$name]] <- NA_character_
          gbif_sf
        }
      )
    }

    # standardize name column for downstream
    if (picked$name != "GADM_name" && picked$name %in% names(joined)) {
      names(joined)[names(joined) == picked$name] <- "GADM_name"
    }
    if (!"iso2c" %in% names(joined)) joined$iso2c <- cc

    out[[cc]] <- joined
  }

  out
}
#new

#' Load, Download, and Cache GADM Polygons per Country and Level
#'
#' Downloads, converts to \code{sf}, and caches GADM polygons as `.rds` files for a set of countries at the specified administrative level.
#' If a cached file exists, it is loaded and normalized. Otherwise, data are downloaded via \code{geodata::gadm()}, converted to \code{sf},
#' augmented with standardized \code{GADM_name} and \code{GID} columns (via \code{.pick_join_cols()}), and saved to the specified cache directory.
#'
#' @param iso2c Character vector of ISO2 country codes (e.g., "FR", "ZA", "US").
#' @param level Integer. GADM administrative level to retrieve (0 = country, 1 = first admin, 2 = second admin, etc.).
#' @param cache_dir Directory to store cached GADM `.rds` files (default: user cache).
#' @param quiet Logical; suppress messages (default: FALSE).
#'
#' @return Named list of \code{sf} objects, one per country, with standardized \code{GADM_name} and \code{GID} columns.
#' @examples
#' \dontrun{
#' gadm_list <- load_gadm(c("FR", "DE"), level = 2)
#' }
#' @export
load_gadm <- function(iso2c, level = 1, cache_dir = gadm_cache_dir, quiet = FALSE) {
  if (!requireNamespace("geodata", quietly = TRUE)) stop("Package 'geodata' is required.")
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")
  if (!requireNamespace("countrycode", quietly = TRUE)) stop("Package 'countrycode' is required.")

  if (is.null(cache_dir)) {
    cache_dir <- rappdirs::user_cache_dir("biofetchR/gadm")
  }
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  gadm_list <- list()
  msg <- function(...) if (!quiet) message(...)

  for (cc in iso2c) {
    msg(glue::glue("ℹ 📥 Loading GADM Level {level} for {cc}"))
    iso3 <- countrycode::countrycode(cc, "iso2c", "iso3c")
    fname <- file.path(cache_dir, paste0("gadm41_", iso3, "_", level, "_pk.rds"))

    if (file.exists(fname)) {
      gadm_sf <- tryCatch(readRDS(fname), error = function(e) NULL)
      if (!is.null(gadm_sf)) {
        # normalize even cached objects
        picked <- .pick_join_cols(gadm_sf, level)
        gadm_sf <- picked$g
        msg(glue::glue("✔ ✅ Loaded cached GADM Level {level} for {cc}"))
      } else {
        msg(glue::glue("⚠️ Cache corrupted for {cc}, re-downloading..."))
        gadm_obj <- geodata::gadm(country = cc, level = level, path = cache_dir)
        gadm_sf <- sf::st_make_valid(sf::st_as_sf(gadm_obj))
        picked <- .pick_join_cols(gadm_sf, level)
        gadm_sf <- picked$g
        saveRDS(gadm_sf, fname)
      }
    } else {
      gadm_obj <- geodata::gadm(country = cc, level = level, path = cache_dir)
      gadm_sf <- sf::st_make_valid(sf::st_as_sf(gadm_obj))
      picked <- .pick_join_cols(gadm_sf, level)
      gadm_sf <- picked$g
      saveRDS(gadm_sf, fname)
      msg(glue::glue("✔ ✅ Downloaded and cached GADM Level {level} for {cc}"))
    }

    # keep essentials
    gadm_sf <- gadm_sf[, c("GID", "GADM_name", "geometry")]
    gadm_list[[cc]] <- gadm_sf
  }

  return(gadm_list)
}


#' Spatially Join GBIF Points to GADM Polygons
#'
#' Joins a list of GBIF point `sf` objects to cached GADM polygons.
#'
#' @param results_list Named list of GBIF `sf` point data by ISO2 code.
#' @param gadm_data Named list of GADM `sf` polygons (from `load_gadm()`).
#' @param level GADM administrative level used (for column naming).
#' @param use_planar Logical; whether to use planar projection (default: TRUE).
#'
#' @return Named list of joined `sf` objects with `GADM_name` column.
#' @export
gadm_join <- function(results_list,
                      gadm_data,
                      level = 1,
                      use_planar = TRUE) {

  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")

  picked <- .pick_join_cols(gadm_sf, level)
  gadm_sf <- picked$g
  join_cols <- c(picked$name, picked$gid)
  joined_list <- list()

  if (use_planar) sf::sf_use_s2(FALSE) else sf::sf_use_s2(TRUE)

  for (cc in names(results_list)) {
    gbif_sf <- results_list[[cc]]
    if (is.null(gbif_sf) || nrow(gbif_sf) == 0) {
      gbif_sf$GADM_name <- NA_character_
      joined_list[[cc]] <- gbif_sf
      next
    }

    gadm_sf <- gadm_data[[cc]]
    if (is.null(gadm_sf)) {
      message(sprintf("❌ No GADM data available for %s", cc))
      gbif_sf$GADM_name <- NA_character_
      joined_list[[cc]] <- gbif_sf
      next
    }

    crs_target <- if (use_planar) 3857 else sf::st_crs(gadm_sf)
    gbif_sf <- sf::st_transform(gbif_sf, crs_target)
    gadm_sf <- sf::st_transform(gadm_sf, crs_target)

    joined <- tryCatch(
      sf::st_join(gbif_sf, gadm_sf[, c(name_col, gid_col)], join = sf::st_within),
      error = function(e) {
        message(sprintf("❌ st_within join failed for %s: %s", cc, e$message))
        gbif_sf
      }
    )

    if (!(name_col %in% names(joined)) || all(is.na(joined[[name_col]]))) {
      message(sprintf("⚠️ No %s matches for %s — using nearest feature fallback", name_col, cc))
      joined <- tryCatch(
        sf::st_join(gbif_sf, gadm_sf[, c(name_col, gid_col)], join = sf::st_nearest_feature),
        error = function(e) {
          message(sprintf("❌ st_nearest_feature failed for %s: %s", cc, e$message))
          gbif_sf[[name_col]] <- NA_character_
          gbif_sf
        }
      )
    }

    # Rename to standardized column name for downstream compatibility
    names(joined)[names(joined) == name_col] <- "GADM_name"

    if (!"iso2c" %in% names(joined)) joined$iso2c <- cc

    joined_list[[cc]] <- joined
  }

  return(joined_list)
}

#' Spatially Join GBIF Points to EEZ Polygons
#'
#' Joins a list of GBIF point `sf` objects to Exclusive Economic Zones (EEZ) polygons.
#' Downloads EEZ polygons via the `mregions2` package if not cached locally.
#'
#' **Encoding:** EEZ names (`geoname`) are normalized to valid UTF-8 internally
#' to prevent Windows/locale-related issues during messaging and export.
#'
#' @param results_list Named list of `sf` point objects, keyed by species or identifiers.
#'   If a single `sf` object is provided, it will be wrapped into a length-1 list.
#'   If the list is unnamed, synthetic names are assigned.
#' @param use_planar Logical; if TRUE disables spherical geometry for spatial joins (default: TRUE).
#' @param eez_cache_file Optional path to a local `.gpkg` file where EEZ polygons will be cached and reused.
#'   Defaults to `"biofetchR_eez_cache.gpkg"` in the user cache directory.
#'
#' @return Named list of `sf` objects with EEZ name (`geoname`) and other metadata columns appended.
#'   The list has the same keys as the input `results_list`.
#' @export
#'
#' @importFrom sf st_as_sf st_transform st_crs st_join st_within st_nearest_feature st_make_valid read_sf write_sf st_use_s2
#' @importFrom mregions2 mrp_get
eez_join <- function(results_list,
                     use_planar = TRUE,
                     eez_cache_file = NULL) {

  # -- helper: normalize character vectors to valid UTF-8 (best effort)
  utf8_safe <- function(x) {
    if (is.null(x)) return(x)
    y <- suppressWarnings(iconv(x, from = "",        to = "UTF-8", sub = ""))
    bad <- is.na(y)
    if (any(bad)) y[bad] <- suppressWarnings(iconv(x[bad], from = "latin1", to = "UTF-8", sub = ""))
    bad <- is.na(y)
    if (any(bad)) y[bad] <- suppressWarnings(iconv(x[bad], from = "CP1252", to = "UTF-8", sub = ""))
    y[is.na(y)] <- x[is.na(y)]
    y
  }

  # Ensure required package is available
  if (!requireNamespace("mregions2", quietly = TRUE)) {
    stop("❌ Package 'mregions2' is required. Please install it with: remotes::install_github('iobis/mregions2', subdir = 'R')")
  }

  # Normalize input shape: allow a single sf, ensure names
  if (inherits(results_list, "sf")) {
    results_list <- list(results_list)
  }
  if (is.null(names(results_list))) {
    names(results_list) <- paste0("item_", seq_along(results_list))
  }

  # Set EEZ cache file path
  if (is.null(eez_cache_file)) {
    eez_cache_file <- file.path(rappdirs::user_cache_dir("biofetchR"), "biofetchR_eez_cache.gpkg")
  }

  # Create cache dir if needed
  if (!dir.exists(dirname(eez_cache_file))) {
    dir.create(dirname(eez_cache_file), recursive = TRUE)
  }

  # Load or download EEZ polygons
  if (file.exists(eez_cache_file)) {
    message("🌐 Using cached EEZ polygons: ", eez_cache_file)
    eez_sf <- sf::read_sf(eez_cache_file)
  } else {
    message("🌐 Downloading EEZ polygons via 'mregions2::mrp_get()'...")
    eez_sf <- tryCatch({
      mregions2::mrp_get("eez")
    }, error = function(e) {
      stop("❌ Failed to download EEZ polygons: ", e$message)
    })
    # Cache the EEZ polygons locally
    sf::write_sf(eez_sf, eez_cache_file)
  }

  # S2 mode per user choice
  if (isTRUE(use_planar)) sf::sf_use_s2(FALSE) else sf::sf_use_s2(TRUE)

  # Ensure valid geometries
  eez_sf <- sf::st_make_valid(eez_sf)

  # Normalize EEZ names to UTF-8 for safe display/logging
  if ("geoname" %in% names(eez_sf)) {
    eez_sf$geoname <- utf8_safe(eez_sf$geoname)
  } else {
    stop("EEZ polygons are missing the 'geoname' attribute.")
  }

  # Precompute target CRS and a transformed EEZ once
  crs_target <- if (isTRUE(use_planar)) 3857 else sf::st_crs(eez_sf)
  eez_sf_tr  <- sf::st_transform(eez_sf, crs_target)

  # Prepare output
  joined_list <- list()

  for (key in names(results_list)) {
    gbif_sf <- results_list[[key]]

    if (is.null(gbif_sf) || nrow(gbif_sf) == 0) {
      # keep structure; ensure geoname exists for downstream split
      gbif_sf$geoname <- NA_character_
      joined_list[[key]] <- gbif_sf
      next
    }

    # Coerce plain data.frame to sf if needed
    if (!inherits(gbif_sf, "sf")) {
      if (!all(c("decimalLongitude","decimalLatitude") %in% names(gbif_sf))) {
        warning(sprintf("Input '%s' is not sf and lacks lon/lat; skipping.", key))
        next
      }
      gbif_sf <- sf::st_as_sf(gbif_sf, coords = c("decimalLongitude","decimalLatitude"),
                              crs = 4326, remove = FALSE)
    }

    gbif_sf <- sf::st_transform(gbif_sf, crs_target)

    joined <- tryCatch(
      sf::st_join(gbif_sf, eez_sf_tr[, c("geoname", "mrgid")], join = sf::st_within),
      error = function(e) {
        msg_key <- utf8_safe(key)
        message(sprintf("❌ st_within join failed for %s: %s", msg_key, e$message))
        gbif_sf
      }
    )

    # If within-EEZ yields no matches, try nearest EEZ
    if (!("geoname" %in% names(joined)) || all(is.na(joined$geoname))) {
      msg_key <- utf8_safe(key)
      message(sprintf("⚠️ No EEZ matches for '%s' — trying nearest feature...", msg_key))
      joined <- tryCatch(
        sf::st_join(gbif_sf, eez_sf_tr[, c("geoname", "mrgid")], join = sf::st_nearest_feature),
        error = function(e) {
          message(sprintf("❌ st_nearest_feature failed for %s: %s", msg_key, e$message))
          gbif_sf$geoname <- NA_character_
          gbif_sf
        }
      )
    }

    # Ensure UTF-8 labels on the result as well
    if ("geoname" %in% names(joined)) {
      joined$geoname <- utf8_safe(joined$geoname)
    }

    joined_list[[key]] <- joined
  }

  return(joined_list)
}

#' Spatially Join GBIF Points to GADM Polygons
#'
#' Joins a list of GBIF point `sf` objects to cached GADM polygons.
#'
#' @param results_list Named list of GBIF `sf` point data by ISO2 code.
#' @param gadm_data Named list of GADM `sf` polygons (from `load_gadm()`).
#' @param level GADM administrative level used (for column naming).
#' @param use_planar Logical; whether to use planar projection (default: TRUE).
#'
#' @return Named list of joined `sf` objects with `GADM_name` column.
#' @export
gadm_join <- function(results_list,
                      gadm_data,
                      level = 1,
                      use_planar = TRUE) {

  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")

  if (use_planar) {
    sf::sf_use_s2(FALSE)
  } else {
    sf::sf_use_s2(TRUE)
  }

  joined_list <- list()

  for (cc in names(results_list)) {
    gbif_sf <- results_list[[cc]]

    # Skip empty results
    if (is.null(gbif_sf) || nrow(gbif_sf) == 0) {
      gbif_sf$GADM_name <- NA_character_
      if (!"iso2c" %in% names(gbif_sf)) gbif_sf$iso2c <- cc
      joined_list[[cc]] <- gbif_sf
      next
    }

    gadm_sf <- gadm_data[[cc]]
    if (is.null(gadm_sf)) {
      message(sprintf("❌ No GADM data available for %s", cc))
      gbif_sf$GADM_name <- NA_character_
      if (!"iso2c" %in% names(gbif_sf)) gbif_sf$iso2c <- cc
      joined_list[[cc]] <- gbif_sf
      next
    }

    # standardize join cols for this country
    picked   <- .pick_join_cols(gadm_sf, level)
    gadm_sf  <- picked$g
    name_col <- picked$name
    gid_col  <- picked$gid

    # harmonize CRS
    crs_target <- if (use_planar) 3857 else sf::st_crs(gadm_sf)
    gbif_sf <- sf::st_transform(gbif_sf, crs_target)
    gadm_sf <- sf::st_transform(gadm_sf, crs_target)

    # primary join: st_within
    joined <- tryCatch(
      sf::st_join(gbif_sf, gadm_sf[, c(name_col, gid_col)], join = sf::st_within),
      error = function(e) {
        message(sprintf("❌ st_within join failed for %s: %s", cc, e$message))
        gbif_sf
      }
    )

    # fallback: nearest feature if join failed or NA
    if (!(name_col %in% names(joined)) || all(is.na(joined[[name_col]]))) {
      message(sprintf("⚠️ No %s matches for %s — using nearest feature fallback", name_col, cc))
      joined <- tryCatch(
        sf::st_join(gbif_sf, gadm_sf[, c(name_col, gid_col)], join = sf::st_nearest_feature),
        error = function(e) {
          message(sprintf("❌ st_nearest_feature failed for %s: %s", cc, e$message))
          gbif_sf[[name_col]] <- NA_character_
          gbif_sf
        }
      )
    }

    # standardize to "GADM_name"
    if (name_col %in% names(joined)) {
      names(joined)[names(joined) == name_col] <- "GADM_name"
    } else {
      joined$GADM_name <- NA_character_
    }

    if (!"iso2c" %in% names(joined)) joined$iso2c <- cc

    joined_list[[cc]] <- joined
  }

  return(joined_list)
}

#' Spatially Join GBIF Points to EEZ Polygons (canonical ISO fields)
#'
#' Joins a list of GBIF point `sf` objects to Exclusive Economic Zones (EEZ)
#' polygons and appends a consistent set of attributes used downstream for
#' GRIIS/GBIF/WoRMS status resolution.
#'
#' The function **normalizes and guarantees** the following canonical columns on
#' the joined points (creating them if absent in the source EEZ layer):
#' - `geoname`  : Human-readable EEZ name (UTF-8 normalized)
#' - `mrgid`    : Marine Regions ID
#' - `SOVEREIGN1`: Sovereign name (UTF-8 normalized)
#' - `TERRITORY1`: Territory name (UTF-8 normalized)
#' - `ISO_SOV1` : ISO3 code for sovereign (when provided in EEZ)
#' - `ISO_Ter1` : ISO3 code for territory (when provided in EEZ)
#' - `iso3_best`: **Derived ISO3**, preferring `ISO_Ter1`, then `ISO_SOV1`,
#'                and finally name-based inference via \pkg{countrycode}
#'                (if available).
#'
#' EEZ attributes often vary by release (e.g., `ISO_TERR1`, `ISO_TERR`, case
#' variants). This function soft-maps those variants into the canonical names
#' above so downstream code can always rely on consistent fields.
#'
#' **Encoding:** All text attributes (`geoname`, `SOVEREIGN1`, `TERRITORY1`) are
#' normalized to valid UTF-8 internally to avoid Windows/locale export issues.
#'
#' @param results_list Named list of `sf` point objects, keyed by species or IDs.
#'   If a single `sf` is provided, it is wrapped into a length-1 list. If the list
#'   is unnamed, synthetic names are assigned.
#' @param use_planar Logical; if `TRUE` (default) disables spherical geometry
#'   (`sf_use_s2(FALSE)`) and joins in EPSG:3857; if `FALSE`, uses the layer CRS.
#' @param eez_cache_file Optional path to a local `.gpkg` to cache EEZ polygons
#'   retrieved via \pkg{mregions2}. Defaults to
#'   `file.path(rappdirs::user_cache_dir('biofetchR'),'biofetchR_eez_cache.gpkg')`.
#'
#' @return A named list of `sf` objects (same keys as input), with the canonical
#'   columns appended: `geoname`, `mrgid`, `SOVEREIGN1`, `TERRITORY1`,
#'   `ISO_SOV1`, `ISO_Ter1`, and `iso3_best`.
#'
#' @details
#' If an occurrence is not strictly within any EEZ polygon, a nearest-feature
#' join is attempted as a fallback. The canonical columns are guaranteed to exist
#' on the output (filled with `NA` when unavailable).
#'
#' @examples
#' \dontrun{
#' # Minimal example with one species data.frame -> sf
#' pts <- data.frame(
#'   decimalLongitude = c(-5, 10),
#'   decimalLatitude  = c(50, 35)
#' )
#' sf_pts <- sf::st_as_sf(pts, coords = c("decimalLongitude","decimalLatitude"), crs = 4326)
#'
#' out <- eez_join(list("Carcinus maenas" = sf_pts), use_planar = TRUE)
#' names(out[["Carcinus maenas"]])
#' # [1] ... "geoname" "mrgid" "SOVEREIGN1" "TERRITORY1" "ISO_SOV1" "ISO_Ter1" "iso3_best"
#' }
#'
#' @export
#' @importFrom sf st_as_sf st_transform st_crs st_join st_within st_nearest_feature
#' @importFrom sf st_make_valid read_sf write_sf sf_use_s2
#' @importFrom mregions2 mrp_get
#' @importFrom rappdirs user_cache_dir
eez_join <- function(results_list,
                     use_planar = TRUE,
                     eez_cache_file = NULL) {

  # -- helper: normalize character vectors to valid UTF-8 (best effort)
  utf8_safe <- function(x) {
    if (is.null(x)) return(x)
    y <- suppressWarnings(iconv(x, from = "",        to = "UTF-8", sub = ""))
    bad <- is.na(y)
    if (any(bad)) y[bad] <- suppressWarnings(iconv(x[bad], from = "latin1", to = "UTF-8", sub = ""))
    bad <- is.na(y)
    if (any(bad)) y[bad] <- suppressWarnings(iconv(x[bad], from = "CP1252", to = "UTF-8", sub = ""))
    y[is.na(y)] <- x[is.na(y)]
    y
  }

  # -- helper: canonicalize EEZ attribute names and compute iso3_best
  canonicalize_eez_attrs <- function(x) {
    nms <- names(x)
    first_present <- function(cands) cands[match(TRUE, cands %in% nms, nomatch = 0L)]

    src_SOVEREIGN1 <- first_present(c("SOVEREIGN1","Sovereign1","sovereign1"))
    src_TERRITORY1 <- first_present(c("TERRITORY1","Territory1","territory1"))
    src_ISO_SOV1   <- first_present(c("ISO_SOV1","ISO_Sov1","iso_sov1"))
    src_ISO_Ter1   <- first_present(c("ISO_Ter1","ISO_TER1","ISO_TERR1","ISO_TERR","iso_ter1","iso_terr1","iso_terr"))

    if (!"SOVEREIGN1" %in% nms) x$SOVEREIGN1 <- NA_character_
    if (!"TERRITORY1" %in% nms) x$TERRITORY1 <- NA_character_
    if (!"ISO_SOV1"   %in% nms) x$ISO_SOV1   <- NA_character_
    if (!"ISO_Ter1"   %in% nms) x$ISO_Ter1   <- NA_character_

    if (!is.null(src_SOVEREIGN1)) x$SOVEREIGN1 <- as.character(x[[src_SOVEREIGN1]])
    if (!is.null(src_TERRITORY1)) x$TERRITORY1 <- as.character(x[[src_TERRITORY1]])
    if (!is.null(src_ISO_SOV1))   x$ISO_SOV1   <- toupper(trimws(as.character(x[[src_ISO_SOV1]])))
    if (!is.null(src_ISO_Ter1))   x$ISO_Ter1   <- toupper(trimws(as.character(x[[src_ISO_Ter1]])))

    x$geoname    <- utf8_safe(x$geoname)
    x$SOVEREIGN1 <- utf8_safe(x$SOVEREIGN1)
    x$TERRITORY1 <- utf8_safe(x$TERRITORY1)

    iso3_best <- x$ISO_Ter1
    need <- is.na(iso3_best) | !nzchar(iso3_best)
    iso3_best[need] <- x$ISO_SOV1[need]

    need <- is.na(iso3_best) | !nzchar(iso3_best)
    if (any(need) && requireNamespace("countrycode", quietly = TRUE)) {
      cc <- suppressWarnings(countrycode::countrycode(x$TERRITORY1[need], origin = "country.name", destination = "iso3c"))
      iso3_best[need] <- toupper(cc)
      need <- is.na(iso3_best) | !nzchar(iso3_best)
      if (any(need)) {
        cc2 <- suppressWarnings(countrycode::countrycode(x$SOVEREIGN1[need], origin = "country.name", destination = "iso3c"))
        iso3_best[need] <- toupper(cc2)
      }
    }

    x$iso3_best <- toupper(iso3_best)
    x
  }

  # -- requirements & inputs
  if (!requireNamespace("mregions2", quietly = TRUE)) {
    stop("❌ Package 'mregions2' is required. Install: remotes::install_github('iobis/mregions2', subdir = 'R')")
  }

  if (inherits(results_list, "sf")) results_list <- list(results_list)
  if (is.null(names(results_list))) names(results_list) <- paste0("item_", seq_along(results_list))

  if (is.null(eez_cache_file)) {
    eez_cache_file <- file.path(rappdirs::user_cache_dir("biofetchR"), "biofetchR_eez_cache.gpkg")
  }
  if (!dir.exists(dirname(eez_cache_file))) dir.create(dirname(eez_cache_file), recursive = TRUE)

  # -- load or download EEZ polygons
  if (file.exists(eez_cache_file)) {
    message("🌐 Using cached EEZ polygons: ", eez_cache_file)
    eez_sf <- sf::read_sf(eez_cache_file)
  } else {
    message("🌐 Downloading EEZ polygons via 'mregions2::mrp_get()'...")
    eez_sf <- tryCatch(mregions2::mrp_get("eez"),
                       error = function(e) stop("❌ Failed to download EEZ polygons: ", e$message))
    sf::write_sf(eez_sf, eez_cache_file)
  }

  if (isTRUE(use_planar)) sf::sf_use_s2(FALSE) else sf::sf_use_s2(TRUE)

  eez_sf <- sf::st_make_valid(eez_sf)
  if (!("geoname" %in% names(eez_sf))) stop("EEZ polygons are missing the 'geoname' attribute.")

  # canonicalize attributes before join
  eez_sf <- canonicalize_eez_attrs(eez_sf)
  keep_cols <- intersect(c("geoname","mrgid","SOVEREIGN1","TERRITORY1","ISO_SOV1","ISO_Ter1","iso3_best"),
                         names(eez_sf))

  crs_target <- if (isTRUE(use_planar)) 3857 else sf::st_crs(eez_sf)
  eez_sf_tr  <- sf::st_transform(eez_sf, crs_target)

  # -- join loop
  joined_list <- vector("list", length(results_list)); names(joined_list) <- names(results_list)

  for (key in names(results_list)) {
    gbif_sf <- results_list[[key]]

    if (is.null(gbif_sf) || nrow(gbif_sf) == 0) {
      gbif_sf$geoname <- NA_character_
      joined_list[[key]] <- gbif_sf
      next
    }

    if (!inherits(gbif_sf, "sf")) {
      if (!all(c("decimalLongitude","decimalLatitude") %in% names(gbif_sf))) {
        warning(sprintf("Input '%s' is not sf and lacks lon/lat; skipping.", key))
        next
      }
      gbif_sf <- sf::st_as_sf(gbif_sf, coords = c("decimalLongitude","decimalLatitude"),
                              crs = 4326, remove = FALSE)
    }

    gbif_sf <- sf::st_transform(gbif_sf, crs_target)

    joined <- tryCatch(
      sf::st_join(gbif_sf, eez_sf_tr[, keep_cols, drop = FALSE], join = sf::st_within),
      error = function(e) { message(sprintf("❌ st_within join failed for %s: %s", utf8_safe(key), e$message)); gbif_sf }
    )

    if (!("geoname" %in% names(joined)) || all(is.na(joined$geoname))) {
      message(sprintf("⚠️ No EEZ matches for '%s' — trying nearest feature...", utf8_safe(key)))
      joined <- tryCatch(
        sf::st_join(gbif_sf, eez_sf_tr[, keep_cols, drop = FALSE], join = sf::st_nearest_feature),
        error = function(e) {
          message(sprintf("❌ st_nearest_feature failed for %s: %s", utf8_safe(key), e$message))
          gbif_sf$geoname <- NA_character_; gbif_sf
        }
      )
    }

    # normalize text fields
    for (cc in intersect(c("geoname","SOVEREIGN1","TERRITORY1"), names(joined))) {
      joined[[cc]] <- utf8_safe(joined[[cc]])
    }

    # guarantee canonical columns exist
    for (cc in c("geoname","mrgid","SOVEREIGN1","TERRITORY1","ISO_SOV1","ISO_Ter1","iso3_best")) {
      if (!(cc %in% names(joined))) joined[[cc]] <- NA
    }

    joined_list[[key]] <- joined
  }

  joined_list
}

eez_join <- function(results_list,
                     use_planar = TRUE,
                     eez_cache_file = NULL) {

  utf8_safe <- function(x) {
    if (is.null(x)) return(x)
    y <- suppressWarnings(iconv(x, from = "",        to = "UTF-8", sub = ""))
    bad <- is.na(y)
    if (any(bad)) y[bad] <- suppressWarnings(iconv(x[bad], from = "latin1", to = "UTF-8", sub = ""))
    bad <- is.na(y)
    if (any(bad)) y[bad] <- suppressWarnings(iconv(x[bad], from = "CP1252", to = "UTF-8", sub = ""))
    y[is.na(y)] <- x[is.na(y)]
    y
  }

  # allow a single sf
  if (inherits(results_list, "sf")) results_list <- list(results_list)
  if (is.null(names(results_list))) names(results_list) <- paste0("item_", seq_along(results_list))

  if (!requireNamespace("mregions2", quietly = TRUE)) {
    stop("❌ Package 'mregions2' is required. Install with remotes::install_github('iobis/mregions2', subdir = 'R')")
  }

  if (is.null(eez_cache_file)) {
    eez_cache_file <- file.path(rappdirs::user_cache_dir("biofetchR"), "biofetchR_eez_cache.gpkg")
  }
  if (!dir.exists(dirname(eez_cache_file))) dir.create(dirname(eez_cache_file), recursive = TRUE)

  # load / fetch
  if (file.exists(eez_cache_file)) {
    message("🌐 Using cached EEZ polygons: ", eez_cache_file)
    eez_sf <- sf::read_sf(eez_cache_file)
  } else {
    message("🌐 Downloading EEZ polygons via 'mregions2::mrp_get()'...")
    eez_sf <- tryCatch(mregions2::mrp_get("eez"), error = function(e) {
      stop("❌ Failed to download EEZ polygons: ", e$message)
    })
    sf::write_sf(eez_sf, eez_cache_file)
  }

  if (isTRUE(use_planar)) sf::sf_use_s2(FALSE) else sf::sf_use_s2(TRUE)
  eez_sf <- sf::st_make_valid(eez_sf)

  if (!("geoname" %in% names(eez_sf))) stop("EEZ polygons are missing the 'geoname' attribute.")
  eez_sf$geoname <- utf8_safe(eez_sf$geoname)

  # ---- Rename variants to canonical -----------------------------------------
  variants <- list(
    SOVEREIGN1 = c("SOVEREIGN1","Sovereign1","sovereign1"),
    TERRITORY1 = c("TERRITORY1","Territory1","territory1"),
    ISO_SOV1   = c("ISO_SOV1","ISO_Sov1","iso_sov1"),
    ISO_Ter1   = c("ISO_Ter1","ISO_TER1","iso_ter1","ISO_TERR1","ISO_TERR")
  )
  for (canon in names(variants)) {
    hits <- intersect(variants[[canon]], names(eez_sf))
    if (length(hits)) {
      if (!(canon %in% names(eez_sf))) eez_sf[[canon]] <- eez_sf[[hits[1]]]
      # if the canonical already exists, leave it; otherwise set from first variant
    } else {
      eez_sf[[canon]] <- NA_character_
    }
  }

  # Provide "best" ISO for downstream matching
  iso3_best <- eez_sf$ISO_SOV1
  need <- is.na(iso3_best) | !nzchar(iso3_best)
  iso3_best[need] <- eez_sf$ISO_Ter1[need]
  eez_sf$iso3_best <- toupper(iso3_best)

  if (requireNamespace("countrycode", quietly = TRUE)) {
    eez_sf$iso2_best <- suppressWarnings(countrycode::countrycode(eez_sf$iso3_best, "iso3c", "iso2c"))
    eez_sf$iso2_best <- toupper(eez_sf$iso2_best)
  } else {
    eez_sf$iso2_best <- NA_character_
  }

  # only keep needed columns + geometry
  keep <- c("geoname","mrgid","SOVEREIGN1","TERRITORY1","ISO_SOV1","ISO_Ter1","iso3_best","iso2_best")
  keep <- intersect(keep, names(eez_sf))
  crs_target <- if (isTRUE(use_planar)) 3857 else sf::st_crs(eez_sf)
  eez_sf_tr  <- sf::st_transform(eez_sf[, keep, drop = FALSE], crs_target)

  joined_list <- list()

  for (key in names(results_list)) {
    gbif_sf <- results_list[[key]]
    if (is.null(gbif_sf) || nrow(gbif_sf) == 0) {
      gbif_sf$geoname <- NA_character_
      joined_list[[key]] <- gbif_sf
      next
    }
    if (!inherits(gbif_sf, "sf")) {
      if (!all(c("decimalLongitude","decimalLatitude") %in% names(gbif_sf))) {
        warning(sprintf("Input '%s' is not sf and lacks lon/lat; skipping.", key))
        next
      }
      gbif_sf <- sf::st_as_sf(gbif_sf, coords = c("decimalLongitude","decimalLatitude"),
                              crs = 4326, remove = FALSE)
    }
    gbif_sf <- sf::st_transform(gbif_sf, crs_target)

    joined <- tryCatch(
      sf::st_join(gbif_sf, eez_sf_tr, join = sf::st_within),
      error = function(e) {
        message(sprintf("❌ st_within join failed for %s: %s", utf8_safe(key), e$message))
        gbif_sf
      }
    )
    if (!("geoname" %in% names(joined)) || all(is.na(joined$geoname))) {
      message(sprintf("⚠️ No EEZ matches for '%s' — trying nearest feature...", utf8_safe(key)))
      joined <- tryCatch(
        sf::st_join(gbif_sf, eez_sf_tr, join = sf::st_nearest_feature),
        error = function(e) {
          gbif_sf$geoname <- NA_character_
          gbif_sf
        }
      )
    }

    # UTF-8 cleanup
    for (cc in intersect(c("geoname","SOVEREIGN1","TERRITORY1"), names(joined))) {
      joined[[cc]] <- utf8_safe(joined[[cc]])
    }

    # guarantee canonical columns exist
    for (cc in setdiff(c("geoname","mrgid","SOVEREIGN1","TERRITORY1","ISO_SOV1","ISO_Ter1","iso3_best","iso2_best"),
                       names(joined))) {
      joined[[cc]] <- NA
    }

    joined_list[[key]] <- joined
  }

  joined_list
}


#' Harmonize column types across a list of `sf`/data frames
#'
#' Ensures a common schema so `dplyr::bind_rows()` won’t fail on
#' type clashes (e.g., `logical` vs `character`). The function:
#' (1) takes the union of all columns, (2) adds any missing columns
#' as `NA_character_`, (3) resolves conflicting classes by coercing
#' to **character** (non-geometry columns only), and (4) assigns a
#' common CRS to `sf` objects using the CRS of the largest element
#' (via `sf::st_set_crs()`; it does **not** reproject).
#'
#' @param sf_list A list of `sf` objects or plain data frames to be harmonized.
#'
#' @return A list of objects (same length/order as input) with:
#' - the same set of columns in the same order, and
#' - non-geometry columns coerced to a compatible type (character when in doubt),
#' - `sf` objects assigned a common CRS (if any was available).
#'
#' @details
#' - Geometry columns are left untouched (identified via `sf::st_geometry_name()`).
#' - When multiple non-`sfc` classes are observed for a column, the target
#'   class is set to `"character"`.
#' - Missing columns are created as `NA_character_` (safe default).
#' - CRS is **set** (not transformed) on `sf` objects to the CRS of the
#'   largest input (by row count), if available.
#'
#' @examples
#' \dontrun{
#' # Suppose you have a list of per-species sf frames with slight schema differences:
#' clean_list <- harmonize_column_types(per_species_sf)
#' combined   <- dplyr::bind_rows(clean_list, .id = "species_eez")
#' }
#'
#' @seealso [dplyr::bind_rows()], [sf::st_set_crs()], [sf::st_crs()]
#'
#' @importFrom sf st_crs st_geometry_name st_set_crs
#' @keywords internal
#' @noRd
harmonize_column_types <- function(sf_list){
  # drop NULLs early
  sf_list <- sf_list[!vapply(sf_list, is.null, logical(1))]
  if (!length(sf_list)) return(sf_list)

  # union of columns across all frames
  all_cols <- unique(unlist(lapply(sf_list, names), use.names = FALSE))

  # pick a reference CRS (first non-NULL by largest nrow)
  ref_idx <- which.max(vapply(sf_list, nrow, integer(1)))
  ref_crs <- tryCatch(sf::st_crs(sf_list[[ref_idx]]), error = function(e) NULL)

  # helper: class of a column in an sf (treat geometry as 'sfc')
  col_class <- function(x, col) {
    if (!col %in% names(x)) return(NA_character_)
    if (inherits(x, "sf") && col == sf::st_geometry_name(x)) return("sfc")
    class(x[[col]])[1]
  }

  # build target class per column:
  # - if multiple non-NA, non-'sfc' classes appear, target == "character"
  # - otherwise keep the single observed class (default to character if unknown)
  class_matrix <- lapply(all_cols, function(col) {
    vapply(sf_list, col_class, character(1), col = col)
  })
  names(class_matrix) <- all_cols

  target_class <- vapply(class_matrix, function(v) {
    v2 <- unique(na.omit(v[v != "sfc"]))
    if (!length(v2)) "character" else if (length(v2) == 1) v2 else "character"
  }, character(1))
  # geometry stays sfc; we will skip coercion for the geometry column of each object

  # coerce each element to the target schema
  out <- lapply(sf_list, function(x) {
    # add missing columns as NA_character_ (safe default)
    miss <- setdiff(all_cols, names(x))
    for (m in miss) x[[m]] <- NA_character_

    # coerce non-geometry columns
    gcol <- if (inherits(x, "sf")) sf::st_geometry_name(x) else NA_character_
    for (nm in setdiff(all_cols, gcol)) {
      tgt <- target_class[[nm]]
      if (tgt == "character") {
        x[[nm]] <- as.character(x[[nm]])
      } else if (tgt == "integer") {
        x[[nm]] <- suppressWarnings(as.integer(x[[nm]]))
      } else if (tgt == "numeric") {
        x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
      } else if (tgt == "logical") {
        x[[nm]] <- as.logical(x[[nm]])
      } else if (tgt == "factor") {
        x[[nm]] <- as.factor(x[[nm]])
      } else {
        x[[nm]] <- as.character(x[[nm]])  # safe fallback
      }
    }

    # unify CRS for sf objects (set, not transform)
    if (inherits(x, "sf") && !is.null(ref_crs)) {
      try(x <- sf::st_set_crs(x, ref_crs), silent = TRUE)
    }
    x
  })

  out
}

eez_join <- function(results_list, use_planar = TRUE, eez_cache_file = NULL) {
  if (!requireNamespace("mregions2", quietly = TRUE)) {
    stop("❌ Package 'mregions2' is required. Install with: remotes::install_github('iobis/mregions2', subdir='R')")
  }

  # Cache path
  if (is.null(eez_cache_file)) {
    eez_cache_file <- file.path(rappdirs::user_cache_dir("biofetchR"), "biofetchR_eez_cache.gpkg")
  }
  if (!dir.exists(dirname(eez_cache_file))) dir.create(dirname(eez_cache_file), recursive = TRUE)

  # Load/download EEZ polygons
  if (file.exists(eez_cache_file)) {
    message("🌐 Using cached EEZ polygons: ", eez_cache_file)
    eez_sf <- sf::read_sf(eez_cache_file)
  } else {
    message("🌐 Downloading EEZ polygons via mregions2::mrp_get('eez')...")
    eez_sf <- mregions2::mrp_get("eez")
    sf::write_sf(eez_sf, eez_cache_file)
  }

  if (use_planar) sf::sf_use_s2(FALSE) else sf::sf_use_s2(TRUE)
  eez_sf <- sf::st_make_valid(eez_sf)

  # Keep useful attributes for ISO matching
  wanted <- c(
    "geoname","mrgid",
    "SOVEREIGN1","TERRITORY1",
    "ISO_SOV1","ISO_Ter1","ISO_SOV2","ISO_Ter2",
    "ISO2_SOV1","ISO2_Ter1","ISO2_SOV2","ISO2_Ter2"
  )
  have <- intersect(wanted, names(eez_sf))
  eez_keep <- eez_sf[, have, drop = FALSE]

  # Pre-compute best ISO3/ISO2 on polygons (so they copy onto points)
  coalesce_char <- function(...) {
    args <- list(...)
    out <- rep(NA_character_, nrow(eez_keep))
    for (v in args) {
      if (is.null(v)) next
      v <- as.character(v)
      fill <- is.na(out) | !nzchar(out)
      out[fill] <- v[fill]
    }
    out
  }

  # Best ISO3: prefer Territory, then Sovereign; derive from name if needed
  eez_keep$iso3_best <- toupper(coalesce_char(eez_keep$ISO_Ter1, eez_keep$ISO_SOV1))
  if (any(!nzchar(eez_keep$iso3_best)) && requireNamespace("countrycode", quietly = TRUE)) {
    idx <- !nzchar(eez_keep$iso3_best)
    nm  <- coalesce_char(eez_keep$TERRITORY1, eez_keep$SOVEREIGN1)
    eez_keep$iso3_best[idx] <- toupper(suppressWarnings(
      countrycode::countrycode(nm[idx], "country.name", "iso3c")
    ))
  }

  # Best ISO2 (if present; otherwise derive from iso3_best if countrycode available)
  if (any(c("ISO2_Ter1","ISO2_SOV1") %in% names(eez_keep))) {
    eez_keep$iso2_best <- toupper(coalesce_char(eez_keep$ISO2_Ter1, eez_keep$ISO2_SOV1))
  } else if (requireNamespace("countrycode", quietly = TRUE)) {
    eez_keep$iso2_best <- toupper(suppressWarnings(
      countrycode::countrycode(eez_keep$iso3_best, "iso3c", "iso2c")
    ))
  } else {
    eez_keep$iso2_best <- NA_character_
  }

  # Join each species' points to EEZ polygons (carry all fields above)
  joined_list <- list()
  for (key in names(results_list)) {
    gbif_sf <- results_list[[key]]
    if (is.null(gbif_sf) || nrow(gbif_sf) == 0) {
      gbif_sf$geoname <- NA_character_
      joined_list[[key]] <- gbif_sf
      next
    }

    crs_target <- if (use_planar) 3857 else sf::st_crs(eez_keep)
    gbif_sf <- sf::st_transform(gbif_sf, crs_target)
    eez_use <- sf::st_transform(eez_keep, crs_target)

    joined <- tryCatch(
      sf::st_join(gbif_sf, eez_use, join = sf::st_within, left = TRUE),
      error = function(e) gbif_sf
    )

    # Rescue with nearest if within-join finds nothing
    if (!("geoname" %in% names(joined)) || all(is.na(joined$geoname))) {
      joined <- tryCatch(
        sf::st_join(gbif_sf, eez_use, join = sf::st_nearest_feature, left = TRUE),
        error = function(e) { gbif_sf$geoname <- NA_character_; gbif_sf }
      )
    }

    joined_list[[key]] <- joined
  }

  joined_list
}

#' Join GBIF points to EEZ polygons (with robust ISO fields)
#'
#' Loads (or caches) EEZ polygons via {mregions2}, joins a list of GBIF
#' point `sf` objects to EEZs, and returns the same list with EEZ attributes
#' added. The function normalizes common EEZ attribute names across
#' multiple releases and computes `iso3_best` and `iso2_best` for downstream
#' status mapping (GRIIS/GBIF distributions).
#'
#' @param results_list Named `list` of `sf` data frames (points). Names are kept.
#' @param use_planar Logical; if `TRUE` (default) use planar ops (s2 off).
#' @param eez_cache_file Optional path to a `.gpkg` cache for EEZ polygons.
#'   Defaults to a user cache directory.
#'
#' @return A `list` of `sf` objects with added columns:
#'   `geoname`, `mrgid`, `territory1`, `sovereign1`, `iso3_best`, `iso2_best`.
#' @export
eez_join <- function(results_list, use_planar = TRUE, eez_cache_file = NULL) {
  if (!requireNamespace("mregions2", quietly = TRUE)) {
    stop("Package 'mregions2' is required. Install with: remotes::install_github('iobis/mregions2', subdir='R')")
  }

  # Cache path
  if (is.null(eez_cache_file)) {
    eez_cache_file <- file.path(rappdirs::user_cache_dir("biofetchR"), "biofetchR_eez_cache.gpkg")
  }
  if (!dir.exists(dirname(eez_cache_file))) dir.create(dirname(eez_cache_file), recursive = TRUE)

  # Load/download EEZ polygons
  if (file.exists(eez_cache_file)) {
    message("Using cached EEZ polygons: ", eez_cache_file)
    eez_sf <- sf::read_sf(eez_cache_file)
  } else {
    message("Downloading EEZ polygons via mregions2::mrp_get('eez')...")
    eez_sf <- mregions2::mrp_get("eez")
    sf::write_sf(eez_sf, eez_cache_file)
  }

  if (use_planar) sf::sf_use_s2(FALSE) else sf::sf_use_s2(TRUE)
  eez_sf <- sf::st_make_valid(eez_sf)

  # ---- helpers (case-insensitive pick & string clean) -----------------------
  pick_ci <- function(x, candidates) {
    nx <- names(x)
    for (cand in candidates) {
      hit <- which(tolower(nx) == tolower(cand))
      if (length(hit)) return(nx[hit[1]])
    }
    NULL
  }
  coalesce_char <- function(...) {
    args <- list(...)
    n <- nrow(eez_sf)
    out <- rep(NA_character_, n)
    for (v in args) {
      if (is.null(v)) next
      v <- as.character(v)
      fill <- is.na(out) | !nzchar(out)
      out[fill] <- v[fill]
    }
    out
  }
  clean_geoname_to_country <- function(s) {
    s <- as.character(s)
    s <- gsub("Exclusive Economic Zone", "", s, ignore.case = TRUE)
    s <- gsub("Joint regime area", "", s, ignore.case = TRUE)
    s <- gsub("Overlapping claim area", "", s, ignore.case = TRUE)
    s <- gsub("\\(.*?\\)", "", s)
    s <- gsub("-{1,2}\\s*EEZ", "", s, ignore.case = TRUE)
    s <- gsub("\\s+EEZ", "", s, ignore.case = TRUE)
    s <- gsub("^\\s+|\\s+$", "", s)
    s
  }

  # ---- normalize column handles (case-insensitive) --------------------------
  geoname_col   <- pick_ci(eez_sf, c("geoname","name","geoname_lng","geoname_long"))
  mrgid_col     <- pick_ci(eez_sf, c("mrgid","MRGID"))
  terr1_col     <- pick_ci(eez_sf, c("TERRITORY1","territory1","territory"))
  sov1_col      <- pick_ci(eez_sf, c("SOVEREIGN1","sovereign1","sovereign"))
  iso_ter1_col  <- pick_ci(eez_sf, c("ISO_Ter1","iso_ter1"))
  iso_sov1_col  <- pick_ci(eez_sf, c("ISO_SOV1","iso_sov1"))
  iso2_ter1_col <- pick_ci(eez_sf, c("ISO2_Ter1","iso2_ter1"))
  iso2_sov1_col <- pick_ci(eez_sf, c("ISO2_SOV1","iso2_sov1"))

  keep_cols <- unique(na.omit(c(
    geoname_col, mrgid_col, terr1_col, sov1_col,
    iso_ter1_col, iso_sov1_col, iso2_ter1_col, iso2_sov1_col
  )))
  eez_keep <- eez_sf[, keep_cols, drop = FALSE]

  # Standardize names in eez_keep
  if (!is.null(geoname_col))    names(eez_keep)[names(eez_keep) == geoname_col]    <- "geoname"
  if (!is.null(mrgid_col))      names(eez_keep)[names(eez_keep) == mrgid_col]      <- "mrgid"
  if (!is.null(terr1_col))      names(eez_keep)[names(eez_keep) == terr1_col]      <- "territory1"
  if (!is.null(sov1_col))       names(eez_keep)[names(eez_keep) == sov1_col]       <- "sovereign1"
  if (!is.null(iso_ter1_col))   names(eez_keep)[names(eez_keep) == iso_ter1_col]   <- "iso_ter1"
  if (!is.null(iso_sov1_col))   names(eez_keep)[names(eez_keep) == iso_sov1_col]   <- "iso_sov1"
  if (!is.null(iso2_ter1_col))  names(eez_keep)[names(eez_keep) == iso2_ter1_col]  <- "iso2_ter1"
  if (!is.null(iso2_sov1_col))  names(eez_keep)[names(eez_keep) == iso2_sov1_col]  <- "iso2_sov1"

  # Some datasets only have 'name'—promote to geoname
  if (!"geoname" %in% names(eez_keep)) {
    eez_keep$geoname <- if ("name" %in% names(eez_keep)) eez_keep$name else NA_character_
  }

  # ---- compute iso3_best / iso2_best on the polygons ------------------------
  iso3_direct <- coalesce_char(eez_keep$iso_ter1, eez_keep$iso_sov1)

  if (requireNamespace("countrycode", quietly = TRUE)) {
    need <- is.na(iso3_direct) | !nzchar(iso3_direct)
    if (any(need)) {
      nm_guess <- clean_geoname_to_country(
        coalesce_char(eez_keep$territory1, eez_keep$sovereign1, eez_keep$geoname)
      )
      iso3_from_name <- suppressWarnings(countrycode::countrycode(nm_guess, "country.name", "iso3c"))
      iso3_direct[need] <- toupper(iso3_from_name[need])
    }
  }
  eez_keep$iso3_best <- toupper(iso3_direct)

  iso2_direct <- coalesce_char(eez_keep$iso2_ter1, eez_keep$iso2_sov1)
  if ((is.null(iso2_direct) || all(is.na(iso2_direct))) &&
      requireNamespace("countrycode", quietly = TRUE)) {
    iso2_direct <- suppressWarnings(countrycode::countrycode(eez_keep$iso3_best, "iso3c", "iso2c"))
  }
  eez_keep$iso2_best <- toupper(iso2_direct)

  # ---- join points → EEZ -----------------------------------------------------
  joined_list <- list()
  for (key in names(results_list)) {
    gbif_sf <- results_list[[key]]
    if (is.null(gbif_sf) || nrow(gbif_sf) == 0) {
      gbif_sf$geoname <- NA_character_
      joined_list[[key]] <- gbif_sf
      next
    }

    crs_target <- if (use_planar) 3857 else sf::st_crs(eez_keep)
    gbif_sf <- sf::st_transform(gbif_sf, crs_target)
    eez_use <- sf::st_transform(eez_keep, crs_target)

    joined <- tryCatch(
      sf::st_join(
        gbif_sf,
        eez_use[, c("geoname","mrgid","territory1","sovereign1","iso3_best","iso2_best")],
        join = sf::st_within,
        left = TRUE
      ),
      error = function(e) gbif_sf
    )

    if (!("geoname" %in% names(joined)) || all(is.na(joined$geoname))) {
      joined <- tryCatch(
        sf::st_join(
          gbif_sf,
          eez_use[, c("geoname","mrgid","territory1","sovereign1","iso3_best","iso2_best")],
          join = sf::st_nearest_feature,
          left = TRUE
        ),
        error = function(e) { gbif_sf$geoname <- NA_character_; gbif_sf }
      )
    }

    joined_list[[key]] <- joined
  }

  return(joined_list)
}

#' Join GBIF points to EEZ polygons (robust, with s2 fallback)
#'
#' Reads/caches the EEZ layer, fixes invalid geometries, harmonizes ISO columns,
#' and joins each sf in `results_list` to EEZ polygons. If any s2 error occurs,
#' the join transparently falls back to planar (sf_use_s2(FALSE)).
#'
#' @param results_list named list of sf point data (e.g., species -> sf)
#' @param use_planar logical; if TRUE, run joins with s2 disabled from the start
#' @param eez_cache_file optional path to a cached GPKG of the EEZ layer
#' @return list of sf with EEZ attributes (geoname, mrgid, iso3_best, iso2_best)
eez_join <- function(results_list, use_planar = TRUE, eez_cache_file = NULL) {
  if (!requireNamespace("mregions2", quietly = TRUE)) {
    stop("Package 'mregions2' is required. Install with: remotes::install_github('iobis/mregions2', subdir='R')")
  }

  # cache path
  if (is.null(eez_cache_file)) {
    eez_cache_file <- file.path(rappdirs::user_cache_dir("biofetchR"), "biofetchR_eez_cache.gpkg")
  }
  if (!dir.exists(dirname(eez_cache_file))) dir.create(dirname(eez_cache_file), recursive = TRUE)

  # load or download
  if (file.exists(eez_cache_file)) {
    message("Using cached EEZ polygons: ", eez_cache_file)
    eez_sf <- sf::read_sf(eez_cache_file)
  } else {
    message("Downloading EEZ polygons via mregions2::mrp_get('eez')...")
    eez_sf <- mregions2::mrp_get("eez")
    sf::write_sf(eez_sf, eez_cache_file)
  }

  # fix invalid geoms (prefer lwgeom if present; else buffer(0) trick)
  eez_sf <- suppressWarnings(sf::st_make_valid(eez_sf))
  bad <- !sf::st_is_valid(eez_sf)
  if (any(bad)) {
    if (requireNamespace("lwgeom", quietly = TRUE)) {
      eez_sf[bad, ] <- lwgeom::st_make_valid(eez_sf[bad, ])
    } else {
      eez_sf <- sf::st_buffer(eez_sf, 0)
    }
  }

  # helper: case-insensitive picker + safe column getter (no warnings)
  pick_ci <- function(x, candidates) {
    nx <- names(x); for (cand in candidates) {
      hit <- which(tolower(nx) == tolower(cand)); if (length(hit)) return(nx[hit[1]])
    }; NULL
  }
  gcol <- function(df, nm) if (nm %in% names(df)) df[[nm]] else NULL

  # normalize handles
  geoname_col   <- pick_ci(eez_sf, c("geoname","name","geoname_lng","geoname_long"))
  mrgid_col     <- pick_ci(eez_sf, c("mrgid","MRGID"))
  terr1_col     <- pick_ci(eez_sf, c("TERRITORY1","territory1","territory"))
  sov1_col      <- pick_ci(eez_sf, c("SOVEREIGN1","sovereign1","sovereign"))
  iso_ter1_col  <- pick_ci(eez_sf, c("ISO_Ter1","iso_ter1"))
  iso_sov1_col  <- pick_ci(eez_sf, c("ISO_SOV1","iso_sov1"))
  iso2_ter1_col <- pick_ci(eez_sf, c("ISO2_Ter1","iso2_ter1"))
  iso2_sov1_col <- pick_ci(eez_sf, c("ISO2_SOV1","iso2_sov1"))

  keep_cols <- unique(na.omit(c(geoname_col, mrgid_col,
                                terr1_col,  sov1_col,
                                iso_ter1_col, iso_sov1_col,
                                iso2_ter1_col, iso2_sov1_col)))
  eez_keep <- eez_sf[, keep_cols, drop = FALSE]

  # standardize names
  ren <- function(old, new) if (!is.null(old) && old %in% names(eez_keep)) names(eez_keep)[names(eez_keep)==old] <<- new
  ren(geoname_col, "geoname"); ren(mrgid_col, "mrgid"); ren(terr1_col, "territory1")
  ren(sov1_col, "sovereign1"); ren(iso_ter1_col, "iso_ter1"); ren(iso_sov1_col, "iso_sov1")
  ren(iso2_ter1_col, "iso2_ter1"); ren(iso2_sov1_col, "iso2_sov1")
  if (!"geoname" %in% names(eez_keep) && "name" %in% names(eez_keep)) eez_keep$geoname <- eez_keep$name

  # compute iso3_best / iso2_best without warning on missing cols
  coalesce_char <- function(...) {
    args <- list(...); n <- nrow(eez_keep); out <- rep(NA_character_, n)
    for (v in args) {
      if (is.null(v)) next
      v <- as.character(v); fill <- is.na(out) | !nzchar(out); out[fill] <- v[fill]
    }; out
  }
  iso3_direct <- coalesce_char(gcol(eez_keep, "iso_ter1"), gcol(eez_keep, "iso_sov1"))
  if (requireNamespace("countrycode", quietly = TRUE)) {
    need <- is.na(iso3_direct) | !nzchar(iso3_direct)
    if (any(need)) {
      nm_guess <- coalesce_char(gcol(eez_keep, "territory1"), gcol(eez_keep, "sovereign1"), gcol(eez_keep, "geoname"))
      nm_guess <- gsub("Exclusive Economic Zone|Joint regime area|Overlapping claim area|\\(.*?\\)|-{1,2}\\s*EEZ|\\s+EEZ", "", nm_guess, ignore.case = TRUE)
      nm_guess <- trimws(nm_guess)
      iso3_from_name <- suppressWarnings(countrycode::countrycode(nm_guess, "country.name", "iso3c"))
      iso3_direct[need] <- toupper(iso3_from_name[need])
    }
  }
  eez_keep$iso3_best <- toupper(iso3_direct)

  iso2_direct <- coalesce_char(gcol(eez_keep, "iso2_ter1"), gcol(eez_keep, "iso2_sov1"))
  if ((is.null(iso2_direct) || all(is.na(iso2_direct))) && requireNamespace("countrycode", quietly = TRUE)) {
    iso2_direct <- suppressWarnings(countrycode::countrycode(eez_keep$iso3_best, "iso3c", "iso2c"))
  }
  eez_keep$iso2_best <- toupper(iso2_direct)

  # s2 mode handling + auto-fallback on error
  old_s2 <- sf::sf_use_s2()
  on.exit(try(sf::sf_use_s2(old_s2), silent = TRUE), add = TRUE)
  if (isTRUE(use_planar)) sf::sf_use_s2(FALSE)

  joined_list <- list()
  for (key in names(results_list)) {
    gbif_sf <- results_list[[key]]
    if (is.null(gbif_sf) || !nrow(gbif_sf)) {
      gbif_sf$geoname <- NA_character_
      joined_list[[key]] <- gbif_sf
      next
    }

    # both layers in the same CRS; try s2 if enabled, else planar; fallback if needed
    try_join <- function(force_planar = FALSE) {
      if (force_planar) sf::sf_use_s2(FALSE)
      crs_target <- if (!force_planar && !isTRUE(use_planar)) sf::st_crs(eez_keep) else 3857
      gb  <- try(suppressWarnings(sf::st_transform(gbif_sf, crs_target)), silent = TRUE)
      ply <- try(suppressWarnings(sf::st_transform(eez_keep, crs_target)), silent = TRUE)
      if (inherits(gb, "try-error") || inherits(ply, "try-error")) stop("transform_failed")
      out <- try(suppressWarnings(
        sf::st_join(gb, ply[, c("geoname","mrgid","territory1","sovereign1","iso3_best","iso2_best")],
                    join = sf::st_within, left = TRUE)
      ), silent = TRUE)
      if (inherits(out, "try-error") || (!("geoname" %in% names(out)) || all(is.na(out$geoname)))) {
        out <- try(suppressWarnings(
          sf::st_join(gb, ply[, c("geoname","mrgid","territory1","sovereign1","iso3_best","iso2_best")],
                      join = sf::st_nearest_feature, left = TRUE)
        ), silent = TRUE)
      }
      if (inherits(out, "try-error")) stop(attr(out, "condition")$message)
      out
    }

    res <- try(try_join(force_planar = FALSE), silent = TRUE)
    if (inherits(res, "try-error")) {
      # any s2/wk failure -> hard fallback to planar
      res <- try(try_join(force_planar = TRUE), silent = TRUE)
      if (inherits(res, "try-error")) {
        # as a last resort, return input with empty columns
        gbif_sf$geoname <- NA_character_
        gbif_sf$mrgid <- NA_integer_
        gbif_sf$iso3_best <- NA_character_
        gbif_sf$iso2_best <- NA_character_
        res <- gbif_sf
      }
    }
    joined_list[[key]] <- res
  }

  joined_list
}

#' Join GBIF points to Marine Regions overlays (robust, with s2 fallback)
#'
#' Reads/caches a Marine Regions polygon layer (EEZ/LME/IHO/MEOW/FAO, etc.),
#' fixes invalid geometries, harmonizes name/id columns, and joins each sf in
#' `results_list` to the polygons. If any s2 error occurs, the join
#' transparently falls back to planar (sf_use_s2(FALSE)).
#'
#' @param results_list Named list of sf point data (e.g., species -> sf).
#' @param overlay Character; which overlay family to use. One of
#'   c("eez","lme","iho","meow","fao","custom"). Only affects defaults if
#'   `overlay_layer` is not supplied.
#' @param overlay_layer Character; exact layer id for Marine Regions. If given,
#'   it is used directly. If NULL, a canonical id is chosen for the `overlay`.
#'   Examples: "eez", "lme", "iho", "meow", "fao".
#' @param overlay_sf Optional `sf` polygon object (preloaded overlay). If
#'   provided, no download/caching is attempted.
#' @param use_planar Logical; if TRUE, run joins with s2 disabled from the
#'   start (planar ops). Regardless, the function auto-falls back to planar
#'   if a spherical join fails.
#' @param cache_dir Directory to cache the overlay as a GPKG. Defaults to a
#'   user cache directory via {rappdirs}.
#'
#' @return Named list of sf with overlay attributes standardized as:
#'   - `geoname` (string): human-readable region label
#'   - `mrgid`   (integer): region id when present on overlay (e.g., EEZ)
#'   - `iso3_best`, `iso2_best` (character): best-effort ISO codes (EEZ only or NA)
#'
#' @details
#' - Uses **mregions2** if available (preferred: `mrp_list`/`mrp_get`), falls
#'   back to **mregions** (`mr_layers`/`mr_as_sf`) otherwise.
#' - For non-EEZ overlays, `iso3_best`/`iso2_best` will likely be `NA`; your
#'   pipeline’s `ensure_iso3_on_points()` will still fill ISO by overlaying
#'   a world-country layer later.
#'
#' @examples
#' # joined <- overlay_join(list("Carcinus maenas" = pts_sf), overlay = "eez")
overlay_join <- function(results_list,
                         overlay = c("eez","lme","iho","meow","fao","custom"),
                         overlay_layer = NULL,
                         overlay_sf = NULL,
                         use_planar = TRUE,
                         cache_dir = rappdirs::user_cache_dir("biofetchR")) {
  overlay <- match.arg(overlay)

  `%||%` <- function(a, b) if (is.null(a)) b else a

  # ---- helpers to discover & fetch overlays ---------------------------------
  list_overlays <- function(pattern = NULL) {
    has_mr2 <- requireNamespace("mregions2", quietly = TRUE)
    has_mr  <- requireNamespace("mregions",  quietly = TRUE)

    if (has_mr2 && "mrp_list" %in% getNamespaceExports("mregions2")) {
      L <- mregions2::mrp_list()
      L <- as.data.frame(L)
      if (!"layername" %in% names(L) && "layer" %in% names(L)) L$layername <- L$layer
    } else if (has_mr && "mr_layers" %in% getNamespaceExports("mregions")) {
      L <- mregions::mr_layers()
      if (!"layername" %in% names(L) && "name" %in% names(L)) L$layername <- L$name
    } else {
      stop("Need mregions2 (mrp_list/mrp_get) or mregions (mr_layers/mr_as_sf).")
    }

    if (!is.null(pattern)) {
      keep <- grepl(pattern, L$title, ignore.case = TRUE) |
        ("layername" %in% names(L) && grepl(pattern, L$layername, ignore.case = TRUE))
      L <- L[keep, , drop = FALSE]
    }
    L
  }

  get_overlay_sf <- function(id_or_family, which = 1, cache_dir = cache_dir, quiet = FALSE) {
    # If a direct id works, use it; else search catalog by family/id.
    L <- try(list_overlays(), silent = TRUE)
    layer <- NA_character_

    if (!inherits(L, "try-error")) {
      # choose canonical id for family if none given
      canon <- switch(tolower(id_or_family),
                      eez  = "eez",
                      lme  = "lme",
                      iho  = "iho",
                      meow = "meow",
                      fao  = "fao",
                      id_or_family)
      if (canon %in% L$layername) {
        layer <- canon
      } else {
        hits <- list_overlays(id_or_family)
        if (!nrow(hits)) stop("No layers match: ", id_or_family)
        if (which < 1 || which > nrow(hits)) stop("`which` out of range (", nrow(hits), " hits).")
        layer <- hits$layername[which]
      }
    } else {
      layer <- id_or_family
    }

    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    safe <- gsub("[^A-Za-z0-9_-]+", "_", layer)
    cache_path <- file.path(cache_dir, paste0("overlay_", safe, ".gpkg"))

    if (file.exists(cache_path)) {
      if (!quiet) message("Using cached overlay: ", cache_path)
      ov <- try(sf::read_sf(cache_path), silent = TRUE)
      if (!inherits(ov, "try-error")) return(ov)
    }

    if (requireNamespace("mregions2", quietly = TRUE) &&
        "mrp_get" %in% getNamespaceExports("mregions2")) {
      ov <- mregions2::mrp_get(layer)
    } else if (requireNamespace("mregions", quietly = TRUE) &&
               "mr_as_sf" %in% getNamespaceExports("mregions")) {
      ov <- mregions::mr_as_sf(layer)
    } else {
      stop("Need mregions2::mrp_get or mregions::mr_as_sf to fetch overlay.")
    }

    ov <- suppressWarnings(sf::st_make_valid(ov))
    # cache (best effort)
    try(sf::write_sf(ov, cache_path), silent = TRUE)
    ov
  }

  pick_ci <- function(x, candidates) {
    nx <- names(x)
    for (cand in candidates) {
      hit <- which(tolower(nx) == tolower(cand))
      if (length(hit)) return(nx[hit[1]])
    }
    NULL
  }

  # ---- fetch overlay (preloaded, or download+cache) -------------------------
  ov <- overlay_sf
  if (is.null(ov)) {
    layer_id <- overlay_layer %||% overlay
    ov <- get_overlay_sf(layer_id, which = 1, cache_dir = cache_dir, quiet = FALSE)
  }

  # ---- try to standardize a label and id column across overlays --------------
  # label candidates by family (first hit wins)
  name_col <- pick_ci(ov, c(
    # EEZ
    "geoname","GEONAME","GeoName","name",
    # LME
    "LME_NAME","lme_name","LME_NAME_LME66",
    # IHO
    "NAME","name","IHO_SEA",
    # MEOW
    "ECOREGION","ecoregion","ECO_NAME","eco_name",
    # FAO
    "F_AREA","F_LEVEL","OCEAN","title"
  ))

  # id (MRGID is only guaranteed on some layers, e.g., EEZ)
  mrgid_col <- pick_ci(ov, c("mrgid","MRGID","Mrgid","MRGID_EEZ","MRGID_POLY"))

  # Try to carry useful country-ish fields for EEZ ISO derivation
  terr1_col     <- pick_ci(ov, c("TERRITORY1","territory1","territory"))
  sov1_col      <- pick_ci(ov, c("SOVEREIGN1","sovereign1","sovereign"))
  iso_ter1_col  <- pick_ci(ov, c("ISO_Ter1","iso_ter1","ISO_TER1"))
  iso_sov1_col  <- pick_ci(ov, c("ISO_SOV1","iso_sov1","ISO_SOV"))
  iso2_ter1_col <- pick_ci(ov, c("ISO2_Ter1","iso2_ter1","ISO2_TER1"))
  iso2_sov1_col <- pick_ci(ov, c("ISO2_SOV1","iso2_sov1","ISO2_SOV"))

  keep_cols <- unique(na.omit(c(name_col, mrgid_col,
                                terr1_col,  sov1_col,
                                iso_ter1_col, iso_sov1_col,
                                iso2_ter1_col, iso2_sov1_col)))
  ov_keep <- ov[, keep_cols, drop = FALSE]

  # standardize names used downstream
  if (!is.null(name_col)  && name_col  %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==name_col]   <- "geoname"
  if (!is.null(mrgid_col) && mrgid_col %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==mrgid_col]  <- "mrgid"
  if (!is.null(terr1_col) && terr1_col %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==terr1_col]  <- "territory1"
  if (!is.null(sov1_col)  && sov1_col  %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==sov1_col]   <- "sovereign1"
  if (!is.null(iso_ter1_col)  && iso_ter1_col  %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==iso_ter1_col]  <- "iso_ter1"
  if (!is.null(iso_sov1_col)  && iso_sov1_col  %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==iso_sov1_col]  <- "iso_sov1"
  if (!is.null(iso2_ter1_col) && iso2_ter1_col %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==iso2_ter1_col] <- "iso2_ter1"
  if (!is.null(iso2_sov1_col) && iso2_sov1_col %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==iso2_sov1_col] <- "iso2_sov1"

  if (!"geoname" %in% names(ov_keep)) {
    # last resort: take any non-geometry column
    non_geom <- setdiff(names(ov_keep), attr(ov_keep, "sf_column"))
    if (length(non_geom)) ov_keep$geoname <- as.character(ov_keep[[non_geom[1]]]) else ov_keep$geoname <- NA_character_
  }

  # ---- compute iso3_best / iso2_best (robust for EEZ; NA for other overlays) -
  gcol <- function(df, nm) if (nm %in% names(df)) df[[nm]] else NULL
  coalesce_char <- function(...) {
    args <- list(...); n <- nrow(ov_keep); out <- rep(NA_character_, n)
    for (v in args) {
      if (is.null(v)) next
      v <- as.character(v); fill <- is.na(out) | !nzchar(out); out[fill] <- v[fill]
    }
    out
  }

  iso3_best <- iso2_best <- rep(NA_character_, nrow(ov_keep))
  if (tolower(overlay) == "eez" || (!is.null(overlay_layer) && grepl("^eez$", overlay_layer, ignore.case = TRUE))) {
    iso3_direct <- coalesce_char(gcol(ov_keep, "iso_ter1"), gcol(ov_keep, "iso_sov1"))
    if (requireNamespace("countrycode", quietly = TRUE)) {
      need <- is.na(iso3_direct) | !nzchar(iso3_direct)
      if (any(need)) {
        nm_guess <- coalesce_char(gcol(ov_keep, "territory1"), gcol(ov_keep, "sovereign1"), gcol(ov_keep, "geoname"))
        nm_guess <- gsub("Exclusive Economic Zone|Joint regime area|Overlapping claim area|\\(.*?\\)|-{1,2}\\s*EEZ|\\s+EEZ",
                         "", nm_guess, ignore.case = TRUE)
        nm_guess <- trimws(nm_guess)
        iso3_from_name <- suppressWarnings(countrycode::countrycode(nm_guess, "country.name", "iso3c"))
        iso3_direct[need] <- toupper(iso3_from_name[need])
      }
    }
    iso3_best <- toupper(iso3_direct)

    iso2_direct <- coalesce_char(gcol(ov_keep, "iso2_ter1"), gcol(ov_keep, "iso2_sov1"))
    if ((is.null(iso2_direct) || all(is.na(iso2_direct))) && requireNamespace("countrycode", quietly = TRUE)) {
      iso2_direct <- suppressWarnings(countrycode::countrycode(iso3_best, "iso3c", "iso2c"))
    }
    iso2_best <- toupper(iso2_direct)
  }

  ov_keep$iso3_best <- iso3_best
  ov_keep$iso2_best <- iso2_best

  # ---- s2 handling and robust fallback --------------------------------------
  old_s2 <- sf::sf_use_s2()
  on.exit(try(sf::sf_use_s2(old_s2), silent = TRUE), add = TRUE)
  if (isTRUE(use_planar)) sf::sf_use_s2(FALSE)

  out_list <- list()
  cols_to_join <- intersect(c("geoname","mrgid","territory1","sovereign1","iso3_best","iso2_best"), names(ov_keep))

  for (key in names(results_list)) {
    pts <- results_list[[key]]
    if (is.null(pts) || !inherits(pts, "sf") || !nrow(pts)) {
      if (inherits(pts, "sf")) {
        if (!"geoname"   %in% names(pts)) pts$geoname   <- NA_character_
        if (!"mrgid"     %in% names(pts)) pts$mrgid     <- NA_integer_
        if (!"iso3_best" %in% names(pts)) pts$iso3_best <- NA_character_
        if (!"iso2_best" %in% names(pts)) pts$iso2_best <- NA_character_
      }
      out_list[[key]] <- pts
      next
    }

    try_join <- function(force_planar = FALSE) {
      if (force_planar) sf::sf_use_s2(FALSE)
      crs_target <- if (!force_planar && !isTRUE(use_planar)) sf::st_crs(ov_keep) else 3857
      pts_tr <- try(suppressWarnings(sf::st_transform(pts, crs_target)), silent = TRUE)
      ply_tr <- try(suppressWarnings(sf::st_transform(ov_keep, crs_target)), silent = TRUE)
      if (inherits(pts_tr, "try-error") || inherits(ply_tr, "try-error")) stop("transform_failed")

      out <- try(suppressWarnings(
        sf::st_join(pts_tr, ply_tr[, cols_to_join, drop = FALSE], join = sf::st_within, left = TRUE)
      ), silent = TRUE)

      if (inherits(out, "try-error") || (!("geoname" %in% names(out)) || all(is.na(out$geoname)))) {
        # nearest as a last-ditch effort (kept because it salvages geometries that barely miss boundaries)
        out <- try(suppressWarnings(
          sf::st_join(pts_tr, ply_tr[, cols_to_join, drop = FALSE], join = sf::st_nearest_feature, left = TRUE)
        ), silent = TRUE)
      }
      if (inherits(out, "try-error")) stop(attr(out, "condition")$message)
      out
    }

    res <- try(try_join(FALSE), silent = TRUE)
    if (inherits(res, "try-error")) {
      res <- try(try_join(TRUE),  silent = TRUE)
      if (inherits(res, "try-error")) {
        # absolute fallback: return input plus empty columns
        if (!"geoname"   %in% names(pts)) pts$geoname   <- NA_character_
        if (!"mrgid"     %in% names(pts)) pts$mrgid     <- NA_integer_
        if (!"iso3_best" %in% names(pts)) pts$iso3_best <- NA_character_
        if (!"iso2_best" %in% names(pts)) pts$iso2_best <- NA_character_
        res <- pts
      }
    }
    out_list[[key]] <- res
  }

  out_list
}

#' Load & cache WWF Terrestrial Ecoregions (TEOW)
#'
#' Downloads (once) and caches the WWF terrestrial ecoregions as an sf object.
#' Prefer `method = "mapme"` if the \pkg{mapme.biodiversity} package is installed;
#' otherwise, use `method = "direct"` with a `source_url` to the official TEOW ZIP.
#'
#' The cached, normalized layer is written to `<cache_dir>/teow/teow.gpkg`
#' with the ecoregion name exposed as `geoname` and the geometry column
#' standardized to `"geometry"`.
#'
#' @export
bf_load_teow <- function(
    cache_dir = tools::R_user_dir("biofetchR", "data"),
    method = c("auto", "mapme", "direct"),
    source_url = NULL,
    aoi = NULL,
    bbox = NULL,
    force_refresh = FALSE,
    simplify = FALSE,
    simplify_tolerance = 0.01,
    quiet = FALSE
) {
  method    <- match.arg(method)
  teow_dir  <- file.path(cache_dir, "teow")
  teow_gpkg <- file.path(teow_dir, "teow.gpkg")
  if (!dir.exists(teow_dir)) dir.create(teow_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------- helpers ----------
  .global_bbox_sf <- function() {
    coords <- matrix(
      c(-179.999, -89.999,
        179.999, -89.999,
        179.999,  89.999,
        -179.999,  89.999,
        -179.999, -89.999),
      ncol = 2, byrow = TRUE
    )
    sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(coords)), crs = 4326))
  }
  .bbox_to_sf <- function(bx) {
    stopifnot(is.numeric(bx), length(bx) == 4)
    coords <- matrix(
      c(bx[1], bx[2],
        bx[3], bx[2],
        bx[3], bx[4],
        bx[1], bx[4],
        bx[1], bx[2]),
      ncol = 2, byrow = TRUE
    )
    sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(coords)), crs = 4326))
  }
  .as_sf <- function(g) {
    if (inherits(g, "sf"))  return(g)
    if (inherits(g, "sfc")) return(sf::st_sf(geometry = g))
    stop("AOI must be 'sf' or 'sfc'.")
  }
  .bf_make_valid <- function(x) {
    if ("st_make_valid" %in% getNamespaceExports("sf")) {
      return(suppressWarnings(sf::st_make_valid(x)))
    }
    suppressWarnings(sf::st_buffer(x, 0))
  }
  # --- NEW: works on older sf (no st_geometry_name needed)
  .standardize_geometry_name <- function(x) {
    stopifnot(inherits(x, "sf"))
    gcol <- attr(x, "sf_column", exact = TRUE)
    if (is.null(gcol) || !gcol %in% names(x)) {
      # Try common geometry column names
      cand <- intersect(c("geometry","geom","wkb_geometry","the_geom"), names(x))
      if (length(cand)) {
        gcol <- cand[1]
      } else {
        # As a last resort, find the first sfc column by class
        sfc_cols <- names(x)[vapply(x, function(col) inherits(col, "sfc"), logical(1))]
        if (length(sfc_cols)) {
          gcol <- sfc_cols[1]
        } else {
          stop("No geometry column found in TEOW layer.")
        }
      }
    }
    # Rename to "geometry" if needed and set sf_column attribute
    if (!identical(gcol, "geometry")) {
      names(x)[names(x) == gcol] <- "geometry"
    }
    attr(x, "sf_column") <- "geometry"
    x
  }
  .read_best_teow <- function(path, quiet = FALSE) {
    if (grepl("\\.gpkg$", path, ignore.case = TRUE)) {
      lay <- tryCatch(sf::st_layers(path), error = function(e) NULL)
      if (!is.null(lay)) {
        best <- NULL; best_n <- -Inf
        for (ln in lay$name) {
          x <- tryCatch(sf::read_sf(path, layer = ln, quiet = quiet), error = function(e) NULL)
          if (is.null(x)) next
          has_cols <- any(c("ECO_NAME","eco_name","ECO_ID","ECOREGION","ecoregion","ECO_CODE") %in% names(x))
          if (has_cols && nrow(x) > best_n) { best <- x; best_n <- nrow(x) }
        }
        if (!is.null(best)) return(best)
      }
    }
    sf::read_sf(path, quiet = quiet)
  }

  aoi_clip  <- aoi
  bbox_clip <- bbox

  # ---------- fast path: cached ----------
  if (!force_refresh && file.exists(teow_gpkg)) {
    out <- .read_best_teow(teow_gpkg, quiet = quiet)
    out <- .standardize_geometry_name(out)
    if (nrow(out) < 100 && !force_refresh) {
      if (!quiet) message("Cached TEOW has only ", nrow(out), " features; refetching via 'direct'.")
      force_refresh <- TRUE
    } else {
      return(.bf_teow_clip(out, aoi_clip, bbox_clip, quiet))
    }
  }

  # ---------- choose method ----------
  if (method == "auto") {
    method <- if (requireNamespace("mapme.biodiversity", quietly = TRUE)) "mapme" else "direct"
  }

  # ---------- fetch ----------
  teow <- NULL

  if (identical(method, "mapme")) {
    if (!requireNamespace("mapme.biodiversity", quietly = TRUE)) {
      stop("mapme.biodiversity not installed; use method='direct' with a source_url.")
    }
    aoi_dl <- if (is.null(aoi) && is.null(bbox)) .global_bbox_sf()
    else if (!is.null(bbox)) .bbox_to_sf(bbox)
    else sf::st_transform(.as_sf(aoi), 4326)

    mapme.biodiversity::mapme_options(outdir = teow_dir, verbose = FALSE)
    invisible(mapme.biodiversity::get_resources(aoi_dl, mapme.biodiversity::get_teow()))

    cand <- c(
      list.files(teow_dir, pattern = "wwf_terr_ecos\\.gpkg$", full.names = TRUE, recursive = TRUE),
      list.files(teow_dir, pattern = "\\.gpkg$",              full.names = TRUE, recursive = TRUE),
      list.files(teow_dir, pattern = "wwf_terr_ecos\\.shp$",  full.names = TRUE, recursive = TRUE),
      list.files(teow_dir, pattern = "\\.shp$",               full.names = TRUE, recursive = TRUE)
    )
    cand <- unique(cand)
    if (length(cand)) {
      info <- file.info(cand)
      path <- cand[order(info$mtime, decreasing = TRUE)][1]
      teow <- .read_best_teow(path, quiet = quiet)
    } else {
      if (is.null(source_url)) {
        source_url <- "https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip"
      }
      if (!quiet) message("Mapme did not materialize TEOW; falling back to direct.")
      method <- "direct"
    }
  }

  if (identical(method, "direct")) {
    if (is.null(source_url)) stop("For method='direct', provide source_url to the official TEOW ZIP.")
    zip_path  <- file.path(teow_dir, "teow.zip")
    unzip_dir <- file.path(teow_dir, "unzipped")
    utils::download.file(source_url, zip_path, mode = "wb", quiet = quiet)
    if (dir.exists(unzip_dir)) unlink(unzip_dir, recursive = TRUE, force = TRUE)
    dir.create(unzip_dir, showWarnings = FALSE)
    utils::unzip(zip_path, exdir = unzip_dir)
    shp <- list.files(unzip_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)
    if (!length(shp)) stop("No .shp found inside the TEOW ZIP.")
    teow <- sf::read_sf(shp[1], quiet = quiet)
  }

  if (!inherits(teow, "sf")) stop("TEOW layer could not be read as an sf object.")

  # Ensure WGS84
  if (!is.na(sf::st_crs(teow)) && sf::st_crs(teow)$epsg != 4326) {
    teow <- sf::st_transform(teow, 4326)
  }

  # Standardize geometry column name to "geometry" (no st_geometry_name used)
  teow <- .standardize_geometry_name(teow)

  # Sanity check: if suspiciously small and not already 'direct', refetch via 'direct'
  if (nrow(teow) < 100 && !identical(method, "direct")) {
    if (!quiet) message("TEOW layer has only ", nrow(teow), " features; refetching via 'direct'.")
    zip_path  <- file.path(teow_dir, "teow.zip")
    unzip_dir <- file.path(teow_dir, "unzipped")
    if (is.null(source_url)) {
      source_url <- "https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip"
    }
    utils::download.file(source_url, zip_path, mode = "wb", quiet = quiet)
    if (dir.exists(unzip_dir)) unlink(unzip_dir, recursive = TRUE, force = TRUE)
    dir.create(unzip_dir, showWarnings = FALSE)
    utils::unzip(zip_path, exdir = unzip_dir)
    shp <- list.files(unzip_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)
    if (!length(shp)) stop("No .shp found inside the TEOW ZIP (refetch).")
    teow <- sf::read_sf(shp[1], quiet = quiet)
    if (!is.na(sf::st_crs(teow)) && sf::st_crs(teow)$epsg != 4326) {
      teow <- sf::st_transform(teow, 4326)
    }
    teow <- .standardize_geometry_name(teow)
  }

  # Normalize fields
  name_candidates <- c("ECO_NAME","eco_name","ECOREGION","ecoregion","ECO_NAME_E","ECO_ID","ECO_CODE")
  nm <- intersect(name_candidates, names(teow))
  teow$geoname <- if (length(nm)) as.character(teow[[nm[1]]]) else paste0("ECO_", seq_len(nrow(teow)))

  # Keep useful attrs + geometry
  keep <- intersect(c("geoname","ECO_NAME","ECO_ID","ECO_CODE","BIOME_NAME","BIOME","REALM"), names(teow))
  teow <- teow[, unique(c(keep, "geometry")), drop = FALSE]

  # Validity + optional simplify
  teow <- .bf_make_valid(teow)
  if (isTRUE(simplify)) {
    teow <- sf::st_simplify(teow, dTolerance = simplify_tolerance, preserveTopology = TRUE)
  }

  # Optional clip (use original user AOI/bbox)
  teow <- .bf_teow_clip(teow, aoi_clip, bbox_clip, quiet)

  # Write normalized cache
  if (file.exists(teow_gpkg)) unlink(teow_gpkg)
  sf::write_sf(teow, teow_gpkg, quiet = quiet)
  if (!quiet) message("TEOW cached at: ", teow_gpkg)

  teow
}

# helper: clip to AOI or bbox
.bf_teow_clip <- function(x, aoi = NULL, bbox = NULL, quiet = TRUE) {
  stopifnot(inherits(x, "sf"))

  # nothing to do
  if (is.null(aoi) && is.null(bbox)) return(x)

  # --- helpers (local, no lwgeom) ---
  .as_sf <- function(g) if (inherits(g, "sf")) g else if (inherits(g, "sfc")) sf::st_sf(geometry = g) else stop("AOI must be sf/sfc.")
  .bbox_to_sf <- function(bx) {
    stopifnot(is.numeric(bx), length(bx) == 4)
    coords <- matrix(c(bx[1],bx[2],  bx[3],bx[2],  bx[3],bx[4],  bx[1],bx[4],  bx[1],bx[2]), ncol = 2, byrow = TRUE)
    sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(coords)), crs = 4326))
  }
  .bf_make_valid <- function(g) {
    if ("st_make_valid" %in% getNamespaceExports("sf")) return(suppressWarnings(sf::st_make_valid(g)))
    suppressWarnings(sf::st_buffer(g, 0))
  }

  # build clip geometry (in WGS84)
  clipper <- if (!is.null(bbox)) .bbox_to_sf(bbox) else sf::st_transform(.as_sf(aoi), 4326)

  # validity (sf only)
  x       <- .bf_make_valid(x)
  clipper <- .bf_make_valid(clipper)

  # fast prefilter by intersection (reduces nodes before precise clip)
  x_pref <- tryCatch(
    sf::st_filter(x, clipper, .predicate = sf::st_intersects, drop = FALSE),
    error = function(e) x
  )

  # precise clip; if it fails, return prefiltered
  out <- tryCatch(
    suppressWarnings(sf::st_intersection(x_pref, clipper)),
    error = function(e) x_pref
  )
  if (!quiet) message("TEOW clipped to AOI/BBOX (", nrow(out), " features).")
  out
}

#' Return cache path/info for TEOW
#' @inheritParams bf_load_teow
#' @return list(path=..., exists=TRUE/FALSE, size_bytes=..., modified=POSIXct)
#' @export
bf_teow_cache_info <- function(cache_dir = tools::R_user_dir("biofetchR", "data")) {
  teow_gpkg <- file.path(cache_dir, "teow", "teow.gpkg")
  info <- if (file.exists(teow_gpkg)) file.info(teow_gpkg) else NULL
  list(
    path = teow_gpkg,
    exists = file.exists(teow_gpkg),
    size_bytes = if (!is.null(info)) unname(info$size) else NA_real_,
    modified = if (!is.null(info)) info$mtime else NA
  )
}

#' Clear the cached TEOW dataset
#' @inheritParams bf_load_teow
#' @param ask If TRUE, prompt before deleting. Set FALSE for non-interactive use.
#' @export
bf_teow_clear_cache <- function(cache_dir = tools::R_user_dir("biofetchR", "data"), ask = interactive()) {
  teow_dir <- file.path(cache_dir, "teow")
  if (!dir.exists(teow_dir)) return(invisible(TRUE))
  if (ask) {
    ans <- utils::menu(c("No", "Yes"), title = sprintf("Delete TEOW cache at %s ?", teow_dir))
    if (ans != 2) return(invisible(FALSE))
  }
  unlink(teow_dir, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

#' Load & cache WWF Freshwater Ecoregions (FEOW)
#'
#' Downloads (once) and caches the WWF FEOW layer as an sf object.
#' No dependency on mapme.biodiversity (there is no get_feow()).
#'
#' The normalized layer is cached at `<cache_dir>/feow/feow.gpkg`
#' with the ecoregion name in `geoname`.
#'
#' @param cache_dir Directory to store cache. Default: user data dir
#'   (`tools::R_user_dir("biofetchR","data")`).
#' @param source_url A URL **or local path** to a FEOW vector dataset:
#'   - a .zip (containing .gpkg or .shp), or
#'   - a .gpkg or .shp file directly.
#'   If cache already exists and `force_refresh = FALSE`, this is optional.
#' @param aoi Optional sf; clip to this AOI (assumed/forced WGS84).
#' @param bbox Optional numeric c(xmin,ymin,xmax,ymax) in lon/lat; alternative to `aoi`.
#' @param force_refresh Logical; re-build cache even if present.
#' @param simplify Logical; simplify geometry before caching.
#' @param simplify_tolerance Numeric tolerance for simplification (degrees; WGS84).
#' @param quiet Logical; suppress messages.
#'
#' @return sf with at least `geoname` and geometry.
#' @export
bf_load_feow <- function(
    cache_dir = tools::R_user_dir("biofetchR","data"),
    source_url = NULL,
    aoi = NULL,
    bbox = NULL,
    force_refresh = FALSE,
    simplify = FALSE,
    simplify_tolerance = 0.01,
    quiet = FALSE
) {
  feow_dir  <- file.path(cache_dir, "feow")
  feow_gpkg <- file.path(feow_dir, "feow.gpkg")
  if (!dir.exists(feow_dir)) dir.create(feow_dir, recursive = TRUE, showWarnings = FALSE)

  # return cached if available
  if (!force_refresh && file.exists(feow_gpkg)) {
    out <- sf::read_sf(feow_gpkg, quiet = quiet)
    return(.feow_clip_internal(out, aoi, bbox, quiet))
  }

  # helpers ------------
  .is_http <- function(x) is.character(x) && grepl("^https?://", x, ignore.case = TRUE)
  .mkvalid <- function(x) {
    if ("st_make_valid" %in% getNamespaceExports("sf")) {
      return(suppressWarnings(sf::st_make_valid(x)))
    } else {
      return(suppressWarnings(sf::st_buffer(x, 0)))
    }
  }
  .std_geom_name <- function(x) {
    stopifnot(inherits(x, "sf"))
    # Try attribute first
    gcol <- tryCatch(attr(x, "sf_column"), error = function(e) NULL)
    if (!is.null(gcol) && gcol %in% names(x)) {
      if (!identical(gcol, "geometry")) {
        names(x)[names(x) == gcol] <- "geometry"
        sf::st_geometry(x) <- "geometry"
      }
      return(x)
    }
    # Fallbacks
    cand <- intersect(c("geometry","geom","wkb_geometry"), names(x))
    if (length(cand)) {
      if (!identical(cand[1], "geometry")) {
        names(x)[names(x) == cand[1]] <- "geometry"
      }
      sf::st_geometry(x) <- "geometry"
      return(x)
    }
    stop("Could not find geometry column in FEOW layer.")
  }
  .read_from_source <- function(src, dest_dir, q = quiet) {
    if (is.null(src) || is.na(src)) {
      stop("FEOW cache not present and no `source_url` (URL or local path) was provided.")
    }
    if (.is_http(src)) {
      zip_path  <- file.path(dest_dir, "feow.zip")
      if (!q) message("Downloading FEOW from URL…")
      utils::download.file(src, zip_path, mode = "wb", quiet = q)
      unzip_dir <- file.path(dest_dir, "unzipped")
      if (dir.exists(unzip_dir)) unlink(unzip_dir, recursive = TRUE, force = TRUE)
      dir.create(unzip_dir, showWarnings = FALSE)
      utils::unzip(zip_path, exdir = unzip_dir)
      cand <- c(
        list.files(unzip_dir, pattern = "\\.gpkg$", full.names = TRUE, recursive = TRUE),
        list.files(unzip_dir, pattern = "\\.shp$",  full.names = TRUE, recursive = TRUE)
      )
      if (!length(cand)) stop("No .gpkg or .shp found inside FEOW zip.")
      return(sf::read_sf(cand[1], quiet = q))
    } else {
      # Local path: can be a file or a directory
      if (dir.exists(src)) {
        cand <- c(
          list.files(src, pattern = "\\.gpkg$", full.names = TRUE, recursive = TRUE),
          list.files(src, pattern = "\\.shp$",  full.names = TRUE, recursive = TRUE)
        )
        if (!length(cand)) stop("No .gpkg or .shp found under local FEOW directory: ", src)
        return(sf::read_sf(cand[1], quiet = q))
      } else {
        stopifnot(file.exists(src))
        if (grepl("\\.zip$", src, ignore.case = TRUE)) {
          unzip_dir <- file.path(dest_dir, "unzipped_local")
          if (dir.exists(unzip_dir)) unlink(unzip_dir, recursive = TRUE, force = TRUE)
          dir.create(unzip_dir, showWarnings = FALSE)
          utils::unzip(src, exdir = unzip_dir)
          cand <- c(
            list.files(unzip_dir, pattern = "\\.gpkg$", full.names = TRUE, recursive = TRUE),
            list.files(unzip_dir, pattern = "\\.shp$",  full.names = TRUE, recursive = TRUE)
          )
          if (!length(cand)) stop("No .gpkg or .shp found inside local FEOW zip: ", src)
          return(sf::read_sf(cand[1], quiet = q))
        } else {
          # .gpkg or .shp directly
          return(sf::read_sf(src, quiet = q))
        }
      }
    }
  }
  .name_col <- function(x) {
    # common FEOW name fields
    nm <- intersect(c("FEOW_NAME","ECO_NAME","ECOREGION","ecoregion","NAME","name"), names(x))
    if (length(nm)) return(nm[1])
    return(NULL)
  }
  .feow_clip_internal <- function(x, aoi, bbox, q) {
    if (is.null(aoi) && is.null(bbox)) return(x)
    if (!is.null(bbox)) {
      stopifnot(is.numeric(bbox), length(bbox) == 4)
      aoi <- sf::st_as_sfc(sf::st_bbox(c(xmin=bbox[1], ymin=bbox[2], xmax=bbox[3], ymax=bbox[4]), crs = sf::st_crs(4326)))
    } else {
      # force AOI to WGS84
      aoi <- if (inherits(aoi, "sf")) aoi else if (inherits(aoi, "sfc")) sf::st_sf(geometry = aoi) else stop("AOI must be sf/sfc.")
      if (is.na(sf::st_crs(aoi))) sf::st_crs(aoi) <- 4326
      aoi <- sf::st_transform(aoi, 4326)
    }
    # Transform AOI if FEOW isn't WGS84
    if (!is.na(sf::st_crs(x)) && sf::st_crs(x) != sf::st_crs(4326)) {
      aoi_use <- tryCatch(sf::st_transform(aoi, sf::st_crs(x)), error = function(e) aoi)
    } else {
      aoi_use <- aoi
    }
    suppressWarnings(sf::st_crop(x, sf::st_bbox(aoi_use)))
  }

  # ------------- read -------------
  feow <- .read_from_source(source_url, feow_dir, q = quiet)
  if (!inherits(feow, "sf")) stop("Could not read FEOW as sf.")

  # Ensure/standardize CRS + geometry column
  if (is.na(sf::st_crs(feow))) sf::st_crs(feow) <- 4326
  if (!is.na(sf::st_crs(feow)) && sf::st_crs(feow)$epsg != 4326) feow <- sf::st_transform(feow, 4326)
  feow <- .std_geom_name(feow)

  # Normalize name
  nm <- .name_col(feow)
  feow$geoname <- if (!is.null(nm)) as.character(feow[[nm]]) else paste0("FEOW_", seq_len(nrow(feow)))

  # Keep a compact set of attributes
  keep <- intersect(c("geoname","FEOW_NAME","ECO_NAME","ECOREGION","REALM","BIOME"), names(feow))
  feow <- feow[, keep, drop = FALSE]

  # Validity + optional simplify
  feow <- .mkvalid(feow)
  if (isTRUE(simplify)) feow <- sf::st_simplify(feow, dTolerance = simplify_tolerance, preserveTopology = TRUE)

  # Optional clip
  feow <- .feow_clip_internal(feow, aoi, bbox, quiet)

  # Cache
  if (file.exists(feow_gpkg)) unlink(feow_gpkg)
  sf::write_sf(feow, feow_gpkg, quiet = quiet)
  if (!quiet) message("FEOW cached at: ", feow_gpkg)

  # quick sanity hint
  if (!quiet && nrow(feow) < 100) {
    message("⚠ FEOW has only ", nrow(feow), " features — this seems low. Check the source file.")
  }

  feow
}

#' Load Freshwater Ecoregions of the World (FEOW) polygons
#'
#' Downloads FEOW directly from The Nature Conservancy's ArcGIS FeatureServer
#' (no local files needed), standardizes key fields, and caches to a GeoPackage.
#'
#' @param force_refresh logical; re-download even if cache exists.
#' @param quiet logical; suppress messages.
#' @return An sf object with columns feow_id, ecoregion, biome, realm plus originals.
#' @export
bf_load_feow <- local({

  .FEOW_ARCGIS_LAYER <- "https://services.arcgis.com/uUvqNMGPm7axC2dD/arcgis/rest/services/FEOWv1_TNC/FeatureServer/0"

  .bf_read_and_standardise_feow <- function(gpkg_path, quiet = TRUE) {
    x  <- sf::st_read(gpkg_path, quiet = quiet)
    nm <- names(x)
    x$feow_id   <- if ("ECO_ID_U" %in% nm) x$ECO_ID_U else if ("ECO_ID" %in% nm) x$ECO_ID else if ("FID" %in% nm) x$FID else NA
    x$ecoregion <- if ("ECOREGION" %in% nm) x$ECOREGION else NA
    x$biome     <- if ("MHT_TXT" %in% nm) x$MHT_TXT else if ("BIOME" %in% nm) x$BIOME else NA
    x$realm     <- if ("REALM" %in% nm) x$REALM else NA
    x <- sf::st_make_valid(x)
    sf::st_set_agr(x, "constant")
    x
  }

  .bf_download_file <- function(url, destfile, quiet = TRUE) {
    mode <- "wb"
    if (requireNamespace("curl", quietly = TRUE)) {
      h <- curl::new_handle(followlocation = 1L, timeout = 180L)
      curl::curl_download(url, destfile, handle = h, mode = mode, quiet = quiet)
    } else {
      utils::download.file(url, destfile, mode = mode, quiet = quiet)
    }
    invisible(destfile)
  }

  .bf_arcgis_get_page <- function(layer_url, offset = 0L, n = 2000L, quiet = TRUE) {
    # Build a GeoJSON query URL (WGS84), with paging
    q <- paste0(
      layer_url,
      "/query?",
      "where=1%3D1",
      "&outFields=*",
      "&returnGeometry=true",
      "&outSR=4326",
      "&f=geojson",
      "&resultRecordCount=", n,
      "&resultOffset=", offset,
      "&orderByFields=FID"
    )
    tf <- tempfile(fileext = ".geojson")
    .bf_download_file(q, tf, quiet = quiet)
    # st_read returns 0-row sf if page empty
    g <- try(suppressWarnings(sf::st_read(tf, quiet = quiet)), silent = TRUE)
    if (inherits(g, "try-error")) return(NULL)
    if (!inherits(g, "sf")) return(NULL)
    if (nrow(g) == 0L) return(NULL)
    g
  }

  .bf_arcgis_fetch_all <- function(layer_url, quiet = TRUE) {
    pgsize <- 2000L
    out <- list()
    off <- 0L
    i <- 1L
    repeat {
      g <- .bf_arcgis_get_page(layer_url, offset = off, n = pgsize, quiet = quiet)
      if (is.null(g)) break
      out[[i]] <- g
      if (nrow(g) < pgsize) break
      off <- off + pgsize
      i <- i + 1L
    }
    if (!length(out)) return(NULL)
    # Base rbind works for sf with consistent crs/geometry
    x <- do.call(rbind, out)
    # Ensure CRS is WGS84
    if (is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326
    x
  }

  function(force_refresh = FALSE, quiet = TRUE) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required for bf_load_feow(). Please install it.")
    }

    cache_root <- tools::R_user_dir("biofetchR", which = "cache")
    feow_dir   <- file.path(cache_root, "feow")
    if (!dir.exists(feow_dir)) dir.create(feow_dir, recursive = TRUE, showWarnings = FALSE)
    feow_gpkg <- file.path(feow_dir, "feow.gpkg")

    # Use cache if present
    if (file.exists(feow_gpkg) && !force_refresh) {
      if (!quiet) message("biofetchR: FEOW cache found: ", feow_gpkg)
      return(.bf_read_and_standardise_feow(feow_gpkg, quiet = quiet))
    }

    if (!quiet) message("biofetchR: Fetching FEOW from ArcGIS FeatureServer...")

    # Primary: pull directly from ArcGIS FeatureServer (TNC)
    feow <- .bf_arcgis_fetch_all(.FEOW_ARCGIS_LAYER, quiet = quiet)

    # Fallback: {feowR} if the service is down
    if (is.null(feow)) {
      if (!quiet) message("biofetchR: ArcGIS fetch failed; attempting feowR::read_feow()")
      if (!requireNamespace("feowR", quietly = TRUE)) {
        stop("Could not fetch FEOW from ArcGIS, and package 'feowR' is not installed for fallback.")
      }
      feow <- feowR::read_feow()
      # Ensure WGS84
      if (is.na(sf::st_crs(feow))) sf::st_crs(feow) <- 4326
    }

    # Cache to gpkg and return standardized
    sf::st_write(feow, feow_gpkg, delete_dsn = TRUE, quiet = TRUE)
    .bf_read_and_standardise_feow(feow_gpkg, quiet = quiet)
  }
})


#' Load Freshwater Ecoregions of the World (FEOW) polygons
#'
#' Directly fetches FEOW polygons from a live source (no user files needed),
#' standardizes fields, and caches to a single GeoPackage. Backward-compatible
#' with older calls that pass \code{cache_dir}, \code{method}, and \code{source_url}.
#'
#' @param force_refresh logical; re-download even if cache exists.
#' @param quiet logical; suppress messages.
#' @param cache_dir optional character; override cache root directory used for the gpkg.
#' @param method character; one of "auto","arcgis","download","feowR" (case-insensitive).
#'   "auto" tries arcgis -> download -> feowR.
#' @param source_url optional character; explicit download URL (zip or vector). Used when
#'   \code{method="download"}. If NULL, the FEOW downloads page is scraped.
#' @param ... ignored (forward/backward compatibility).
#'
#' @return An sf object with standardized columns: feow_id, ecoregion, biome, realm.
#' @importFrom tools R_user_dir
#' @importFrom utils unzip readLines download.file
#' @export
bf_load_feow <- function(force_refresh = FALSE,
                         quiet = TRUE,
                         cache_dir = NULL,
                         method = "auto",
                         source_url = NULL,
                         ...) {
  FEOW_ARCGIS_LAYER <- "https://services.arcgis.com/uUvqNMGPm7axC2dD/arcgis/rest/services/FEOWv1_TNC/FeatureServer/0"

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for bf_load_feow(). Please install it.")
  }

  method <- tolower(method)
  if (!method %in% c("auto","arcgis","download","feowr")) {
    stop("bf_load_feow(): unknown method: ", method)
  }

  cache_root <- if (!is.null(cache_dir)) cache_dir else tools::R_user_dir("biofetchR", which = "cache")
  feow_dir   <- file.path(cache_root, "feow")
  if (!dir.exists(feow_dir)) dir.create(feow_dir, recursive = TRUE, showWarnings = FALSE)
  feow_gpkg  <- file.path(feow_dir, "feow.gpkg")

  # Use cache if present
  if (file.exists(feow_gpkg) && !isTRUE(force_refresh)) {
    if (!quiet) message("biofetchR: FEOW cache found: ", feow_gpkg)
    return(.bf_read_and_standardise_feow(feow_gpkg, quiet = quiet))
  }

  # Fetchers
  fetch_arcgis <- function() {
    if (!quiet) message("biofetchR: Fetching FEOW from ArcGIS FeatureServer...")
    x <- .bf_arcgis_fetch_all(FEOW_ARCGIS_LAYER, quiet = quiet)
    if (!is.null(x) && is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326
    x
  }

  fetch_download <- function() {
    url <- source_url
    if (is.null(url)) {
      if (!quiet) message("biofetchR: Discovering FEOW download URL...")
      url <- .bf_feow_find_zip_url("https://www.feow.org/download")
    }
    if (is.null(url)) return(NULL)
    if (!quiet) message("biofetchR: Downloading FEOW archive...")

    tf <- tempfile(fileext = if (grepl("\\.zip($|\\?)", url, ignore.case = TRUE)) ".zip" else "")
    .bf_download_file(url, tf, quiet = quiet)

    root <- if (grepl("\\.zip($|\\?)", tolower(tf))) {
      td <- tempfile("feow_unzip_"); dir.create(td, recursive = TRUE, showWarnings = FALSE)
      utils::unzip(tf, exdir = td)
      .bf_unzip_nested(td)
      td
    } else {
      dirname(tf)
    }

    .bf_read_any_vector(root, quiet = quiet)
  }

  fetch_feowr <- function() {
    if (!requireNamespace("feowR", quietly = TRUE)) return(NULL)
    if (!quiet) message("biofetchR: Fetching FEOW via feowR::read_feow()...")
    x <- feowR::read_feow()
    if (!is.null(x) && is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326
    x
  }

  # Orchestrate per method
  feow <- switch(method,
                 arcgis   = fetch_arcgis(),
                 download = fetch_download(),
                 feowr    = fetch_feowr(),
                 auto     = { x <- fetch_arcgis(); if (is.null(x)) x <- fetch_download(); if (is.null(x)) x <- fetch_feowr(); x }
  )

  if (is.null(feow)) {
    stop("bf_load_feow(): failed to obtain FEOW via method='", method,
         "'. Tried ArcGIS, downloadable archive, and feowR (if installed).")
  }

  suppressWarnings(sf::st_write(feow, feow_gpkg, delete_dsn = TRUE, quiet = TRUE))
  .bf_read_and_standardise_feow(feow_gpkg, quiet = quiet)
}

# ---------- helpers (internal) ----------

.bf_read_and_standardise_feow <- function(gpkg_path, quiet = TRUE) {
  x  <- sf::st_read(gpkg_path, quiet = quiet)
  nm <- names(x)
  x$feow_id   <- if ("ECO_ID_U" %in% nm) x$ECO_ID_U else if ("ECO_ID" %in% nm) x$ECO_ID else if ("FEOW_ID" %in% nm) x$FEOW_ID else if ("ECO_CODE" %in% nm) x$ECO_CODE else NA
  x$ecoregion <- if ("ECOREGION" %in% nm) x$ECOREGION else if ("ECO_NAME" %in% nm) x$ECO_NAME else NA
  x$biome     <- if ("MHT_TXT"  %in% nm) x$MHT_TXT  else if ("BIOME" %in% nm) x$BIOME else NA
  x$realm     <- if ("REALM"    %in% nm) x$REALM    else NA
  x <- sf::st_make_valid(x)
  sf::st_set_agr(x, "constant")
  x
}

.bf_download_file <- function(url, destfile, quiet = TRUE) {
  mode <- "wb"
  if (requireNamespace("curl", quietly = TRUE)) {
    h <- curl::new_handle(followlocation = 1L, timeout = 180L)
    curl::curl_download(url, destfile, handle = h, mode = mode, quiet = quiet)
  } else {
    utils::download.file(url, destfile, mode = mode, quiet = quiet)
  }
  invisible(destfile)
}

.bf_feow_find_zip_url <- function(feow_downloads_page) {
  raw <- tryCatch(
    suppressWarnings(utils::readLines(feow_downloads_page, warn = FALSE)),
    error = function(e) NULL
  )
  if (is.null(raw)) return(NULL)
  m     <- regmatches(raw, gregexpr('href\\s*=\\s*"([^"]+)"', raw, perl = TRUE))
  hrefs <- unique(unlist(m))
  if (!length(hrefs)) return(NULL)
  urls <- sub('^href\\s*=\\s*"', "", hrefs)
  urls <- sub('"$', "", urls)
  urls <- urls[grepl("\\.zip(\\?.*)?$", urls, ignore.case = TRUE)]
  if (!length(urls)) return(NULL)
  urls <- ifelse(
    grepl("^https?://", urls, ignore.case = TRUE),
    urls,
    paste0("https://www.feow.org/", sub("^/", "", urls))
  )
  urls[[1]]
}

.bf_unzip_nested <- function(root_dir) {
  inner_zips <- list.files(root_dir, pattern = "\\.zip$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (!length(inner_zips)) return(invisible())
  for (z in inner_zips) {
    target <- file.path(dirname(z), tools::file_path_sans_ext(basename(z)))
    dir.create(target, showWarnings = FALSE, recursive = TRUE)
    utils::unzip(z, exdir = target)
  }
  invisible()
}

.bf_read_any_vector <- function(dirpath, quiet = TRUE) {
  gpkg <- list.files(dirpath, pattern = "\\.gpkg$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(gpkg)) return(sf::st_read(gpkg[1], quiet = quiet))

  gdbs <- list.dirs(dirpath, full.names = TRUE, recursive = TRUE)
  gdbs <- gdbs[grepl("\\.gdb$", gdbs, ignore.case = TRUE)]
  if (length(gdbs)) {
    lay <- try(sf::st_layers(gdbs[1]), silent = TRUE)
    if (!inherits(lay, "try-error")) {
      nms <- lay$name
      gtypes <- NULL
      if (!is.null(lay$geomtype)) gtypes <- lay$geomtype else if (!is.null(lay$geometry_type)) gtypes <- lay$geometry_type
      score <- tolower(paste(nms, gtypes))
      idx <- grep("poly", score)
      lyr <- if (length(idx)) nms[idx[1]] else nms[1]
      return(sf::st_read(gdbs[1], layer = lyr, quiet = quiet))
    }
  }

  shp <- list.files(dirpath, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(shp)) return(sf::st_read(shp[1], quiet = quiet))

  jsn <- list.files(dirpath, pattern = "\\.(geo)?json$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(jsn)) return(sf::st_read(jsn[1], quiet = quiet))

  kml <- list.files(dirpath, pattern = "\\.kml$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(kml)) return(sf::st_read(kml[1], quiet = quiet))

  kmz <- list.files(dirpath, pattern = "\\.kmz$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(kmz)) {
    td <- tempfile("kmz_"); dir.create(td, showWarnings = FALSE)
    utils::unzip(kmz[1], exdir = td)
    kml2 <- list.files(td, pattern = "\\.kml$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    if (length(kml2)) return(sf::st_read(kml2[1], quiet = quiet))
  }

  NULL
}

.bf_arcgis_get_page <- function(layer_url, offset = 0L, n = 2000L, quiet = TRUE) {
  q <- paste0(
    layer_url,
    "/query?",
    "where=1%3D1",
    "&outFields=*",
    "&returnGeometry=true",
    "&outSR=4326",
    "&f=geojson",
    "&resultRecordCount=", n,
    "&resultOffset=", offset,
    "&orderByFields=FID"
  )
  tf <- tempfile(fileext = ".geojson")
  .bf_download_file(q, tf, quiet = quiet)
  g <- try(suppressWarnings(sf::st_read(tf, quiet = quiet)), silent = TRUE)
  if (inherits(g, "try-error") || !inherits(g, "sf") || nrow(g) == 0L) return(NULL)
  g
}

.bf_arcgis_fetch_all <- function(layer_url, quiet = TRUE) {
  pgsize <- 2000L
  out <- list()
  off <- 0L
  i <- 1L
  repeat {
    g <- .bf_arcgis_get_page(layer_url, offset = off, n = pgsize, quiet = quiet)
    if (is.null(g)) break
    out[[i]] <- g
    if (nrow(g) < pgsize) break
    off <- off + pgsize
    i <- i + 1L
  }
  if (!length(out)) return(NULL)
  x <- do.call(rbind, out)
  if (is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326
  x
}

#' Load global lake polygons (HydroLAKES)
#'
#' Downloads HydroLAKES once (RAW), standardizes key fields, then applies
#' optional filters in-memory on each call. Cached as a single GeoPackage.
#'
#' @param force_refresh logical; re-download even if cache exists.
#' @param quiet logical; suppress messages.
#' @param cache_dir optional character; override cache root.
#' @param lakes_only logical; keep only natural lakes (exclude reservoirs/controls).
#' @param min_area_km2 optional numeric; drop lakes smaller than this area.
#' @return sf with standardized columns: lake_id, name, is_reservoir, area_km2.
#' @importFrom tools R_user_dir
#' @importFrom utils unzip download.file
#' @export
bf_load_lakes <- function(force_refresh = FALSE,
                          quiet = TRUE,
                          cache_dir = NULL,
                          lakes_only = TRUE,
                          min_area_km2 = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for bf_load_lakes(). Please install it.")
  }

  HYDROL_SHP_ZIP <- "https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip"

  cache_root <- if (!is.null(cache_dir)) cache_dir else tools::R_user_dir("biofetchR", which = "cache")
  lakes_dir  <- file.path(cache_root, "lakes")
  if (!dir.exists(lakes_dir)) dir.create(lakes_dir, recursive = TRUE, showWarnings = FALSE)

  gpkg_raw <- file.path(lakes_dir, "hydrolakes_raw.gpkg")

  # --- ensure RAW cache ----------------------------------------------------
  if (!file.exists(gpkg_raw) || isTRUE(force_refresh)) {
    if (!quiet) message("biofetchR: Downloading HydroLAKES (first run; large file)…")
    tmp_zip <- tempfile(fileext = ".zip")
    .bf_download_file(HYDROL_SHP_ZIP, tmp_zip, quiet = quiet)

    tmp_dir <- tempfile("hydrolakes_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(tmp_zip, exdir = tmp_dir)

    shp <- list.files(tmp_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    if (!length(shp)) stop("HydroLAKES archive did not contain a shapefile.")
    x <- sf::st_read(shp[1], quiet = quiet)

    # Standardize key fields (per HydroLAKES schema)
    nm <- names(x)
    x$lake_id      <- if ("Hylak_id"  %in% nm) x$Hylak_id  else NA_integer_
    x$name         <- if ("Lake_name" %in% nm) x$Lake_name else NA_character_
    x$is_reservoir <- if ("Lake_type" %in% nm) x$Lake_type != 1L else NA
    x$area_km2     <- if ("Lake_area" %in% nm) x$Lake_area else NA_real_

    x <- sf::st_make_valid(x)
    if (is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326

    suppressWarnings(sf::st_write(x, gpkg_raw, delete_dsn = TRUE, quiet = TRUE))
    if (!quiet) message("biofetchR: HydroLAKES cached at: ", gpkg_raw)
  }

  # --- read RAW cache, then filter in-memory --------------------------------
  x <- sf::st_read(gpkg_raw, quiet = quiet)

  # apply filters without altering the RAW cache
  nm <- names(x)
  if (isTRUE(lakes_only) && "Lake_type" %in% nm) {
    x <- x[is.na(x$Lake_type) | x$Lake_type == 1L, , drop = FALSE]
  }
  if (!is.null(min_area_km2) && "Lake_area" %in% nm) {
    x <- x[is.na(x$Lake_area) | x$Lake_area >= min_area_km2, , drop = FALSE]
  }

  # return standardized + original attrs
  std_cols <- c("lake_id","name","is_reservoir","area_km2")
  std_cols <- std_cols[std_cols %in% names(x)]
  x <- x[, unique(c(std_cols, setdiff(names(x), std_cols))), drop = FALSE]
  sf::st_set_agr(x, "constant")
  x
}

# internal helper (shared style with your other loaders)
.bf_download_file <- function(url, destfile, quiet = TRUE) {
  mode <- "wb"
  if (requireNamespace("curl", quietly = TRUE)) {
    h <- curl::new_handle(followlocation = 1L, timeout = 300L)
    curl::curl_download(url, destfile, handle = h, mode = mode, quiet = quiet)
  } else {
    utils::download.file(url, destfile, mode = mode, quiet = quiet)
  }
  invisible(destfile)
}

#' Load global river reaches (HydroRIVERS) with optional caching
#'
#' Downloads + standardizes HydroRIVERS. If `cache="disk"`, writes an RDS cache
#' (and optionally a GPKG). If `cache="memory"`, returns in-memory only.
#' Filters (`min_strahler`, `min_discharge_cms`) are applied on read, not baked into cache.
#'
#' @param force_refresh logical; re-download even if cache exists (disk mode).
#' @param quiet logical; suppress messages.
#' @param cache_dir optional character; cache root (ignored if cache="memory").
#' @param cache "disk" or "memory". Disk persists a cache; memory does not write files.
#' @param cache_format "rds","gpkg","both" (used only when cache="disk"; default "rds").
#' @param regions character in {"global","af","ar","as","au","eu","gr","na","sa","si"}; NULL -> "global".
#' @param min_strahler optional integer; keep ORD_STRA >= this (applied on read).
#' @param min_discharge_cms optional numeric; keep DIS_AV_CMS >= this (applied on read).
#' @return sf LINESTRING with canonical fields:
#'   hyriv_id, ord_stra, ord_clas, ord_flow, dis_av_cms, length_km, hybas_l12
#' @export
bf_load_rivers <- function(force_refresh = FALSE,
                           quiet = TRUE,
                           cache_dir = NULL,
                           cache = c("disk","memory"),
                           cache_format = c("rds","gpkg","both"),
                           regions = NULL,
                           min_strahler = NULL,
                           min_discharge_cms = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for bf_load_rivers(). Please install it.")
  }
  cache <- match.arg(cache)
  cache_format <- match.arg(cache_format)

  base <- "https://data.hydrosheds.org/file/hydrorivers"
  pick_urls <- function(regs) {
    if (is.null(regs) || identical(regs, "global"))
      return(file.path(base, "HydroRIVERS_v10_shp.zip"))
    regs <- match.arg(regs, c("af","ar","as","au","eu","gr","na","sa","si"), several.ok = TRUE)
    file.path(base, paste0("HydroRIVERS_v10_", regs, "_shp.zip"))
  }

  has_gpkg <- function() {
    drv <- try(sf::st_drivers(), silent = TRUE)
    !inherits(drv, "try-error") && any(drv$name == "GPKG" & drv$write)
  }

  # build one big sf in memory
  build_in_memory <- function(regs) {
    urls <- pick_urls(regs)
    if (!quiet) message("biofetchR: Downloading HydroRIVERS (", length(urls), " file(s)) …")
    layers <- vector("list", length(urls))
    for (i in seq_along(urls)) {
      u <- urls[i]
      tmp_zip <- tempfile(fileext = ".zip")
      .bf_download_file(u, tmp_zip, quiet = quiet)
      tmp_dir <- tempfile("hydrorivers_"); dir.create(tmp_dir, TRUE, FALSE)
      utils::unzip(tmp_zip, exdir = tmp_dir)
      shp <- list.files(tmp_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
      if (!length(shp)) stop("HydroRIVERS archive did not contain a shapefile: ", u)
      layers[[i]] <- sf::st_read(shp[1], quiet = quiet)
    }
    x <- do.call(rbind, layers)

    # ---- FIXED: pure/local renamer (no <<-) ----
    ren <- function(obj, old, new) {
      if (old %in% names(obj)) {
        nm <- names(obj); nm[nm == old] <- new; names(obj) <- nm
      }
      obj
    }
    x <- ren(x, "HYRIV_ID","hyriv_id")
    x <- ren(x, "ORD_STRA","ord_stra")
    x <- ren(x, "ORD_CLAS","ord_clas")
    x <- ren(x, "ORD_FLOW","ord_flow")
    x <- ren(x, "DIS_AV_CMS","dis_av_cms")
    x <- ren(x, "LENGTH_KM","length_km")
    x <- ren(x, "HYBAS_L12","hybas_l12")

    for (nm in c("hyriv_id","ord_stra","ord_clas","ord_flow","dis_av_cms","length_km","hybas_l12")) {
      if (!nm %in% names(x)) x[[nm]] <- NA
    }
    x <- x[, !duplicated(tolower(names(x))), drop = FALSE]

    x <- suppressWarnings(sf::st_make_valid(x))
    if (is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326
    sf::st_set_agr(x, "constant")
    x
  }

  if (cache == "memory") {
    x <- build_in_memory(regions)
    if (!is.null(min_strahler) && "ord_stra" %in% names(x))
      x <- x[is.na(x$ord_stra) | x$ord_stra >= min_strahler, , drop = FALSE]
    if (!is.null(min_discharge_cms) && "dis_av_cms" %in% names(x))
      x <- x[is.na(x$dis_av_cms) | x$dis_av_cms >= min_discharge_cms, , drop = FALSE]
    return(x)
  }

  # disk caching
  cache_root <- if (!is.null(cache_dir)) cache_dir else tools::R_user_dir("biofetchR", which = "cache")
  rivers_dir <- file.path(cache_root, "rivers"); if (!dir.exists(rivers_dir)) dir.create(rivers_dir, TRUE, FALSE)
  rds_raw    <- file.path(rivers_dir, "hydrorivers_raw.rds")
  gpkg_raw   <- file.path(rivers_dir, "hydrorivers_raw.gpkg")

  write_cache <- function(x) {
    if (cache_format %in% c("rds","both")) {
      tmp_rds <- tempfile(fileext = ".rds"); saveRDS(x, tmp_rds)
      if (file.exists(rds_raw)) unlink(rds_raw)
      file.rename(tmp_rds, rds_raw)
    }
    if (cache_format %in% c("gpkg","both") && has_gpkg()) {
      tmp_gpkg <- tempfile(fileext = ".gpkg")
      ok <- try({
        suppressWarnings(sf::st_write(x, tmp_gpkg, driver = "GPKG", delete_dsn = TRUE, quiet = TRUE))
        TRUE
      }, silent = TRUE)
      if (isTRUE(ok)) {
        if (file.exists(gpkg_raw)) unlink(gpkg_raw)
        file.rename(tmp_gpkg, gpkg_raw)
      } else {
        if (!quiet) message("biofetchR: GPKG write failed; continuing with RDS cache.")
        if (file.exists(tmp_gpkg)) unlink(tmp_gpkg)
      }
    }
  }

  read_cached <- function() {
    if (file.exists(rds_raw)) {
      x <- try(readRDS(rds_raw), silent = TRUE)
      if (!inherits(x, "try-error") && inherits(x, "sf")) return(x)
    }
    if (file.exists(gpkg_raw)) {
      x <- try(sf::st_read(gpkg_raw, quiet = quiet), silent = TRUE)
      if (!inherits(x, "try-error") && inherits(x, "sf")) return(x)
    }
    NULL
  }

  corrupt_file <- function(p) file.exists(p) && isTRUE(tryCatch(file.info(p)$size < 1024, error = function(e) FALSE))
  need_build <- isTRUE(force_refresh) ||
    (is.null(read_cached()) && !(file.exists(rds_raw) || file.exists(gpkg_raw))) ||
    corrupt_file(gpkg_raw)

  if (need_build) {
    x <- build_in_memory(regions)
    write_cache(x)
  } else {
    x <- read_cached()
    if (is.null(x)) {
      if (!quiet) message("biofetchR: Cache unreadable; rebuilding …")
      x <- build_in_memory(regions); write_cache(x)
    }
  }

  if (!is.null(min_strahler) && "ord_stra" %in% names(x))
    x <- x[is.na(x$ord_stra) | x$ord_stra >= min_strahler, , drop = FALSE]
  if (!is.null(min_discharge_cms) && "dis_av_cms" %in% names(x))
    x <- x[is.na(x$dis_av_cms) | x$dis_av_cms >= min_discharge_cms, , drop = FALSE]

  x
}

#' Annotate HydroRIVERS reaches with OSM river names
#' @param rivers sf LINESTRING with field `hyriv_id` (from bf_load_rivers()).
#' @param iso2c character ISO2 country codes used to build a bbox query.
#' @param max_snap_m numeric, max distance (m) to accept a name match.
#' @param quiet logical.
#' @return same sf with new column `river_name` (character).
#' @export
bf_name_rivers_osm <- function(rivers, iso2c, max_snap_m = 500, quiet = TRUE) {
  stopifnot(inherits(rivers, "sf"), "hyriv_id" %in% names(rivers))
  if (!requireNamespace("osmdata", quietly = TRUE) ||
      !requireNamespace("rnaturalearth", quietly = TRUE)) {
    stop("Packages 'osmdata' and 'rnaturalearth' are required.")
  }

  # 1) Country bbox (WGS84)
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  sel   <- world[world$iso_a2 %in% toupper(iso2c), , drop = FALSE]
  if (!nrow(sel)) stop("No countries found for: ", paste(iso2c, collapse=", "))
  bb    <- sf::st_bbox(sf::st_transform(sel, 4326))

  # 2) OSM query for named waterways
  q <- osmdata::opq(bbox = bb) |>
    osmdata::add_osm_feature(key = "waterway",
                             value = c("river","stream","canal","drain")) |>
    osmdata::add_osm_feature(key = "name")
  osm <- osmdata::osmdata_sf(q)
  lines <- osm$osm_lines
  if (is.null(lines) || !nrow(lines)) {
    if (!quiet) message("No named OSM waterways returned for bbox.")
    rivers$river_name <- NA_character__
    return(rivers)
  }
  lines <- lines[!is.na(lines$name) & nzchar(lines$name), c("name","geometry")]
  lines <- sf::st_transform(lines, 4326)

  # 3) Snap by nearest feature with distance screen (meters)
  riv4326 <- sf::st_transform(rivers, 4326)
  idx <- sf::st_nearest_feature(riv4326, lines)

  # distance check in a metric CRS
  riv_m  <- tryCatch(sf::st_transform(riv4326, 3857), error = function(e) riv4326)
  line_m <- tryCatch(sf::st_transform(lines[idx, ], 3857), error = function(e) lines[idx, ])
  d <- suppressWarnings(sf::st_distance(sf::st_geometry(riv_m),
                                        sf::st_geometry(line_m),
                                        by_element = TRUE))
  nm <- lines$name[idx]
  nm[as.numeric(d) > max_snap_m] <- NA_character_

  riv4326$river_name <- nm
  riv4326
}


#' Load nested basins (HydroBASINS; Pfafstetter)
#'
#' Downloads selected regional tiles at a chosen Pfafstetter level (L1–L12),
#' caches once, and returns polygons with canonical fields:
#' hybas_id, next_down, next_sink, main_bas, sub_area, up_area, endo, order.
#'
#' @param level integer 1..12 (Pfafstetter level).
#' @param with_lakes logical; use the customized “with lakes” variant.
#' @param regions character in {"af","ar","as","au","eu","gr","na","sa","si"} (default all).
#' @param force_refresh,quiet,cache_dir see bf_load_rivers.
#' @return sf POLYGON
#' @export
bf_load_basins <- function(level = 12, with_lakes = FALSE,
                           cache_dir = tools::R_user_dir("biofetchR","data"),
                           quiet = TRUE) {
  stopifnot(requireNamespace("sf", quietly = TRUE))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  hb_dir <- file.path(cache_dir, "hydrobasins")
  dir.create(hb_dir, recursive = TRUE, showWarnings = FALSE)

  regions <- getOption("biofetchR.basins.regions",
                       c("eu","na","sa","af","as","au","ar","an"))

  url_roots <- c(
    "https://data.hydrosheds.org/file/hydrobasins",
    "https://data.hydrosheds.org/file/hydrobasins/standard"
  )

  # candidate filenames (some mirrors use gdb.zip)
  fname_patterns <- c(
    sprintf("hybas_%%s_lev%02d_v1c.zip", level),
    sprintf("hybas_%%s_lev%02d_v1c.gdb.zip", level),
    sprintf("hybas_%%s_lev%02d_v1c_gdb.zip", level)
  )

  # downloader (quietly tries several URL variants)
  fetch_zip <- function(region) {
    dest <- file.path(hb_dir, sprintf("hybas_%s_lev%02d_v1c.zip", region, level))
    if (file.exists(dest)) return(dest)

    for (root in url_roots) {
      for (pat in fname_patterns) {
        fname <- sprintf(pat, region)
        url <- file.path(root, fname)
        tmp <- tempfile(fileext = ".zip")
        ok <- tryCatch({
          curl::curl_download(url, tmp, quiet = quiet)
          TRUE
        }, error = function(e) FALSE)
        if (ok) {
          file.rename(tmp, dest)
          if (!quiet) message("HydroBASINS: downloaded ", url)
          return(dest)
        }
      }
    }
    if (!quiet) message("HydroBASINS: FAILED to fetch region '", region,
                        "' for level ", level, " (tried ", paste(url_roots, collapse=" | "), ")")
    return(NA_character_)
  }

  zips <- vapply(regions, fetch_zip, character(1))
  zips <- zips[!is.na(zips) & file.exists(zips)]
  if (!length(zips)) stop("HydroBASINS download failed for all regions.")

  # unzip + read the target level shapefile (or gdb)
  read_one <- function(zip) {
    exdir <- tempfile("hybas_")
    utils::unzip(zip, exdir = exdir)
    # prefer levXX shp
    shp <- list.files(exdir, pattern = sprintf("lev%02d_v1c\\.shp$", level),
                      recursive = TRUE, full.names = TRUE)
    if (!length(shp)) {
      # any shapefile as fallback
      shp <- list.files(exdir, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
    }
    if (length(shp)) {
      sf <- suppressMessages(sf::st_read(shp[1], quiet = TRUE))
    } else {
      # try geodatabase
      gdb <- list.files(exdir, pattern = "\\.gdb$", recursive = TRUE, full.names = TRUE)
      if (!length(gdb)) return(NULL)
      # layer name usually includes level; read first layer if unknown
      layers <- suppressWarnings(sf::st_layers(gdb[1])$name)
      lyr <- if (length(grep(sprintf("lev%02d", level), layers))) {
        layers[grep(sprintf("lev%02d", level), layers)][1]
      } else layers[1]
      sf <- suppressMessages(sf::st_read(gdb[1], layer = lyr, quiet = TRUE))
    }
    # normalize id column
    if (!"hybas_id" %in% names(sf)) {
      if ("HYBAS_ID" %in% names(sf)) sf$hybas_id <- sf$HYBAS_ID
      else if ("HYBAS_ID" %in% toupper(names(sf))) {
        nm <- names(sf)[toupper(names(sf)) == "HYBAS_ID"][1]
        sf$hybas_id <- sf[[nm]]
      } else {
        sf$hybas_id <- NA_integer_
      }
    }
    # keep only id + geometry to stay light
    geom_col <- attr(sf, "sf_column") %||% "geometry"
    sf[, c("hybas_id", geom_col)]
  }

  sfs <- Filter(Negate(is.null), lapply(zips, read_one))
  if (!length(sfs)) stop("HydroBASINS: unzip/read produced no layers.")
  out <- do.call(rbind, sfs)
  # make valid just in case
  out <- tryCatch(suppressWarnings(sf::st_make_valid(out)), error = function(e) out)
  out
}

#' GMBA hierarchy-aware point-in-polygon join (0–360° safe)
#'
#' For each input point, finds all intersecting GMBA Mountains v2 polygons and
#' returns the **single most specific** polygon label per point. “Specificity”
#' is chosen by either the deepest hierarchy level (preferred) or the smallest
#' polygon area (fallback).
#'
#' @param pts_sf  An `sf` POINT layer (any CRS).
#' @param gmba_sf An `sf` POLYGON/MULTIPOLYGON layer for GMBA v2 (any CRS). Must
#'   contain a label column named in `name_col` (default `"geoname_gmba"`).
#' @param name_col Character scalar; column in `gmba_sf` to return as the label
#'   (e.g., `"geoname_gmba"`, or a native GMBA name/ID column you prepared).
#' @param level_cols Character vector of candidate columns encoding hierarchical
#'   depth (e.g., `c("LEVEL","level","HIERARCHY","hierarchy")`). The first
#'   present column that can be coerced to numeric is used. Higher numbers are
#'   treated as deeper (more specific) levels.
#' @param prefer One of `c("deepest_level","smallest_area")`. If
#'   `"deepest_level"` is requested but no usable level column exists, the
#'   function falls back to `"smallest_area"`.
#' @param quiet Logical; suppress progress messages.
#'
#' @return A character vector of length `nrow(pts_sf)` with the chosen GMBA
#'   label for each point, or `NA` when a point intersects no GMBA polygon.
#'
#' @details
#' - **0–360° safety:** If `gmba_sf` spans 0–360° longitudes (xmin ≥ 0 & xmax > 180),
#'   points with negative longitudes are shifted by +360° for the spatial join,
#'   mirroring the behavior used elsewhere in the pipeline.
#' - **Tie-breaking:** When multiple polygons intersect a point:
#'   1) Use the polygon with the **deepest hierarchy level** (largest numeric
#'      value in the chosen `level_cols`).
#'   2) If levels are unavailable or tied, pick the **smallest area** polygon
#'      (computed in a projected CRS when possible).
#'
#' @examples
#' \dontrun{
#' pts  <- sf::st_as_sf(data.frame(x=c(8, 86), y=c(46, 28)),
#'                      coords = c("x","y"), crs = 4326)
#' gmba <- load_gmba_v2()              # your loader; ensure it has `geoname_gmba`
#' labs <- .join_gmba_labels(pts, gmba) # character vector of labels
#' }
#'
#' @keywords internal
#' @noRd
.join_gmba_labels <- function(
    pts_sf, gmba_sf,
    name_col   = "geoname_gmba",
    level_cols = c("LEVEL","level","HIERARCHY","hierarchy"),
    prefer     = c("deepest_level","smallest_area"),
    quiet      = FALSE
) {
  stopifnot(inherits(pts_sf, "sf"), inherits(gmba_sf, "sf"), name_col %in% names(gmba_sf))
  prefer <- match.arg(prefer)

  # Work in WGS84
  pts4326  <- if (is.na(sf::st_crs(pts_sf))  || sf::st_crs(pts_sf)  != sf::st_crs(4326)) sf::st_transform(pts_sf,  4326) else pts_sf
  gmba4326 <- if (is.na(sf::st_crs(gmba_sf)) || sf::st_crs(gmba_sf) != sf::st_crs(4326)) sf::st_transform(gmba_sf, 4326) else gmba_sf

  .is_0360 <- function(x) {
    bb <- sf::st_bbox(x)
    isTRUE(!is.na(bb["xmin"]) && !is.na(bb["xmax"]) && bb["xmin"] >= 0 && bb["xmax"] > 180)
  }

  if (!.is_0360(gmba4326)) {
    pts_for_join <- pts4326
  } else {
    cc <- sf::st_coordinates(pts4326)
    lon360 <- ifelse(cc[,1] < 0, cc[,1] + 360, cc[,1])
    lat    <- cc[,2]
    tmp <- cbind(sf::st_drop_geometry(pts4326), ..lon360 = lon360, ..lat = lat)
    pts_for_join <- sf::st_as_sf(tmp, coords = c("..lon360","..lat"), crs = 4326, remove = TRUE)
    if (!quiet) message("GMBA join using 0–360° point shift")
  }

  # Specificity metrics
  lvl_col <- intersect(level_cols, names(gmba4326))
  lvl_vec <- if (length(lvl_col)) suppressWarnings(as.numeric(gmba4326[[lvl_col[1]]])) else rep(NA_real_, nrow(gmba4326))

  gmba_m <- tryCatch(sf::st_transform(gmba4326, 3857), error = function(e) gmba4326)
  areas  <- tryCatch(as.numeric(sf::st_area(gmba_m)),    error = function(e) rep(NA_real_, nrow(gmba4326)))

  hits <- sf::st_intersects(pts_for_join, gmba4326)

  pick_one <- function(cands) {
    if (!length(cands)) return(NA_integer_)
    if (prefer == "deepest_level" && any(!is.na(lvl_vec[cands]))) {
      cands <- cands[order(-lvl_vec[cands], areas[cands], na.last = TRUE)]
      return(cands[1])
    }
    cands[order(areas[cands], na.last = TRUE)][1]
  }

  sel_idx <- vapply(hits, pick_one, integer(1))
  out <- rep(NA_character_, nrow(pts_sf))
  ok  <- !is.na(sel_idx)
  out[ok] <- as.character(gmba4326[[name_col]][sel_idx[ok]])

  if (!quiet) message("GMBA join (most-specific): ", sum(ok), " / ", length(out), " labeled")
  out
}
