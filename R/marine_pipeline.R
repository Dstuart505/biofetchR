#' Harmonize column types across a list of sf/data.frames
#'
#' Ensures that columns with the same name share a common type across all
#' elements, so they can be safely combined with `dplyr::bind_rows()`.
#' The common type is chosen by simple precedence:
#' character/factor > POSIXct > Date > numeric > logical.
#' The geometry column is preserved and all `sf` objects are coerced to EPSG:4326.
#'
#' @param lst A list of `sf` or `data.frame` objects.
#' @return A list of objects with harmonized column types.
#' @keywords internal
harmonize_column_types <- function(lst) {
  if (length(lst) == 0) return(lst)

  # Standardize geometry column name and CRS on sf objects
  lst <- lapply(lst, function(x) {
    if (inherits(x, "sf")) {
      # CRS → EPSG:4326 (lon/lat)
      if (!identical(sf::st_crs(x), sf::st_crs(4326))) {
        x <- sf::st_transform(x, 4326)
      }
      # Ensure sf column is named "geometry"
      sf_col <- attr(x, "sf_column", exact = TRUE)
      if (!identical(sf_col, "geometry")) {
        names(x)[names(x) == sf_col] <- "geometry"
        attr(x, "sf_column") <- "geometry"
      }
    }
    x
  })

  all_cols <- Reduce(union, lapply(lst, names))
  all_cols <- setdiff(all_cols, "geometry")  # handle geometry implicitly

  # Determine a target type per column, then cast/add missing columns
  for (col in all_cols) {
    # collect classes present
    classes <- unique(unlist(lapply(lst, function(x) if (col %in% names(x)) class(x[[col]])[1] else NULL)))

    target <- if (any(classes %in% c("character","factor"))) {
      "character"
    } else if (any(classes %in% c("POSIXct","POSIXt"))) {
      "POSIXct"
    } else if (any(classes %in% c("Date"))) {
      "Date"
    } else if (any(classes %in% c("numeric","double","integer"))) {
      "numeric"
    } else if (any(classes %in% c("logical"))) {
      "logical"
    } else {
      "character"  # safe fallback
    }

    for (i in seq_along(lst)) {
      n <- nrow(lst[[i]])
      if (!(col %in% names(lst[[i]]))) {
        # add missing column of target type
        lst[[i]][[col]] <- switch(
          target,
          "character" = rep(NA_character_, n),
          "POSIXct"  = rep(as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC"), n),
          "Date"     = rep(as.Date(NA_real_, origin = "1970-01-01"), n),
          "numeric"  = rep(NA_real_, n),
          "logical"  = rep(NA, n),
          rep(NA_character_, n)
        )
      } else {
        v <- lst[[i]][[col]]
        lst[[i]][[col]] <- switch(
          target,
          "character" = as.character(v),
          "POSIXct"  = {
            if (inherits(v, "POSIXct")) v else suppressWarnings(as.POSIXct(v, tz = "UTC", origin = "1970-01-01"))
          },
          "Date"     = if (inherits(v, "Date")) v else suppressWarnings(as.Date(v)),
          "numeric"  = suppressWarnings(as.numeric(v)),
          "logical"  = as.logical(v),
          as.character(v)
        )
      }
    }
  }

  lst
}

#' Normalize a character vector to valid UTF-8 (best effort)
#' @param x character
#' @return character, UTF-8 where possible
#' @keywords internal
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

#' Process GBIF points against a chosen marine boundary overlay (EEZ by default)
#' with optional alien/native tagging from DWC, GRIIS, WoRMS, GISD, and GBIF distributions
#'
#' Submits one GBIF download per batch of species (global scope), joins the
#' occurrences to a boundary overlay (EEZ by default), optionally performs
#' spatial thinning, and (optionally) annotates each record with an
#' `alien/native/unknown` status using several providers.
#'
#' @section Provider availability & fallbacks:
#' If a provider package is not installed, that provider is skipped
#' (a warning is printed) and the pipeline continues. Use
#' [install_optional_deps()] to install providers. Mapping order can be
#' controlled with `*_first` flags. Typical order when `griis_first = TRUE` is:
#' GRIIS → DWC → WoRMS → GBIF → lookup/custom. You can also enable GISD:
#' GISD can be tried early (set `gisd_first = TRUE`) or as a fill-in later.
#'
#' @section ISO mapping:
#' Providers are matched by ISO country codes. The overlay join ideally carries
#' ISO-3 (preferred) and/or ISO-2 codes; the pipeline also honors `iso3_best`/`iso2_best`
#' (if the join function created them) and can convert ISO2↔ISO3 via {countrycode}.
#' When ISO3 is missing, the function falls back to a world-country overlay
#' to populate `iso3_best` on the points.
#'
#' @section Debugging:
#' Set `options(biofetchR.debug = TRUE)` to print per-provider hit counts.
#'
#' @param df Data frame with at least a column `species` (binomial name).
#' @param output_dir Directory for CSV outputs (created if missing).
#' @param user,pwd,email GBIF credentials for download requests.
#' @param batch_size Integer; number of species per GBIF download batch.
#' @param apply_thinning Logical; apply spatial thinning (per-region) to points?
#' @param dist_km Numeric; minimum distance (km) for thinning.
#' @param return_all_results Logical; return combined results in memory?
#' @param export_summary Logical; write a CSV summary table to `output_dir`.
#' @param store_in_memory Logical; keep per-species data in memory (vs. only CSV).
#' @param use_planar Logical; if TRUE, disable s2 for joins/thinning (planar ops).
#'
#' @param add_status Logical; if TRUE, add per-record `alien/native/unknown`
#'   column based on selected providers.
#' @param status_method Character; one of `"dwc"`, `"lookup"`, `"custom"`,
#'   `"dwc_then_lookup"`.
#' @param status_lookup Optional data.frame for `"lookup"`/`"dwc_then_lookup"`.
#'   Columns: `species` (binomial) and one of `mrgid` (integer), `region_id`
#'   (overlay region name; UTF-8), or `region_label_safe` (ASCII-safe). Must include
#'   `status` in `c("alien","native","unknown")`.
#' @param status_fun Optional function for `"custom"`; signature
#'   `function(species, region_human, sf)` returning a scalar or `nrow(sf)` vector
#'   with values in `c("alien","native","unknown")`.
#' @param status_col Name of the new status column (default `"alien_status"`).
#'
#' @param use_worms Logical; if TRUE, fill from WoRMS distributions by MRGID
#'   (requires {worrms}). Missing package ⇒ skipped with a warning.
#' @param worms_first Logical; if TRUE, prefer WoRMS before DWC.
#'
#' @param use_griis Logical; if TRUE, fill from GRIIS country/territory
#'   checklists via {originr}. Missing package ⇒ skipped with a warning.
#' @param griis_first Logical; if TRUE, prefer GRIIS before DWC.
#'
#' @param use_gisd Logical; if TRUE, fill using GISD (Global Invasive Species
#'   Database) via {originr}. Missing package ⇒ skipped with a warning.
#' @param gisd_first Logical; if TRUE, prefer GISD before DWC.
#'
#' @param use_gbif_distributions Logical; if TRUE, fill using GBIF species-level
#'   distributions via {rgbif}. Requires ISO3 on joined records.
#' @param gbif_first Logical; if TRUE, prefer GBIF distributions before DWC.
#'
#' @param overlay Character; which boundary overlay to use. Defaults to `"eez"`.
#'   If `overlay != "eez"`, you must supply `overlay_sf`.
#'   Common {mregions2} overlays include `"lme"` (Large Marine Ecosystems),
#'   `"meow"` (Marine Ecoregions of the World), `"iho"` (seas & oceans),
#'   and `"fao"` (FAO fishing areas).
#' @param overlay_sf sf POLYGON/MULTIPOLYGON layer to join when `overlay != "eez"`.
#'   Fetch with {mregions2}, for example:
#'   \preformatted{
#'   # install.packages("mregions2")
#'   lme_sf  <- mregions2::mrp_get("lme")   # Large Marine Ecosystems
#'   meow_sf <- mregions2::mrp_get("meow")  # Marine Ecoregions of the World
#'   iho_sf  <- mregions2::mrp_get("iho")   # IHO seas & oceans
#'   fao_sf  <- mregions2::mrp_get("fao")   # FAO fishing areas
#'   }
#'
#' @return If `return_all_results = TRUE` and `store_in_memory = TRUE`, returns
#'   a combined `tibble` of all per-species×region results (GBIF fields, overlay attrs,
#'   and—if `add_status=TRUE`—the `status_col`). Otherwise, returns `invisible(NULL)`.
#'
#' @examples
#' \dontrun{
#' # EEZ (no overlay_sf needed)
#' res_eez <- process_gbif_eez_pipeline(
#'   df = tibble::tibble(species = "Carcinus maenas"),
#'   output_dir = "out_eez",
#'   user = user, pwd = pwd, email = email,
#'   overlay = "eez",
#'   add_status = TRUE, use_griis = TRUE, griis_first = TRUE
#' )
#'
#' # LME overlay (must fetch via mregions2)
#' lme_sf <- mregions2::mrp_get("lme")
#' res_lme <- process_gbif_eez_pipeline(
#'   df = tibble::tibble(species = "Carcinus maenas"),
#'   output_dir = "out_lme",
#'   user = user, pwd = pwd, email = email,
#'   overlay = "lme",
#'   overlay_sf = lme_sf,
#'   add_status = TRUE, use_griis = TRUE, griis_first = TRUE
#' )
#' }
#'
#' @export
process_gbif_eez_pipeline <- function(
    df,
    output_dir,
    user, pwd, email,
    batch_size = 5,
    apply_thinning = TRUE,
    dist_km = 5,
    return_all_results = TRUE,
    export_summary = TRUE,
    store_in_memory = TRUE,
    use_planar = FALSE,
    add_status = FALSE,
    status_method = c("dwc","lookup","custom","dwc_then_lookup"),
    status_lookup = NULL,
    status_fun = NULL,
    status_col = "alien_status",
    use_worms = FALSE,
    worms_first = FALSE,
    use_griis = FALSE,
    griis_first = FALSE,
    use_gisd = FALSE,
    gisd_first = FALSE,
    use_gbif_distributions = FALSE,
    gbif_first = FALSE,
    overlay = "eez",
    overlay_sf = NULL
) {
  # ---------------- helpers ----------------

  # Generic overlay join (EEZ routes to existing eez_join(); otherwise join overlay_sf)
  overlay_join <- function(named_points_list, overlay, overlay_sf, use_planar = FALSE) {
    if (tolower(overlay) == "eez") {
      return(eez_join(named_points_list, use_planar = use_planar))
    }

    if (is.null(overlay_sf) || !inherits(overlay_sf, "sf")) {
      stop(
        "overlay != 'eez' but `overlay_sf` is NULL or not an sf object.\n",
        "Fetch polygons via {mregions2} and pass them to `overlay_sf`, e.g.:\n\n",
        "  # install.packages('mregions2')\n",
        "  lme_sf  <- mregions2::mrp_get('lme')   # Large Marine Ecosystems\n",
        "  meow_sf <- mregions2::mrp_get('meow')  # Marine Ecoregions of the World\n",
        "  iho_sf  <- mregions2::mrp_get('iho')   # IHO seas & oceans\n",
        "  fao_sf  <- mregions2::mrp_get('fao')   # FAO fishing areas\n\n",
        "Then call, for example:\n",
        "  process_gbif_eez_pipeline(..., overlay = 'lme', overlay_sf = lme_sf)\n"
      )
    } else {
      cli::cli_alert_info(paste0(
        "Using user-supplied overlay '", toupper(overlay), "' (", nrow(overlay_sf), " polygons)."
      ))
    }

    # Normalize polygon attrs so downstream code can rely on `geoname` + `mrgid`
    poly <- overlay_sf
    name_candidates <- c(
      "geoname","GeoName","geoname_en","name","NAME","name_en","NAME_EN",
      "lme_name","LME_NAME","ecoregion","ECOREGION","province","PROVINCE",
      "realm","REALM","ocean","OCEAN","sea_name","SEA_NAME","ECO_NAME",
      "F_AREA","SUB_AREA","FAO_NAME","FAO"
    )
    nm <- intersect(name_candidates, names(poly))
    if (length(nm)) {
      poly$geoname <- as.character(poly[[nm[1]]])
    } else {
      poly$geoname <- paste0(toupper(overlay), "_", seq_len(nrow(poly)))
    }

    id_candidates <- c("mrgid","MRGID","id","ID","OBJECTID","ECO_CODE","LME_NUMBER","LME_ID","F_CODE","POLY_ID")
    im <- intersect(id_candidates, names(poly))
    if (length(im)) {
      suppressWarnings(poly$mrgid <- as.integer(poly[[im[1]]]))
    } else if (!"mrgid" %in% names(poly)) {
      poly$mrgid <- NA_integer_
    }

    old_s2 <- sf::sf_use_s2()
    on.exit(try(sf::sf_use_s2(old_s2), silent = TRUE), add = TRUE)
    if (isTRUE(use_planar)) sf::sf_use_s2(FALSE)

    poly <- suppressWarnings(try(sf::st_make_valid(poly), silent = TRUE))
    if (inherits(poly, "try-error")) poly <- overlay_sf
    bad <- !sf::st_is_valid(poly)
    if (any(bad)) poly <- sf::st_buffer(poly, 0)

    out <- list()
    for (nm in names(named_points_list)) {
      pts <- named_points_list[[nm]]
      if (!inherits(pts, "sf")) next

      poly2 <- try(sf::st_transform(poly, sf::st_crs(pts)), silent = TRUE)
      if (inherits(poly2, "try-error")) poly2 <- poly

      jj <- suppressWarnings(try(sf::st_join(pts, poly2, join = sf::st_within, left = TRUE), silent = TRUE))
      if (inherits(jj, "try-error") || all(is.na(jj$geoname))) {
        jj <- suppressWarnings(try(sf::st_join(pts, poly2, join = sf::st_intersects, left = TRUE), silent = TRUE))
      }
      if (!"geoname" %in% names(jj)) jj$geoname <- NA_character_
      if (!"mrgid"   %in% names(jj)) jj$mrgid   <- NA_integer_
      out[[nm]] <- jj
    }
    out
  }

  # Fill missing iso3_best on point data by joining to country polygons (planar-first; no lwgeom)
  ensure_iso3_on_points <- function(x) {
    if (!inherits(x, "sf")) return(x)
    if (!requireNamespace("rnaturalearth", quietly = TRUE)) return(x)
    if (!requireNamespace("countrycode",   quietly = TRUE)) return(x)

    if (!"iso3_best" %in% names(x)) x$iso3_best <- NA_character_
    need <- is.na(x$iso3_best) | x$iso3_best == ""
    if (!any(need, na.rm = TRUE)) return(x)

    world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")
    iso3_col <- intersect(c("adm0_a3","iso_a3_eh","iso_a3","ADM0_A3","ISO_A3_EH","ISO_A3"), names(world))
    if (!length(iso3_col)) return(x)
    world <- world[, c(iso3_col[1], "geometry")]
    names(world)[1] <- "iso3"

    old_s2 <- sf::sf_use_s2()
    on.exit(try(sf::sf_use_s2(old_s2), silent = TRUE), add = TRUE)

    sf::sf_use_s2(FALSE)
    world <- suppressWarnings(try(sf::st_make_valid(world), silent = TRUE))
    if (inherits(world, "try-error")) {
      world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")[, c(iso3_col[1], "geometry")]
      names(world)[1] <- "iso3"
    }
    if (any(!sf::st_is_valid(world))) world <- suppressWarnings(sf::st_buffer(world, 0))

    try_planar_join <- function() {
      pts <- try(sf::st_transform(x[need, ], 3857),  silent = TRUE)
      ply <- try(sf::st_transform(world,   3857),    silent = TRUE)
      if (inherits(pts, "try-error") || inherits(ply, "try-error")) stop("transform_failed")
      suppressWarnings(sf::st_join(pts, ply, join = sf::st_within, left = TRUE))
    }
    j <- try(try_planar_join(), silent = TRUE)

    if (inherits(j, "try-error")) {
      sf::sf_use_s2(TRUE)
      ply <- try(sf::st_transform(world, sf::st_crs(x)), silent = TRUE)
      if (!inherits(ply, "try-error")) {
        j <- suppressWarnings(try(sf::st_join(x[need, ], ply, join = sf::st_within, left = TRUE), silent = TRUE))
      }
      if (inherits(j, "try-error")) return(x)
    }

    fill <- suppressWarnings(countrycode::countrycode(j$iso3, "iso3c", "iso3c"))
    good <- !is.na(fill) & nzchar(fill)
    if (any(good)) {
      idx <- which(need)
      x$iso3_best[idx[good]] <- toupper(fill[good])
    }
    x
  }

  make_safe_label <- function(x) {
    if (length(x) == 0 || is.na(x) || !nzchar(x)) return("EEZ_UNKNOWN")
    x <- as.character(x[[1L]])
    to_utf8 <- function(s) {
      s1 <- suppressWarnings(iconv(s, from = "",        to = "UTF-8", sub = ""))
      if (is.na(s1)) s1 <- suppressWarnings(iconv(s, from = "latin1", to = "UTF-8", sub = ""))
      if (is.na(s1)) s1 <- suppressWarnings(iconv(s, from = "CP1252", to = "UTF-8", sub = ""))
      if (is.na(s1)) s1 <- enc2utf8(s)
      s1
    }
    y <- to_utf8(x); if (is.na(y)) y <- x
    y2 <- suppressWarnings(iconv(y, from = "UTF-8", to = "ASCII//TRANSLIT", sub = ""))
    if (is.na(y2)) {
      if (requireNamespace("stringi", quietly = TRUE)) {
        y2 <- stringi::stri_trans_general(y, "Latin-ASCII")
      } else y2 <- y
    }
    y2 <- gsub("[^[:alnum:]]+", "_", y2, perl = TRUE, useBytes = TRUE)
    y2 <- gsub("_+", "_", y2, perl = TRUE, useBytes = TRUE)
    y2 <- gsub("^_+|_+$", "", y2, perl = TRUE, useBytes = TRUE)
    if (!nzchar(y2)) y2 <- "EEZ_UNKNOWN"
    substr(y2, 1, 150)
  }

  utf8_safe <- function(x) {
    if (is.null(x)) return(x)
    y <- suppressWarnings(iconv(x, from = "", to = "UTF-8", sub = ""))
    bad <- is.na(y)
    if (any(bad)) y[bad] <- suppressWarnings(iconv(x[bad], from = "latin1", to = "UTF-8", sub = ""))
    bad <- is.na(y)
    if (any(bad)) y[bad] <- suppressWarnings(iconv(x[bad], from = "CP1252", to = "UTF-8", sub = ""))
    y[is.na(y)] <- x[is.na(y)]
    y
  }

  map_dwc_status <- function(establishmentMeans, degreeOfEstablishment, occurrenceStatus = NULL) {
    z <- rep(NA_character_, length(establishmentMeans))
    norm <- function(v) tolower(trimws(as.character(v)))
    em <- norm(establishmentMeans); de <- norm(degreeOfEstablishment); os <- norm(occurrenceStatus)
    native_terms <- c("native","endemic","domesticated")
    alien_terms  <- c("introduced","invasive","naturalised","naturalized","alien","exotic","non-native","managed")
    deg_alien    <- c("invasive","naturalised","naturalized","casual","casual-introduced","established non-native","colonising","widespread non-native")
    z[em %in% native_terms] <- "native"
    z[em %in% alien_terms]  <- "alien"
    z[is.na(z) & de %in% deg_alien] <- "alien"
    z[is.na(z)] <- "unknown"
    z
  }

  lookup_status_scalar <- function(species, region_human, raw_sf, lk) {
    if (is.null(lk) || !nrow(lk)) return(NA_character_)
    s_match <- tolower(trimws(species))
    lk_s <- lk[tolower(trimws(lk$species)) == s_match, , drop = FALSE]
    if (!nrow(lk_s)) return(NA_character_)
    if ("mrgid" %in% names(lk_s) && "mrgid" %in% names(raw_sf)) {
      mrgid_vals <- unique(raw_sf$mrgid); mrgid_vals <- mrgid_vals[!is.na(mrgid_vals)]
      if (length(mrgid_vals)) {
        lk_m <- lk_s[lk_s$mrgid %in% mrgid_vals, , drop = FALSE]
        if (nrow(lk_m)) return(as.character(lk_m$status[1]))
      }
    }
    if ("region_id" %in% names(lk_s)) {
      reg_h <- utf8_safe(region_human)
      lk_r <- lk_s[lk_s$region_id == reg_h, , drop = FALSE]
      if (nrow(lk_r)) return(as.character(lk_r$status[1]))
    }
    if ("region_label_safe" %in% names(lk_s)) {
      safe <- make_safe_label(region_human)
      lk_a <- lk_s[lk_s$region_label_safe == safe, , drop = FALSE]
      if (nrow(lk_a)) return(as.character(lk_a$status[1]))
    }
    NA_character_
  }

  .dbg <- isTRUE(getOption("biofetchR.debug", FALSE))
  .say <- function(...) if (.dbg) cli::cli_alert_info(paste0(...))

  # ---------------- provider mappers (WoRMS / GRIIS / GISD / GBIF) ----------------
  worms_status_cache <- new.env(parent = emptyenv())
  get_worms_map_for_species <- function(species) {
    if (!requireNamespace("worrms", quietly = TRUE)) return(NULL)
    key <- tolower(trimws(species))
    if (exists(key, envir = worms_status_cache, inherits = FALSE)) return(get(key, envir = worms_status_cache))
    recs <- try(worrms::wm_records_name(name = species, marine_only = TRUE), silent = TRUE)
    if (inherits(recs, "try-error") || is.null(recs) || length(recs) == 0) { assign(key, NULL, worms_status_cache); return(NULL) }
    recs <- as.data.frame(recs, stringsAsFactors = FALSE)
    pick <- which(tolower(recs$status) == "accepted")[1]; if (is.na(pick)) pick <- 1L
    aphia <- recs$AphiaID[pick]
    dist <- try(worrms::wm_distribution(aphia), silent = TRUE)
    if (inherits(dist, "try-error") || is.null(dist) || length(dist) == 0) { assign(key, NULL, worms_status_cache); return(NULL) }
    dist <- as.data.frame(dist, stringsAsFactors = FALSE)
    norm_status <- function(s) {
      s <- tolower(trimws(as.character(s)))
      if (s %in% c("native","endemic")) "native"
      else if (s %in% c("alien","introduced","invasive","naturalised","naturalized","casual","established non-native")) "alien"
      else "unknown"
    }
    if (!"status" %in% names(dist)) dist$status <- NA_character_
    dist$status_std <- vapply(dist$status, norm_status, character(1))
    if (!"locationID" %in% names(dist)) dist$locationID <- NA
    out <- unique(dist[, c("locationID","status_std")])
    names(out) <- c("mrgid","status")
    suppressWarnings(out$mrgid <- as.integer(out$mrgid))
    assign(key, out, worms_status_cache)
    out
  }

  griis_status_cache <- new.env(parent = emptyenv())
  get_griis_map_for_species <- function(species) {
    if (!requireNamespace("originr", quietly = TRUE)) return(NULL)
    key <- tolower(trimws(species))
    if (exists(key, envir = griis_status_cache, inherits = FALSE)) return(get(key, envir = griis_status_cache))
    res <- NULL
    fun1 <- get0("griis", envir = asNamespace("originr"), inherits = FALSE)
    if (is.function(fun1)) res <- suppressWarnings(try(fun1(species = species), silent = TRUE))
    if (inherits(res, "try-error") || is.null(res)) res <- NULL
    if (is.null(res)) {
      fun2 <- get0("griis_taxa", envir = asNamespace("originr"), inherits = FALSE)
      if (is.function(fun2)) res <- suppressWarnings(try(fun2(species), silent = TRUE))
      if (inherits(res, "try-error")) res <- NULL
    }
    if (is.null(res) || length(res) == 0) { assign(key, NULL, griis_status_cache); return(NULL) }
    df <- try(as.data.frame(res, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(df, "try-error") || !nrow(df)) { assign(key, NULL, griis_status_cache); return(NULL) }

    nms <- names(df)
    pick_col <- function(choices) { ch <- intersect(choices, nms); if (length(ch)) ch[1] else NA_character_ }
    col_status <- pick_col(c("status","establishmentMeans","origin","alien_native","alien.status","Alien_status"))
    if (is.na(col_status)) { assign(key, NULL, griis_status_cache); return(NULL) }

    col_country  <- pick_col(c("country","Country","jurisdiction","territory","location"))
    col_iso3_in  <- pick_col(c("iso3","ISO3","countryCode","country_code","iso","countryCode_alpha3","alpha3"))
    col_iso2_in  <- pick_col(c("countryCode_alpha2","alpha2","ISO2","iso2","country_code_alpha2"))

    stat <- tolower(trimws(as.character(df[[col_status]])))
    stat[stat %in% c("alien","introduced","invasive","naturalised","naturalized")] <- "alien"
    stat[stat %in% c("native","endemic")] <- "native"
    stat[!(stat %in% c("alien","native")) | is.na(stat)] <- "unknown"

    iso3 <- if (!is.na(col_iso3_in)) toupper(trimws(as.character(df[[col_iso3_in]]))) else rep(NA_character_, nrow(df))
    iso2 <- if (!is.na(col_iso2_in)) toupper(trimws(as.character(df[[col_iso2_in]]))) else rep(NA_character_, nrow(df))

    if (requireNamespace("countrycode", quietly = TRUE)) {
      need3 <- is.na(iso3) | !nzchar(iso3)
      if (any(need3) && !all(is.na(iso2))) {
        iso3[need3] <- suppressWarnings(countrycode::countrycode(iso2, "iso2c", "iso3c"))[need3]
      }
      need3 <- is.na(iso3) | !nzchar(iso3)
      if (any(need3) && !is.na(col_country)) {
        iso3[need3] <- suppressWarnings(countrycode::countrycode(df[[col_country]], "country.name", "iso3c"))[need3]
      }
    }

    out <- unique(data.frame(
      iso3   = toupper(ifelse(is.na(iso3), "", iso3)),
      iso2   = toupper(ifelse(is.na(iso2), "", iso2)),
      status = stat, stringsAsFactors = FALSE
    ))
    out <- out[(nzchar(out$iso3) | nzchar(out$iso2)), , drop = FALSE]
    if (!nrow(out)) { assign(key, NULL, griis_status_cache); return(NULL) }

    assign(key, out, griis_status_cache)
    out
  }

  gisd_status_cache <- new.env(parent = emptyenv())
  get_gisd_map_for_species <- function(species) {
    if (!requireNamespace("originr", quietly = TRUE)) return(NULL)
    key <- paste0("gisd::", tolower(trimws(species)))
    if (exists(key, envir = gisd_status_cache, inherits = FALSE)) return(get(key, envir = gisd_status_cache))

    res <- suppressWarnings(try(originr::gisd(species), silent = TRUE))
    if (inherits(res, "try-error") || is.null(res)) { assign(key, NULL, gisd_status_cache); return(NULL) }
    df <- try(as.data.frame(res, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(df, "try-error") || !nrow(df)) { assign(key, NULL, gisd_status_cache); return(NULL) }

    nms <- names(df)
    pick <- function(choices) { ch <- intersect(choices, nms); if (length(ch)) ch[1] else NA_character_ }
    col_country <- pick(c("country","Country","location","where"))
    col_status  <- pick(c("status","Status","occurrence","invasiveness","invasive"))
    if (is.na(col_country) || is.na(col_status)) { assign(key, NULL, gisd_status_cache); return(NULL) }

    stat <- tolower(trimws(as.character(df[[col_status]])))
    alien_terms <- c("invasive","alien","introduced","established","established non-native",
                     "naturalised","naturalized","casual")
    nat_terms   <- c("native","endemic")
    out_status  <- ifelse(stat %in% alien_terms, "alien",
                          ifelse(stat %in% nat_terms, "native", "unknown"))

    iso3 <- rep(NA_character_, nrow(df))
    if (requireNamespace("countrycode", quietly = TRUE)) {
      iso3 <- suppressWarnings(countrycode::countrycode(df[[col_country]], "country.name", "iso3c"))
    }
    out <- unique(data.frame(iso3 = toupper(iso3), status = out_status, stringsAsFactors = FALSE))
    out <- out[!is.na(out$iso3) & nzchar(out$iso3), , drop = FALSE]
    if (!nrow(out)) { assign(key, NULL, gisd_status_cache); return(NULL) }

    assign(key, out, gisd_status_cache)
    out
  }

  gbif_status_cache <- new.env(parent = emptyenv())
  get_gbif_map_for_species <- function(species) {
    if (!requireNamespace("rgbif", quietly = TRUE)) return(NULL)
    key <- tolower(trimws(species))
    if (exists(key, envir = gbif_status_cache, inherits = FALSE)) return(get(key, envir = gbif_status_cache))

    bb <- try(rgbif::name_backbone(name = species), silent = TRUE)
    if (inherits(bb, "try-error") || is.null(bb) || is.na(bb$usageKey)) { assign(key, NULL, gbif_status_cache); return(NULL) }

    dis <- try(rgbif::name_usage(key = bb$usageKey, data = "distributions", limit = 1000), silent = TRUE)
    if (inherits(dis, "try-error") || is.null(dis) || !NROW(dis$data)) { assign(key, NULL, gbif_status_cache); return(NULL) }
    df <- dis$data

    norm_status <- function(s) {
      s <- tolower(trimws(as.character(s)))
      if (s %in% c("native","endemic")) "native"
      else if (s %in% c("introduced","invasive","alien","naturalised","naturalized")) "alien"
      else "unknown"
    }
    if (!"establishmentMeans" %in% names(df)) df$establishmentMeans <- NA_character_
    stat <- vapply(df$establishmentMeans, norm_status, character(1))

    iso3 <- rep(NA_character_, nrow(df))
    if ("countryCode" %in% names(df) && requireNamespace("countrycode", quietly = TRUE)) {
      cc2 <- toupper(trimws(as.character(df$countryCode)))
      iso3 <- suppressWarnings(countrycode::countrycode(cc2, "iso2c", "iso3c"))
    }
    if (all(is.na(iso3)) && "country" %in% names(df) && requireNamespace("countrycode", quietly = TRUE)) {
      iso3 <- suppressWarnings(countrycode::countrycode(df$country, "country.name", "iso3c"))
    }

    out <- unique(data.frame(iso3 = toupper(iso3), status = stat, stringsAsFactors = FALSE))
    out <- out[!is.na(out$iso3) & nzchar(out$iso3), , drop = FALSE]
    if (!nrow(out)) { assign(key, NULL, gbif_status_cache); return(NULL) }

    assign(key, out, gbif_status_cache)
    out
  }

  coalesce2 <- function(a, b) {
    if (is.null(a) && is.null(b)) return(NULL)
    a <- if (is.null(a)) rep(NA_character_, length(b)) else as.character(a)
    b <- if (is.null(b)) rep(NA_character_, length(a)) else as.character(b)
    i <- is.na(a) | a == ""
    a[i] <- b[i]
    a
  }

  iso_cols_from_sf <- function(x) {
    pick <- function(cands) { nm <- intersect(cands, names(x)); if (length(nm)) nm[1] else NULL }
    iso3 <- if ("iso3_best" %in% names(x)) toupper(trimws(as.character(x$iso3_best))) else NULL
    iso2 <- if ("iso2_best" %in% names(x)) toupper(trimws(as.character(x$iso2_best))) else NULL
    if (is.null(iso3)) { nm3 <- pick(c("ISO_SOV1","ISO_Ter1","ISO_Ter2","ISO_SOV2","SOV_A3","ISO_A3")); if (!is.null(nm3)) iso3 <- toupper(trimws(as.character(x[[nm3]]))) }
    if (is.null(iso2)) { nm2 <- pick(c("ISO2","ISO_A2","ISO2_Ter1","ISO2_SOV1")); if (!is.null(nm2)) iso2 <- toupper(trimws(as.character(x[[nm2]]))) }
    if (is.null(iso2) && "countryCode" %in% names(x)) iso2 <- toupper(trimws(as.character(x$countryCode)))
    if (is.null(iso2) && !is.null(iso3) && requireNamespace("countrycode", quietly = TRUE)) {
      iso2 <- toupper(suppressWarnings(countrycode::countrycode(iso3, "iso3c", "iso2c")))
    }
    list(iso2 = if (!is.null(iso2)) iso2 else NULL,
         iso3 = if (!is.null(iso3)) iso3 else NULL)
  }

  harmonize_column_types <- function(sf_list) {
    sf_list <- sf_list[!vapply(sf_list, is.null, logical(1))]
    if (!length(sf_list)) return(sf_list)
    all_cols <- unique(unlist(lapply(sf_list, names), use.names = FALSE))
    ref_idx  <- which.max(vapply(sf_list, nrow, integer(1)))
    ref_crs  <- tryCatch(sf::st_crs(sf_list[[ref_idx]]), error = function(e) NULL)
    get_gcol <- function(x) if (inherits(x, "sf")) attr(x, "sf_column") else NA_character_
    col_class <- function(x, col) {
      if (!col %in% names(x)) return(NA_character_)
      gcol <- get_gcol(x)
      if (!is.na(gcol) && col == gcol) return("sfc")
      class(x[[col]])[1]
    }
    class_matrix <- lapply(all_cols, function(col) { vapply(sf_list, col_class, character(1), col = col) })
    names(class_matrix) <- all_cols
    target_class <- vapply(class_matrix, function(v) {
      v2 <- unique(na.omit(v[v != "sfc"]))
      if (!length(v2)) "character" else if (length(v2) == 1) v2 else "character"
    }, character(1))
    out <- lapply(sf_list, function(x) {
      miss <- setdiff(all_cols, names(x))
      for (m in miss) x[[m]] <- NA_character_
      gcol <- get_gcol(x)
      for (nm in setdiff(all_cols, gcol)) {
        tgt <- target_class[[nm]]
        x[[nm]] <- switch(tgt,
                          character = as.character(x[[nm]]),
                          integer   = suppressWarnings(as.integer(x[[nm]])),
                          numeric   = suppressWarnings(as.numeric(x[[nm]])),
                          logical   = as.logical(x[[nm]]),
                          factor    = as.factor(x[[nm]]),
                          as.character(x[[nm]]))
      }
      if (inherits(x, "sf") && !is.null(ref_crs)) { try(x <- sf::st_set_crs(x, ref_crs), silent = TRUE) }
      x
    })
    out
  }

  # ---------------- checks & setup ----------------

  status_method <- match.arg(status_method)
  if (!"species" %in% names(df)) stop("Marine dataframe must contain a 'species' column.")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (use_planar) sf::sf_use_s2(FALSE)

  has_originr <- requireNamespace("originr", quietly = TRUE)
  has_worrms  <- requireNamespace("worrms",  quietly = TRUE)
  has_rgbif   <- requireNamespace("rgbif",   quietly = TRUE)

  if (use_griis && !has_originr) {
    cli::cli_alert_warning("use_griis=TRUE but package 'originr' is not available. GRIIS fill will be skipped.\nRun: install_optional_deps(griis=TRUE, prefer='runiverse')")
    use_griis <- FALSE
  }
  if (use_gisd && !has_originr) {
    cli::cli_alert_warning("use_gisd=TRUE but package 'originr' is not available. GISD fill will be skipped.\nRun: install_optional_deps(griis=TRUE, prefer='runiverse')")
    use_gisd <- FALSE
  }
  if (use_worms && !has_worrms) {
    cli::cli_alert_warning("use_worms=TRUE but package 'worrms' is not available. WoRMS fill will be skipped.")
    use_worms <- FALSE
  }
  if (use_gbif_distributions && !has_rgbif) {
    cli::cli_alert_warning("use_gbif_distributions=TRUE but package 'rgbif' is not available. GBIF distributions fill will be skipped.")
    use_gbif_distributions <- FALSE
  }

  cli::cli_inform(c("ℹ" = paste0("🌊 Processing ", nrow(df), " marine species — ", toupper(overlay), " overlay")))

  batches      <- split(df, ceiling(seq_len(nrow(df)) / batch_size))
  summary_tbl  <- initialize_summary()
  all_results  <- list()

  # ---------------- main loop ----------------

  for (batch in batches) {
    cli::cli_alert_info(paste0("🌊 Processing batch: ", paste(batch$species, collapse = ", "),
                               " — ", toupper(overlay), " overlay"))
    download_keys <- download_gbif_batch_eez(batch, user, pwd, email)
    if (is.null(download_keys)) next

    imported <- wait_and_import_gbif(download_keys)
    if (length(imported) == 0) next

    if (is.list(imported) && !inherits(imported, "data.frame")) {
      if (is.null(names(imported))) {
        names(imported) <- vapply(imported, function(x) if ("species" %in% names(x)) as.character(x$species[1]) else "item", character(1))
      }
      split_by_species <- imported
    } else {
      if (!"species" %in% names(imported)) stop("wait_and_import_gbif() returned a data.frame without a 'species' column.")
      split_by_species <- split(imported, imported$species)
    }

    for (species_name in names(split_by_species)) {
      sp_df <- split_by_species[[species_name]]
      if (is.null(sp_df) || nrow(sp_df) == 0) next

      if (!inherits(sp_df, "sf")) {
        if (!all(c("decimalLongitude","decimalLatitude") %in% names(sp_df))) next
        sp_df <- sf::st_as_sf(sp_df, coords = c("decimalLongitude","decimalLatitude"), crs = 4326, remove = FALSE)
      }

      # preserve any publisher-provided MRGID before the join (avoid collision)
      if ("mrgid" %in% names(sp_df) && !"mrgid_src" %in% names(sp_df)) {
        sp_df$mrgid_src <- suppressWarnings(as.integer(sp_df$mrgid))
        sp_df$mrgid <- NULL
      }
      if ("locationID" %in% names(sp_df)) {
        lid <- tolower(as.character(sp_df$locationID))
        if (any(grepl("mrgid", lid, useBytes = TRUE, ignore.case = TRUE), na.rm = TRUE)) {
          mrgid_from_loc <- rep(NA_integer_, length(lid))
          sel <- grepl("mrgid", lid, useBytes = TRUE, ignore.case = TRUE)
          mrgid_from_loc[sel] <- suppressWarnings(as.integer(sub(".*?mrgid\\D*(\\d+).*", "\\1", lid[sel])))
          if (!"mrgid_src" %in% names(sp_df)) sp_df$mrgid_src <- NA_integer_
          fill <- is.na(sp_df$mrgid_src) & !is.na(mrgid_from_loc)
          sp_df$mrgid_src[fill] <- mrgid_from_loc[fill]
        }
      }

      named_in <- setNames(list(sp_df), species_name)
      joined_by_species <- overlay_join(named_in, overlay = overlay, overlay_sf = overlay_sf, use_planar = use_planar)

      joined <- joined_by_species[[species_name]]
      if (is.null(joined) || nrow(joined) == 0) {
        summary_tbl <- append_summary_row(summary_tbl, species = species_name,
                                          region_id = NA_character_, region_type = toupper(overlay),
                                          n_total = 0, n_cleaned = 0, n_thinned = 0,
                                          output_file = if (store_in_memory) NA_character_ else NA_character_,
                                          status = "no_data")
        next
      }

      joined$geoname <- ifelse(is.na(joined$geoname), paste0(toupper(overlay), "_UNKNOWN"), joined$geoname)
      joined$geoname <- utf8_safe(joined$geoname)
      split_by_region <- split(joined, joined$geoname, drop = TRUE)

      # name variants (GBIF backbone)
      name_variants_for <- function(name) {
        out <- unique(c(name))
        if (has_rgbif) {
          bb <- try(rgbif::name_backbone(name = name), silent = TRUE)
          if (!inherits(bb, "try-error") && !is.null(bb)) {
            out <- unique(c(out, bb$scientificName, bb$canonicalName))
          }
        }
        unique(out[!is.na(out) & nzchar(out)])
      }

      for (region_human in names(split_by_region)) {
        raw_sf <- split_by_region[[region_human]]
        if (is.null(raw_sf) || nrow(raw_sf) == 0) {
          region_human_disp <- utf8_safe(region_human)
          summary_tbl <- append_summary_row(summary_tbl, species = species_name, region_id = region_human_disp,
                                            region_type = toupper(overlay), n_total = 0, n_cleaned = 0, n_thinned = 0,
                                            output_file = if (store_in_memory) NA_character_ else NA_character_,
                                            status = "no_data")
          next
        }

        n_total    <- nrow(raw_sf)
        thinned_sf <- if (isTRUE(apply_thinning)) thin_spatial_points(raw_sf, dist_km = dist_km, quiet = TRUE) else raw_sf

        # ensure ISO3 exists for provider mapping
        thinned_sf <- ensure_iso3_on_points(thinned_sf)
        n_thinned <- nrow(thinned_sf)

        if (isTRUE(add_status)) {
          status_vec <- rep(NA_character_, nrow(thinned_sf))

          iso_list <- iso_cols_from_sf(thinned_sf)
          iso2_vec <- iso_list$iso2
          iso3_vec <- iso_list$iso3
          iso3_for_matching <- iso3_vec
          if (is.null(iso3_for_matching) && !is.null(iso2_vec) && requireNamespace("countrycode", quietly = TRUE)) {
            iso3_for_matching <- toupper(trimws(as.character(
              suppressWarnings(countrycode::countrycode(iso2_vec, "iso2c", "iso3c"))
            )))
          }

          species_variants <- name_variants_for(species_name)
          get_map_try_variants <- function(get_fun) {
            for (nm in species_variants) {
              mm <- get_fun(nm)
              if (!is.null(mm) && nrow(mm)) return(mm)
            }
            NULL
          }

          if (isTRUE(use_gbif_distributions) && isTRUE(gbif_first) && !is.null(iso3_for_matching)) {
            gm_gbif <- get_map_try_variants(get_gbif_map_for_species)
            if (!is.null(gm_gbif) && nrow(gm_gbif)) {
              idx <- match(iso3_for_matching, gm_gbif$iso3)
              hit <- !is.na(idx); status_vec[hit] <- gm_gbif$status[idx[hit]]
              .say(paste0(species_name, " @ ", region_human, ": GBIF-first matched ", sum(hit), " / ", length(status_vec)))
            }
          }

          if (isTRUE(use_griis) && isTRUE(griis_first) && (!is.null(iso3_vec) || !is.null(iso2_vec))) {
            gm <- get_map_try_variants(get_griis_map_for_species)
            if (!is.null(gm) && nrow(gm)) {
              needs <- is.na(status_vec)
              hit3 <- hit2 <- rep(FALSE, length(status_vec))
              if (!is.null(iso3_vec) && "iso3" %in% names(gm)) {
                idx3 <- match(toupper(iso3_vec), gm$iso3)
                hit3 <- !is.na(idx3) & needs
                status_vec[hit3] <- gm$status[idx3[hit3]]
                needs <- is.na(status_vec)
              }
              if (any(needs) && !is.null(iso2_vec) && "iso2" %in% names(gm)) {
                idx2 <- match(toupper(iso2_vec), toupper(gm$iso2))
                hit2 <- !is.na(idx2) & needs
                status_vec[hit2] <- gm$status[idx2[hit2]]
              }
              .say(paste0(species_name, " @ ", region_human, ": GRIIS-first matched ", sum(hit3)+sum(hit2), " / ", length(status_vec)))
            }
          }

          if (isTRUE(use_gisd) && isTRUE(gisd_first) && !is.null(iso3_for_matching)) {
            gm_gisd <- get_map_try_variants(get_gisd_map_for_species)
            if (!is.null(gm_gisd) && nrow(gm_gisd)) {
              idx <- match(iso3_for_matching, gm_gisd$iso3)
              needs <- is.na(status_vec)
              hit <- !is.na(idx) & needs
              status_vec[hit] <- gm_gisd$status[idx[hit]]
              .say(paste0(species_name, " @ ", region_human, ": GISD-first matched ", sum(hit), " / ", length(status_vec)))
            }
          }

          if (status_method %in% c("dwc","dwc_then_lookup")) {
            em <- if ("establishmentMeans" %in% names(thinned_sf)) thinned_sf$establishmentMeans else NA_character_
            de <- if ("degreeOfEstablishment" %in% names(thinned_sf)) thinned_sf$degreeOfEstablishment else NA_character_
            os <- if ("occurrenceStatus" %in% names(thinned_sf)) thinned_sf$occurrenceStatus else NA_character_
            dwc_vec <- map_dwc_status(em, de, os)
            fill <- is.na(status_vec); status_vec[fill] <- dwc_vec[fill]
          }

          if (isTRUE(use_worms)) {
            wm <- get_map_try_variants(get_worms_map_for_species)
            if (!is.null(wm) && nrow(wm)) {
              mrgid_for_match <- rep(NA_integer_, nrow(thinned_sf))
              if ("mrgid" %in% names(thinned_sf)) {
                mrgid_for_match <- suppressWarnings(as.integer(thinned_sf$mrgid))
              }
              if ("mrgid_src" %in% names(thinned_sf)) {
                src <- suppressWarnings(as.integer(thinned_sf$mrgid_src))
                mrgid_for_match[is.na(mrgid_for_match) & !is.na(src)] <- src
              }
              if ("locationID" %in% names(thinned_sf)) {
                lid <- tolower(as.character(thinned_sf$locationID))
                loc <- rep(NA_integer_, length(lid))
                sel <- grepl("mrgid", lid, useBytes = TRUE, ignore.case = TRUE)
                if (any(sel, na.rm = TRUE)) {
                  loc[sel] <- suppressWarnings(as.integer(sub(".*?mrgid\\D*(\\d+).*", "\\1", lid[sel])))
                  need_fill <- is.na(mrgid_for_match) & !is.na(loc)
                  mrgid_for_match[need_fill] <- loc[need_fill]
                }
              }
              needs <- is.na(status_vec) | status_vec == "unknown"
              idx   <- match(mrgid_for_match, wm$mrgid)
              hit   <- !is.na(idx) & needs
              status_vec[hit] <- wm$status[idx[hit]]
              .say(paste0(
                species_name, " @ ", region_human,
                ": WoRMS-fill matched ", sum(hit), " / ", length(status_vec),
                " (", sum(!is.na(mrgid_for_match)), " rows had a usable MRGID)"
              ))
            }
          }

          if (isTRUE(use_griis) && (!isTRUE(griis_first)) && (!is.null(iso3_vec) || !is.null(iso2_vec))) {
            gm <- get_map_try_variants(get_griis_map_for_species)
            if (!is.null(gm) && nrow(gm)) {
              needs <- is.na(status_vec) | status_vec == "unknown"
              hit3 <- hit2 <- rep(FALSE, length(status_vec))
              if (any(needs) && !is.null(iso3_vec) && "iso3" %in% names(gm)) {
                idx3 <- match(toupper(iso3_vec), gm$iso3)
                hit3 <- !is.na(idx3) & needs
                status_vec[hit3] <- gm$status[idx[hit3]]
                needs <- is.na(status_vec) | status_vec == "unknown"
              }
              if (any(needs) && !is.null(iso2_vec) && "iso2" %in% names(gm)) {
                idx2 <- match(toupper(iso2_vec), toupper(gm$iso2))
                hit2 <- !is.na(idx2) & needs
                status_vec[hit2] <- gm$status[idx2[hit2]]
              }
              .say(paste0(species_name, " @ ", region_human, ": GRIIS-fill matched ", sum(hit3)+sum(hit2), " / ", length(status_vec)))
            }
          }

          if (isTRUE(use_gisd) && !isTRUE(gisd_first) && !is.null(iso3_for_matching)) {
            gm_gisd <- get_map_try_variants(get_gisd_map_for_species)
            if (!is.null(gm_gisd) && nrow(gm_gisd)) {
              idx <- match(iso3_for_matching, gm_gisd$iso3)
              needs <- is.na(status_vec) | status_vec == "unknown"
              hit <- !is.na(idx) & needs
              status_vec[hit] <- gm_gisd$status[idx[hit]]
              .say(paste0(species_name, " @ ", region_human, ": GISD-fill matched ", sum(hit), " / ", length(status_vec)))
            }
          }

          if (isTRUE(use_gbif_distributions) && !isTRUE(gbif_first) && !is.null(iso3_for_matching)) {
            gm_gbif <- get_map_try_variants(get_gbif_map_for_species)
            if (!is.null(gm_gbif) && nrow(gm_gbif)) {
              idx  <- match(iso3_for_matching, gm_gbif$iso3)
              needs <- is.na(status_vec) | status_vec == "unknown"
              hit  <- !is.na(idx) & needs
              status_vec[hit] <- gm_gbif$status[idx[hit]]
              .say(paste0(species_name, " @ ", region_human, ": GBIF-fill matched ", sum(hit), " / ", length(status_vec)))
            }
          }

          if (status_method %in% c("lookup","dwc_then_lookup")) {
            s_scalar <- lookup_status_scalar(species_name, region_human, raw_sf, status_lookup)
            if (!is.na(s_scalar)) {
              needs <- is.na(status_vec) | status_vec == "unknown"
              status_vec[needs] <- s_scalar
            }
          }

          if (identical(status_method, "custom")) {
            if (is.null(status_fun) || !is.function(status_fun)) {
              stop("`status_method='custom'` requires a valid `status_fun`.")
            }
            s_custom <- status_fun(species_name, region_human, thinned_sf)
            if (length(s_custom) == 1L)      status_vec[] <- as.character(s_custom)
            else if (length(s_custom) == nrow(thinned_sf)) status_vec <- as.character(s_custom)
            else stop("`status_fun` must return a scalar or a vector of length nrow(sf).")
          }

          status_vec <- tolower(trimws(as.character(status_vec)))
          ok <- c("alien","native","unknown")
          status_vec[!(status_vec %in% ok) | is.na(status_vec)] <- "unknown"
          thinned_sf[[status_col]] <- status_vec
        }

        region_safe <- tryCatch(
          make_safe_label(region_human),
          error = function(e) {
            if ("mrgid" %in% names(raw_sf) && length(raw_sf$mrgid) && !is.na(raw_sf$mrgid[1])) paste0(toupper(overlay), "_", as.character(raw_sf$mrgid[1])) else paste0(toupper(overlay), "_UNKNOWN")
          }
        )

        out_path_or_obj <- export_gbif_csv(
          sf_df = thinned_sf, region_label = region_safe, species_name = species_name,
          output_dir = output_dir, thinning = apply_thinning, dist_km = dist_km,
          store_in_memory = store_in_memory, workflow = "eez"  # keep legacy label to avoid breaking downstream code
        )

        region_human_disp <- utf8_safe(region_human)
        summary_tbl <- append_summary_row(
          summary_tbl, species = species_name, region_id = region_human_disp, region_type = toupper(overlay),
          n_total = n_total, n_cleaned = n_total, n_thinned = n_thinned,
          output_file = if (isTRUE(store_in_memory)) NA_character_ else out_path_or_obj,
          status = "success"
        )

        if (isTRUE(return_all_results) && isTRUE(store_in_memory)) {
          key <- paste(species_name, region_safe, sep = "_")
          all_results[[key]] <- thinned_sf
        }
      }
    }
  }

  # ---------------- outputs ----------------

  if (isTRUE(export_summary)) {
    readr::write_csv(summary_tbl, file.path(output_dir, "gbif_summary.csv"))
  }

  if (isTRUE(return_all_results) && isTRUE(store_in_memory)) {
    if (length(all_results) == 0) {
      cli::cli_alert_info("No overlay-joined results to return (all_results is empty).")
      return(tibble::tibble())
    }
    harmonized <- harmonize_column_types(all_results)
    return(dplyr::bind_rows(harmonized, .id = "species_region"))
  } else {
    return(invisible(NULL))
  }
}


