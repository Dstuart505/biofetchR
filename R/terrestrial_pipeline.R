#' Terrestrial GBIF Processing Pipeline
#' (GADM + TEOW/FEOW/LAKES/RIVERS/BASINS/GMBA; auto-cached, 0–360°-safe)
#'
#' Join GBIF occurrences to administrative/ecoregion/hydro/mountain layers:
#' - **GADM** (admin polygons)
#' - **TEOW** (Terrestrial Ecoregions of the World)
#' - **FEOW** (Freshwater Ecoregions of the World)
#' - **LAKES** (HydroLAKES; natural lakes by default)
#' - **RIVERS** (HydroRIVERS; points snapped to nearest reach → `hyriv_id` only)
#' - **BASINS** (HydroBASINS; Pfafstetter levels)
#' - **GMBA** (GMBA Mountain Inventory v2 polygons; IDs & hierarchy)
#'
#' Polygon joins are robust to 0–360° overlays: if an overlay is in 0–360°,
#' points are transparently shifted for the spatial join. Rivers are line
#' features; points are snapped to the nearest reach within a distance
#' threshold (meters) using a memory-safe buffer window and global fallback.
#' If the initial tolerance yields zero matches, the tolerance is relaxed and
#' the snap is retried (also with the uncropped rivers layer).
#'
#' ## GMBA Mountains v2 (no-URL, auto-download)
#' The GMBA Mountain Inventory v2 is fetched automatically from EarthEnv
#' (DOI: 10.48601/earthenv-t9k2-1407) and cached. Defaults: `gmba_layer="all"`,
#' `gmba_extent="standard"`. You can switch to `gmba_layer="basic"` (non-overlapping
#' units) or `gmba_extent="broad"` if desired. No URL input is required.
#'
#' @param df data.frame with columns `species` and `iso2c` (ISO2).
#' @param output_dir Directory for exported CSVs and summary.
#' @param user,pwd,email GBIF credentials.
#' @param region_source Character vector from
#'   `c("gadm","teow","feow","lakes","rivers","basins","gmba")`.
#'   The first entry is used when `overlay_mode="single"`.
#' @param overlay_mode "single"|"dual_separate"|"dual_intersection".
#' @param gadm_unit Integer 0/1/2 (ignored unless "gadm" used).
#' @param cache_dir_gadm Cache dir for GADM.
#' @param teow_cache_dir Cache dir for TEOW.
#' @param teow_method "auto"|"mapme"|"direct".
#' @param teow_url URL when `teow_method="direct"`.
#' @param teow_force_refresh Logical.
#' @param teow_clip_to_countries Logical; clip TEOW to input-country bbox.
#' @param feow_cache_dir Cache dir for FEOW.
#' @param feow_method "auto"|"mapme"|"direct" (mapped internally to bf_load_feow()).
#' @param feow_url URL when `feow_method="direct"`.
#' @param feow_force_refresh Logical.
#' @param feow_clip_to_countries Logical; clip FEOW to input-country bbox.
#' @param lakes_cache_dir Cache dir for HydroLAKES.
#' @param lakes_force_refresh Logical.
#' @param lakes_only Logical; keep natural lakes only (drop reservoirs/controls).
#' @param lakes_min_area_km2 Numeric; drop lakes smaller than this (speed-up).
#' @param rivers_cache_dir Cache dir for HydroRIVERS.
#' @param rivers_cache "disk" or "memory"; passed to [bf_load_rivers()].
#' @param rivers_cache_format "rds","gpkg","both"; on-disk format for rivers cache.
#' @param rivers_force_refresh Logical; force rebuild of rivers cache.
#' @param rivers_regions character subset of HydroSHEDS regions (e.g. `c("eu","na")`), or `NULL` for global.
#' @param rivers_min_strahler integer; keep reaches with Strahler order >= this (optional).
#' @param rivers_min_discharge_cms numeric; keep reaches with mean discharge >= this (optional).
#' @param rivers_max_snap_m numeric; max snap distance (meters) from point to nearest reach.
#' @param basins_cache_dir Cache dir for HydroBASINS.
#' @param basins_level integer 1..12; Pfafstetter detail level.
#' @param basins_with_lakes logical; use HydroBASINS “with lakes” variant.
#' @param basins_clip_to_countries logical; crop basins to input-country bbox.
#' @param gmba_cache_dir Cache dir for GMBA Mountains v2.
#' @param gmba_layer "all"|"basic" (use "all" for full hierarchy; "basic" for non-overlapping units).
#' @param gmba_extent "standard"|"broad" (GMBA definitions; see EarthEnv docs).
#' @param gmba_force_refresh Logical; force (re)download of GMBA archive.
#' @param gmba_clip_to_countries Logical; crop GMBA polygons to input-country bbox.
#' @param batch_size Integer; species per batch.
#' @param dist_km Thinning distance (km).
#' @param apply_thinning Logical; apply thinning before export.
#' @param apply_cleaning Logical; pass-through to your thinning helper.
#' @param return_all_results Logical.
#' @param export_summary Logical; write summary CSV.
#' @param store_in_memory Logical.
#' @param use_planar Logical; if TRUE disables s2 during ops.
#' @param quiet Logical.
#'
#' @return If `return_all_results && store_in_memory`, returns combined `sf`;
#'   else `invisible(NULL)`.
#' @export
process_gbif_terrestrial_pipeline <- function(
    df, output_dir, user, pwd, email,
    region_source          = c("gadm","teow","feow","lakes","rivers","basins","gmba"),
    overlay_mode           = c("single","dual_separate","dual_intersection"),
    gadm_unit              = 1,
    cache_dir_gadm         = "gadm_cache",
    teow_cache_dir         = tools::R_user_dir("biofetchR","data"),
    teow_method            = c("auto","mapme","direct"),
    teow_url               = NULL,
    teow_force_refresh     = FALSE,
    teow_clip_to_countries = TRUE,
    feow_cache_dir         = tools::R_user_dir("biofetchR","data"),
    feow_method            = c("auto","mapme","direct"),
    feow_url               = NULL,
    feow_force_refresh     = FALSE,
    feow_clip_to_countries = TRUE,
    lakes_cache_dir        = tools::R_user_dir("biofetchR","data"),
    lakes_force_refresh    = FALSE,
    lakes_only             = TRUE,
    lakes_min_area_km2     = 1,
    rivers_cache_dir       = tools::R_user_dir("biofetchR","data"),
    rivers_cache           = c("disk","memory"),
    rivers_cache_format    = c("rds","gpkg","both"),
    rivers_force_refresh   = FALSE,
    rivers_regions         = NULL,
    rivers_min_strahler    = NULL,
    rivers_min_discharge_cms = NULL,
    rivers_max_snap_m      = 1000,
    basins_cache_dir       = tools::R_user_dir("biofetchR","data"),
    basins_level           = 12,
    basins_with_lakes      = FALSE,
    basins_clip_to_countries = TRUE,
    gmba_cache_dir         = tools::R_user_dir("biofetchR","data"),
    gmba_layer             = c("all","basic"),
    gmba_extent            = c("standard","broad"),
    gmba_force_refresh     = FALSE,
    gmba_clip_to_countries = TRUE,
    batch_size             = 5,
    dist_km                = 5,
    apply_thinning         = FALSE,
    apply_cleaning         = TRUE,
    return_all_results     = TRUE,
    export_summary         = TRUE,
    store_in_memory        = TRUE,
    use_planar             = TRUE,
    quiet                  = FALSE
) {

  # ---- options & basic checks ---------------------------------------------
  sources      <- unique(match.arg(region_source, c("gadm","teow","feow","lakes","rivers","basins","gmba"), several.ok = TRUE))
  overlay_mode <- match.arg(overlay_mode)
  teow_method  <- match.arg(teow_method)
  feow_method  <- match.arg(feow_method)
  rivers_cache        <- match.arg(rivers_cache)
  rivers_cache_format <- match.arg(rivers_cache_format)
  gmba_layer  <- match.arg(gmba_layer)
  gmba_extent <- match.arg(gmba_extent)

  # FEOW method mapping for your loader API
  feow_method_resolved <- switch(
    feow_method,
    auto   = "auto",
    mapme  = "arcgis",
    direct = "download",
    feow_method
  )

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  old_s2 <- sf::sf_use_s2(); on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  if (use_planar) sf::sf_use_s2(FALSE) else sf::sf_use_s2(TRUE)

  stopifnot(all(c("species","iso2c") %in% names(df)))
  if ("gadm" %in% sources && !gadm_unit %in% c(0,1,2)) stop("`gadm_unit` must be 0,1,2.")
  df$iso2c <- toupper(df$iso2c)

  species_list <- unique(df$species)
  batches      <- split(species_list, ceiling(seq_along(species_list) / batch_size))
  summary_tbl  <- initialize_summary()
  all_results  <- list()
  countries    <- unique(df$iso2c)

  # ---- helpers -------------------------------------------------------------
  .bf_make_valid <- function(x) {
    if ("st_make_valid" %in% getNamespaceExports("sf")) return(suppressWarnings(sf::st_make_valid(x)))
    suppressWarnings(sf::st_buffer(x, 0))
  }
  .geom_name <- function(x) {
    nm <- tryCatch(attr(x, "sf_column"), error = function(e) NULL)
    if (!is.null(nm) && nm %in% names(x)) return(nm)
    cand <- intersect(c("geometry","geom","wkb_geometry"), names(x))
    if (length(cand)) { sf::st_geometry(x) <- cand[1]; return(cand[1]) }
    stop("No geometry column found.")
  }
  .is_0360 <- function(x) {
    bb <- sf::st_bbox(x)
    isTRUE(!is.na(bb["xmin"]) && !is.na(bb["xmax"]) && bb["xmin"] >= 0 && bb["xmax"] > 180)
  }
  .join_overlay_labels <- function(pts_sf, overlay_sf, label_col, tag = "overlay", quiet = FALSE) {
    stopifnot(inherits(pts_sf, "sf"), inherits(overlay_sf, "sf"), label_col %in% names(overlay_sf))
    pts <- if (is.na(sf::st_crs(pts_sf)) || sf::st_crs(pts_sf) != sf::st_crs(4326)) sf::st_transform(pts_sf, 4326) else pts_sf
    ovl <- if (is.na(sf::st_crs(overlay_sf)) || sf::st_crs(overlay_sf) != sf::st_crs(4326)) sf::st_transform(overlay_sf, 4326) else overlay_sf
    if (!.is_0360(ovl)) {
      jt <- suppressWarnings(sf::st_join(pts, ovl, join = sf::st_intersects, left = TRUE))
      if (!quiet) message(tag, " join (−180…180): ", sum(!is.na(jt[[label_col]])), " / ", nrow(jt), " labeled")
      return(jt[[label_col]])
    }
    coords <- sf::st_coordinates(pts)
    if (!nrow(coords)) return(rep(NA_character_, nrow(pts)))
    lon360 <- ifelse(coords[,1] < 0, coords[,1] + 360, coords[,1])
    lat    <- coords[,2]
    tmp <- cbind(sf::st_drop_geometry(pts), ..lon360 = lon360, ..lat = lat)
    tmp <- sf::st_as_sf(tmp, coords = c("..lon360","..lat"), crs = 4326, remove = TRUE)
    jt  <- suppressWarnings(sf::st_join(tmp, ovl, join = sf::st_intersects, left = TRUE))
    if (!quiet) message(tag, " join (0…360 with point shift): ", sum(!is.na(jt[[label_col]])), " / ", nrow(jt), " labeled")
    jt[[label_col]]
  }

  # Snap points to HydroRIVERS (ID-only)
  .join_to_rivers <- function(pts_sf, rivers_sf, snap_m = 1000, quiet = FALSE) {
    stopifnot(inherits(pts_sf, "sf"), inherits(rivers_sf, "sf"))
    pts4326 <- if (is.na(sf::st_crs(pts_sf))) sf::st_set_crs(pts_sf, 4326) else sf::st_transform(pts_sf, 4326)
    riv4326 <- if (is.na(sf::st_crs(rivers_sf))) sf::st_set_crs(rivers_sf, 4326) else sf::st_transform(rivers_sf, 4326)
    pts_m <- tryCatch(sf::st_transform(pts4326, 3857), error = function(e) pts4326)
    riv_m <- tryCatch(sf::st_transform(riv4326, 3857),  error = function(e) riv4326)

    use_global <- FALSE
    cand_idx <- tryCatch({
      buf_union <- suppressWarnings(sf::st_union(sf::st_buffer(sf::st_geometry(pts_m), dist = snap_m)))
      hits_log  <- lengths(sf::st_intersects(riv_m, buf_union)) > 0
      which(hits_log)
    }, error = function(e) integer(0))
    if (!length(cand_idx)) use_global <- TRUE

    if (!use_global) {
      riv_cand_4326 <- riv4326[cand_idx, , drop = FALSE]
      idx <- sf::st_nearest_feature(pts4326, riv_cand_4326)
      riv_cand_m <- tryCatch(sf::st_transform(riv_cand_4326[idx, ], 3857), error = function(e) riv_cand_4326[idx, ])
      d  <- suppressWarnings(sf::st_distance(sf::st_geometry(pts_m), sf::st_geometry(riv_cand_m), by_element = TRUE))
      id <- riv_cand_4326$hyriv_id[idx]
      id[as.numeric(d) > snap_m] <- NA_integer_
    } else {
      idx <- sf::st_nearest_feature(pts4326, riv4326)
      riv_m_idx <- tryCatch(sf::st_transform(riv4326[idx, ], 3857), error = function(e) riv4326[idx, ])
      d  <- suppressWarnings(sf::st_distance(sf::st_geometry(pts_m), sf::st_geometry(riv_m_idx), by_element = TRUE))
      id <- riv4326$hyriv_id[idx]
      id[as.numeric(d) > snap_m] <- NA_integer_
    }

    if (!quiet) {
      msg <- if (use_global) "global-nearest" else "buffer-windowed"
      message("RIVERS snap (", msg, "): ", sum(!is.na(id)), " / ", length(id), " within ", snap_m, " m")
    }
    id
  }

  .crop_overlay_to_countries <- function(overlay, iso2, quiet = FALSE, tag = "overlay") {
    if (!length(iso2)) return(overlay)
    bb_sfc_4326 <- tryCatch({
      world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")
      sel   <- world[world$iso_a2 %in% toupper(iso2), , drop = FALSE]
      if (!nrow(sel)) sf::st_as_sfc(sf::st_bbox(c(xmin=-180,ymin=-90,xmax=180,ymax=90), crs = sf::st_crs(4326)))
      else sf::st_as_sfc(sf::st_bbox(sf::st_transform(sel, 4326)))
    }, error = function(e) sf::st_as_sfc(sf::st_bbox(c(xmin=-180,ymin=-90,xmax=180,ymax=90), crs = sf::st_crs(4326))))
    bb_for_crop <- bb_sfc_4326
    if (!is.na(sf::st_crs(overlay)) && sf::st_crs(overlay) != sf::st_crs(4326)) {
      bb_for_crop <- tryCatch(sf::st_transform(bb_sfc_4326, sf::st_crs(overlay)), error = function(e) bb_sfc_4326)
    }
    if (!quiet) {
      bbx <- sf::st_bbox(bb_for_crop)
      message(tag, " cropping to bbox [", sprintf("%.2f, %.2f, %.2f, %.2f", bbx["xmin"],bbx["ymin"],bbx["xmax"],bbx["ymax"]), "]")
    }
    pre_n <- nrow(overlay)
    cropped <- suppressWarnings(sf::st_crop(overlay, sf::st_bbox(bb_for_crop)))
    if (!quiet) message(tag, " after crop: ", nrow(cropped), " / ", pre_n, " features.")
    if (nrow(cropped) > 0) cropped else {
      if (!quiet) message("Crop produced 0 features; using uncropped ", tag, " for joins.")
      overlay
    }
  }

  # --- GMBA loader (auto-download; no URL needed) --------------------------
  .load_gmba_inventory <- function(cache_dir = gmba_cache_dir,
                                   layer = gmba_layer,         # "all"|"basic"
                                   extent = gmba_extent,       # "standard"|"broad"
                                   force_refresh = gmba_force_refresh,
                                   quiet = TRUE) {
    stopifnot(layer %in% c("all","basic"), extent %in% c("standard","broad"))
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    # Canonical EarthEnv endpoints (DOI: 10.48601/earthenv-t9k2-1407)
    # Standard filenames used on the official downloads page.
    base <- "https://data.earthenv.org/mountains"
    fname <- switch(paste(layer, extent, sep = "_"),
                    "all_standard"   = "GMBA_Inventory_v2.0_standard.zip",
                    "all_broad"      = "GMBA_Inventory_v2.0_broad.zip",
                    "basic_standard" = "GMBA_Basic_v2.0_standard.zip",
                    "basic_broad"    = "GMBA_Basic_v2.0_broad.zip"
    )
    # Path layout is /{extent}/<zip>
    url   <- file.path(base, extent, fname)
    tgt   <- file.path(cache_dir, fname)
    unz   <- file.path(cache_dir, sub("\\.zip$", "", fname))

    need_dl <- !file.exists(tgt) || force_refresh
    if (need_dl) {
      if (!quiet) message("GMBA: downloading ", url)
      utils::download.file(url, tgt, mode = "wb", quiet = quiet)
    }
    if (!dir.exists(unz) || force_refresh) {
      if (!quiet) message("GMBA: unpacking → ", unz)
      utils::unzip(tgt, exdir = unz)
    }
    # Find a shapefile (there is typically a single .shp)
    shp <- list.files(unz, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
    if (!length(shp)) stop("GMBA: no .shp found after unzip (archive layout may have changed).")
    gmba <- sf::st_read(shp[1], quiet = quiet)

    # Normalize ID column → gmba_id (robust to field name variants)
    id_cand <- intersect(tolower(names(gmba)),
                         tolower(c("GMBA_ID","GMBA_V2_ID","ID","RANGE_ID","ID_","id_gmba","id")))
    if (length(id_cand)) {
      orig <- names(gmba)[match(id_cand[1], tolower(names(gmba)))]
      gmba$gmba_id <- as.character(gmba[[orig]])
    } else {
      gmba$gmba_id <- as.character(seq_len(nrow(gmba)))
    }

    gcol <- .geom_name(gmba)
    gmba <- gmba[, c("gmba_id", gcol), drop = FALSE]
    .bf_make_valid(gmba)
  }

  # Label columns lookup
  label_cols <- list(
    gadm   = "geoname_gadm",
    teow   = "geoname_teow",
    feow   = "geoname_feow",
    lakes  = "geoname_lakes",
    rivers = "hyriv_id",
    basins = "hybas_id",
    gmba   = "gmba_id"
  )

  # ---- load overlays (once) -----------------------------------------------
  overlay_gadm <- overlay_teow <- overlay_feow <- overlay_lakes <- overlay_rivers <- overlay_rivers_full <- overlay_basins <- overlay_gmba <- NULL

  if ("gadm" %in% sources) {
    overlay_gadm <- load_gadm(iso2c = countries, level = gadm_unit, cache_dir = cache_dir_gadm, quiet = TRUE)
    nm <- intersect(c("NAME_0","NAME_1","NAME_2","name","NAME","GID_0","GID_1","GID_2"), names(overlay_gadm))
    overlay_gadm$geoname_gadm <- if (length(nm)) as.character(overlay_gadm[[nm[1]]]) else paste0("GADM_", gadm_unit, "_", seq_len(nrow(overlay_gadm)))
    gcol <- .geom_name(overlay_gadm)
    overlay_gadm <- overlay_gadm[, c("geoname_gadm", gcol), drop = FALSE]
    overlay_gadm <- .bf_make_valid(overlay_gadm)
  }

  if ("teow" %in% sources) {
    overlay_teow <- bf_load_teow(cache_dir = teow_cache_dir, method = teow_method,
                                 source_url = teow_url, force_refresh = teow_force_refresh, quiet = TRUE)
    if (!"geoname" %in% names(overlay_teow)) {
      nm <- intersect(c("ECO_NAME","eco_name","ECOREGION","ecoregion"), names(overlay_teow))
      overlay_teow$geoname <- if (length(nm)) as.character(overlay_teow[[nm[1]]]) else paste0("ECO_", seq_len(nrow(overlay_teow)))
    }
    names(overlay_teow)[names(overlay_teow) == "geoname"] <- "geoname_teow"
    gcol <- .geom_name(overlay_teow)
    overlay_teow <- overlay_teow[, c("geoname_teow", gcol), drop = FALSE]
    overlay_teow <- .bf_make_valid(overlay_teow)
    if (isTRUE(teow_clip_to_countries)) overlay_teow <- .crop_overlay_to_countries(overlay_teow, countries, quiet, tag = "TEOW")
  }

  if ("feow" %in% sources) {
    overlay_feow <- bf_load_feow(cache_dir = feow_cache_dir, method = feow_method_resolved,
                                 source_url = feow_url, force_refresh = feow_force_refresh, quiet = TRUE)
    if (!"geoname" %in% names(overlay_feow)) {
      nm <- intersect(c("ecoregion","ECOREGION","ECO_NAME","FEOW_NAME","NAME"), names(overlay_feow))
      overlay_feow$geoname <- if (length(nm)) as.character(overlay_feow[[nm[1]]]) else paste0("FEOW_", seq_len(nrow(overlay_feow)))
    }
    names(overlay_feow)[names(overlay_feow) == "geoname"] <- "geoname_feow"
    gcol <- .geom_name(overlay_feow)
    overlay_feow <- overlay_feow[, c("geoname_feow", gcol), drop = FALSE]
    overlay_feow <- .bf_make_valid(overlay_feow)
    if (isTRUE(feow_clip_to_countries)) overlay_feow <- .crop_overlay_to_countries(overlay_feow, countries, quiet, tag = "FEOW")
  }

  if ("lakes" %in% sources) {
    overlay_lakes <- bf_load_lakes(
      cache_dir     = lakes_cache_dir,
      force_refresh = lakes_force_refresh,
      lakes_only    = lakes_only,
      min_area_km2  = lakes_min_area_km2,
      quiet         = TRUE
    )
    nm <- intersect(c("name","Lake_name","lake_name"), names(overlay_lakes))
    if (!length(nm)) {
      if ("lake_id" %in% names(overlay_lakes)) {
        overlay_lakes$geoname_lakes <- paste0("LAKE_", overlay_lakes$lake_id)
      } else {
        overlay_lakes$geoname_lakes <- paste0("LAKE_", seq_len(nrow(overlay_lakes)))
      }
    } else {
      overlay_lakes$geoname_lakes <- as.character(overlay_lakes[[nm[1]]])
    }
    gcol <- .geom_name(overlay_lakes)
    overlay_lakes <- overlay_lakes[, c("geoname_lakes", gcol), drop = FALSE]
    overlay_lakes <- .bf_make_valid(overlay_lakes)
    overlay_lakes <- .crop_overlay_to_countries(overlay_lakes, countries, quiet, tag = "LAKES")
  }

  if ("rivers" %in% sources) {
    overlay_rivers <- bf_load_rivers(
      cache_dir           = rivers_cache_dir,
      cache               = rivers_cache,
      cache_format        = rivers_cache_format,
      force_refresh       = rivers_force_refresh,
      regions             = rivers_regions,
      min_strahler        = rivers_min_strahler,
      min_discharge_cms   = rivers_min_discharge_cms,
      quiet               = TRUE
    )
    if (!"hyriv_id" %in% names(overlay_rivers) && "HYRIV_ID" %in% names(overlay_rivers)) {
      overlay_rivers$hyriv_id <- overlay_rivers$HYRIV_ID
    }
    gcol <- .geom_name(overlay_rivers)
    overlay_rivers <- overlay_rivers[, c("hyriv_id", gcol), drop = FALSE]
    overlay_rivers <- .bf_make_valid(overlay_rivers)
    overlay_rivers_full <- overlay_rivers
    overlay_rivers      <- .crop_overlay_to_countries(overlay_rivers, countries, quiet, tag = "RIVERS")
  }

  if ("basins" %in% sources) {
    overlay_basins <- bf_load_basins(
      level        = basins_level,
      with_lakes   = basins_with_lakes,
      cache_dir    = basins_cache_dir,
      quiet        = TRUE
    )
    if (!"hybas_id" %in% names(overlay_basins) && "HYBAS_ID" %in% names(overlay_basins)) {
      overlay_basins$hybas_id <- overlay_basins$HYBAS_ID
    }
    gcol <- .geom_name(overlay_basins)
    overlay_basins <- overlay_basins[, c("hybas_id", gcol), drop = FALSE]
    overlay_basins <- .bf_make_valid(overlay_basins)
    if (isTRUE(basins_clip_to_countries)) overlay_basins <- .crop_overlay_to_countries(overlay_basins, countries, quiet, tag = "BASINS")
  }

  if ("gmba" %in% sources) {
    overlay_gmba <- .load_gmba_inventory(cache_dir = gmba_cache_dir,
                                         layer = gmba_layer,
                                         extent = gmba_extent,
                                         force_refresh = gmba_force_refresh,
                                         quiet = TRUE)
    if (isTRUE(gmba_clip_to_countries)) {
      overlay_gmba <- .crop_overlay_to_countries(overlay_gmba, countries, quiet, tag = "GMBA")
    }
  }

  # ---- main loop -----------------------------------------------------------
  for (batch in batches) {
    for (species_name in batch) {

      iso2_codes <- df |>
        dplyr::filter(species == species_name) |>
        dplyr::pull(iso2c) |>
        unique()

      valid_countries <- iso2_codes[unlist(lapply(iso2_codes, check_gbif_presence, species = species_name))]
      if (!length(valid_countries)) next

      download_keys <- download_gbif_batch_gadm(
        species    = species_name,
        iso2_codes = valid_countries,
        user = user, pwd = pwd, email = email
      )
      if (is.null(download_keys)) next

      results_list <- wait_and_import_gbif(download_keys)
      if (!length(results_list)) next

      for (cc in names(results_list)) {
        raw_sf <- results_list[[cc]]

        if (is.null(raw_sf) || !nrow(raw_sf)) {
          summary_tbl <- append_summary_row(summary_tbl, species_name, cc,
                                            region_type = paste(sources, collapse = "+"),
                                            n_total = 0, n_cleaned = 0, n_thinned = 0,
                                            output_file = NA, status = "no_data")
          next
        }

        # ensure sf in WGS84
        if (!inherits(raw_sf, "sf") &&
            all(c("decimalLongitude","decimalLatitude") %in% names(raw_sf))) {
          raw_sf <- sf::st_as_sf(raw_sf, coords = c("decimalLongitude","decimalLatitude"),
                                 crs = 4326, remove = FALSE)
        } else if (!is.na(sf::st_crs(raw_sf)) && sf::st_crs(raw_sf)$epsg != 4326) {
          raw_sf <- sf::st_transform(raw_sf, 4326)
        }

        # spatial joins
        joined <- raw_sf
        if (!is.null(overlay_gadm)) {
          jg <- suppressWarnings(sf::st_join(raw_sf, overlay_gadm, join = sf::st_within, left = TRUE))
          stopifnot(nrow(jg) == nrow(raw_sf))
          joined$geoname_gadm <- jg$geoname_gadm
        }
        if (!is.null(overlay_teow)) {
          joined$geoname_teow <- .join_overlay_labels(raw_sf, overlay_teow, "geoname_teow", tag = "TEOW", quiet = quiet)
        }
        if (!is.null(overlay_feow)) {
          joined$geoname_feow <- .join_overlay_labels(raw_sf, overlay_feow, "geoname_feow", tag = "FEOW", quiet = quiet)
        }
        if (!is.null(overlay_lakes)) {
          jl <- suppressWarnings(sf::st_join(raw_sf, overlay_lakes, join = sf::st_within, left = TRUE))
          stopifnot(nrow(jl) == nrow(raw_sf))
          joined$geoname_lakes <- jl$geoname_lakes
        }
        if (!is.null(overlay_rivers)) {
          # primary snap (cropped, requested tolerance)
          joined$hyriv_id <- .join_to_rivers(raw_sf, overlay_rivers, snap_m = rivers_max_snap_m, quiet = quiet)

          # if zero matches, relax tolerance (up to 5 km) using cropped
          if (all(is.na(joined$hyriv_id))) {
            relax <- max(2000, min(5000, rivers_max_snap_m * 5))
            if (!quiet) message("biofetchR: 0 snaps at ", rivers_max_snap_m, " m → trying ", relax, " m (cropped)")
            id_relax <- .join_to_rivers(raw_sf, overlay_rivers, snap_m = relax, quiet = quiet)
            joined$hyriv_id <- id_relax
          }

          # if still zero, retry with uncropped rivers (relaxed tolerance)
          if (all(is.na(joined$hyriv_id)) && !is.null(overlay_rivers_full) && nrow(overlay_rivers_full) > 0) {
            relax <- max(2000, min(5000, rivers_max_snap_m * 5))
            if (!quiet) message("biofetchR: still 0 snaps → retrying with uncropped rivers at ", relax, " m")
            id_unc <- .join_to_rivers(raw_sf, overlay_rivers_full, snap_m = relax, quiet = quiet)
            joined$hyriv_id <- id_unc
          }
        }
        if (!is.null(overlay_basins)) {
          jb <- suppressWarnings(sf::st_join(raw_sf, overlay_basins, join = sf::st_within, left = TRUE))
          stopifnot(nrow(jb) == nrow(raw_sf))
          joined$hybas_id <- jb$hybas_id
        }
        if (!is.null(overlay_gmba)) {
          joined$gmba_id <- .join_overlay_labels(raw_sf, overlay_gmba, "gmba_id", tag = "GMBA", quiet = quiet)
        }

        export_one_group <- function(sf_group, region_id, region_type, workflow_tag) {
          n_total <- nrow(sf_group)
          thinned_sf <- thin_spatial_points(sf_group, dist_km = dist_km,
                                            filter_uncertain = apply_cleaning, quiet = TRUE)
          n_cleaned <- if (apply_cleaning) n_total else NA
          n_thinned <- nrow(thinned_sf)
          out_data <- export_gbif_csv(
            sf_df = thinned_sf,
            region_label = region_id,
            species_name = species_name,
            output_dir = output_dir,
            thinning = apply_thinning,
            dist_km = dist_km,
            store_in_memory = store_in_memory,
            workflow = workflow_tag
          )
          summary_tbl <<- append_summary_row(
            summary_tbl, species = species_name, region_id = region_id,
            region_type = region_type, n_total = n_total, n_cleaned = n_cleaned,
            n_thinned = n_thinned, output_file = if (store_in_memory) NA else out_data,
            status = "success"
          )
          if (return_all_results && store_in_memory) {
            all_results[[paste(species_name, region_id, workflow_tag, sep = "_")]] <<- out_data
          }
        }

        # export logic
        if (overlay_mode == "single" || length(sources) == 1) {
          src <- sources[1]
          if      (src == "gadm")   { x <- joined[!is.na(joined$geoname_gadm),  , drop = FALSE]; rt <- paste0("GADM_level_", gadm_unit); wf <- "gadm";   col <- "geoname_gadm"
          } else if (src == "teow") { x <- joined[!is.na(joined$geoname_teow),  , drop = FALSE]; rt <- "TEOW";                      wf <- "teow";   col <- "geoname_teow"
          } else if (src == "feow") { x <- joined[!is.na(joined$geoname_feow),  , drop = FALSE]; rt <- "FEOW";                      wf <- "feow";   col <- "geoname_feow"
          } else if (src == "lakes"){ x <- joined[!is.na(joined$geoname_lakes), , drop = FALSE]; rt <- "LAKES";                     wf <- "lakes";  col <- "geoname_lakes"
          } else if (src == "rivers"){x <- joined[!is.na(joined$hyriv_id),      , drop = FALSE]; rt <- "RIVERS";                    wf <- "rivers"; col <- "hyriv_id"
          } else if (src == "gmba") { x <- joined[!is.na(joined$gmba_id),       , drop = FALSE]; rt <- paste0("GMBA_", gmba_layer, "_", gmba_extent); wf <- "gmba";   col <- "gmba_id"
          } else {                     x <- joined[!is.na(joined$hybas_id),      , drop = FALSE]; rt <- paste0("BASINS_L", basins_level); wf <- "basins"; col <- "hybas_id" }
          if (nrow(x)) for (nm in names(split(x, x[[col]])))
            export_one_group(split(x, x[[col]])[[nm]], nm, rt, wf)

        } else if (overlay_mode == "dual_separate") {
          for (src in intersect(c("gadm","teow","feow","lakes","rivers","basins","gmba"), sources)) {
            col <- label_cols[[src]]; if (is.null(col) || !col %in% names(joined)) next
            x <- joined[!is.na(joined[[col]]), , drop = FALSE]; if (!nrow(x)) next
            rt <- switch(src,
                         gadm   = paste0("GADM_level_", gadm_unit),
                         teow   = "TEOW",
                         feow   = "FEOW",
                         lakes  = "LAKES",
                         rivers = "RIVERS",
                         basins = paste0("BASINS_L", basins_level),
                         gmba   = paste0("GMBA_", gmba_layer, "_", gmba_extent))
            for (nm in names(split(x, x[[col]])))
              export_one_group(split(x, x[[col]])[[nm]], nm, rt, src)
          }

        } else if (overlay_mode == "dual_intersection") {
          if (length(sources) != 2) {
            message("overlay_mode='dual_intersection' requires exactly two sources; falling back to dual_separate.")
            for (src in intersect(c("gadm","teow","feow","lakes","rivers","basins","gmba"), sources)) {
              col <- label_cols[[src]]; if (is.null(col) || !col %in% names(joined)) next
              x <- joined[!is.na(joined[[col]]), , drop = FALSE]; if (!nrow(x)) next
              rt <- switch(src,
                           gadm   = paste0("GADM_level_", gadm_unit),
                           teow   = "TEOW",
                           feow   = "FEOW",
                           lakes  = "LAKES",
                           rivers = "RIVERS",
                           basins = paste0("BASINS_L", basins_level),
                           gmba   = paste0("GMBA_", gmba_layer, "_", gmba_extent))
              for (nm in names(split(x, x[[col]])))
                export_one_group(split(x, x[[col]])[[nm]], nm, rt, src)
            }
          } else {
            s1 <- sources[1]; s2 <- sources[2]
            c1 <- label_cols[[s1]]; c2 <- label_cols[[s2]]
            stopifnot(!is.null(c1), !is.null(c2))
            xi <- joined[!is.na(joined[[c1]]) & !is.na(joined[[c2]]), , drop = FALSE]
            if (nrow(xi)) {
              xi$geoname_combo <- paste0(xi[[c1]], " :: ", xi[[c2]])
              pair_label <- paste0(toupper(s1), "×", toupper(s2))
              for (nm in names(split(xi, xi$geoname_combo)))
                export_one_group(split(xi, xi$geoname_combo)[[nm]], nm, pair_label, paste0(s1, "_", s2))
            }
          }
        }
      }
    }
  }

  if (export_summary) readr::write_csv(summary_tbl, file.path(output_dir, "gbif_summary.csv"))

  if (return_all_results && store_in_memory) {
    results_list <- harmonize_column_types(all_results)
    return(dplyr::bind_rows(results_list))
  } else {
    invisible(NULL)
  }
}

