#' Thin spatial points by distance, density, and coordinate uncertainty
#'
#' Reduces spatial clustering by enforcing a minimum distance between points or a density threshold,
#' while optionally filtering or flagging by coordinate uncertainty (from GBIF).
#'
#' @param sf_obj A spatial `sf` object with point geometry and GBIF columns.
#' @param dist_km Minimum allowed distance between retained points (in kilometers). Default is 5.
#' @param priority_col Column to prioritize point retention (e.g., "individualCount"). Default is `"individualCount"`.
#' @param return_all Logical. If `TRUE`, return all rows with a `.thinned` column. Default is `FALSE`.
#' @param plot Logical. If `TRUE`, generates a before/after thinning map. Default is `FALSE`.
#' @param max_density Optional numeric. If set, limits retained points to this many per km².
#' @param max_uncertainty_m Numeric. Remove points with `coordinateUncertaintyInMeters` > this value. Default is 10000.
#' @param warn_uncertainty_m Numeric. Flag retained points with uncertainty > this value. Default is 5000.
#' @param filter_uncertain Logical. If `TRUE`, apply uncertainty filtering. Default is `TRUE`.
#' @param quiet Logical. If `TRUE`, suppress messages. Default is `FALSE`.
#'
#' @return A thinned `sf` object or full object with `.thinned` (and optionally `.high_uncertainty`) columns.
#' @export
thin_spatial_points <- function(sf_obj,
                                dist_km = 5,
                                priority_col = "individualCount",
                                return_all = FALSE,
                                plot = FALSE,
                                max_density = NULL,
                                max_uncertainty_m = 10000,
                                warn_uncertainty_m = 5000,
                                filter_uncertain = TRUE,
                                quiet = FALSE) {
  if (!requireNamespace("geosphere", quietly = TRUE)) {
    cli::cli_abort("❌ The {.pkg geosphere} package is required.")
  }

  if (nrow(sf_obj) == 0) {
    if (!quiet) cli::cli_alert_warning("No input points provided; returning empty result.")
    return(sf_obj)
  }

  coords <- sf::st_coordinates(sf_obj)
  n_pts <- nrow(coords)

  if (n_pts < 2) {
    if (!quiet) cli::cli_alert_info("Only one point provided; skipping thinning.")
    if (return_all) {
      sf_obj$.thinned <- TRUE
      if ("coordinateUncertaintyInMeters" %in% names(sf_obj)) {
        sf_obj$.high_uncertainty <- FALSE
      }
      return(sf_obj)
    } else {
      return(sf_obj)
    }
  }

  # ── Step 1: Remove high-uncertainty points ──────────────────────────────
  if (filter_uncertain && "coordinateUncertaintyInMeters" %in% names(sf_obj)) {
    unc <- sf_obj$coordinateUncertaintyInMeters
    keep_uncertainty <- is.na(unc) | unc <= max_uncertainty_m
    n_removed_unc <- sum(!keep_uncertainty, na.rm = TRUE)
    if (!quiet && n_removed_unc > 0) {
      cli::cli_alert_info("Filtered {.val {n_removed_unc}} records with uncertainty > {max_uncertainty_m} m.")
    }
    sf_obj <- sf_obj[keep_uncertainty, ]
    coords <- sf::st_coordinates(sf_obj)  # update after subsetting
  }

  # ── Step 2: Sort by priority column ─────────────────────────────────────
  count_vals <- if (!is.null(priority_col) && priority_col %in% names(sf_obj)) {
    val <- sf_obj[[priority_col]]
    ifelse(is.na(val), -Inf, val)
  } else {
    rep(0, nrow(sf_obj))
  }

  # Penalize uncertainty (if present)
  if ("coordinateUncertaintyInMeters" %in% names(sf_obj)) {
    unc_vals <- sf_obj$coordinateUncertaintyInMeters
    penalty <- ifelse(is.na(unc_vals), 0, unc_vals / 1000)
    count_vals <- count_vals - penalty
  }

  order_priority <- order(-count_vals)
  coords <- coords[order_priority, , drop = FALSE]
  sf_obj <- sf_obj[order_priority, ]

  # ── Step 3: Optional density-based thinning ─────────────────────────────
  if (!is.null(max_density)) {
    bbox <- sf::st_bbox(sf_obj)
    area_km2 <- geosphere::areaPolygon(matrix(c(
      bbox["xmin"], bbox["ymin"],
      bbox["xmin"], bbox["ymax"],
      bbox["xmax"], bbox["ymax"],
      bbox["xmax"], bbox["ymin"],
      bbox["xmin"], bbox["ymin"]
    ), ncol = 2, byrow = TRUE)) / 1e6

    target_n <- ceiling(area_km2 * max_density)
    if (target_n < nrow(sf_obj)) {
      if (!quiet) {
        cli::cli_alert_info("Thinning to ~{target_n} points based on {max_density} pts/km² over ~{round(area_km2)} km².")
      }
      keep_idx <- seq_len(target_n)
      if (return_all) {
        sf_obj$.thinned <- FALSE
        sf_obj$.thinned[keep_idx] <- TRUE
        return(sf_obj)
      } else {
        return(sf_obj[keep_idx, ])
      }
    } else {
      if (!quiet) {
        cli::cli_alert_info("Max density not exceeded; returning all {nrow(sf_obj)} points.")
      }
      if (return_all) sf_obj$.thinned <- TRUE
      return(sf_obj)
    }
  }

  # ── Step 4: Distance-based thinning ─────────────────────────────────────
  keep <- logical(nrow(coords))
  remaining <- seq_len(nrow(coords))

  while (length(remaining) > 0) {
    idx <- remaining[1]
    keep[idx] <- TRUE
    dists <- geosphere::distHaversine(coords[idx, , drop = FALSE], coords[remaining, , drop = FALSE]) / 1000
    remaining <- remaining[dists > dist_km]
  }

  # ── Step 5: Output format and summary ───────────────────────────────────
  n_retained <- sum(keep)
  n_removed <- length(keep) - n_retained

  if (!quiet) {
    cli::cli_alert_info("Distance-based thinning retained {.val {n_retained}} of {.val {length(keep)}} points.")
    if ("coordinateUncertaintyInMeters" %in% names(sf_obj)) {
      unc_vals <- sf_obj$coordinateUncertaintyInMeters
      high_unc <- !is.na(unc_vals) & unc_vals > warn_uncertainty_m
      n_flagged <- sum(high_unc & keep, na.rm = TRUE)
      if (n_flagged > 0) {
        cli::cli_alert_warning("{.val {n_flagged}} retained points have coordinate uncertainty > {warn_uncertainty_m} m.")
      }
    }
  }

  out <- if (return_all) {
    sf_obj$.thinned <- keep
    out <- sf_obj
  } else {
    out <- sf_obj[keep, ]
  }

  # Optional flag for retained uncertainty
  if ("coordinateUncertaintyInMeters" %in% names(out) && warn_uncertainty_m > 0) {
    unc_vals <- out$coordinateUncertaintyInMeters
    out$.high_uncertainty <- !is.na(unc_vals) & unc_vals > warn_uncertainty_m
  }

  # ── Step 6: Optional visualization ──────────────────────────────────────
  if (plot && requireNamespace("ggplot2", quietly = TRUE)) {
    coords_df <- as.data.frame(sf::st_coordinates(sf_obj))
    coords_df$.thinned <- if (return_all) sf_obj$.thinned else keep
    coords_df$category <- ifelse(coords_df$.thinned, "Retained", "Removed")

    p <- ggplot2::ggplot(coords_df, ggplot2::aes(X, Y, color = category)) +
      ggplot2::geom_point(alpha = 0.7, size = 1.5) +
      ggplot2::coord_fixed() +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Spatial thinning results", x = "Longitude", y = "Latitude")
    print(p)
  }

  return(out)
}
