# ---- summary helpers ---------------------------------------------------------
initialize_summary <- function() {
  tibble::tibble(
    species      = character(),
    region_id    = character(),
    region_type  = character(),
    n_total      = integer(),
    n_cleaned    = integer(),
    n_thinned    = integer(),
    output_file  = character(),
    status       = character()
  )
}

append_summary_row <- function(summary_tbl,
                               species,
                               region_id,
                               region_type,
                               n_total,
                               n_cleaned,
                               n_thinned,
                               output_file,
                               status) {
  dplyr::bind_rows(
    summary_tbl,
    tibble::tibble(
      species     = as.character(species),
      region_id   = if (is.null(region_id)) NA_character_ else as.character(region_id),
      region_type = as.character(region_type),
      n_total     = as.integer(n_total),
      n_cleaned   = as.integer(n_cleaned),
      n_thinned   = as.integer(n_thinned),
      output_file = if (is.null(output_file)) NA_character_ else as.character(output_file),
      status      = as.character(status)
    )
  )
}

# Make all list elements share the same columns (lightweight harmonizer)
harmonize_column_types <- function(lst) {
  if (!length(lst)) return(lst)
  all_names <- Reduce(union, lapply(lst, names))
  lapply(lst, function(x) {
    # add missing columns as NA (character NA by default; that’s fine for tibble bind)
    miss <- setdiff(all_names, names(x))
    for (nm in miss) x[[nm]] <- NA
    # drop sf geometry to ensure smooth bind_rows later
    if (inherits(x, "sf")) x <- sf::st_drop_geometry(x)
    x
  })
}

# ---- exporter (writes CSV when not storing in memory) ------------------------
export_gbif_csv <- function(sf_df,
                            region_label,
                            species_name,
                            output_dir,
                            thinning,
                            dist_km,
                            store_in_memory,
                            workflow = "eez") {
  # local safe label
  make_safe_label_local <- function(x) {
    x <- as.character(x %||% "")
    y <- suppressWarnings(iconv(x, from = "", to = "UTF-8", sub = ""))
    y[is.na(y)] <- x[is.na(y)]
    y2 <- suppressWarnings(iconv(y, from = "UTF-8", to = "ASCII//TRANSLIT", sub = ""))
    y2[is.na(y2)] <- y[is.na(y2)]
    y2 <- gsub("[^[:alnum:]]+", "_", y2, perl = TRUE)
    y2 <- gsub("_+", "_", y2, perl = TRUE)
    y2 <- gsub("^_+|_+$", "", y2, perl = TRUE)
    if (!nzchar(y2)) y2 <- "UNK"
    substr(y2, 1, 150)
  }
  `%||%` <- function(a,b) if (is.null(a)) b else a
  
  if (isTRUE(store_in_memory)) {
    # caller ignores returned path when storing in memory
    return(sf_df)
  }
  
  sp_safe  <- make_safe_label_local(species_name)
  eez_safe <- make_safe_label_local(region_label)
  subdir   <- file.path(output_dir, sp_safe)
  if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  
  thin_tag <- if (isTRUE(thinning)) paste0(dist_km, "km") else "noThin"
  file_nm  <- sprintf("%s__%s__%s__%s.csv", sp_safe, eez_safe, thin_tag, workflow)
  path     <- file.path(subdir, file_nm)
  
  df_out <- if (inherits(sf_df, "sf")) sf::st_drop_geometry(sf_df) else sf_df
  readr::write_csv(df_out, path)
  path
}
