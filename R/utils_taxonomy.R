#' Check for NCBI Entrez API Key
#'
#' Ensures the user has set an Entrez API key for NCBI-based taxonomic functions.
#' Silently prevents repeated taxize warnings by setting an empty key if not present.
#'
#' @return NULL. Alerts the user if no key is found.
#' @export
check_entrez_key <- function() {
  key <- Sys.getenv("ENTREZ_KEY")
  if (!nzchar(key)) {
    cli::cli_alert_warning("⚠️  No ENTREZ API key found. Some taxonomic functions may fail or be rate-limited.")
    cli::cli_alert_info("👉  Get one at https://www.ncbi.nlm.nih.gov/account and run:")
    cli::cli_code('taxize::use_entrez("YOUR_API_KEY")')
    Sys.setenv(ENTREZ_KEY = "")  # Prevent repeated warnings
  } else {
    cli::cli_alert_success("✅  Entrez API key detected.")
  }
}

#' Resolve and Standardize Species Names
#'
#' Attempts to correct misspelled, outdated, or synonymous species names using multiple taxonomic databases.
#'
#' @param names Character vector of species names to resolve.
#' @param sources Character vector of sources for taxize::gnr_resolve (default = c("GBIF", "ITIS")).
#' @param use_taxize Logical. Use `taxize::gnr_resolve()` as a secondary match engine (default = TRUE).
#' @param verbose Logical. If TRUE, prints matching progress.
#'
#' @return A tibble with columns: `original_name`, `resolved_name`, `match_type`, `source`, `score`.
#' @export
#'
#' @examples
#' \dontrun{
#' resolve_species_names(c("Panthera leo", "Homo sapiens", "Lantana camra", "Ailuropoda melanoleuca"))
#' }
resolve_species_names <- function(names,
                                  sources = c("GBIF", "ITIS"),
                                  use_taxize = TRUE,
                                  verbose = TRUE) {
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("tibble", quietly = TRUE))
  stopifnot(requireNamespace("rgbif", quietly = TRUE))

  results <- list()

  for (name in names) {
    if (verbose) cli::cli_alert_info("🔍 Resolving: {name}")

    # --- First attempt: GBIF backbone ---
    match_gbif <- tryCatch({
      res <- rgbif::name_backbone(name = name)
      if (!is.null(res$scientificName)) {
        tibble::tibble(
          original_name = name,
          resolved_name = res$scientificName,
          match_type = res$matchType,
          source = "GBIF",
          score = NA_real_
        )
      } else {
        NULL
      }
    }, error = function(e) NULL)

    if (!is.null(match_gbif)) {
      results[[name]] <- match_gbif
      next
    }

    # --- Second attempt: Taxize resolver ---
    if (use_taxize && requireNamespace("taxize", quietly = TRUE)) {
      match_taxize <- tryCatch({
        tax_res <- taxize::gnr_resolve(name, preferred_data_sources = sources, best_match_only = TRUE)
        if (nrow(tax_res) > 0) {
          tibble::tibble(
            original_name = name,
            resolved_name = tax_res$matched_name2[1],
            match_type = "resolved",
            source = tax_res$data_source_title[1],
            score = tax_res$score[1]
          )
        } else {
          NULL
        }
      }, error = function(e) NULL)

      if (!is.null(match_taxize)) {
        results[[name]] <- match_taxize
        next
      }
    }

    # --- Final fallback: unresolved ---
    results[[name]] <- tibble::tibble(
      original_name = name,
      resolved_name = NA_character_,
      match_type = "unresolved",
      source = NA_character_,
      score = NA_real_
    )
  }

  return(dplyr::bind_rows(results))
}

#' Determine if a species is marine
#'
#' Uses taxonomic hints (via `taxize::classification`) and optionally known marine groups or name patterns.
#'
#' @param species_name Scientific name of the species.
#' @param verbose Logical. Print classification info.
#' @return Logical. TRUE if the species is likely marine.
#' @export
is_marine_species <- function(species_name, verbose = FALSE) {
  if (!requireNamespace("taxize", quietly = TRUE)) stop("Package 'taxize' is required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' is required.")

  known_marine_clades <- c(
    "Echinodermata", "Cnidaria", "Porifera", "Mollusca", "Crustacea",
    "Chondrichthyes", "Actinopterygii", "Elasmobranchii",
    "Cephalopoda", "Bivalvia", "Gastropoda", "Decapoda",
    "Bryozoa", "Ascidiacea", "Hydrozoa", "Anthozoa"
  )

  classification <- tryCatch({
    out <- taxize::classification(species_name, db = "ncbi")
    if (!is.null(out) && is.list(out) && is.data.frame(out[[1]])) {
      out[[1]]
    } else {
      out <- taxize::classification(species_name, db = "gbif")
      if (!is.null(out) && is.list(out) && is.data.frame(out[[1]])) out[[1]] else NULL
    }
  }, error = function(e) NULL)

  clades <- if (!is.null(classification)) classification$name else character(0)
  is_marine_clade <- any(clades %in% known_marine_clades)

  # Optional fallback: check genus name
  marine_genus_patterns <- c("Mytilus", "Scomber", "Gadus", "Pomacentrus", "Sepia", "Engraulis")
  is_marine_by_name <- any(stringr::str_detect(species_name, paste(marine_genus_patterns, collapse = "|")))

  is_marine <- is_marine_clade || is_marine_by_name

  if (verbose) {
    cli::cli_alert_info("Taxonomic path for {.val {species_name}}: {paste(clades, collapse = ' > ')}")
    if (is_marine_clade) cli::cli_alert_success("✔ Identified as marine from taxonomy.")
    if (is_marine_by_name && !is_marine_clade) cli::cli_alert_success("✔ Fallback: identified as marine from name pattern.")
    if (!is_marine) cli::cli_alert_warning("✖ Not identified as marine.")
  }

  return(is_marine)
}
