#' Install optional data-provider dependencies
#'
#' @param griis  Use GRIIS via {originr}
#' @param worms  Use WoRMS via {worrms}
#' @param gbif   Use GBIF distributions via {rgbif}
#' @param prefer "runiverse" (default) or "github"
#' @return (invisible) TRUE on success, FALSE if any package failed
#' @export
install_optional_deps <- function(griis = TRUE, worms = TRUE, gbif = TRUE,
                                  prefer = c("runiverse","github")) {
  prefer <- match.arg(prefer)

  want <- unique(c(
    if (griis)  c("originr","countrycode") else character(0),
    if (worms)  "worrms"                    else character(0),
    if (gbif)   c("rgbif","countrycode")    else character(0)
  ))

  need <- want[!vapply(want, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(need)) { message("All optional deps already installed."); return(invisible(TRUE)) }

  # R-universe repos (include archive for originr) + CRAN
  repos_runiverse <- c(
    "https://ropensci-archive.r-universe.dev",
    "https://ropensci.r-universe.dev",
    "https://cloud.r-project.org"
  )

  # Correct GitHub homes for fallbacks
  gh_map <- c(
    originr     = "ropensci-archive/originr",
    worrms      = "ropensci/worrms",
    rgbif       = "ropensci/rgbif",
    countrycode = "vincentarelbundock/countrycode"
  )

  if (prefer == "runiverse") {
    message("Installing from R-universe/CRAN where available...")
    try(install.packages(need, repos = repos_runiverse), silent = TRUE)

    # fall back to GitHub for anything still missing
    still_need <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
    if (length(still_need)) {
      if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
      for (p in still_need) {
        repo <- unname(gh_map[p]); if (is.na(repo)) next
        message("Falling back to GitHub: ", repo)
        try(remotes::install_github(repo, upgrade = "never"), silent = TRUE)
      }
    }
  } else {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    for (p in need) {
      repo <- unname(gh_map[p]); if (is.na(repo)) next
      message("Installing from GitHub: ", repo)
      try(remotes::install_github(repo, upgrade = "never"), silent = TRUE)
    }
  }

  missing <- want[!vapply(want, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    warning("Some optional deps failed to install: ", paste(missing, collapse = ", "))
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

#' Install optional data-provider dependencies
#'
#' Convenience installer for optional data sources used by biofetchR.
#'
#' @param griis  Install packages for GRIIS via {originr}.
#' @param worms  Install packages for WoRMS via {worrms}.
#' @param gbif   Install packages for GBIF via {rgbif}.
#' @param osm    Install packages for OSM-based river naming ({osmdata}, {rnaturalearth}, {rnaturalearthdata}).
#' @param prefer "runiverse" (default) or "github".
#'
#' @return (invisible) TRUE on success, FALSE if any package failed.
#' @export
install_optional_deps <- function(
    griis  = TRUE,
    worms  = TRUE,
    gbif   = TRUE,
    osm    = FALSE,
    prefer = c("runiverse","github")
) {
  prefer <- match.arg(prefer)

  want <- unique(c(
    if (griis) c("originr","countrycode") else character(0),
    if (worms) "worrms"                   else character(0),
    if (gbif)  c("rgbif","countrycode")   else character(0),
    if (osm)   c("osmdata","rnaturalearth","rnaturalearthdata") else character(0)
  ))

  need <- want[!vapply(want, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(need)) { message("All requested optional deps already installed."); return(invisible(TRUE)) }

  # R-universe repos + CRAN (covers most packages quickly)
  repos_runiverse <- c(
    "https://ropensci-archive.r-universe.dev",  # originr (archived)
    "https://ropensci.r-universe.dev",          # worrms, rgbif, osmdata, etc.
    "https://cloud.r-project.org"               # CRAN (rnaturalearth, rnaturalearthdata, etc.)
  )

  # GitHub fallbacks
  gh_map <- c(
    originr            = "ropensci-archive/originr",
    worrms             = "ropensci/worrms",
    rgbif              = "ropensci/rgbif",
    countrycode        = "vincentarelbundock/countrycode",
    osmdata            = "ropensci/osmdata",
    rnaturalearth      = "ropensci/rnaturalearth",
    rnaturalearthdata  = "ropensci/rnaturalearthdata"
  )

  if (prefer == "runiverse") {
    message("Installing from R-universe/CRAN where available...")
    try(install.packages(need, repos = repos_runiverse), silent = TRUE)

    still_need <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
    if (length(still_need)) {
      if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
      for (p in still_need) {
        repo <- unname(gh_map[p]); if (is.na(repo)) next
        message("Falling back to GitHub: ", repo)
        try(remotes::install_github(repo, upgrade = "never"), silent = TRUE)
      }
    }
  } else {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    for (p in need) {
      repo <- unname(gh_map[p]); if (is.na(repo)) next
      message("Installing from GitHub: ", repo)
      try(remotes::install_github(repo, upgrade = "never"), silent = TRUE)
    }
  }

  missing <- want[!vapply(want, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    warning("Some optional deps failed to install: ", paste(missing, collapse = ", "))
    return(invisible(FALSE))
  }
  invisible(TRUE)
}
