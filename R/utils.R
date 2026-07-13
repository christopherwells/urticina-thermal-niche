# utils.R
# Shared utilities sourced by analysis.qmd, tpc.qmd, and make_figures.R:
#   source(here::here("R", "utils.R"))
#
# In this reproducibility repository the analysis-ready data live in data/clean/,
# so load_output()/save_output() read and write there (the full pipeline used a
# separate output/ directory; only the read/write root differs here).

library(here)

# --- Coordinate Reference Systems ----
crs_project <- "EPSG:32619"   # UTM Zone 19N (Gulf of Maine)
crs_lonlat  <- "EPSG:4326"    # WGS84 geographic

# --- Study Extent Bounding Box (approximate) ----
# Massachusetts through Newfoundland, with buffer for context
study_bbox <- c(xmin = -71, xmax = -52, ymin = 41.5, ymax = 51)

# Note: effort-tracker taxonomic filtering is handled pre-download via
# the iNaturalist project settings. All observations in the download
# are marine intertidal taxa and serve as effort trackers.

# --- Output Helpers ----
# Analysis-ready inputs and derived outputs both live in data/clean/.
save_output <- function(obj, name) {
  path <- here("data", "clean", paste0(name, ".rds"))
  saveRDS(obj, path)
  message("Saved: ", path)
}

load_output <- function(name) {
  path <- here("data", "clean", paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop(
      "Data file not found: ", path, "\n",
      "Render analysis.qmd (and tpc.qmd) first to generate derived outputs.",
      call. = FALSE
    )
  }
  readRDS(path)
}

# --- ggplot2 Themes ----

# Non-map panels (TPC, regression diagnostics, etc.)
theme_anemone <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold")
    )
}

# Map theme (shared across all spatial figures)
theme_map <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(color = "black", size = base_size),
      axis.title = ggplot2::element_text(color = "black", size = base_size),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.ticks.length = ggplot2::unit(2, "pt"),
      legend.text = ggplot2::element_text(color = "black", size = base_size),
      legend.title = ggplot2::element_text(color = "black", size = base_size),
      legend.position = "bottom",
      plot.tag = ggplot2::element_text(color = "black", size = base_size)
    )
}

# Standard colorbar guide (black frame, inward ticks)
guide_cb <- ggplot2::guide_colorbar(
  frame.colour = "black", frame.linewidth = 0.5,
  ticks.colour = "black", ticks.linewidth = 0.5,
  ticks.length = ggplot2::unit(-2, "pt"),
  title.position = "left", title.vjust = 0.5
)

# Standard map scales (longitude/latitude breaks)
map_scales <- list(
  ggplot2::scale_x_continuous(breaks = seq(-80, -50, by = 5)),
  ggplot2::scale_y_continuous(breaks = seq(40, 55, by = 3))
)

# --- Display Bounding Box ----
# Same as study_bbox for consistency between analyses and figures
display_bbox <- study_bbox

# --- Base Map ----
base_map <- function(bbox = display_bbox) {
  x1 <- unname(bbox["xmin"]); x2 <- unname(bbox["xmax"])
  y1 <- unname(bbox["ymin"]); y2 <- unname(bbox["ymax"])
  crop_box <- sf::st_bbox(
    c(xmin = x1 - 3, ymin = y1 - 0.5, xmax = x2 + 0.5, ymax = y2 + 3),
    crs = sf::st_crs(crs_lonlat)
  )

  land <- sf::st_make_valid(rnaturalearthhires::countries10)
  land_crop <- suppressWarnings(sf::st_crop(land, crop_box))

  list(
    land = land_crop,
    coord = ggplot2::coord_sf(xlim = c(x1, x2), ylim = c(y1, y2),
                               crs = sf::st_crs(crs_lonlat), expand = FALSE)
  )
}

# Convenience: build a complete base map ggplot
base_map_plot <- function(bbox = display_bbox) {
  bm <- base_map(bbox)
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = bm$land, fill = "grey90", color = "grey40", linewidth = 0.2) +
    bm$coord +
    map_scales +
    theme_map() +
    ggplot2::labs(x = "Longitude", y = "Latitude")
}

# --- Suitability extraction to H3 hexes ---
# The MaxEnt SDM is fit at hex level, so 05-sdm-prediction.rds is an sf
# with one row per coastal hex (h3_address + suitability). This helper joins
# the suitability values back onto a target hex sf by h3_address, returning a
# vector aligned to the target's row order. Hexes not present in the source
# (offshore, missing predictors) come back as NA.
#
# A SpatRaster fallback path is preserved for older RDS outputs that stored a
# pixel-level suitability raster; that path does pixel-to-hex zonal-mean
# extraction with 1-ring neighbor-fill for land-dominated hexes.
extract_suitability <- function(hex_sf, suit_data) {
  if (inherits(suit_data, c("sf", "data.frame"))) {
    df <- if (inherits(suit_data, "sf")) sf::st_drop_geometry(suit_data) else suit_data
    if (!all(c("h3_address", "suitability") %in% names(df))) {
      stop("suit_data must have columns 'h3_address' and 'suitability'")
    }
    return(df$suitability[match(hex_sf$h3_address, df$h3_address)])
  }

  hex_proj <- sf::st_transform(hex_sf, terra::crs(suit_data))
  raw <- terra::extract(suit_data, terra::vect(hex_proj), fun = mean, na.rm = TRUE)
  suit_col <- setdiff(names(raw), "ID")[1]
  vals <- raw[[suit_col]]

  na_idx <- which(is.na(vals))
  if (length(na_idx) > 0 && requireNamespace("h3jsr", quietly = TRUE)) {
    h3_ids <- hex_sf$h3_address
    for (i in na_idx) {
      neighbors <- h3jsr::get_disk(h3_ids[i], ring_size = 1)[[1]]
      nb_idx <- match(neighbors, h3_ids)
      nb_idx <- nb_idx[!is.na(nb_idx)]
      nb_vals <- vals[nb_idx]
      nb_vals <- nb_vals[!is.na(nb_vals)]
      if (length(nb_vals) > 0) vals[i] <- mean(nb_vals)
    }
  }
  vals
}
