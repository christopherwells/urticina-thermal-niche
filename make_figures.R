# make_figures.R
# -----------------------------------------------------------------------------
# Reproduces manuscript Figures 1-5 for the Urticina thermal-niche paper.
#
#   Figure 1  Study area with Urticina occurrences (fig1-study-area.png)
#   Figure 2  NB-predicted encounter rate + MaxEnt suitability maps
#             (fig2-sdm-cpue-comparison.png)
#   Figure 3  NB encounter rate vs MaxEnt suitability scatter
#             (fig3-nb-maxent-scatter.png)
#   Figure 4  Best-fitting Sharpe-Schoolfield multiplicative TPC
#             (fig4-tpc-bestmodel.png)
#   Figure 5  MaxEnt response curves for non-zero predictors
#             (fig5-maxent-response-curves.png)
#
# All inputs are read from data/clean/ via load_output(). Figures 1-3 and 5
# depend on the fitted NB and MaxEnt objects, so render analysis.qmd FIRST
# (it writes 04-*/05-* into data/clean/). Figure 4 reads the fitted TPC objects
# in data/clean/06-tpc-fits.rds (shipped), so it needs no re-sampling.
#
# Output PNGs are written to figures/ (created here; gitignored).
#
# Run:  Rscript make_figures.R   (from the repository root)
# -----------------------------------------------------------------------------

source(here::here("R", "utils.R"))

library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyterra)
library(patchwork)
library(cowplot)
library(ENMeval)
library(maxnet)
library(glmmTMB)
library(mgcv)
library(scales)
library(viridisLite)
library(rnaturalearthhires)

dir.create(here::here("figures"), showWarnings = FALSE, recursive = TRUE)

# --- Load upstream outputs ---------------------------------------------------
urticina_sf  <- load_output("01-urticina-obs")
h3_data      <- load_output("02-cpue-hex")
h3_uncorr    <- load_output("03-predictors-uncorr")
suitability  <- load_output("05-sdm-prediction")
sdm_model    <- load_output("05-maxent-model")
nb_model     <- load_output("04-nb-model")
tpc_outputs  <- load_output("06-tpc-fits")
scale_params <- load_output("04-scale-params")
nb_selected  <- load_output("04-nb-selected-predictors")

# High-resolution (1:10m) land polygons for map figures
land <- countries10

# --- Shared styling ----------------------------------------------------------
shared_coord <- coord_sf(
  xlim = c(display_bbox["xmin"], display_bbox["xmax"]),
  ylim = c(display_bbox["ymin"], display_bbox["ymax"]),
  crs = st_crs(crs_lonlat), expand = FALSE
)

shared_theme <- theme_minimal(base_size = 10) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(2, "pt"),
    legend.text = element_text(color = "black", size = 10),
    legend.title = element_text(color = "black", size = 10),
    legend.position = "right",
    plot.tag = element_text(color = "black", size = 10)
  )

shared_scales <- list(
  scale_x_continuous(breaks = seq(-80, -50, by = 5)),
  scale_y_continuous(breaks = seq(40, 55, by = 3))
)

tpc_theme <- theme_minimal(base_size = 10) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(2, "pt"),
    legend.text = element_text(color = "black", size = 10),
    legend.title = element_text(color = "black", size = 10),
    legend.position = "right",
    plot.tag = element_text(color = "black", size = 10)
  )

# Sentence-case plain-language labels for predictor axes
pred_label_plain <- c(
  sst_summer_max    = "summer maximum SST",
  sst_winter_min    = "winter minimum SST",
  salinity          = "salinity",
  chlorophyll       = "chlorophyll-a concentration",
  ph                = "pH",
  sea_ice           = "sea ice concentration",
  bottom_temp       = "bottom temperature",
  cloud_cover       = "cloud cover",
  wind_speed_summer = "wind speed",
  air_temp          = "air temperature",
  precip_summer     = "precipitation",
  vpd_summer        = "vapor pressure deficit",
  solar_rad_summer  = "solar radiation",
  shoreline_km      = "shoreline length",
  pct_developed     = "developed coastline",
  park_fraction     = "park fraction"
)
plain <- function(v) dplyr::if_else(v %in% names(pred_label_plain),
                                    unname(pred_label_plain[v]), v)
sentence_case <- function(s) {
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

# --- Variable importance (non-zero MaxEnt feature coefficients) --------------
vimp_df <- NULL
if (!is.null(sdm_model$betas)) {
  vcounts <- sapply(nb_selected, function(v) {
    sum(grepl(paste0("^", v), names(sdm_model$betas)) & sdm_model$betas != 0)
  })
  vimp_df <- data.frame(variable = nb_selected, n_features = as.integer(vcounts),
                        stringsAsFactors = FALSE) |>
    dplyr::arrange(dplyr::desc(n_features))
}

# --- Concordance data: MaxEnt suitability vs NB-predicted encounter rate -----
# Coastal tiles used in Figures 2, 3, and 5 (shoreline > 0, complete predictors).
concordance_df <- local({
  hx <- h3_uncorr |> filter(shoreline_km > 0)
  keep <- complete.cases(st_drop_geometry(hx)[, scale_params$variable])
  hx <- hx[keep, ]
  suit_vals_c <- extract_suitability(hx, suitability)
  suit_vals_c[is.na(suit_vals_c)] <- 0
  nb_pred_c <- sf::st_drop_geometry(hx)
  for (v in scale_params$variable) {
    sp <- scale_params[scale_params$variable == v, ]
    if (v %in% names(nb_pred_c)) {
      nb_pred_c[[v]] <- (nb_pred_c[[v]] - sp$center) / sp$scale
    }
  }
  med_log_eff <- median(log(h3_data$effort_n[h3_data$effort_n > 0]))
  nb_pred_c$log_effort <- med_log_eff
  pred_counts_c <- predict(nb_model, newdata = nb_pred_c, type = "response")
  cpue_c <- pred_counts_c / exp(med_log_eff) * 1000
  data.frame(
    suitability    = suit_vals_c,
    predicted_cpue = cpue_c,
    has_urticina   = hx$urticina_n > 0
  )
})
concordance_rho <- {
  ok <- !is.na(concordance_df$suitability) & !is.na(concordance_df$predicted_cpue)
  round(cor(concordance_df$suitability[ok], concordance_df$predicted_cpue[ok],
            method = "spearman"), 2)
}

# =============================================================================
# FIGURE 1 -- Study area with Urticina occurrences
# =============================================================================
message("Building Figure 1 (study area) ...")

# Snap Urticina points to nearest coastline for display
sf_use_s2(FALSE)
land_crop <- land |>
  st_make_valid() |>
  st_crop(st_bbox(c(
    xmin = unname(display_bbox["xmin"]),
    ymin = unname(display_bbox["ymin"]),
    xmax = unname(display_bbox["xmax"]),
    ymax = unname(display_bbox["ymax"])
  ), crs = crs_lonlat)) |>
  st_union() |>
  st_make_valid()
coastline <- st_boundary(land_crop)

snapped_coords <- matrix(NA, nrow = nrow(urticina_sf), ncol = 2)
for (i in seq_len(nrow(urticina_sf))) {
  np <- st_nearest_points(urticina_sf[i, ], coastline)
  snapped_coords[i, ] <- st_coordinates(np)[2, 1:2]
}
urticina_snapped <- st_as_sf(
  data.frame(x = snapped_coords[, 1], y = snapped_coords[, 2]),
  coords = c("x", "y"), crs = crs_lonlat
)
sf_use_s2(TRUE)

p_study <- ggplot() +
  geom_sf(data = land, fill = "#F5E6C8", color = "grey40", linewidth = 0.2) +
  geom_sf(data = urticina_snapped, fill = "darkorange", color = "black",
          shape = 21, size = 2.5, alpha = 0.7, stroke = 0.3) +
  shared_scales +
  shared_coord +
  shared_theme +
  theme(panel.background = element_rect(fill = "aliceblue", color = NA)) +
  labs(x = "Longitude", y = "Latitude")

# Optional inset photograph of the study organism, bottom-right of the panel.
# The photo is not archived in this repository; drop a copy at the path below
# to include the inset. Without it, Figure 1 is the base study-area map.
inset_path <- here::here("figures", "fig1-organism-inset.png")
if (file.exists(inset_path)) {
  inset_png  <- png::readPNG(inset_path)
  inset_grob <- grid::grobTree(
    grid::rasterGrob(inset_png, interpolate = TRUE,
                     width = grid::unit(1, "npc"), height = grid::unit(1, "npc")),
    grid::rectGrob(gp = grid::gpar(fill = NA, col = "black", lwd = 0.7))
  )
  inset_size <- 0.33; margin_r <- 0.025; margin_b <- 0.035; panel_aspect <- 1.39
  p_study_inset <- p_study +
    patchwork::inset_element(
      inset_grob,
      left  = 1 - margin_r - inset_size / panel_aspect, bottom = margin_b,
      right = 1 - margin_r,                             top    = margin_b + inset_size,
      align_to = "panel"
    )
} else {
  message("  (organism inset photo not found at ", inset_path,
          "; saving study-area map without inset)")
  p_study_inset <- p_study
}

ggsave(here::here("figures", "fig1-study-area.png"),
       p_study_inset, width = 7, height = 5.1, dpi = 600)

# =============================================================================
# FIGURE 2 -- NB-predicted encounter rate (A) + MaxEnt suitability (B)
# =============================================================================
message("Building Figure 2 (encounter rate + suitability maps) ...")

# Common hex set: all coastal cells with complete predictor data and shoreline
display_hexes <- h3_uncorr |>
  filter(shoreline_km > 0)
complete_mask <- complete.cases(st_drop_geometry(display_hexes)[, scale_params$variable])
display_hexes <- display_hexes[complete_mask, ]

# Panel B: MaxEnt suitability (zonal mean + neighbor-fill, consistent with pipeline)
h3_grid_suit <- display_hexes |>
  mutate(suitability = extract_suitability(display_hexes, suitability),
         suitability = if_else(is.na(suitability), 0, suitability))

# Helper: standalone vertical colorbar with right-side-only ticks
make_colorbar <- function(palette, title, limits, breaks, labels = breaks,
                          trans = "identity") {
  if (trans == "log1p") {
    trans_limits <- log1p(limits)
    trans_breaks <- log1p(breaks)
    grad_vals <- seq(trans_limits[1], trans_limits[2], length.out = 300)
  } else {
    trans_breaks <- breaks
    grad_vals <- seq(limits[1], limits[2], length.out = 300)
  }
  grad_df <- data.frame(x = 1, y = grad_vals)

  ggplot(grad_df, aes(x = x, y = y, fill = y)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_viridis_c(option = palette, guide = "none") +
    scale_y_continuous(breaks = trans_breaks, labels = labels,
                       expand = c(0, 0), position = "right") +
    scale_x_continuous(expand = c(0, 0)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = NA, color = "black", linewidth = 0.3) +
    ggtitle(title) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      axis.text.y.right = element_text(color = "black", size = 8,
                                       margin = margin(l = 2)),
      axis.ticks.y.right = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length.y.right = unit(2, "pt"),
      plot.title = element_text(color = "black", size = 8, hjust = 0.5,
                                margin = margin(b = 3)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 2, 5, 2)
    )
}

p_sdm <- ggplot() +
  geom_sf(data = land, fill = "#F5E6C8", color = "grey40", linewidth = 0.2) +
  geom_sf(data = h3_grid_suit, aes(fill = suitability), color = "grey40",
          linewidth = 0.02, alpha = 0.8) +
  scale_fill_viridis_c(
    option = "viridis", name = "Suitability",
    limits = c(0, 1), breaks = c(0, 0.5, 1.0),
    na.value = "transparent", guide = "none"
  ) +
  shared_scales +
  shared_coord +
  shared_theme +
  theme(panel.background = element_rect(fill = "aliceblue", color = NA)) +
  labs(x = "Longitude", y = "Latitude", tag = "B")

# Panel A: NB-predicted encounter rate
nb_pred_data <- display_hexes |> st_drop_geometry()
for (v in scale_params$variable) {
  sp <- scale_params[scale_params$variable == v, ]
  if (v %in% names(nb_pred_data)) {
    nb_pred_data[[v]] <- (nb_pred_data[[v]] - sp$center) / sp$scale
  }
}
median_log_effort <- median(log(h3_data$effort_n[h3_data$effort_n > 0]))
nb_pred_data$log_effort <- median_log_effort

pred_counts <- predict(nb_model, newdata = nb_pred_data, type = "response")
cpue_per <- 1000
nb_pred_sf <- display_hexes |>
  mutate(predicted_cpue = pred_counts / exp(median_log_effort) * cpue_per)

p_cpue <- ggplot() +
  geom_sf(data = land, fill = "#F5E6C8", color = "grey40", linewidth = 0.2) +
  geom_sf(data = nb_pred_sf,
          aes(fill = predicted_cpue), color = "grey40", linewidth = 0.02, alpha = 0.8) +
  scale_fill_viridis_c(
    option = "magma", name = "Predicted\nencounter rate\n(per 1000 obs.)",
    trans = "log1p", breaks = c(1, 10, 100),
    na.value = "transparent", guide = "none"
  ) +
  shared_scales +
  shared_coord +
  shared_theme +
  theme(panel.background = element_rect(fill = "aliceblue", color = NA)) +
  labs(x = "Longitude", y = "Latitude", tag = "A")

leg_cpue <- make_colorbar("magma", "Predicted\nencounter rate\n(per 1000 obs.)",
                          limits = c(0, 103), breaks = c(1, 10, 100),
                          trans = "log1p")
leg_suit <- make_colorbar("viridis", "Suitability",
                          limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1.0))

maps_col <- plot_grid(
  p_cpue + labs(x = NULL),
  p_sdm,
  ncol = 1, align = "hv", axis = "lr"
)

leg_cpue_padded <- plot_grid(NULL, leg_cpue, NULL, ncol = 1,
                             rel_heights = c(0.2, 0.6, 0.2))
leg_suit_padded <- plot_grid(NULL, leg_suit, NULL, ncol = 1,
                             rel_heights = c(0.2, 0.6, 0.2))
legs_col <- plot_grid(leg_cpue_padded, leg_suit_padded, ncol = 1)

p_sdm_cpue <- plot_grid(maps_col, legs_col,
                         ncol = 2, rel_widths = c(1, 0.09)) +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("figures", "fig2-sdm-cpue-comparison.png"),
       p_sdm_cpue, width = 6.5, height = 8, dpi = 600)

# =============================================================================
# FIGURE 3 -- NB encounter rate vs MaxEnt suitability scatter
# =============================================================================
message("Building Figure 3 (concordance scatter) ...")

scatter_df <- concordance_df |>
  dplyr::filter(!is.na(suitability), !is.na(predicted_cpue)) |>
  dplyr::mutate(status = factor(
    ifelse(has_urticina, "Present", "Absent"),
    levels = c("Present", "Absent")
  ))

absent_df  <- dplyr::filter(scatter_df, status == "Absent")
present_df <- dplyr::filter(scatter_df, status == "Present")

fit_full <- mgcv::gam(log1p(predicted_cpue) ~ s(suitability, k = 3, bs = "tp"),
                       data = scatter_df)
x_grid <- seq(min(scatter_df$suitability), max(scatter_df$suitability),
              length.out = 200)
pred_se <- predict(fit_full,
                   newdata = data.frame(suitability = x_grid),
                   se.fit = TRUE)
band_df <- data.frame(
  suitability = x_grid,
  fit = pmax(0, expm1(as.numeric(pred_se$fit))),
  lwr = pmax(0, expm1(as.numeric(pred_se$fit - 1.96 * pred_se$se.fit))),
  upr = pmax(0, expm1(as.numeric(pred_se$fit + 1.96 * pred_se$se.fit)))
)

p_scatter <- ggplot() +
  geom_ribbon(data = band_df,
              aes(x = suitability, ymin = lwr, ymax = upr),
              fill = "grey70", alpha = 0.8) +
  geom_line(data = band_df, aes(x = suitability, y = fit),
            linewidth = 0.7, color = "black") +
  geom_point(data = absent_df,
             aes(x = suitability, y = predicted_cpue, fill = status),
             shape = 21, color = "black", size = 1.4, stroke = 0.3) +
  geom_point(data = present_df,
             aes(x = suitability, y = predicted_cpue, fill = status),
             shape = 21, color = "black", size = 2.2, stroke = 0.3) +
  scale_fill_manual(
    values = c("Present" = "darkorange", "Absent" = NA),
    name = NULL,
    drop = FALSE,
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        color = "black",
        size = c(2.2, 1.4),
        stroke = 0.3,
        fill = c("darkorange", NA)
      )
    )
  ) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 10, 100)) +
  labs(x = "MaxEnt habitat suitability",
       y = "Predicted encounter rate (per 1000 obs.)") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(2, "pt"),
    legend.position = c(0.02, 0.92),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 10, color = "black"),
    legend.key.height = unit(12, "pt"),
    legend.margin = margin(0, 4, 0, 0),
    legend.spacing.y = unit(0, "pt")
  )

# Overlay a three-piece title so only rho renders italic. Text widths are
# measured in the output coordinate system via a ragg device.
rho_value  <- format(concordance_rho, nsmall = 2)
title_y    <- 0.95
x_spearman <- 0.116

meas_dev <- function(expr) {
  ragg::agg_png(tempfile(fileext = ".png"),
                width = 6.5, height = 5, units = "in", res = 600)
  on.exit(dev.off())
  par(family = "sans", ps = 10)
  force(expr)
}
meas_dev({
  w_sp  <- strwidth("Spearman ", units = "inches")
  w_rho <- strwidth("ρ ", units = "inches", font = 3)
})
plot_w_in <- 6.5
x_rho <- x_spearman + w_sp  / plot_w_in
x_eq  <- x_rho      + w_rho / plot_w_in

p_final <- cowplot::ggdraw(p_scatter) +
  cowplot::draw_label("Spearman", x = x_spearman, y = title_y,
                      hjust = 0, size = 10, fontface = "plain") +
  cowplot::draw_label("ρ", x = x_rho, y = title_y,
                      hjust = 0, size = 10, fontface = "italic") +
  cowplot::draw_label(paste0("= ", rho_value),
                      x = x_eq, y = title_y,
                      hjust = 0, size = 10, fontface = "plain")

ggsave(here::here("figures", "fig3-nb-maxent-scatter.png"),
       p_final, width = 6.5, height = 5, dpi = 600)

# =============================================================================
# FIGURE 4 -- Best-fitting Sharpe-Schoolfield multiplicative TPC
# =============================================================================
message("Building Figure 4 (best-model TPC) ...")

tpc_data     <- tpc_outputs$tpc_data
tpc_pred_exp <- tpc_outputs$tpc_pred_expanded
tpc_pred_clo <- tpc_outputs$tpc_pred_closed

# Per-treatment-temperature mean and 95% bootstrap CI of measured rates
set.seed(2026)
tpc_boot <- tpc_data |>
  dplyr::group_by(temperature = set.temp) |>
  dplyr::reframe({
    x <- rate_mass_specific[!is.na(rate_mass_specific)]
    boots <- replicate(5000, mean(sample(x, replace = TRUE)))
    qs <- quantile(boots, c(0.025, 0.975))
    data.frame(
      mean_rate = mean(x),
      ci_lo = unname(qs[1]),
      ci_hi = unname(qs[2]),
      mean_openness = mean(avg.openness, na.rm = TRUE)
    )
  })

# Mean-state curve: multiplicative model is linear in expansion state
mean_state <- mean(tpc_data$avg.openness, na.rm = TRUE)
pred_mean_state <- data.frame(
  temperature = tpc_pred_clo$temperature,
  rate = tpc_pred_clo$rate + mean_state * (tpc_pred_exp$rate - tpc_pred_clo$rate)
)

ribbon_df <- data.frame(
  temperature = tpc_pred_exp$temperature,
  rate_closed = tpc_pred_clo$rate,
  rate_expanded = tpc_pred_exp$rate
)
rate_ymax <- max(c(ribbon_df$rate_expanded, tpc_boot$ci_hi), na.rm = TRUE)

# Manual horizontal colorbar at data coords, centered below the curves
cbar_xmin <- 21.5
cbar_xmax <- 30.5
cbar_yc <- 60
cbar_hh <- 12
cbar_n <- 200
cbar_state_max <- 0.6
cbar_breaks <- c(0, 0.2, 0.4, 0.6)

cbar_df <- data.frame(
  xmin = cbar_xmin + (seq_len(cbar_n) - 1) / cbar_n * (cbar_xmax - cbar_xmin),
  xmax = cbar_xmin + seq_len(cbar_n)       / cbar_n * (cbar_xmax - cbar_xmin),
  ymin = cbar_yc - cbar_hh,
  ymax = cbar_yc + cbar_hh,
  fill_col = viridisLite::viridis(cbar_n)
)
tick_df <- data.frame(
  x = cbar_xmin + cbar_breaks / cbar_state_max * (cbar_xmax - cbar_xmin),
  label = format(cbar_breaks, nsmall = 1)
)

p_tpc_bestmodel <- ggplot(ribbon_df, aes(x = temperature)) +
  geom_ribbon(aes(ymin = rate_closed, ymax = rate_expanded),
              fill = "grey75", alpha = 0.5) +
  geom_line(aes(y = rate_expanded, linetype = "Fully expanded (state = 1)"),
            color = "black", linewidth = 0.9) +
  geom_line(aes(y = rate_closed, linetype = "Fully closed (state = 0)"),
            color = "black", linewidth = 0.9) +
  geom_line(data = pred_mean_state,
            aes(x = temperature, y = rate,
                linetype = "Mean expansion state"),
            color = "black", linewidth = 0.9,
            inherit.aes = FALSE) +
  geom_errorbar(data = tpc_boot,
                aes(x = temperature, ymin = ci_lo, ymax = ci_hi),
                width = 0.8, color = "black", linewidth = 0.4,
                inherit.aes = FALSE) +
  geom_point(data = tpc_boot,
             aes(x = temperature, y = mean_rate, fill = mean_openness),
             size = 3.5, color = "black", shape = 21, stroke = 0.5,
             inherit.aes = FALSE) +
  scale_fill_viridis_c(
    limits = c(0, cbar_state_max),
    breaks = cbar_breaks,
    oob = scales::squish,
    guide = "none"
  ) +
  geom_rect(data = cbar_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = cbar_df$fill_col, color = NA, inherit.aes = FALSE) +
  annotate("rect",
           xmin = cbar_xmin, xmax = cbar_xmax,
           ymin = cbar_yc - cbar_hh, ymax = cbar_yc + cbar_hh,
           color = "black", fill = NA, linewidth = 0.4) +
  geom_segment(data = tick_df,
               aes(x = x, xend = x,
                   y = cbar_yc - cbar_hh,
                   yend = cbar_yc - cbar_hh - 6),
               color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  geom_text(data = tick_df,
            aes(x = x, y = cbar_yc - cbar_hh - 10, label = label),
            hjust = 0.5, vjust = 1, size = 2.8, inherit.aes = FALSE) +
  annotate("text", x = (cbar_xmin + cbar_xmax) / 2,
           y = cbar_yc + cbar_hh + 22,
           label = "Mean expansion state", size = 3) +
  scale_linetype_manual(
    values = c("Fully expanded (state = 1)" = "dashed",
               "Mean expansion state"       = "solid",
               "Fully closed (state = 0)"   = "dotdash"),
    breaks = c("Fully expanded (state = 1)",
               "Mean expansion state",
               "Fully closed (state = 0)"),
    name = NULL,
    guide = guide_legend(order = 1)
  ) +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(limits = c(0, rate_ymax * 1.08), expand = c(0, 0)) +
  labs(x = expression("Temperature ("*degree*"C)"),
       y = expression(dot(V)[O[2]]~"(mg"~O[2]~g^{-1}~h^{-1}~")")) +
  tpc_theme +
  theme(
    legend.position = c(0.03, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.height = unit(12, "pt"),
    legend.key.width  = unit(28, "pt"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.spacing.y = unit(4, "pt"),
    legend.box.just = "left"
  )

ggsave(here::here("figures", "fig4-tpc-bestmodel.png"),
       p_tpc_bestmodel, width = 5, height = 3.8, dpi = 600)

# =============================================================================
# FIGURE 5 -- MaxEnt response curves for non-zero predictors
# =============================================================================
message("Building Figure 5 (MaxEnt response curves) ...")

# Order the non-zero predictors by |beta| so panels follow the Results order
non_zero_vars <- vimp_df$variable[vimp_df$n_features > 0]
beta_abs_order <- order(abs(sapply(non_zero_vars, function(v) {
  bi <- which(grepl(paste0("^", v), names(sdm_model$betas)))
  sum(sdm_model$betas[bi])
})), decreasing = TRUE)
non_zero_vars <- non_zero_vars[beta_abs_order]

# Two views of the coastal hex set: natural units (point x-values) and
# standardized (SDM prediction).
hex_natural <- h3_uncorr |>
  dplyr::filter(shoreline_km > 0) |>
  sf::st_drop_geometry()
keep_natural <- complete.cases(hex_natural[, scale_params$variable])
hex_natural <- hex_natural[keep_natural, ]

hex_std <- hex_natural
for (v in scale_params$variable) {
  sp <- scale_params[scale_params$variable == v, ]
  if (v %in% names(hex_std))
    hex_std[[v]] <- (hex_std[[v]] - sp$center) / sp$scale
}

hex_natural$suitability <- as.numeric(predict(sdm_model, hex_std,
                                              type = "cloglog", clamp = TRUE))
hex_natural$urticina_n <- h3_data$urticina_n[match(hex_natural$h3_address,
                                                   h3_data$h3_address)]
hex_natural$urticina_n[is.na(hex_natural$urticina_n)] <- 0
hex_natural$status <- factor(
  ifelse(hex_natural$urticina_n > 0, "Present", "Absent"),
  levels = c("Present", "Absent")
)

medians <- sapply(hex_std[, scale_params$variable, drop = FALSE],
                  median, na.rm = TRUE)

build_response_panel <- function(v, tag) {
  rng <- range(hex_std[[v]], na.rm = TRUE)
  grid_v <- seq(rng[1], rng[2], length.out = 200)
  newdata <- as.data.frame(matrix(rep(medians, length(grid_v)),
                                  nrow = length(grid_v), byrow = TRUE))
  names(newdata) <- names(medians)
  newdata[[v]] <- grid_v
  pred <- as.numeric(predict(sdm_model, newdata,
                             type = "cloglog", clamp = TRUE))

  sp <- scale_params[scale_params$variable == v, ]
  x_natural <- grid_v * sp$scale + sp$center

  point_data <- data.frame(
    x = hex_natural[[v]],
    y = hex_natural$suitability,
    status = hex_natural$status
  )
  absent_pts  <- dplyr::filter(point_data, status == "Absent")
  present_pts <- dplyr::filter(point_data, status == "Present")

  ggplot() +
    geom_point(data = absent_pts, aes(x = x, y = y),
               shape = 21, fill = NA, color = "black",
               size = 0.6, stroke = 0.2, alpha = 0.55) +
    geom_point(data = present_pts, aes(x = x, y = y),
               shape = 21, fill = "darkorange", color = "black",
               size = 1.2, stroke = 0.3, alpha = 0.95) +
    geom_line(data = data.frame(x = x_natural, y = pred),
              aes(x = x, y = y),
              linewidth = 0.7, color = "black") +
    scale_x_continuous(breaks = scales::breaks_extended(n = 6),
                       expand = expansion(mult = c(0.02, 0.06))) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0.02),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = sentence_case(pred_label_plain[v]),
         y = "Suitability",
         tag = tag) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text  = element_text(color = "black", size = 9),
      axis.title = element_text(color = "black", size = 9),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(2, "pt"),
      plot.tag = element_text(color = "black", size = 9),
      plot.tag.position = c(0.06, 0.95)
    )
}

panel_tags <- LETTERS[seq_along(non_zero_vars)]
panels <- lapply(seq_along(non_zero_vars), function(i)
  build_response_panel(non_zero_vars[i], panel_tags[i]))

fig5 <- cowplot::plot_grid(plotlist = c(panels, list(NULL)),
                           ncol = 3, align = "hv", axis = "lrtb")

ggsave(here::here("figures", "fig5-maxent-response-curves.png"),
       fig5, width = 6.5, height = 4, dpi = 600)

message("\nAll figures written to ", here::here("figures"))
