# =====================================================================
# Figure 1: the spatial lattices. Panel (a) the full lattice of 3,152
# localities with contiguity, panel (b) the estimation lattice of 666 with
# the 5 nearest neighbours. Greyscale, for a black and white manuscript.
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(sf); library(spdep); library(ggplot2)
})
source("00_config.R")

P   <- as.data.table(readRDS(file.path(OUT, "panel_uat.rds")))
S   <- readRDS(file.path(OUT, "spatial_weights.rds"))
lau <- st_read(file.path(GEO, "lau2021.gpkg"), quiet = TRUE)
lau <- lau[lau$CNTR_CODE == "RO", ]
lau$siruta <- suppressWarnings(as.integer(lau$LAU_ID))
nuts <- st_read(file.path(GEO, "nuts3.geojson"), quiet = TRUE)
nuts <- nuts[substr(nuts$NUTS_ID, 1, 2) == "RO", ]
nuts <- st_transform(nuts, st_crs(lau))

# ---- the edges of a neighbourhood graph, as segments -------------------
edges <- function(nb, pts) {
  ii <- unlist(lapply(seq_along(nb), function(k)
    if (nb[[k]][1] != 0) rep(k, length(nb[[k]])) else integer(0)))
  jj <- unlist(lapply(nb, function(v) if (v[1] != 0) v else integer(0)))
  keep <- ii < jj                      # each edge only once
  ii <- ii[keep]; jj <- jj[keep]
  data.table(x = pts[ii, 1], y = pts[ii, 2], xend = pts[jj, 1], yend = pts[jj, 2])
}

# ---- panel (a): the full lattice --------------------------------------
lf_a <- lau[lau$siruta %in% S$siruta, ]
lf_a <- lf_a[order(lf_a$siruta), ]
pt_a <- suppressWarnings(st_coordinates(st_centroid(st_geometry(lf_a))))
e_a  <- edges(S$nb_uat, pt_a)
d_a  <- data.table(x = pt_a[, 1], y = pt_a[, 2])

# ---- panel (b): the estimation lattice --------------------------------
keep <- P[, .(m = min(firme)), by = siruta][m >= 50, siruta]
lf_b <- lau[lau$siruta %in% keep, ]
lf_b <- lf_b[order(lf_b$siruta), ]
pt_b <- suppressWarnings(st_coordinates(st_centroid(st_geometry(lf_b))))
nb_b <- make.sym.nb(knn2nb(knearneigh(pt_b, k = 5)))
e_b  <- edges(nb_b, pt_b)
d_b  <- data.table(x = pt_b[, 1], y = pt_b[, 2])

cat("panel a:", nrow(d_a), "localities,", nrow(e_a), "edges\n")
cat("panel b:", nrow(d_b), "localities,", nrow(e_b), "edges\n")

GREY_BORDER <- "#C8C8C8"; GREY_EDGE <- "#8C8C8C"; INK <- "#1A1A1A"

panel <- function(dd, ee, lab, ptsize, ewidth, ealpha) {
  ggplot() +
    geom_sf(data = nuts, fill = NA, colour = GREY_BORDER, linewidth = 0.22) +
    geom_segment(data = ee, aes(x = x, y = y, xend = xend, yend = yend),
                 colour = GREY_EDGE, linewidth = ewidth, alpha = ealpha) +
    geom_point(data = dd, aes(x = x, y = y), colour = INK,
               size = ptsize, shape = 16) +
    labs(subtitle = lab) +
    theme_void(base_family = "serif", base_size = 11) +
    theme(plot.subtitle = element_text(size = 10, colour = INK, hjust = 0,
                                       margin = margin(b = 4)),
          plot.margin = margin(4, 6, 4, 6))
}

pa <- panel(d_a, e_a, "(a) Full lattice: 3,152 localities, queen contiguity",
            0.18, 0.10, 0.35)
pb <- panel(d_b, e_b, "(b) Estimation lattice: 666 localities, 5 nearest neighbours",
            0.42, 0.16, 0.50)

ok <- requireNamespace("patchwork", quietly = TRUE)
if (ok) {
  library(patchwork)
  fig <- pa | pb
  ggsave(file.path(OUT, "Figure1_lattice.png"), fig,
         width = 11, height = 4.2, dpi = 400, bg = "white")
  ggsave(file.path(OUT, "Figure1a_full_lattice.png"), pa,
         width = 6.0, height = 4.2, dpi = 400, bg = "white")
  ggsave(file.path(OUT, "Figure1b_estimation_lattice.png"), pb,
         width = 6.0, height = 4.2, dpi = 400, bg = "white")
} else {
  ggsave(file.path(OUT, "Figure1a_full_lattice.png"), pa,
         width = 5.2, height = 5.2, dpi = 400, bg = "white")
  ggsave(file.path(OUT, "Figure1b_estimation_lattice.png"), pb,
         width = 5.2, height = 5.2, dpi = 400, bg = "white")
  cat("patchwork missing, two separate files were saved instead\n")
}
cat("saved in", OUT, "\n")
