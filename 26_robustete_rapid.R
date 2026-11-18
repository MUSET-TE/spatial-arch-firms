# =====================================================================
# 26_robustete_rapid.R
# Estimation of the model under several definitions of the neighbourhood:
# Delaunay and Gabriel (parameter free), distance bands (which give the radius
# of transmission in km) and k-nearest with k = 4 and 8, the classic check.
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(sf); library(spdep); library(Matrix)
  library(spGARCH); library(parallel)
})
options(datatable.showProgress = FALSE)
source("00_config.R")
FIS <- file.path(OUT, "robustete_W.csv")
INC <- file.path(OUT, "robustete_W_partial.csv")
if (file.exists(INC)) file.remove(INC)
if (file.exists(FIS)) file.remove(FIS)

P   <- as.data.table(readRDS(file.path(OUT, "panel_uat.rds")))
lau <- st_read(file.path(GEO, "lau2021.gpkg"), quiet = TRUE)
lau <- lau[lau$CNTR_CODE == "RO", ]
lau$siruta <- suppressWarnings(as.integer(lau$LAU_ID))
keep <- P[, .(m = min(firme)), by = siruta][m >= 50, siruta]
lf   <- lau[lau$siruta %in% keep, ]; lf <- lf[order(lf$siruta), ]
pts  <- suppressWarnings(st_coordinates(st_centroid(st_geometry(lf))))
n    <- nrow(pts)
cat("n =", n, "| CRS in metres:", st_is_longlat(lf) == FALSE, "\n")

# distance to the nearest neighbour, used to choose the bands
k1 <- knn2nb(knearneigh(pts, k = 1))
d1 <- unlist(nbdists(k1, pts)) / 1000
cat(sprintf("distance to the nearest neighbour (km): median %.1f, p95 %.1f, max %.1f\n\n",
            median(d1), quantile(d1, .95), max(d1)))

# ---- the neighbourhood definitions -----------------------------------
mk <- list()
mk[["Delaunay"]]  <- tri2nb(pts)
g <- gabrielneigh(pts); mk[["Gabriel"]] <- graph2nb(g, sym = TRUE)
if (FALSE) for (dkm in c(40, 60, 80, 100)) {
  mk[[paste0("Distance ", dkm, " km")]] <- dnearneigh(pts, 0, dkm * 1000)
}
mk[["k = 4"]] <- make.sym.nb(knn2nb(knearneigh(pts, k = 4)))
mk[["k = 5"]] <- make.sym.nb(knn2nb(knearneigh(pts, k = 5)))   # baseline
mk[["k = 8"]] <- make.sym.nb(knn2nb(knearneigh(pts, k = 8)))

info <- rbindlist(lapply(names(mk), function(nm) {
  nb <- mk[[nm]]
  data.table(W = nm, vecini_medii = round(mean(card(nb)), 2),
             izolate = sum(card(nb) == 0),
             muchii = sum(card(nb)) / 2)
}))
print(info)
cat("\n")

# the matrices, row standardised
WL <- lapply(mk, function(nb) as(nb2mat(nb, style = "W", zero.policy = TRUE), "dgCMatrix"))
LW <- lapply(mk, function(nb) nb2listw(nb, style = "W", zero.policy = TRUE))

# ---- fixed effect residuals, identical to those in the paper ---------
D <- P[siruta %in% lf$siruta]
D[, ly := log(ca / pmax(sal, 1))]
q <- quantile(D$ly, c(.005, .995), na.rm = TRUE)
D[, ly := pmin(pmax(ly, q[1]), q[2])]
D[, e := ly - mean(ly), by = siruta]
D[, e := e - mean(e), by = an]
ANI <- sort(unique(D$an))
YS <- lapply(ANI, function(a) D[an == a][data.table(siruta = lf$siruta), on = "siruta"]$e)
names(YS) <- as.character(ANI)

job <- CJ(W = names(mk), an = ANI, sorted = FALSE)
cat("to estimate:", nrow(job), "models\n\n")

unul <- function(r) {
  nm <- job$W[r]; a <- job$an[r]
  y <- YS[[as.character(a)]]; Wm <- WL[[nm]]
  t0 <- Sys.time()
  m <- try(qml.spARCH(y ~ 1, W = Wm, type = "spARCH"), silent = TRUE)
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  mi <- try(moran.test(y^2, LW[[nm]], zero.policy = TRUE), silent = TRUE)
  if (inherits(m, "try-error")) {
    return(data.table(W = nm, an = a, rho = NA_real_, rho_t = NA_real_,
                      moran = NA_real_, moran_p = NA_real_, sec = round(sec, 1)))
  }
  co <- summary(m)$coef; rr <- grep("rho", rownames(co))
  data.table(W = nm, an = a, rho = co[rr, 1], rho_t = co[rr, 3],
             moran = if (inherits(mi, "try-error")) NA_real_ else unname(mi$estimate[1]),
             moran_p = if (inherits(mi, "try-error")) NA_real_ else mi$p.value,
             sec = round(sec, 1))
}

cl <- makeCluster(16)
on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({
  library(spGARCH); library(Matrix); library(spdep); library(data.table)}))
clusterExport(cl, c("job", "YS", "WL", "LW"), envir = environment())
t0 <- Sys.time()
res <- rbindlist(lapply(split(seq_len(nrow(job)), ceiling(seq_len(nrow(job))/16)), function(ix) {
  r <- rbindlist(parLapply(cl, ix, unul), fill = TRUE)
  fwrite(r, INC, append = file.exists(INC))
  cat("  ", nrow(r), "estimations done\n"); flush.console(); r
}), fill = TRUE)
cat("total time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minutes\n\n")
fwrite(res, FIS)

sint <- res[, .(ani_estimati = sum(is.finite(rho)),
                rho_mediu = round(mean(rho, na.rm = TRUE), 3),
                ani_semnif_moran = sum(moran_p < .05, na.rm = TRUE),
                rho_in_anii_moran = round(mean(rho[moran_p < .05], na.rm = TRUE), 3)), by = W]
sint <- merge(info, sint, by = "W", sort = FALSE)
print(sint)
fwrite(sint, file.path(OUT, "robustete_W_sinteza.csv"))
cat("\nSaved:", FIS, "\n")
