# =====================================================================
# Checking the estimator: we simulate data with a known rho and see what
# qml.spARCH recovers. The key question: on data WITHOUT a spatial
# ARCH effect (rho = 0), does the estimator still report a high rho?
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(sf); library(spdep); library(Matrix); library(spGARCH)
})
set.seed(20260820)
source("00_config.R")

P   <- readRDS(file.path(OUT, "panel_uat.rds"))
lau <- st_read(file.path(GEO, "lau2021.gpkg"), quiet = TRUE)
lau <- lau[lau$CNTR_CODE == "RO", ]
lau$siruta <- suppressWarnings(as.integer(lau$LAU_ID))
keep <- P[, .(m = min(firme)), by = siruta][m >= 200, siruta]
lf   <- lau[lau$siruta %in% keep, ]; lf <- lf[order(lf$siruta), ]
cent <- suppressWarnings(st_centroid(st_geometry(lf)))
nb   <- make.sym.nb(knn2nb(knearneigh(st_coordinates(cent), k = 5)))
W    <- as(nb2mat(nb, style = "W"), "dgCMatrix")
lw   <- nb2listw(nb, style = "W")
n    <- nrow(W)
cat("Same lattice as in the real estimations: n =", n, "\n\n")

REP <- 30
scenarii <- list(
  list(nume = "iid normal noise (true rho = 0)", gen = "iid",   tip = NA,          rho = 0),
  list(nume = "simulated spARCH, rho = 0.3",           gen = "sim",    tip = "spARCH",    rho = 0.3),
  list(nume = "simulated log-spARCH, rho = 0.5",       gen = "sim",    tip = "log-spARCH", rho = 0.5)
)

rez <- list(); k <- 0
for (sc in scenarii) {
  cat("### ", sc$nume, "\n", sep = "")
  for (r in seq_len(REP)) {
    y <- if (sc$gen == "iid") {
      rnorm(n)
    } else {
      as.numeric(sim.spARCH(n = n, rho = sc$rho, alpha = 1, W = W, type = sc$tip))
    }
    y <- y - mean(y)
    mi2 <- moran.test(y^2, lw, zero.policy = TRUE)
    for (tp in c("spARCH", "log-spARCH")) {
      m <- try(qml.spARCH(y ~ 1, W = W, type = tp), silent = TRUE)
      if (inherits(m, "try-error")) next
      co <- summary(m)$coef
      rr <- grep("rho", rownames(co))
      k <- k + 1
      rez[[k]] <- data.table(scenariu = sc$nume, rho_adevarat = sc$rho, rep = r,
                             estimator = tp, rho_est = co[rr, 1], rho_t = co[rr, 3],
                             moran_y2 = unname(mi2$estimate[1]), moran_p = mi2$p.value)
    }
  }
}
R <- rbindlist(rez)
fwrite(R, file.path(OUT, "verificare_estimator.csv"))

cat("\n================ REZULTATE ================\n")
print(R[, .(repetari = .N,
            rho_mediu = round(mean(rho_est), 3),
            rho_mediana = round(median(rho_est), 3),
            rho_min = round(min(rho_est), 3),
            rho_max = round(max(rho_est), 3),
            pct_semnificativ = round(100 * mean(abs(rho_t) > 1.96), 1),
            moran_y2_mediu = round(mean(moran_y2), 4),
            pct_moran_semnif = round(100 * mean(moran_p < .05), 1)),
        by = .(scenariu, rho_adevarat, estimator)])
cat("\nReading: on the row with a true rho of 0, any large value of\n")
cat("'rho_mediu' or a 'pct_semnificativ' well above 5 per cent means that\n")
cat("estimatorul inventeaza efect spatial acolo unde nu exista.\n")
