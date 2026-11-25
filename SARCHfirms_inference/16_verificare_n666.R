# =====================================================================
# Size and power of the estimator on the lattice ACTUALLY
# used in the estimations (n = 666), not on a smaller one.
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(sf); library(spdep); library(Matrix)
  library(spGARCH); library(parallel)
})
set.seed(20260929)
source(if (file.exists("00_config.R")) "00_config.R" else file.path("..", "00_config.R"))
REP   <- 50L
NCORE <- 16L

P   <- readRDS(file.path(OUT, "panel_uat.rds"))
lau <- st_read(file.path(GEO, "lau2021.gpkg"), quiet = TRUE)
lau <- lau[lau$CNTR_CODE == "RO", ]
lau$siruta <- suppressWarnings(as.integer(lau$LAU_ID))
keep <- P[, .(m = min(firme)), by = siruta][m >= 50, siruta]
lf   <- lau[lau$siruta %in% keep, ]; lf <- lf[order(lf$siruta), ]
cent <- suppressWarnings(st_centroid(st_geometry(lf)))
nb   <- make.sym.nb(knn2nb(knearneigh(st_coordinates(cent), k = 5)))
W    <- as(nb2mat(nb, style = "W"), "dgCMatrix")
n    <- nrow(W)
cat("Lattice identical to the estimation one: n =", n, "\n\n")

# the empirical distribution of the real residuals, for a realistic null
D <- P[siruta %in% lf$siruta]
D[, ly := log(ca / pmax(sal, 1))]
qq <- quantile(D$ly, c(.005, .995), na.rm = TRUE)
D[, ly := pmin(pmax(ly, qq[1]), qq[2])]
D[, e := ly - mean(ly), by = siruta]
D[, e := e - mean(e), by = an]
E_REAL <- D$e

fit_rho <- function(y, W) {
  m <- try(qml.spARCH(y ~ 1, W = W, type = "spARCH"), silent = TRUE)
  if (inherits(m, "try-error")) return(NA_real_)
  co <- summary(m)$coef; r <- grep("rho", rownames(co))
  if (!length(r)) return(c(NA_real_, NA_real_))
  c(co[r[1], 1], co[r[1], 3])
}

cl <- makeCluster(NCORE)
on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({library(spGARCH); library(Matrix)}))
clusterExport(cl, c("fit_rho"), envir = environment())
clusterSetRNGStream(cl, 424242)

scen <- list(
  list(nm = "null, normal noise (rho = 0)",            gen = "norm", rho = 0),
  list(nm = "nul, reziduuri reale amestecate (rho = 0)", gen = "real", rho = 0),
  list(nm = "alternative, simulated spARCH with rho = 0.20", gen = "sim",  rho = 0.20),
  list(nm = "alternative, simulated spARCH with rho = 0.27", gen = "sim",  rho = 0.27)
)

rez <- list()
for (i in seq_along(scen)) {
  s <- scen[[i]]
  t0 <- Sys.time()
  out <- parLapply(cl, seq_len(REP), function(b, s, W, n, E_REAL) {
    y <- switch(s$gen,
      norm = rnorm(n),
      real = sample(E_REAL, n),
      sim  = as.numeric(sim.spARCH(n = n, rho = s$rho, alpha = 1,
                                   W = W, type = "spARCH")))
    y <- y - mean(y)
    fit_rho(y, W)
  }, s = s, W = W, n = n, E_REAL = E_REAL)
  M <- do.call(rbind, out)
  M <- M[is.finite(M[, 1]), , drop = FALSE]
  rez[[i]] <- data.table(scenariu = s$nm, rho_adevarat = s$rho, n = n,
                         repetari = nrow(M),
                         rho_mediu = mean(M[, 1]), rho_q95 = quantile(M[, 1], .95),
                         rho_max = max(M[, 1]),
                         pct_semnificativ = 100 * mean(abs(M[, 2]) > 1.96, na.rm = TRUE),
                         minute = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
  cat(sprintf("  %-42s mean rho %.4f | q95 %.4f | max %.4f | signif %5.1f%%  [%.1f min]\n",
              s$nm, mean(M[, 1]), quantile(M[, 1], .95), max(M[, 1]),
              100 * mean(abs(M[, 2]) > 1.96, na.rm = TRUE),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
R <- rbindlist(rez)
fwrite(R, file.path(OUT, "verificare_n666.csv"))
cat("\nSaved:", file.path(OUT, "verificare_n666.csv"), "\n")
