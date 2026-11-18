# =====================================================================
# log-spARCH (E-spARCH) on the fixed effect residuals
#   log h = alpha * 1 + rho * W * log(y^(2))
# The main specification of Otto, Schmid & Garthoff (2018),
# Spatial Statistics 26: 125-145. The non-dynamic, purely spatial case.
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(sf); library(spdep); library(Matrix)
  library(spGARCH); library(parallel)
})
source("00_config.R")
NCORE <- 4L

P   <- readRDS(file.path(OUT, "panel_uat.rds"))
lau <- st_read(file.path(GEO, "lau2021.gpkg"), quiet = TRUE)
lau <- lau[lau$CNTR_CODE == "RO", ]
lau$siruta <- suppressWarnings(as.integer(lau$LAU_ID))
keep <- P[, .(m = min(firme)), by = siruta][m >= 50, siruta]
lf   <- lau[lau$siruta %in% keep, ]; lf <- lf[order(lf$siruta), ]
cent <- suppressWarnings(st_centroid(st_geometry(lf)))
nb   <- make.sym.nb(knn2nb(knearneigh(st_coordinates(cent), k = 5)))
W    <- as(nb2mat(nb, style = "W"), "dgCMatrix")
lw   <- nb2listw(nb, style = "W")

D <- P[siruta %in% lf$siruta]
D[, ly := log(ca / pmax(sal, 1))]
qq <- quantile(D$ly, c(.005, .995), na.rm = TRUE)
D[, ly := pmin(pmax(ly, qq[1]), qq[2])]
D[, e := ly - mean(ly), by = siruta]
D[, e := e - mean(e), by = an]

ani <- sort(unique(D$an))
YS <- lapply(ani, function(a) D[an == a][data.table(siruta = lf$siruta), on = "siruta"]$e)
names(YS) <- as.character(ani)

unul <- function(a) {
  y <- YS[[as.character(a)]]
  out <- list(an = a)
  for (tp in c("log-spARCH", "spARCH")) {
    t0 <- Sys.time()
    m <- try(qml.spARCH(y ~ 1, W = W, type = tp), silent = TRUE)
    sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(m, "try-error")) {
      out[[paste0(tp, "_rho")]] <- NA_real_; next
    }
    co <- summary(m)$coef
    rr <- grep("rho", rownames(co)); aa <- grep("alpha", rownames(co))
    out[[paste0(tp, "_rho")]]   <- co[rr, 1]
    out[[paste0(tp, "_se")]]    <- co[rr, 2]
    out[[paste0(tp, "_t")]]     <- co[rr, 3]
    out[[paste0(tp, "_alpha")]] <- co[aa, 1]
    out[[paste0(tp, "_AIC")]]   <- AIC(m)
    out[[paste0(tp, "_sec")]]   <- round(sec, 1)
    # sign of non-convergence: the estimator stays at the starting values
    out[[paste0(tp, "_susp")]]  <- as.integer(abs(co[rr, 2] - 1) < 1e-6 |
                                              abs(co[rr, 3] - 1) < 1e-6 |
                                              abs(co[rr, 1] - 0.4) < 1e-4)
    r2 <- try(moran.test(as.numeric(residuals(m))^2, lw, zero.policy = TRUE), silent = TRUE)
    out[[paste0(tp, "_Irez2")]] <- if (inherits(r2, "try-error")) NA_real_ else unname(r2$estimate[1])
    out[[paste0(tp, "_prez2")]] <- if (inherits(r2, "try-error")) NA_real_ else r2$p.value
  }
  as.data.table(out)
}

cl <- makeCluster(NCORE)
on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({
  library(spGARCH); library(Matrix); library(spdep); library(data.table)}))
clusterExport(cl, c("W", "lw", "YS"), envir = environment())
R <- rbindlist(parLapply(cl, ani, unul), fill = TRUE)
setorder(R, an)
fwrite(R, file.path(OUT, "logsparch_fe.csv"))

cat("=== log-spARCH on fixed effect residuals, n =", nrow(W), "===\n")
print(R[, .(an,
            log_rho = round(`log-spARCH_rho`, 3), log_t = round(`log-spARCH_t`, 2),
            log_nec = `log-spARCH_susp`,
            lin_rho = round(spARCH_rho, 3), lin_t = round(spARCH_t, 2),
            dAIC = round(`log-spARCH_AIC` - spARCH_AIC, 1),
            I_rez2 = round(`log-spARCH_Irez2`, 4), p = round(`log-spARCH_prez2`, 3))])
cat("\nlog-spARCH significant and converged in",
    R[`log-spARCH_susp` == 0 & `log-spARCH_t` > 1.96, .N], "years of", nrow(R), "\n")
cat("linear form significant in", R[spARCH_t > 1.96, .N], "years\n")
cat("neconvergente log:", R[`log-spARCH_susp` == 1, .N], "\n")
cat("\nSaved:", file.path(OUT, "logsparch_fe.csv"), "\n")
