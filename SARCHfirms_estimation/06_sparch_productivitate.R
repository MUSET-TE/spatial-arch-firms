# =====================================================================
# spARCH and SARspARCH on labour productivity at locality level
# Annual cross sections, 2008-2023 (non-dynamic model)
#
#   spARCH   : y = diag(h)^(1/2) eps ,  h = alpha + rho * W y^(2)
#   SARspARCH: y = lambda * B y + X beta + e , e = diag(h)^(1/2) eps
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(sf); library(spdep); library(Matrix); library(spGARCH)
})
options(datatable.showProgress = FALSE)
source(if (file.exists("00_config.R")) "00_config.R" else file.path("..", "00_config.R"))
a_    <- commandArgs(trailingOnly = TRUE)
PRAG  <- if (length(a_) >= 1) as.integer(a_[1]) else 50L
KNN   <- if (length(a_) >= 2) as.integer(a_[2]) else 5L
FIS   <- file.path(OUT, sprintf("sparch_prod_p%d_k%d.csv", PRAG, KNN))

P   <- readRDS(file.path(OUT, "panel_uat.rds"))
lau <- st_read(file.path(GEO, "lau2021.gpkg"), quiet = TRUE)
lau <- lau[lau$CNTR_CODE == "RO", ]
lau$siruta <- suppressWarnings(as.integer(lau$LAU_ID))
keep <- P[, .(m = min(firme)), by = siruta][m >= PRAG, siruta]
lf   <- lau[lau$siruta %in% keep, ]
lf   <- lf[order(lf$siruta), ]
cent <- suppressWarnings(st_centroid(st_geometry(lf)))
nb   <- make.sym.nb(knn2nb(knearneigh(st_coordinates(cent), k = KNN)))
lw   <- nb2listw(nb, style = "W")
W    <- as(nb2mat(nb, style = "W"), "dgCMatrix")
n    <- nrow(W)
cat(sprintf("Lattice: n = %d localities (threshold %d firms), k = %d neighbours\n\n", n, PRAG, KNN))

scrie <- function(dt) {
  fwrite(dt, FIS, append = file.exists(FIS))
}
if (file.exists(FIS)) file.remove(FIS)

extrage <- function(m, nume) {
  co <- summary(m)$coef
  g  <- function(pat) { r <- grep(pat, rownames(co)); if (!length(r)) return(rep(NA_real_, 3)); co[r[1], 1:3] }
  al <- g("alpha"); rh <- g("rho"); lm_ <- g("lambda")
  list(alpha = al[1], alpha_se = al[2],
       rho = rh[1], rho_se = rh[2], rho_t = rh[3],
       lambda = lm_[1], lambda_se = lm_[2], lambda_t = lm_[3],
       logLik = as.numeric(logLik(m)), AIC = AIC(m))
}

for (a in sort(unique(P$an))) {
  d <- P[an == a][data.table(siruta = lf$siruta), on = "siruta"]
  y <- log(d$ca / pmax(d$sal, 1))
  q <- quantile(y, c(.01, .99), na.rm = TRUE)
  y <- pmin(pmax(y, q[1]), q[2])
  y0 <- y - mean(y)

  mi1 <- moran.test(y0, lw, zero.policy = TRUE)
  mi2 <- moran.test(y0^2, lw, zero.policy = TRUE)

  for (tp in c("spARCH", "log-spARCH")) {
    # --- 1. plain spARCH on the centred variable
    t0 <- Sys.time()
    m1 <- try(qml.spARCH(y0 ~ 1, W = W, type = tp), silent = TRUE)
    s1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    # --- 2. SARspARCH: spatial autoregressive mean plus spARCH variance
    t0 <- Sys.time()
    m2 <- try(qml.SARspARCH(y0 ~ 1, B = W, W = W, type = tp), silent = TRUE)
    s2 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    for (nm in c("spARCH", "SARspARCH")) {
      m   <- if (nm == "spARCH") m1 else m2
      sec <- if (nm == "spARCH") s1 else s2
      if (inherits(m, "try-error")) {
        cat(sprintf("  %d %-11s %-10s ERROR [%.0f s]\n", a, tp, nm, sec)); next
      }
      e  <- extrage(m, nm)
      mr <- try(moran.test(as.numeric(residuals(m))^2, lw, zero.policy = TRUE), silent = TRUE)
      row <- data.table(an = a, n = n, prag = PRAG, k = KNN, varianta = tp, model = nm,
                        alpha = e$alpha, alpha_se = e$alpha_se,
                        rho = e$rho, rho_se = e$rho_se, rho_t = e$rho_t,
                        lambda = e$lambda, lambda_se = e$lambda_se, lambda_t = e$lambda_t,
                        logLik = e$logLik, AIC = e$AIC,
                        moran_y = unname(mi1$estimate[1]), moran_y_p = mi1$p.value,
                        moran_y2 = unname(mi2$estimate[1]), moran_y2_p = mi2$p.value,
                        moran_rez2 = if (inherits(mr, "try-error")) NA_real_ else unname(mr$estimate[1]),
                        moran_rez2_p = if (inherits(mr, "try-error")) NA_real_ else mr$p.value,
                        secunde = round(sec, 1))
      scrie(row)
      cat(sprintf("  %d %-11s %-10s rho = %7.4f (t %6.2f)  lambda = %7.4f (t %6.2f)  AIC = %8.1f  [%.0f s]\n",
                  a, tp, nm, e$rho, e$rho_t,
                  ifelse(is.na(e$lambda), NA_real_, e$lambda),
                  ifelse(is.na(e$lambda_t), NA_real_, e$lambda_t), e$AIC, sec))
    }
  }
}
cat("\nSaved:", FIS, "\n")
