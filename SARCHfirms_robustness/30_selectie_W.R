# =====================================================================
# 30_selectie_W.R
# Selection of the weight matrix on a likelihood criterion, as in
# standard practice: the same y, the same number of parameters, only the
# W specification differs, so the log-likelihoods are comparable.
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(sf); library(spdep); library(Matrix)
  library(spGARCH); library(parallel)
})
options(datatable.showProgress = FALSE)
source(if (file.exists("00_config.R")) "00_config.R" else file.path("..", "00_config.R"))
FIS <- file.path(OUT, "selectie_W.csv")
INC <- file.path(OUT, "selectie_W_partial.csv")
for (f in c(FIS, INC)) if (file.exists(f)) file.remove(f)

P   <- as.data.table(readRDS(file.path(OUT, "panel_uat.rds")))
lau <- st_read(file.path(GEO, "lau2021.gpkg"), quiet = TRUE)
lau <- lau[lau$CNTR_CODE == "RO", ]
lau$siruta <- suppressWarnings(as.integer(lau$LAU_ID))
keep <- P[, .(m = min(firme)), by = siruta][m >= 50, siruta]
lf   <- lau[lau$siruta %in% keep, ]; lf <- lf[order(lf$siruta), ]
pts  <- suppressWarnings(st_coordinates(st_centroid(st_geometry(lf))))

mk <- list()
mk[["Delaunay"]] <- tri2nb(pts)
mk[["Gabriel"]]  <- graph2nb(gabrielneigh(pts), sym = TRUE)
mk[["k = 4"]]    <- make.sym.nb(knn2nb(knearneigh(pts, k = 4)))
mk[["k = 5"]]    <- make.sym.nb(knn2nb(knearneigh(pts, k = 5)))
mk[["k = 8"]]    <- make.sym.nb(knn2nb(knearneigh(pts, k = 8)))
WL <- lapply(mk, function(nb) as(nb2mat(nb, style = "W", zero.policy = TRUE), "dgCMatrix"))
LW <- lapply(mk, function(nb) nb2listw(nb, style = "W", zero.policy = TRUE))

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
  y <- YS[[as.character(a)]]
  m <- try(qml.spARCH(y ~ 1, W = WL[[nm]], type = "spARCH"), silent = TRUE)
  if (inherits(m, "try-error"))
    return(data.table(W = nm, an = a, rho = NA_real_, rho_t = NA_real_,
                      loglik = NA_real_, AIC = NA_real_, BIC = NA_real_,
                      moran_rez2 = NA_real_, p_rez2 = NA_real_))
  co <- summary(m)$coef; rr <- grep("rho", rownames(co))
  s  <- summary(m)
  mr <- try(moran.test(as.numeric(residuals(m))^2, LW[[nm]], zero.policy = TRUE), silent = TRUE)
  data.table(W = nm, an = a, rho = co[rr, 1], rho_t = co[rr, 3],
             loglik = as.numeric(logLik(m)),
             AIC = if (!is.null(s$AIC)) as.numeric(s$AIC) else AIC(m),
             BIC = if (!is.null(s$BIC)) as.numeric(s$BIC) else NA_real_,
             moran_rez2 = if (inherits(mr, "try-error")) NA_real_ else unname(mr$estimate[1]),
             p_rez2     = if (inherits(mr, "try-error")) NA_real_ else mr$p.value)
}

cl <- makeCluster(16); on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({
  library(spGARCH); library(Matrix); library(spdep); library(data.table)}))
clusterExport(cl, c("job", "YS", "WL", "LW"), envir = environment())
t0 <- Sys.time()
res <- rbindlist(lapply(split(seq_len(nrow(job)), ceiling(seq_len(nrow(job)) / 16)),
  function(ix) { r <- rbindlist(parLapply(cl, ix, unul), fill = TRUE)
                 fwrite(r, INC, append = file.exists(INC))
                 cat("  batch done,", nrow(r), "estimations\n"); flush.console(); r }), fill = TRUE)
cat("time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minutes\n\n")
fwrite(res, FIS)

S <- res[, .(loglik_total = round(sum(loglik, na.rm = TRUE), 1),
             AIC_total    = round(sum(AIC, na.rm = TRUE), 1),
             BIC_total    = round(sum(BIC, na.rm = TRUE), 1),
             ani_castigati = NA_integer_,
             rho_mediu    = round(mean(rho, na.rm = TRUE), 3)), by = W]
best <- res[, .SD[which.max(loglik)], by = an][, .N, by = W]
S[, ani_castigati := best$N[match(W, best$W)]]
S[is.na(ani_castigati), ani_castigati := 0L]
setorder(S, -loglik_total)
print(S)
fwrite(S, file.path(OUT, "selectie_W_sinteza.csv"))
cat("\nSaved:", FIS, "\n")
