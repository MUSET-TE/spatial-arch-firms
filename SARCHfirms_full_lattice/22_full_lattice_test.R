# =====================================================================
# Estimation on the full lattice of 3,152 localities.
# Differences from the failed attempt: the sparse matrix is built directly,
# starting values taken from the reduced lattice estimation, writing to
# disk as soon as a year finishes.
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(sf); library(spdep); library(Matrix); library(spGARCH)
})
source(if (file.exists("00_config.R")) "00_config.R" else file.path("..", "00_config.R"))
FIS <- file.path(OUT, "sparch_full_lattice.csv")

P <- as.data.table(readRDS(file.path(OUT, "panel_uat.rds")))
S <- readRDS(file.path(OUT, "spatial_weights.rds"))

# sparse matrix straight from the neighbour list, without the dense form
lw <- nb2listw(S$nb_uat, style = "W", zero.policy = TRUE)
# build the sparse matrix directly from the neighbour list
nb <- S$nb_uat
ii <- unlist(lapply(seq_along(nb), function(k) if (nb[[k]][1] != 0) rep(k, length(nb[[k]])) else integer(0)))
jj <- unlist(lapply(nb, function(v) if (v[1] != 0) v else integer(0)))
xx <- unlist(lapply(nb, function(v) if (v[1] != 0) rep(1/length(v), length(v)) else numeric(0)))
W  <- sparseMatrix(i = ii, j = jj, x = xx, dims = c(length(nb), length(nb)))
W  <- as(W, "dgCMatrix")
n  <- nrow(W)
cat("n =", n, "| elemente nenule:", length(W@x),
    "| density:", round(100*length(W@x)/n^2, 3), "per cent\n")
cat("W class:", class(W)[1], "| memory:", round(object.size(W)/1024^2, 1), "MB\n\n")

D <- copy(P)[siruta %in% S$siruta]
D[, ly := log(ca / pmax(sal, 1))]
q <- quantile(D$ly, c(.005, .995), na.rm = TRUE)
D[, ly := pmin(pmax(ly, q[1]), q[2])]
D[, e := ly - mean(ly), by = siruta]
D[, e := e - mean(e), by = an]

# starting values from the reduced lattice
START <- c(0.06, 0.20, 0.00)   # alpha, rho, intercept
cat("starting values: alpha =", START[1], ", rho =", START[2], "\n\n")

# order of the years: we start with the most promising
ANI <- c(2011, 2020, 2021, 2023, 2008, 2012, 2015, 2019,
         2013, 2014, 2016, 2018, 2022, 2009, 2010, 2017)

for (a in ANI) {
  d <- D[an == a][data.table(siruta = S$siruta), on = "siruta"]
  y <- d$e
  cat(sprintf("--- %d  (n=%d) ... ", a, length(y))); flush.console()
  t0 <- Sys.time()
  m <- try(qml.spARCH(y ~ 1, W = W, type = "spARCH", start = START), silent = TRUE)
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(m, "try-error")) {
    cat(sprintf("ERROR after %.1f min: %s\n", sec/60,
                substr(as.character(m), 1, 120)))
    row <- data.table(an = a, n = n, rho = NA_real_, rho_se = NA_real_,
                      rho_t = NA_real_, alpha = NA_real_, loglik = NA_real_,
                      minute = round(sec/60, 2), stare = "error")
  } else {
    co <- summary(m)$coef
    rr <- grep("rho", rownames(co)); aa <- grep("alpha", rownames(co))
    row <- data.table(an = a, n = n,
                      rho = co[rr,1], rho_se = co[rr,2], rho_t = co[rr,3],
                      alpha = co[aa,1], loglik = as.numeric(logLik(m)),
                      minute = round(sec/60, 2), stare = "ok")
    cat(sprintf("rho = %.4f (t %.2f)  [%.1f min]\n", co[rr,1], co[rr,3], sec/60))
  }
  fwrite(row, FIS, append = file.exists(FIS))
}
cat("\nDone. Results in", FIS, "\n")
