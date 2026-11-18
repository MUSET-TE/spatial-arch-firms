# =====================================================================
# Moran on squared residuals, the FULL lattice of 3,152 localities
# (queen contiguity), for comparison with the lattice of 666
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(sf); library(spdep)})
source("00_config.R")

P <- readRDS(file.path(OUT, "panel_uat.rds"))
S <- readRDS(file.path(OUT, "spatial_weights.rds"))
lw <- nb2listw(S$nb_uat, style = "W", zero.policy = TRUE)

D <- copy(P)
D[, ly := log(ca / pmax(sal, 1))]
q <- quantile(D$ly, c(.005, .995), na.rm = TRUE)
D[, ly := pmin(pmax(ly, q[1]), q[2])]
D[, e := ly - mean(ly), by = siruta]
D[, e := e - mean(e), by = an]

R <- rbindlist(lapply(sort(unique(D$an)), function(a) {
  d <- D[an == a][data.table(siruta = S$siruta), on = "siruta"]
  t <- moran.test(d$e^2, lw, zero.policy = TRUE)
  data.table(an = a, n = length(S$siruta),
             I_full = unname(t$estimate[1]), p_full = t$p.value)
}))
fwrite(R, file.path(OUT, "moran_retea_completa.csv"))
cat("n =", length(S$siruta), "| mean neighbours:", round(mean(card(S$nb_uat)), 2), "\n")
cat("significant years:", R[p_full < .05, .N], "of", nrow(R), "\n")
print(R[, .(an, I = round(I_full, 4), p = round(p_full, 4))])
