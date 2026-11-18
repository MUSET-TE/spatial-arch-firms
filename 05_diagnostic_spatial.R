# =====================================================================
# Is there spatial dependence to exploit? Moran's I on levels and on squares,
# for several candidate variables and several lattices.
# Includes the 2023 cross section with balance sheet items (TopFirme plus register).
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(sf); library(spdep)
})
options(datatable.showProgress = FALSE)
source("00_config.R")
norm <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))

P   <- readRDS(file.path(OUT, "panel_uat.rds"))
lau <- st_read(file.path(GEO, "lau2021.gpkg"), quiet = TRUE)
lau <- lau[lau$CNTR_CODE == "RO", ]
lau$siruta <- suppressWarnings(as.integer(lau$LAU_ID))

retea <- function(prag, k = 5) {
  keep  <- P[, .(minf = min(firme)), by = siruta][minf >= prag, siruta]
  lf    <- lau[lau$siruta %in% keep, ]
  lf    <- lf[order(lf$siruta), ]
  cent  <- suppressWarnings(st_centroid(st_geometry(lf)))
  nb    <- knn2nb(knearneigh(st_coordinates(cent), k = k))
  list(siruta = lf$siruta, lw = nb2listw(nb, style = "W"))
}

mi <- function(x, lw) {
  x <- x[is.finite(x)]
  t1 <- moran.test(x, lw, zero.policy = TRUE)
  t2 <- moran.test((x - mean(x))^2, lw, zero.policy = TRUE)
  c(I = unname(t1$estimate[1]), p = t1$p.value,
    I2 = unname(t2$estimate[1]), p2 = t2$p.value)
}

# =====================================================================
# A. variables from the register, 2009-2023
# =====================================================================
cat("=== A. Register: turnover growth and productivity growth ===\n")
A <- rbindlist(lapply(c(50, 200, 500), function(prag) {
  R <- retea(prag)
  D <- P[siruta %in% R$siruta]
  rbindlist(lapply(sort(unique(D[!is.na(r_ca), an])), function(a) {
    d <- D[an == a][data.table(siruta = R$siruta), on = "siruta"]
    if (any(is.na(d$r_ca))) return(NULL)
    w <- function(v) { q <- quantile(v, c(.01, .99)); pmin(pmax(v, q[1]), q[2]) }
    m1 <- mi(w(d$r_ca), R$lw)
    m2 <- mi(w(d$r_prod), R$lw)
    data.table(prag = prag, n = length(R$siruta), an = a,
               I_rca = m1[["I"]], p_rca = m1[["p"]],
               I2_rca = m1[["I2"]], p2_rca = m1[["p2"]],
               I_rprod = m2[["I"]], I2_rprod = m2[["I2"]], p2_rprod = m2[["p2"]])
  }))
}))
print(A[, .(ani = .N,
            I_rca = round(mean(I_rca), 4), I2_rca = round(mean(I2_rca), 4),
            ani_semnif_I2 = sum(p2_rca < .05),
            I_rprod = round(mean(I_rprod), 4), I2_rprod = round(mean(I2_rprod), 4),
            ani_semnif_I2p = sum(p2_rprod < .05)), by = .(prag, n)])
fwrite(A, file.path(OUT, "moran_registru.csv"))

# =====================================================================
# B. the 2023 cross section with balance sheet items
# =====================================================================
cat("\n=== B. TopFirme 2023 plus geography from the register ===\n")
tf <- fread(TOPFIRME, showProgress = FALSE,
            integer64 = "double")
setnames(tf, norm(names(tf)))
tf <- tf[an == 2023]
rg <- fread(file.path(REG, "date2023.CSV"), sep = ",", header = TRUE,
            showProgress = FALSE, fill = TRUE, integer64 = "double")
setnames(rg, norm(names(rg)))
rg <- rg[, .(cod_fiscal = suppressWarnings(as.integer(cod_fiscal)),
             siruta = suppressWarnings(as.integer(siruta_superior_subunitate)))]
rg <- unique(rg[!is.na(cod_fiscal) & !is.na(siruta)], by = "cod_fiscal")
F <- merge(tf, rg, by.x = "fiscal_code", by.y = "cod_fiscal")
cat("Firms with balance sheet and location:", nrow(F), "\n")

F[, active_totale := active_imobilizate + active_circulante]
G <- F[, .(firme = .N,
           ca = sum(cifra_afaceri, na.rm = TRUE),
           pnet = sum(profit_net, na.rm = TRUE),
           active = sum(active_totale, na.rm = TRUE),
           datorii = sum(datorii, na.rm = TRUE),
           cap = sum(capitaluri_total, na.rm = TRUE),
           sal = sum(numar_angajati, na.rm = TRUE)), by = siruta]
G <- G[ca > 0 & active > 0 & sal > 0]
G[, `:=`(roa = pnet / active,
         marja = pnet / ca,
         lev = datorii / active,
         lprod = log(ca / sal))]
cat("Localities with aggregated data:", nrow(G), "\n")

B <- rbindlist(lapply(c(50, 200, 500), function(prag) {
  R <- retea(prag)
  d <- G[data.table(siruta = R$siruta), on = "siruta"]
  if (any(is.na(d$roa))) { d <- d[!is.na(roa)]; return(NULL) }
  w <- function(v) { q <- quantile(v, c(.01, .99), na.rm = TRUE); pmin(pmax(v, q[1]), q[2]) }
  rbindlist(lapply(c("roa","marja","lev","lprod"), function(v) {
    m <- mi(w(d[[v]]), R$lw)
    data.table(prag = prag, n = nrow(d), variabila = v,
               I = round(m[["I"]], 4), p = round(m[["p"]], 4),
               I_patrate = round(m[["I2"]], 4), p_patrate = round(m[["p2"]], 4))
  }))
}))
print(B)
fwrite(B, file.path(OUT, "moran_topfirme2023.csv"))
cat("\nSaved: moran_registru.csv, moran_topfirme2023.csv\n")
