# =====================================================================
# Data aggregated by locality (SIRUTA) and by county, plus spatial weight matrices
# for the spARCH models (Otto, Schmid & Garthoff)
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(sf); library(spdep); library(Matrix)
})
options(datatable.showProgress = FALSE)

source("00_config.R")

norm <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))

# ---- 1. aggregate the register by locality and year -------------------
reg_files <- sprintf("%s/date%d.CSV", REG, 2008:2023)
uat <- rbindlist(lapply(reg_files, function(f) {
  an <- as.integer(gsub("[^0-9]", "", basename(f)))
  d <- fread(f, sep = ",", header = TRUE, showProgress = FALSE,
             fill = TRUE, integer64 = "double")
  setnames(d, norm(names(d)))
  cn <- c(sir = "siruta_superior_subunitate", jud = "judet_subunitate",
          ca = "cifra_de_afaceri", sal = "nr_salariati")
  stopifnot(all(cn %in% names(d)))
  d <- d[, .(sir = get(cn[["sir"]]), jud = get(cn[["jud"]]),
             ca = get(cn[["ca"]]), sal = get(cn[["sal"]]))]
  suppressWarnings(d[, `:=`(sir = as.integer(sir), jud = as.integer(jud),
                            ca = as.numeric(ca), sal = as.numeric(sal))])
  d <- d[!is.na(sir) & !is.na(ca)]
  d[, .(an = an, ca = sum(ca, na.rm = TRUE), sal = sum(sal, na.rm = TRUE),
        firme = .N, jud = jud[1]), by = .(siruta = sir)]
}))
cat("Aggregated locality-year:", nrow(uat), "rows,", uniqueN(uat$siruta), "localities\n")

# ---- 2. LAU geometry and the contiguity matrix -------------------------
lau <- st_read(file.path(GEO, "lau2021.gpkg"), quiet = TRUE)
lau <- lau[lau$CNTR_CODE == "RO", ]
lau$siruta <- suppressWarnings(as.integer(lau$LAU_ID))
lau <- lau[!is.na(lau$siruta), ]
cat("LAU Romania:", nrow(lau), "poligoane\n")
cat("Register localities found in the geometry:",
    uniqueN(uat[siruta %in% lau$siruta, siruta]), "\n")

# the fixed lattice: localities with geometry AND positive turnover every year
ani <- sort(unique(uat$an))
compl <- uat[ca > 0 & siruta %in% lau$siruta, .N, by = siruta][N == length(ani), siruta]
cat("Localities with turnover > 0 in all", length(ani), "years:", length(compl), "\n")

lau_f <- lau[lau$siruta %in% compl, ]
lau_f <- lau_f[order(lau_f$siruta), ]

nb <- poly2nb(lau_f, queen = TRUE)
izolate <- which(card(nb) == 0)
cat("Localities with no neighbours (islands):", length(izolate), "\n")
if (length(izolate)) {
  # link them to the nearest centroid, so that W has no empty rows
  cent <- suppressWarnings(st_coordinates(st_centroid(st_geometry(lau_f))))
  for (i in izolate) {
    dd <- sqrt(rowSums((cent - cent[rep(i, nrow(cent)), ])^2)); dd[i] <- Inf
    j <- which.min(dd)
    nb[[i]] <- as.integer(j)
    nb[[j]] <- sort(unique(c(setdiff(nb[[j]], 0L), as.integer(i))))
  }
}
W_uat <- as(nb2mat(nb, style = "W", zero.policy = TRUE), "dgCMatrix")
cat("W locality:", nrow(W_uat), "x", ncol(W_uat),
    "| mean neighbours:", round(mean(card(nb)), 2), "\n")

siruta_ord <- lau_f$siruta

# ---- 3. the locality panel on the fixed lattice ------------------------
P <- uat[siruta %in% compl]
setkey(P, siruta, an)
P[, `:=`(prod = ca / pmax(sal, 1),
         lca = log(ca))]
setorder(P, siruta, an)
P[, r_ca := lca - shift(lca), by = siruta]                     # log turnover growth
P[, r_prod := log(prod) - shift(log(prod)), by = siruta]       # log productivity growth
P <- merge(P, data.table(siruta = siruta_ord, idx = seq_along(siruta_ord)),
           by = "siruta")
setorder(P, an, idx)
fwrite(P, file.path(OUT, "panel_uat_2008_2023.csv"))

# ---- 4. the same at county level (NUTS3) --------------------------------
nuts <- st_read(file.path(GEO, "nuts3.geojson"), quiet = TRUE)
nuts <- nuts[substr(nuts$NUTS_ID, 1, 2) == "RO", ]
nuts <- nuts[order(nuts$NUTS_ID), ]
cat("\nNUTS3 Romania:", nrow(nuts), "\n")
nb_j <- poly2nb(nuts, queen = TRUE)
W_jud <- as(nb2mat(nb_j, style = "W", zero.policy = TRUE), "dgCMatrix")
cat("W county:", nrow(W_jud), "x", ncol(W_jud),
    "| mean neighbours:", round(mean(card(nb_j)), 2), "\n")

jud_p <- uat[, .(ca = sum(ca), sal = sum(sal), firme = sum(firme)), by = .(jud, an)]
jud_p <- jud_p[!is.na(jud) & ca > 0]
setorder(jud_p, jud, an)
jud_p[, r_ca := log(ca) - shift(log(ca)), by = jud]
fwrite(jud_p, file.path(OUT, "panel_judet_2008_2023.csv"))

# ---- 5. save -----------------------------------------------------------
saveRDS(list(W_uat = W_uat, siruta = siruta_ord, nb_uat = nb,
             W_jud = W_jud, nuts_id = nuts$NUTS_ID, nuts_name = nuts$NAME_LATN,
             lau = st_drop_geometry(lau_f)[, c("siruta", "LAU_NAME", "POP_2021", "AREA_KM2")]),
        file.path(OUT, "spatial_weights.rds"))
saveRDS(P, file.path(OUT, "panel_uat.rds"))
cat("\nSaved: spatial_weights.rds, panel_uat.rds\n")
