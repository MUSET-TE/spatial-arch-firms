# =====================================================================
# Is the spatial effect in fact a sectoral effect?
#
# Neighbouring localities often share the same industrial structure. A sectoral
# shock would then appear as spatial clustering. We test this in three ways:
#   (a) Moran on squared residuals, with a GEOGRAPHIC neighbourhood
#   (b) the same Moran, but with a SECTORAL neighbourhood (similar CAEN
#       composition, regardless of distance)
#   (c) the geographic Moran after the residual is orthogonalised on
#       the principal components of the sectoral structure
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(sf); library(spdep); library(Matrix)
})
options(datatable.showProgress = FALSE)
source("00_config.R")
norm <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))
NPC <- 10L   # principal components of the sectoral structure

# ---- 1. the geographic lattice, identical to the main estimation ------
P   <- readRDS(file.path(OUT, "panel_uat.rds"))
lau <- st_read(file.path(GEO, "lau2021.gpkg"), quiet = TRUE)
lau <- lau[lau$CNTR_CODE == "RO", ]
lau$siruta <- suppressWarnings(as.integer(lau$LAU_ID))
keep <- P[, .(m = min(firme)), by = siruta][m >= 50, siruta]
lf   <- lau[lau$siruta %in% keep, ]; lf <- lf[order(lf$siruta), ]
cent <- suppressWarnings(st_centroid(st_geometry(lf)))
nb_geo <- make.sym.nb(knn2nb(knearneigh(st_coordinates(cent), k = 5)))
lw_geo <- nb2listw(nb_geo, style = "W")
SIR <- lf$siruta
cat("Geographic lattice: n =", length(SIR), "\n")

# ---- 2. the sectoral structure of each locality ------------------------
cat("Reading the register for the CAEN structure...\n")
sect <- rbindlist(lapply(2008:2023, function(an) {
  f <- sprintf("%s/date%d.CSV", REG, an)
  if (!file.exists(f)) return(NULL)
  d <- fread(f, sep = ",", header = TRUE, showProgress = FALSE,
             fill = TRUE, integer64 = "double")
  setnames(d, norm(names(d)))
  d <- d[, .(sir = suppressWarnings(as.integer(siruta_superior_subunitate)),
             caen = suppressWarnings(as.integer(caen_rev2_mama)),
             sal = suppressWarnings(as.numeric(nr_salariati)))]
  d <- d[!is.na(sir) & !is.na(caen) & !is.na(sal) & sir %in% SIR]
  d[, div := caen %/% 100]
  d[, .(sal = sum(sal, na.rm = TRUE)), by = .(sir, div)]
}))
sect <- sect[, .(sal = sum(sal)), by = .(sir, div)]        # average structure
tot  <- sect[, .(tot = sum(sal)), by = sir]
sect <- merge(sect, tot, by = "sir")[tot > 0]
sect[, pondere := sal / tot]

M <- dcast(sect, sir ~ div, value.var = "pondere", fill = 0)
M <- M[data.table(sir = SIR), on = "sir"]
for (j in names(M)[-1]) set(M, which(is.na(M[[j]])), j, 0)
X <- as.matrix(M[, -1])
cat("Sectoral structure matrix:", nrow(X), "localities x", ncol(X), "CAEN divisions\n")

# ---- 3. SECTORAL neighbourhood: the 5 closest structures ---------------
#    Euclidean distance on the shares, with no reference to geography
nb_sec <- make.sym.nb(knn2nb(knearneigh(X, k = 5)))
lw_sec <- nb2listw(nb_sec, style = "W")
supr <- sum(sapply(seq_along(nb_geo), function(i) length(intersect(nb_geo[[i]], nb_sec[[i]]))))
cat("Sectoral lattice built. Pairs shared with the geographic one:", supr,
    "of", sum(card(nb_geo)), "\n")

# ---- 4. principal components of the sectoral structure -----------------
pc <- prcomp(X, center = TRUE, scale. = FALSE)
PCs <- pc$x[, seq_len(NPC), drop = FALSE]
cat("The first", NPC, "components explain",
    round(100 * sum(pc$sdev[1:NPC]^2) / sum(pc$sdev^2), 1), "per cent of the variation\n\n")

# ---- 5. fixed effect residuals, as in the main estimation --------------
D <- P[siruta %in% SIR]
D[, ly := log(ca / pmax(sal, 1))]
qq <- quantile(D$ly, c(.005, .995), na.rm = TRUE)
D[, ly := pmin(pmax(ly, qq[1]), qq[2])]
D[, e := ly - mean(ly), by = siruta]
D[, e := e - mean(e), by = an]

# ---- 6. the three tests, year by year ----------------------------------
rez <- rbindlist(lapply(sort(unique(D$an)), function(a) {
  d <- D[an == a][data.table(siruta = SIR), on = "siruta"]
  e <- d$e
  e_ort <- residuals(lm(e ~ PCs))          # cleaned of the sectoral structure
  g1 <- moran.test(e^2,     lw_geo, zero.policy = TRUE)
  g2 <- moran.test(e^2,     lw_sec, zero.policy = TRUE)
  g3 <- moran.test(e_ort^2, lw_geo, zero.policy = TRUE)
  data.table(an = a,
             I_geo      = unname(g1$estimate[1]), p_geo      = g1$p.value,
             I_sect     = unname(g2$estimate[1]), p_sect     = g2$p.value,
             I_geo_ort  = unname(g3$estimate[1]), p_geo_ort  = g3$p.value,
             R2_sector  = summary(lm(e ~ PCs))$r.squared)
}))
fwrite(rez, file.path(OUT, "confuzie_sectoriala.csv"))

cat("=== I_geo: geographic neighbourhood | I_sect: sectoral neighbourhood\n")
cat("=== I_geo_ort: geographic, after removing the sectoral structure\n\n")
print(rez[, .(an,
              I_geo = round(I_geo, 4), p_geo = round(p_geo, 3),
              I_sect = round(I_sect, 4), p_sect = round(p_sect, 3),
              I_geo_ort = round(I_geo_ort, 4), p_ort = round(p_geo_ort, 3),
              R2_sect = round(R2_sector, 3))])

cat("\nSignificant years (p below 0.05):\n")
cat("  geographic                        :", rez[p_geo < .05, .N], "of", nrow(rez), "\n")
cat("  sectoral                          :", rez[p_sect < .05, .N], "of", nrow(rez), "\n")
cat("  geographic after orthogonalisation:", rez[p_geo_ort < .05, .N], "of", nrow(rez), "\n")
cat("\nSaved:", file.path(OUT, "confuzie_sectoriala.csv"), "\n")
