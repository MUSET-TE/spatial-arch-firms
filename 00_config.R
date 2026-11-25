# =====================================================================
# Paths. Every script sources this file and nothing else is machine
# specific. Set the three variables below, or export the matching
# environment variables, and the whole pipeline runs unchanged.
#
# A script can be run either from the repository root or from inside its
# own Quantlet folder. ROOT below resolves to the repository root in
# both cases, so the data and output folders are always the same ones.
# =====================================================================

ROOT <- if (file.exists("00_config.R")) "." else ".."

# Raw INS business register, 16 files named date2008.CSV ... date2023.CSV.
# These are administrative microdata and are not distributed with this
# repository. See the Data section of the README.
REG <- Sys.getenv("SARCH_REG", unset = file.path(ROOT, "data", "registru"))

# Geometry: lau2021.gpkg and nuts3.geojson, both from Eurostat GISCO.
GEO <- Sys.getenv("SARCH_GEO", unset = file.path(ROOT, "data", "geo"))

# Everything the scripts write.
OUT <- Sys.getenv("SARCH_OUT", unset = file.path(ROOT, "out"))

# Optional, used only by 05_diagnostic_spatial.R for the 2023 cross
# section with balance sheet items.
TOPFIRME <- Sys.getenv("SARCH_TOPFIRME",
                       unset = file.path(ROOT, "data", "TopFirme2023.csv"))

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
