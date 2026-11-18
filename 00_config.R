# =====================================================================
# Paths. Every script sources this file and nothing else is machine
# specific. Set the three variables below, or export the matching
# environment variables, and the whole pipeline runs unchanged.
# =====================================================================

# Raw INS business register, 16 files named date2008.CSV ... date2023.CSV.
# These are administrative microdata and are not distributed with this
# repository. See the Data section of the README.
REG <- Sys.getenv("SARCH_REG", unset = "data/registru")

# Geometry: lau2021.gpkg and nuts3.geojson, both from Eurostat GISCO.
GEO <- Sys.getenv("SARCH_GEO", unset = "data/geo")

# Everything the scripts write.
OUT <- Sys.getenv("SARCH_OUT", unset = "out")

# Optional, used only by 05_diagnostic_spatial.R for the 2023 cross
# section with balance sheet items.
TOPFIRME <- Sys.getenv("SARCH_TOPFIRME", unset = "data/TopFirme2023.csv")

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
