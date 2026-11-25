# Spatial Effects in Local Productivity Volatility

Replication code for the paper

> Mare, C., Achim, M.-A. and Muresan, G.-M. "Spatial Effects in Local Productivity
> Volatility: A Non-Dynamic Spatial ARCH Analysis of Romanian Localities."
> Submitted to *Empirical Economics*.

The paper applies the non-dynamic spatial ARCH model of Otto, Schmid and Garthoff
(2018) to firm-level administrative data rather than to asset prices, where this
literature has concentrated so far. More than 9 million local unit records from the
Romanian business register are aggregated into a fixed lattice of the 666 localities
with at least 50 firms, observed over 16 annual cross sections from 2008 to 2023.
Spatial dependence in the conditional variance of labour productivity is confirmed
in 7 of those years, concentrated in periods of aggregate disturbance.

## Data

Three inputs are needed and none of them are in this repository.

**INS business register, 2008 to 2023.** Sixteen files named `date2008.CSV` through
`date2023.CSV`, roughly 9 million local unit records in total. These are
administrative microdata provided by the National Institute of Statistics of Romania
and cannot be redistributed by the authors. Access can be requested directly from
INS.

**Administrative geometry.** `lau2021.gpkg` and `nuts3.geojson`, both from Eurostat
GISCO and both public. The SIRUTA locality code in the register matches the Eurostat
LAU identifier, which is what makes the register addressable in space.

**TopFirme2023.csv**, optional. Used only by `05_diagnostic_spatial.R` for the 2023
cross section with balance sheet items. Every other script runs without it.

The aggregated locality-year series produced by `03_spatial_data.R` are derived from
the register and are therefore subject to the same restriction. They are available
from the corresponding author on reasonable request.

## Repository layout

The code is organised as eight Quantlets, one folder each, following the convention
of the [QuantLet](https://quantlet.com) organisation. Every folder holds the code
that produces one part of the paper, the estimation output behind it and a
`Metainfo.txt` describing it.

| Quantlet | Contents | Where it appears in the paper |
| --- | --- | --- |
| `SARCHfirms_data` | `03_spatial_data.R`. Aggregates the register to locality level, matches it to the LAU geometry and builds the weight matrices for both lattices | Data section, construction of the 3,152 and 666 lattices |
| `SARCHfirms_moran_screening` | `05_diagnostic_spatial.R`. Moran screening of candidate outcome variables on levels and on squares | Data section, why labour productivity and not turnover growth |
| `SARCHfirms_estimation` | `08_sparch_fe.R`, the main estimation on two-way fixed effect residuals, with `06_sparch_productivitate.R` before filtering and `11_logsparch_fe.R` for the log form | **Table 1**, plus the functional form exercise |
| `SARCHfirms_sectoral` | `13_confuzie_sectoriala.R`. Geographic against sectoral neighbourhood, plus shocks orthogonalised on sectoral composition | **Table 2**, first three rows |
| `SARCHfirms_full_lattice` | `18_moran_retea_completa.R` on the full 3,152 locality lattice and `22_full_lattice_test.R`, the attempt to estimate on it | **Table 2**, last two rows, and the feasibility statement, no convergence after 3 hours against 84 seconds on the reduced lattice |
| `SARCHfirms_inference` | `07_verificare_estimator.R` and `16_verificare_n666.R`. Monte Carlo on the estimator, on the lattice actually used | Inference section, 18 per cent rejection at a nominal 5 per cent under the empirical distribution against 2 per cent under normality |
| `SARCHfirms_robustness` | `26_robustete_rapid.R`, the model under 5 neighbourhood definitions and 16 years each, and `30_selectie_W.R`, likelihood-based selection among them | **Table 3**, and the likelihood ranking reported as a non-result, the total log-likelihood spans less than 10 points out of roughly 9,300 |
| `SARCHfirms_figure` | `23_figura_retea.R`. Builds the two panel map of the lattices | **Figure 1** |

The numbering of the scripts follows the working project and is not consecutive.
Scripts that produce internal reports or that assemble the manuscript are not
included here.

The CSV files shipped inside each folder are the estimation output behind the tables
and the figure, so the numbers in the paper can be checked without access to the
register. These are model results, not data: 5 to 180 rows of year-level or
weight-matrix-level estimates, with no locality identifiers.

## Running the code

Edit `00_config.R` in the repository root, or export the four environment variables
it reads, so that `REG` points to the register files, `GEO` to the geometry and `OUT`
to a writable output folder. Nothing else in the code is machine specific. Each
script finds the config on its own, so it can be run either from the repository root
or from inside its own folder.

Run `SARCHfirms_data/03_spatial_data.R` first. It writes `panel_uat.rds` and
`spatial_weights.rds`, which every other script reads. After that the Quantlets are
independent of one another and can be run in any order.

Two scripts are slow. `26_robustete_rapid.R` performs 80 estimations and takes about
an hour. `30_selectie_W.R` uses a parallel cluster and takes longer still. Both write
their results incrementally, so a run that is interrupted does not lose the work
already done.

## Requirements

R 4.2.2 or later, with `data.table`, `sf`, `spdep`, `Matrix`, `spGARCH`, `ggplot2`,
`patchwork` and `parallel`.

The spatial ARCH estimator comes from the `spGARCH` package of Otto (2019). Note that
its analytic standard errors are not used in the paper. Inference rests on Moran's I
applied to the squared filtered shocks, for the reason documented by
`07_verificare_estimator.R`.

## Funding

Supported by a grant from the Ministry of Research, Innovation and Digitization,
CNCS-UEFISCDI, project number PN-IV-P2-2.1-TE-2023-1317, within PNCDI IV.

## Licence

MIT, see `LICENSE`.
