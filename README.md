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

## Running the code

Edit `00_config.R`, or export the four environment variables it reads, so that `REG`
points to the register files, `GEO` to the geometry and `OUT` to a writable output
folder. Nothing else in the code is machine specific.

Run `03_spatial_data.R` first. It writes `panel_uat.rds` and `spatial_weights.rds`,
which every other script reads. After that the scripts are independent of one another
and can be run in any order.

Two scripts are slow. `26_robustete_rapid.R` performs 80 estimations and takes about
an hour. `30_selectie_W.R` uses a parallel cluster and takes longer still. Both write
their results incrementally, so a run that is interrupted does not lose the work
already done.

## What each script does

| Script | Purpose | Where it appears in the paper |
| --- | --- | --- |
| `03_spatial_data.R` | Aggregates the register to locality level, matches it to the LAU geometry and builds the weight matrices for both lattices | Data section, construction of the 3,152 and 666 lattices |
| `05_diagnostic_spatial.R` | Moran screening of candidate outcome variables on levels and on squares | Data section, why labour productivity and not turnover growth |
| `06_sparch_productivitate.R` | spARCH and SARspARCH on raw productivity, before fixed effects are removed | Shows the spatial autoregressive term in the mean before filtering |
| `07_verificare_estimator.R` | Monte Carlo on the estimator, known rho recovered under normal and under empirical errors | Inference section, 18 per cent rejection at a nominal 5 per cent under the empirical distribution against 2 per cent under normality |
| `08_sparch_fe.R` | Main estimation. spARCH and SARspARCH on two-way fixed effect residuals, 16 annual cross sections | **Table 1** |
| `11_logsparch_fe.R` | Log form of the variance equation, AIC comparison against the linear form | Functional form exercise |
| `13_confuzie_sectoriala.R` | Geographic against sectoral neighbourhood, plus shocks orthogonalised on sectoral composition | **Table 2**, first three rows |
| `16_verificare_n666.R` | Size and power of the estimator on the lattice actually used, not on a smaller one | Inference section |
| `18_moran_retea_completa.R` | Moran's I on the full 3,152 locality lattice, for comparison with the reduced one | **Table 2**, last two rows |
| `22_full_lattice_test.R` | Attempt to estimate the model on the full lattice, with sparse matrices and warm starts | Feasibility statement, no convergence after 3 hours against 84 seconds on the reduced lattice |
| `23_figura_retea.R` | Builds the two panel map of the lattices | **Figure 1** |
| `26_robustete_rapid.R` | The model under 5 neighbourhood definitions, 16 years each | **Table 3** |
| `30_selectie_W.R` | Likelihood-based selection among the weight definitions | Reported as a non-result, the total log-likelihood spans less than 10 points out of roughly 9,300 |

The numbering follows the working project and is not consecutive. Scripts that
produce internal reports or that assemble the manuscript are not included here.

## results/

The estimation output behind the tables and the figure, so the numbers in the paper
can be checked without access to the register. These are model results, not data.

`sparch_fe_p50_k5.csv` and `tabel_master_articol.csv` back Table 1,
`confuzie_sectoriala.csv` and `moran_retea_completa.csv` back Table 2,
`robustete_W.csv` with its summary backs Table 3 and `Figure1_lattice.png` is the
figure as it appears in the paper.

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
