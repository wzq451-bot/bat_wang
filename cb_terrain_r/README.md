# R reproduction of the CB bat terrain analyses

Version: 0.1.0; prepared 2026-10-30. Status: author review and archive metadata pending.

This project independently refits the current seven-bat autumn terrain analyses in R. Python results are immutable numerical references, not inputs to model fitting. No Python process is called by `run_all.R`.

## Run

1. Install R (tested: 4.3.3 on Linux x86_64) and `jsonlite` (tested: 1.8.8).
2. Obtain the author-held reproducibility bundle and extract its `data/` and `reference/` directories into this project. The public code branch does not contain individual locations, device identifiers or raw tracking records. A public data accession is still pending.
3. In a terminal run:

```sh
Rscript -e 'install.packages("jsonlite", repos="https://cloud.r-project.org")'
Rscript run_all.R
```

The run writes numerical results, CSV tables, audit PDFs and a comparison report. It stops on failed numerical checks. Runtime depends on hardware, particularly the 22-parameter conditional-logistic sensitivity model. No network access is needed during analysis after installing jsonlite and supplying the bundle.

For Arial audit plots, install your licensed Arial fonts, ensure R's Cairo/fontconfig can find them, and set `PLOT_FONT=Arial` before running. Font files are not distributed. Default graphics use the system sans-serif font. These PDFs are statistical audit plots, not replacements for the approved map/3D figure compositions.

## Analysis coverage

| Module | Recomputed analyses |
|---|---|
| `reconstruct.R` | Original-record accumulation durations; 3,062 centered-window straightness, path lengths, spans and step counts |
| `vertical.R` | 3,193 fixes, 21 segments, six bats; within- and between-segment slopes; individual-segment slopes; 2,000 circular shifts and one-sided Monte Carlo probability |
| `selection.R` | 3,360 choice strata with one used + 40 available endpoints; additive and height-interaction conditional logistic models; leave-one-bat-out fits; within-bat height; nonnegative height; valid heights; GPS intervals; time-of-night, before/after midnight; bat-specific terrain slopes; individual fits; candidate-terrain variability; conditional TPI predictions |
| `movement.R` | 46 linear-model fits, including equal-bat weights, segment intercepts, nonoverlap, duration, clock time, GPS quality, count/rate thresholds, leave-one-bat-out and synchronized whole windows; 14 individual slopes; activity interaction; fixed-knot natural cubic curves and uncertainty; exact 128-sign wild-cluster tests |
| `tables.R` | Per-bat counts, candidate terrain support, unrounded coefficient source tables |
| `validate.R` | Key-aligned Python/R comparisons, sample counts, missing-value patterns and documented correction |
| `plots.R` | Statistical audit PDFs derived from R outputs |

## Numerical contract

The absolute tolerance for finite numerical outputs is 1e-6; counts and validity masks are also checked. Row keys must be unique and every Python reference row must be covered. Floating point equality is not claimed. Models use the same analytical samples, predictor scaling (sample SD, ddof=1), likelihood, fixed effects, CR1 small-sample covariance correction and t reference with G-1 degrees of freedom. Conditional-logistic optimization uses an analytic gradient, information matrix, Newton refinement, a maximum-mean-score check and positive-information check. No rounding is used before fitting.

The same integer circular-shift schedule is supplied to both implementations. A seed alone does not synchronize NumPy and R random generators. The schedule archives the original 2,000 replicates for exact experiment replication. Wild-cluster inference enumerates every sign assignment and does not depend on a PRNG.

**Documented correction:** in the original Python coefficient export, `sqrt(max(0, NaN))` turned an undefined one-cluster SE into zero for 28 individual-fit coefficient rows. R correctly reports NA for these cluster SEs, intervals and P values. The comparison explicitly checks and records this exception. Individual point estimates and model-based SEs are still compared. Do not interpret a single-bat clustered SE of zero as certainty. The final Python submission table had already replaced these values with NA; this corrects the intermediate export.

## Reproduction boundary

This project starts from the archived analytical tables. It does not claim to regenerate original GPS quality-control decisions, DEM extraction, terrain metrics or the original random available endpoints from raw rasters. DEM source/version, analytical TPI/TRI neighborhood definitions, original candidate-generation rules and initial cohort exclusions require the original preprocessing code/metadata before a raw-data-to-paper claim can be made. Terrain variables and choice sets are retained, not simulated or substituted.

Times are parsed with `tz='UTC'` to preserve existing naive clock labels and compute durations deterministically; this is not evidence that device timestamps were acquired in UTC. Confirm device time zone in Methods. Movement activity is log(count / original-record duration in minutes), not log(count / GPS gap). Spatial model covariates retain the original GPS interval and window-span adjustments.

`reference/` contains fixed Python exports for validation only. Source tables are machine-readable analysis outputs; journal-formatted tables and legends remain in the manuscript figure package. The repository is an analytical reproduction of the current manuscript version, not every superseded exploratory fit.

## Publication checklist

- Confirm contributor names, ORCIDs, title and software license before a formal release.
- Resolve upstream preprocessing metadata and data-sharing access conditions.
- Archive the approved release using Zenodo or another durable software archive and insert its actual DOI into `CITATION.md` and STAR Methods.
- A GitHub URL/commit is a versioned code reference, not a DOI or a data accession. No DOI has been fabricated here.
