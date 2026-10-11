# Reproducibility

## Environment

| Component | Version used |
| - | - |
| R | 4.6.1 |
| R packages | `survey` 4.5, `digest` 0.6.39, `jsonlite` 2.0.0 |
| Python | 3.14 |
| Python packages | `pandas` 3.0.6, `numpy` 2.5.3 |

Other versions have not been tested. `analysis/panel/panel_models.R` prints a warning if R is not the version above.

## Running

From the repository root, with the data in place (see [DATA.md](DATA.md)):

```bash
bash scripts/run_pipeline.sh
python3 scripts/check_headline_results.py
```

Set `PYTHON` or `RSCRIPT` to use interpreters other than `python3` and `Rscript`. The runner stops at the
first failing step (`set -euo pipefail`). Steps run in this order: UK dataset, design and models; UK
sensitivity analyses; Italy dataset, design, models and sensitivity analyses; firm-size diagnostics;
measurement check; CSLS panel; next-wave analyses; table and figure inputs.

## Checks built into the pipeline

- Source files are hashed and compared with the values in `DATA.md`.
- Analysis sample sizes and survey-design degrees of freedom are asserted before models are fitted.
- Each model script writes an `*_integrity_checks.txt` file (convergence, finite coefficients, positive
  standard errors, sample size match). A failed check stops the run.
- The next-wave state-transition models refuse to run unless their three input files match the expected
  SHA-256 values:
  - `configuration_adjacent_pairs.csv`: `c67e973444186d248783f61c355b092e31c139df8f243ebb8eb2dd780060e725`
  - `configuration_two_wave_gap_pairs.csv`: `3fc47646531167a90d82f9741549602d200cecfe553b3fe08bdb836534872f3d`
  - `configuration_retention_eligibility.csv`: `3ae45c4dfdd2dfefbf28aac585c9a0e9625720a26f887db52c8504e1ae56c1a0`

## Headline values

`scripts/check_headline_results.py` compares these values with the files in `results/`, using the rounding
shown:

| Quantity | UK CSBS | Italy ICT |
| - | - | - |
| Analysis N | 2,859 | 16,809 |
| Digital-breadth coefficient (logit scale) | 0.843 | 1.047 |
| Standard error | 0.174 | 0.074 |
| Average marginal effect | 0.149 | 0.251 |

## What can differ between machines

Manifest files record the paths and timestamps of the run, so they differ. Floating-point values can differ
in the last printed digit (up to about 1e-16 in the runs compared) when the numerical libraries differ. Values
reported to three decimals are unaffected.
