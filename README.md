# Digital breadth and technical cyber controls in SMEs

Analysis code for a study of how the breadth of a small or medium-sized enterprise's observable digital
interface relates to the breadth of the technical cyber controls it adopts. The code covers three data
sources:

- the UK Cyber Security Breaches Survey (CSBS), 2017 to 2021, from the UK Data Service (SN 8971);
- the Istat survey on ICT use in enterprises, 2019, in its public-use microdata release;
- the Cybersecurity Longitudinal Survey (CSLS), waves 1 to 4, from the UK Data Service (SN 9356).

The repository holds code only; the survey microdata are not redistributed. [DATA.md](DATA.md) explains how to obtain them. [REPRODUCIBILITY.md](REPRODUCIBILITY.md) lists the software versions, the
checks that the pipeline runs on itself, and the headline values to compare against.

## Quick start

```bash
pip install -r requirements.txt
Rscript scripts/install_r_packages.R
# place the three data sets under data/raw/ as described in DATA.md
bash scripts/run_pipeline.sh
python3 scripts/check_headline_results.py
```

The full pipeline takes about one minute on a laptop. All paths are relative to the repository root, so run
every command from there. Results are written to `results/`, which is not tracked by Git.

## Layout

| Path | Contents |
| - | - |
| `analysis/uk/` | UK dataset construction, survey design, models M0 to M3 |
| `analysis/uk/sensitivity/` | UK sensitivity analyses S1 to S7 and the attenuation comparisons D1 and D2 |
| `analysis/italy/` | Italy dataset construction, survey design, primary model |
| `analysis/italy/sensitivity/` | Italy sensitivity analyses S1 to S8 |
| `analysis/size_gap/` | Diagnostics for the firm-size difference between the two national surveys |
| `analysis/measurement/` | Item-level comparison of the UK and Italian digital-breadth and control items |
| `analysis/panel/` | CSLS panel export and panel models |
| `analysis/next_wave/` | Pair construction, configuration models and state-transition models for the CSLS |
| `analysis/reporting/` | Tables of values that feed the manuscript tables and figures |
| `scripts/` | Pipeline runner, R package installer, headline-value check |
| `tests/` | Structural tests that run without the data |

## Naming

In the code, the digital breadth score is called `footprint_score` (and `footprint_breadth` in model
output). It is the share of the measured digital-interface items that a firm reports, and corresponds to
digital breadth in the paper. Variable names that start with `i_` are binary item indicators.
`E`, `R` and `R3` in the CSLS scripts are the digital-breadth score, the four-item control score and the
three-item control score.

## Citation

Citation details will be added once the paper is published.

## Licence

Code is released under the MIT licence (see [LICENSE](LICENSE)). The survey data are governed by the terms of
their providers and are not part of this repository.
