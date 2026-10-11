# Data

None of the survey microdata is included. Each file below has to be obtained from its provider and placed
under `data/raw/` with the folder names shown. The pipeline stops with an error if a file is missing.

## UK Cyber Security Breaches Survey

- Source: UK Data Service, study SN 8971, DOI [10.5255/UKDA-SN-8971-1](https://doi.org/10.5255/UKDA-SN-8971-1).
- Access: free registration and acceptance of the end user licence.
- Version used: Stata file `dcms_csbs_combined_2016-2022_banded_cost_data_only_v3_public.dta`.
- Expected path:
  `data/raw/UKDA-8971-stata/stata/stata13/dcms_csbs_combined_2016-2022_banded_cost_data_only_v3_public.dta`
- SHA-256: `e80c2ae670f524076037d61cb0a8cbf779e6c7bd3884e201bc48d52af8febf13`

## Istat ICT in enterprises, 2019

- Source: Istat public-use microdata for the 2019 survey on ICT use in enterprises (distributed as
  `ICT_2019_IT.zip`, with documentation in Italian).
- Expected path after unzipping:
  `data/raw/istat_ict_2019_public_use/MICRODATI/ICT_Microdati_2019.txt`
- SHA-256 of the text file: `01ceee6bd5ae1885affdc5f3202fffe9950758dddbc508d8a3fed1474fffb8ee`

## Cybersecurity Longitudinal Survey

- Source: UK Data Service, study SN 9356, DOI [10.5255/UKDA-SN-9356-1](https://doi.org/10.5255/UKDA-SN-9356-1).
- Access: free registration and acceptance of the end user licence.
- Expected path: `data/raw/UKDA-9356-stata/stata/stata13/csls_waves_1-4.dta`
- SHA-256: `ef4632c26fb1a02c3278bafb707ee11c1a8acb26c7a33d7c2dedcc272c71a306`

The panel scripts compare the file against this hash and stop if it differs. An intermediate CSV exported
from it (`results/panel/csls_raw_selected_columns.csv`) is checked against
`d236650ae3576470b82db13586975f8f28d1d67175864d660013de7b294faa87`.

## Checking your copy

```bash
shasum -a 256 data/raw/UKDA-8971-stata/stata/stata13/*.dta \
               data/raw/UKDA-9356-stata/stata/stata13/*.dta \
               data/raw/istat_ict_2019_public_use/MICRODATI/ICT_Microdati_2019.txt
```

On Linux use `sha256sum` instead of `shasum -a 256`.

## Analysis samples

| Sample | Rows | Restriction |
| - | - | - |
| UK CSBS | 2,859 | Business interviews, small and medium size bands (2 and 3), survey years 2017 to 2021, non-missing region, positive weight, valid values on the eight breadth and control items |
| Italy ICT | 16,809 | Size classes cl1 and cl2 in the Istat file, valid values on the analysis items and design variables |
| CSLS panel | 783 organization-wave rows, 318 organizations | Waves 1 to 4 |

The scripts assert the row counts, so a different data release stops the run.
