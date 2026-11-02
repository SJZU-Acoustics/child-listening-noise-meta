# Classroom noise and listening outcomes in children

Code for data analysis for the paper **Effects of classroom-relevant noise on task accuracy, response time, and perceived effort in school-aged children: A three-level meta-analysis**, by Yunting Cui, Xinhao Yang, Yang Liu and Yuan Zhang.

## Requirements

The release was checked with R 4.5.3 and Python 3.13.7 on macOS. Required R packages are `metafor` (5.0.1), `clubSandwich` (0.7.0), `wildmeta` (0.3.2), `ggplot2` (4.0.2), `patchwork` (1.3.2) and `ggrepel` (0.9.8). Python requires `openpyxl` (3.1.5). Figure typography uses Helvetica; font substitution on another operating system can affect appearance.

```r
install.packages(c("metafor", "clubSandwich", "wildmeta", "ggplot2", "patchwork", "ggrepel"))
```

Install `openpyxl` in your Python environment with `python -m pip install openpyxl`.

## File structure

- `run_all.R` runs the complete reproduction from the deposited workbook.
- `R/load_data.py` reads the workbook and prepares disposable CSV inputs. Public ACC and cross-task labels map to the equivalent internal IS and dual-task identifiers; the latter describes outcome pairing, not experimental design.
- `R/analysis_*.R` fits the models used by the manuscript and its supplementary analyses. The numbered modules retain their dependency order.
- `R/tables_*.R`, `R/fig*.R` and `R/prepare_ohat.py` generate the display items and supporting numerical summaries.
- `R/corpus_bibkeys.csv` connects public report identifiers to the manuscript's bibliography keys.
- `data/` holds the downloaded workbook, which is not duplicated in this repository.
- `output/` receives five main figures, four supplementary figures, three main table fragments, ten supplementary table fragments and numerical summaries. Supplementary Table S1 is the manuscript's descriptive search-strategy table and has no fitted result to regenerate. Table fragments refer to the manuscript's bibliography keys and are intended for inclusion in LaTeX documents.
- `.cache/` holds intermediate model outputs, run logs and the input trace; it can be deleted and regenerated.

## Usage

1. Download `Children_listening_effort_classroom_noise_meta-analysis.xlsx` from [Zenodo](https://doi.org/10.5281/zenodo.23891491) and place it in `data/`.
2. Open a terminal in the repository root and run:

```sh
Rscript run_all.R
```

If Python is not available as `python3`, set the `PYTHON` environment variable to your environment's Python executable before running. The workflow includes cluster wild bootstrap calculations with 1,000 resamples and seed 26; allow time for the full run. It stops on a failed module and reports the relevant log. Outputs are overwritten when rerun.

The workbook contains 269 contrasts from 24 studies in 25 reports, 47 aligned ACC–RT pairs, extracted group statistics, risk-of-bias assessments, search-flow counts and dictionaries. Positive aligned Hedges' g represents a more adverse outcome in noise. The search date is 30 March 2026. The deposited workbook is the complete data input; no source-study PDFs or licensed database exports are needed.

## Licence

Code is released under the [MIT Licence](LICENSE). The dataset is separately released on Zenodo under CC BY 4.0.
