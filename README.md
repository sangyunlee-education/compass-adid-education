# Compass Variable-Based Adjusted Difference-in-Differences in Educational Effect Research

This repository provides R code for the empirical example reported in the paper:

**Adjusted Difference-in-Differences Using Compass Variables: A Practical Guide to Application and Reporting**

The code covers sample construction, construction of pretest/posttest/compass variables, relevance-condition diagnostics, conditional local independence diagnostics using tetrad tests, and OLS, DID, and adjusted DID estimation.

The original Gyeonggi Education Longitudinal Study data are **not included**. Users should obtain the raw data from the official data provider and place the `.sav` files in a local `data/` directory.

## Repository structure

```text
R/functions.R
scripts/01_prepare_data.R
scripts/02_run_analysis.R
output/
README.md
.gitignore
```

## Expected raw data files

```text
data/y8STU.sav
data/y9STU.sav
data/Y10STU_학술대회.sav
```

## How to run

```r
source("scripts/01_prepare_data.R")
source("scripts/02_run_analysis.R")
```
