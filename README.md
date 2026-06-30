# Alpha-Gal Risk Mapping in Illinois
Bayesian conditional autoregressive (CAR) spatial modeling of alpha-gal syndrome (AGS) proxy risk across Illinois counties, 2019–2022.

## Overview
No confirmed AGS cases exist in Illinois. This repository contains the analysis code for a county-level proxy risk score combining:
- Lone star tick abundance (county-level)
- Ehrlichiosis case counts (IDPH)
- Tick establishment status (CDC classification)

using a Leroux CAR prior fit via the `CARBayes` R package, with posterior predictive checks and sensitivity analyses across alternative weighting schemes.

## Manuscript
Hussain A, Mateus-Pinilla N, Smith RL. County-Level Risk Mapping of Alpha-Gal Syndrome Using a Bayesian Proxy Approach, Illinois, 2019–2022. [Under review]

## Data availability
Tick abundance data: [link to your previously published paper's data repository]
Ehrlichiosis data: available from IDPH upon request; not redistributed here due to data use restrictions.

## Repository structure
- `R/` — analysis scripts (data joining, CAR models, sensitivity analyses, figures)
- `output/` — generated figures and tables (not version-controlled raw outputs)

## Contact
Abrar Hussain, University of Illinois Urbana-Champaign — abrarh2@illinois.edu
