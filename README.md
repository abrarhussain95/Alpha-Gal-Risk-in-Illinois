
# Alpha-Gal Syndrome Exposure-Risk Mapping in Illinois

Analysis code for **“A Spatially Validated Machine-Learning Framework for Mapping Alpha-Gal Syndrome Exposure Potential in Illinois.”**

## Overview

This project develops county-level indices of **potential exposure to alpha-gal syndrome (AGS)** across all **102 Illinois counties**. Because geographically resolved, clinically confirmed AGS case data were not available for this analysis, the indices combine ecological predictions and an epidemiological proxy; they do not estimate confirmed AGS incidence or individual risk.

The analysis uses **2022–2024 Illinois tick surveillance and ecological data** to predict the proportion of collected ticks identified as *Amblyomma americanum* (lone star ticks).

## Methods

- **Ecological predictors:** Five bioclimatic variables (BIO2, BIO8, BIO9, BIO15, BIO18), forest cover, wetland cover, and white-tailed deer habitat suitability.

- **Model comparison:** Five baseline models (beta-binomial, elastic net, generalized additive model, random forest, and XGBoost) and three BYM2 spatial extensions (beta-binomial–BYM2, random forest–BYM2, and XGBoost–BYM2).

- **Validation:** Five-fold geographic cross-validation across 102 Illinois counties.

- **Primary AGS indices:** 70% spatially cross-validated predicted lone star tick proportion + 30% normalized ehrlichiosis incidence (2022–2024), calculated separately from random forest and XGBoost.

- **Sensitivity analyses:** Alternative ecological-to-ehrlichiosis weights (50:50, 60:40, 80:20), a wildland–urban interface (WUI) population component (70:20:10), and neighboring-county spatial averaging.

- **External evaluation:** Comparison with CDC suspected-AGS laboratory-testing categories (2017–2022) for 97 Illinois counties with available testing information.

## Analysis Scripts

The R Markdown scripts cover:

1. **Outcome preparation:** Calculate and map the observed lone star tick proportion from surveillance collections.

2. **Feature extraction and preparation:** Summarize climate, land cover, and deer habitat suitability by county.

3. **Model preparation:** Combine outcomes and predictors, screen collinearity, and prepare the modeling dataset.

4. **WUI extraction:** Calculate county-level proportions of residents living in WUI areas.

5. **Machine-learning model development:** Fit and evaluate the eight models; interpret random forest and XGBoost; construct, map, and evaluate AGS exposure-risk indices.

The analyses use R packages including `sf`, `terra`, `dplyr`, `glmmTMB`, `glmnet`, `mgcv`, `ranger`, `xgboost`, `INLA`, `spdep`, and `ordinal`, plus plotting and interpretation packages referenced in the scripts.

**Reproducibility note:** The scripts currently contain local file paths and require the input datasets described below. Update those paths to your own directory structure before running the workflows.

## Data Availability

Analysis code is provided in this repository. Tick surveillance and human ehrlichiosis data were obtained from the **Illinois Department of Public Health (IDPH)**. Climate, land cover, deer habitat, and WUI inputs are sourced from their respective cited sources. CDC-suspected AGS categories were obtained from a published CDC map and are used to assess geographic concordance.

## Manuscript

Hussain A, Mateus-Pinilla N, Smith RL. **A Spatially Validated Machine-Learning Framework for Mapping Alpha-Gal Syndrome Exposure Potential in Illinois.**

## Contact

**Abrar Hussain**  
University of Illinois Urbana-Champaign  
abrarh2@illinois.edu
