# Spectral-Diversity-RNV-Sooretama

Reproducible workflow for assessing spatial and temporal patterns of canopy spectral diversity in the Atlantic Forest using Sentinel-2 imagery.

This repository contains all scripts, metadata, and documentation required to reproduce the analyses presented in:

[Full manuscript citation will be added following publication]

## Study area

The analyses were conducted in two contiguous protected areas of the Atlantic Forest in southeastern Brazil:
Vale Natural Reserve (RNV)
Sooretama Biological Reserve (REBIO Sooretama)
Together, these reserves cover more than 40,000 ha of old-growth tropical forest and encompass strong environmental gradients associated with hydrology, topography, and landscape structure.


## Repository contents

### Scripts

| Script | Purpose |
|----------|----------|
| `01_spectral_diversity_computation.R` | Computes spectral species, alpha-diversity metrics (S, H', FRic, FDiv, FDis, RaoQ), beta-diversity ordinations, and Rao's Q maps from Sentinel-2 imagery. |
| `02_environmental_predictors.R` | Generates terrain, hydrological, and landscape predictors, including elevation, slope, aspect-derived variables, TRI, TPI, ruggedness, rumple, roughness, distance to rivers, and distance to forest edges. |
| `03_beta_diversity_alignment.R` | Performs within-year Procrustes alignment of beta-diversity iterations, generates annual consensus ordinations, aligns annual consensus ordinations among years, and calculates Procrustes statistics. |
| `04_environmental_modelling.R` | Performs environmental modelling, environmentally blocked nested cross-validation, spatial sensitivity analyses, Moran's I analyses, and temporal recurrence assessments. |

### Computational environment

| File | Purpose |
|----------|----------|
| `sessionInfo_scripts01-03.txt` | R version, operating system, and package versions used to run Scripts 01–03. |
| `sessionInfo_modelling.txt` | R version, operating system, and package versions used to run Script 04. |

---

## Input data

The workflow requires:

| Sentinel-2 Level-2A imagery |
| SRTM digital elevation model |
| River network |
| MapBiomas Collection 9 land-cover map |
| Detailed descriptions are provided in docs/input_data_description.csv |

---

## Execution order

Scripts should be executed in the following order:

1. `01_spectral_diversity_computation.R`
2. `02_environmental_predictors.R`
3. `03_beta_diversity_alignment.R`
4. `04_environmental_modelling.R`

---

## Outputs

The workflow reproduces:

- Spectral species maps
- Spectral alpha-diversity metrics
- Spectral beta-diversity maps
- Rao's Q maps
- Environmental predictor tables
- Environmental and spatial validation analyses
- Moran's I statistics
- Temporal recurrence analyses
- Main manuscript figures
- Supplementary figures and tables

---




