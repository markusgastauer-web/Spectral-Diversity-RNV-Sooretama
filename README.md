Spectral-Diversity-RNV-Sooretama

Reproducible workflow for assessing spatial and temporal patterns of canopy spectral diversity in the Atlantic Forest using Sentinel-2 imagery.

This repository contains all scripts, metadata, and documentation required to reproduce the analyses presented in:

[Full manuscript citation will be added following publication]

Study area

The analyses were conducted in two contiguous protected areas of the Atlantic Forest in southeastern Brazil:
Vale Natural Reserve (RNV)
Sooretama Biological Reserve (REBIO Sooretama)
Together, these reserves cover more than 40,000 ha of old-growth tropical forest and encompass strong environmental gradients associated with hydrology, topography, and landscape structure.

Repository contents

Scripts and Purpose

01_spectral_diversity_computation.R	(Spectral species mapping, alpha diversity, beta diversity, and Rao's Q)
02_environmental_predictors.R	(errain, hydrological, and landscape predictors)
03_beta_diversity_alignment.R	 (Within-year and among-year Procrustes alignment)
04_environmental_modelling.R	(Environmental modelling, validation, sensitivity analyses, and temporal recurrence)

Input data

The workflow requires:
Sentinel-2 Level-2A imagery
SRTM digital elevation model
River network
MapBiomas Collection 9 land-cover map
Detailed descriptions are provided in docs/input_data_description.csv.

Workflow

Scripts should be executed in the following order:
Spectral diversity computation
Environmental predictor generation
Beta-diversity alignment
Environmental modelling and validation
Reproducibility
All stochastic procedures use fixed random seeds.

Package versions and software information are automatically recorded during workflow execution.

The computational environment can be reconstructed using the renv.lock file provided in this repository.

Outputs

The repository reproduces:
Spectral diversity rasters
Spectral beta-diversity maps
Rao's Q maps
Environmental predictor tables
Environmental and spatial validation analyses
Moran's I statistics
Temporal recurrence analyses
Main manuscript figures
Supplementary figures and tables

Contact

Markus Gastauer Instituto Tecnológico Vale (ITV) Belém, Pará, Brazil - markus.gastauer [at] itv.org
