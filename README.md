# Woolly Rhinoceros

This repository contains the Python code and datasets used for palaeoenvironmental reconstruction and habitat-suitability modelling of the woolly rhinoceros (*Coelodonta antiquitatis*). It includes palaeoclimate datasets from 15 independent climate-state simulations, the corresponding global biome datasets generated with BIOME4, the MaxEnt workflow using temperature, precipitation, and biome predictors, and the code and data used to reconstruct δ¹⁸Owater (VSMOW).

## Palaeoenvironment

The `paleoenvironment/` directory contains the palaeoclimate and biome datasets used in the habitat-suitability analyses.

The palaeoclimate datasets represent 15 independent climate-state simulations from the Pliocene Modelling Intercomparison Project and include temperature and precipitation variables. The corresponding global biome datasets were generated with BIOME4.

### Repository structure

```text
paleoenvironment/
├── paleoclimate/
│   └── paleoclimate_<model_code>.nc
└── biome/
    └── biome_<model_code>.nc
```

Each palaeoclimate dataset has a corresponding BIOME4 biome dataset identified by the same model code.

## Habitat suitability

The `habitat_suitability/` directory contains the data and Python workflow used to model the potential habitat suitability of the woolly rhinoceros (*Coelodonta antiquitatis*) across 15 Pliocene-like and Pleistocene-like palaeoclimate simulations using temperature, precipitation, and biome predictors.

### Requirements

- ArcGIS Pro with ArcPy and a valid Spatial Analyst licence
- Java
- Jupyter Notebook or JupyterLab using the ArcGIS Pro Python environment
- NumPy, pandas, Matplotlib, and openpyxl

The included `maxent.jar` is MaxEnt 3.4.4, distributed under the MIT License. See the [official MaxEnt repository](https://github.com/mrmaxent/Maxent) for its licence and third-party licence information.

### Repository structure

```text
habitat_suitability/
├── MaxEnt_Modeling.ipynb
├── maxent.jar
└── data/
    ├── fossils/Fossils.xls
    ├── aoi/Eurasia.*
    ├── study_regions/
    └── paleoenvironment/
        ├── model_predictors/
        └── predictor_selection/
```

### Running the notebook

Open `habitat_suitability/MaxEnt_Modeling.ipynb`, select the ArcGIS Pro Python kernel, and run all cells in order.

### Outputs

Outputs are written to:

```text
habitat_suitability/output/
└── 20000-22000_Eurasia/
    ├── input/                    Derived analysis inputs
    └── results/                  Final model and projections
        ├── CV_3fold/             Cross-validation outputs
        ├── MESS/                 Extrapolation diagnostics
        ├── Stats/                Summary tables and figures
        └── Tuning_CV_3fold/      Parameter-tuning runs
```

## Reconstructing δ¹⁸Owater (VSMOW)

The `d18Owater_reconstruction/` directory contains the Python code and dataset used to reconstruct water δ¹⁸O values on the VSMOW scale and visualize the associated dual-clumped-isotope data.

### Requirements

The scripts require NumPy, pandas, Matplotlib, SciPy, openpyxl, uncertainties, D95eq, correldata, and LaTeX.

### Repository structure

```text
d18Owater_reconstruction/
├── Supplementary Dataset 5.xlsx
├── plot_dual_clumped_errorbars.py
└── d18O_water.py
```

### Running the scripts

Run both scripts from `d18Owater_reconstruction/`, with `Supplementary Dataset 5.xlsx` in the same directory.

- `plot_dual_clumped_errorbars.py` visualizes the Δ47–Δ48 data with sample uncertainties and the Fiebig 2024 equilibrium calibration.
- `d18O_water.py` reconstructs water δ¹⁸O values on the VSMOW scale, propagates the associated analytical and temperature uncertainties, and produces an Excel table and figures.

## Citation and licences

The original analysis code authored for this repository is released under the Apache License 2.0.

The input datasets remain subject to the licences and terms specified by their respective original providers.

The included `maxent.jar` remains subject to the MaxEnt licence and any applicable third-party licences. Please cite MaxEnt when using this workflow.
