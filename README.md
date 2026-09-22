# Woolly Rhinoceros

This repository contains the Python code and datasets used for palaeoenvironmental reconstruction and habitat-suitability modelling of the woolly rhinoceros (*Coelodonta antiquitatis*). It includes palaeoclimate datasets from 15 independent climate-state simulations, the corresponding global biome datasets generated with BIOME4, the MaxEnt workflow using temperature, precipitation, and biome predictors, and the code and data used to reconstruct δ¹⁸Owater (VSMOW).

## Palaeoenvironment

The `paleoenvironment/` directory contains palaeoclimate model datasets for 15 independent climate-state simulations and the corresponding global biome datasets used in this study.

The palaeoclimate datasets include temperature and precipitation variables from the Pliocene Modelling Intercomparison Project simulations. The corresponding global biome distributions were generated using the BIOME4 equilibrium global vegetation model, using the implementation available at [jedokaplan/BIOME4](https://github.com/jedokaplan/BIOME4), with the palaeoclimate model data as climatic inputs.

### Repository structure

```text
paleoenvironment/
├── paleoclimate/
│   └── paleoclimate_<model_code>.nc
└── biome/
    └── biome_<model_code>.nc
```

Each palaeoclimate dataset and its corresponding BIOME4 biome dataset are identified by the same model code.

### Palaeoclimate simulations

| Scenario | Model code | Climate state | Experimental settings |
|---:|:---:|---|---|
| 1 | `xozzc` | Pliocene – Warm summer orbit | 400 ppm CO₂ / Warm summer orbit |
| 2 | `xozzb` | Pliocene – Standard | 400 ppm CO₂ / Modern orbit |
| 3 | `xozzf` | Pliocene – Warm winter orbit | 400 ppm CO₂ / Warm winter orbit |
| 4 | `tenvl` | Pliocene – High CO₂ | 450 ppm CO₂ / Modern orbit |
| 5 | `tenvk` | Pliocene – Low CO₂ | 350 ppm CO₂ / Modern orbit |
| 6 | `tenvm` | Pliocene – PI CO₂ | 280 ppm CO₂ / Modern orbit |
| 7 | `tdlvb` | Early Pleistocene Glacial – PI CO₂ | 280 ppm CO₂ / M2 orbit / Medium ice |
| 8 | `tdlvd1` | Early Pleistocene Glacial – PI CO₂ | 280 ppm CO₂ / M2 orbit / Large ice |
| 9 | `tdlvc` | Early Pleistocene Glacial – Low CO₂ | 220 ppm CO₂ / M2 orbit / Medium ice |
| 10 | `tdlve` | Early Pleistocene Glacial – Low CO₂ | 220 ppm CO₂ / M2 orbit / Large ice |
| 11 | `xqlne` | Early Pleistocene Glacial – Low CO₂ | 220 ppm CO₂ / Modern orbit / Medium ice |
| 12 | `xqlnf` | Early Pleistocene Glacial – Low CO₂ | 220 ppm CO₂ / Orbit with low NHSI / Medium ice |
| 13 | `xnzib` | Late Pleistocene Glacial | 180 ppm CO₂ / LGM orbit / Standard ice |
| 14 | `xnziw` | Late Pleistocene Glacial | 180 ppm CO₂ / LGM orbit / Medium ice |
| 15 | `xnzix` | Late Pleistocene Glacial | 180 ppm CO₂ / LGM orbit / Large ice |

These model codes are also used for the corresponding environmental predictors and MaxEnt palaeoclimate projections in `habitat_suitability/`.

## Habitat suitability

The `habitat_suitability/` directory contains the data and Python workflow used to model the potential habitat suitability of the woolly rhinoceros (*Coelodonta antiquitatis*) across the 15 palaeoclimate simulations using temperature, precipitation, and biome predictors.

The final MaxEnt predictors are annual mean temperature, annual total precipitation, and biome. The `xnzib` Late Pleistocene glacial simulation is used as the training climate, and additional temperature and precipitation variables for `xnzib` are included for predictor-correlation analysis.

### Requirements

- ArcGIS Pro with ArcPy and a valid Spatial Analyst licence
- Java
- Jupyter Notebook or JupyterLab using the ArcGIS Pro Python environment
- NumPy, pandas, Matplotlib, xlrd, and openpyxl

### Repository structure

```text
habitat_suitability/
├── MaxEnt_Modeling.ipynb
├── maxent.jar
└── data/
    ├── fossils/
    │   └── Fossils.xls
    ├── aoi/
    │   └── Eurasia.*
    ├── study_regions/
    │   ├── JunggarBasin.*
    │   ├── LoessPlateau.*
    │   └── TibetanPlateau.*
    └── paleoenvironment/
        ├── biome/
        │   └── biome_<model_code>.tif
        ├── temperature/
        │   ├── temp_annual_mean_<model_code>.tif
        │   └── temp_<additional_stat>_xnzib.tif
        └── precipitation/
            ├── precip_annual_total_<model_code>.tif
            └── precip_<additional_stat>_xnzib.tif
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

BIOME4 is used through the implementation available at [jedokaplan/BIOME4](https://github.com/jedokaplan/BIOME4).

The included `maxent.jar` is MaxEnt 3.4.4 and remains subject to the MaxEnt licence and any applicable third-party licences. Please cite MaxEnt when using this workflow.
