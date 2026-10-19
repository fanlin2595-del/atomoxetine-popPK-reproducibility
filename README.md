# Atomoxetine population-PK analysis materials

This repository candidate contains the four final NONMEM model control streams, NONMEM dose-simulation controls, synthetic demonstration inputs, and the R scripts used for diagnostic and simulation figure generation. The synthetic records are fictional and are not derived from patient data.

## Portable paths

All R scripts locate the project from their own file location; they do not require a particular working directory or user-specific computer path. By default, private or study-derived inputs belong under `data/analysis_inputs/`, and generated files are written under `results/`.

You may override an input directory or output directory with the environment variables listed in [`docs/FILE_INVENTORY.md`](docs/FILE_INVENTORY.md). An override may be an absolute path or a path relative to the project root. The `data/analysis_inputs/` directory is ignored by Git to help prevent accidental publication of patient-derived files.

Run a script from any working directory with `Rscript /path/to/this-project/analysis/.../script.R`, or open the project folder in RStudio and source the script. The R packages required by each analysis are listed in the script headers and file inventory.

## Models and synthetic demonstrations

- `models/`: four estimation control streams (Models I–IV).
- `simulation/`: three dose-simulation control streams. Model IV has the same structure as Model II in this cohort and uses the equivalent genotype representation.
- `data/synthetic_demo.dat`: fictional concentration/event-history data for demonstrating model input structure.
- `data/simulation/`: fictional NONMEM dose-simulation inputs.
- `data/simulation_demo_combined/`: a harmonized synthetic input grid for inspection; its `DV=0` values are placeholders, not simulated concentration results.
- `scripts/`: synthetic-data preparation and static validation scripts.

For an estimation demonstration, run from a fresh copy of the project:

```sh
Rscript --vanilla scripts/prepare_demo.R
Rscript --vanilla scripts/check_package.R
```

For the NONMEM dose-simulation input demonstration:

```sh
Rscript --vanilla scripts/prepare_simulation_demo.R
Rscript --vanilla scripts/check_simulation_module.R
```

These scripts prepare isolated run folders and synthetic inputs only. NONMEM/NM-TRAN must be separately installed and licensed; no NONMEM run or manuscript estimate is reproduced by these preparation scripts.

## Analysis scripts and required inputs

The collected analysis scripts are:

- `analysis/diagnostics/`: GOF, pc-VPC, and NPDE diagnostics/figures;
- `analysis/simulation/`: C₂h distribution figures, dose/PTA figures, and Ka sensitivity summaries/figures;
- `analysis/supplementary/`: Figure S2 genotype comparison.

The analysis scripts need model outputs produced by the original NONMEM/PsN workflows. Such empirical inputs are deliberately not included because patient-level clinical/genetic records and derived outputs are restricted. Place authorized inputs in the project-relative folders described in [`data/analysis_inputs/README.md`](data/analysis_inputs/README.md), or set the corresponding environment variable to another local folder. Without those inputs the figure scripts stop with a missing-file message; the synthetic event datasets do not reproduce the paper's diagnostic plots, PTA results, or dose-selection figures.

## Data protection and release scope

No raw clinical dataset, patient-level empirical `.fit`/`.tab`/NPDE output, PsN VPC summary, NONMEM executable, password, or credential is intended for public release. The R scripts and NONMEM control streams are released under the MIT License (see [`LICENSE`](LICENSE)). Synthetic demonstration data are fictional and contain no patient records.

This package is for computational transparency and demonstration only. It is not a clinical dosing tool or clinical dosing recommendation.
