# File inventory and release boundaries

## Included materials

| Folder/file | Contents | Notes |
| --- | --- | --- |
| `models/model_I_BW.mod` through `model_IV_CYP2D6.mod` | Four final estimation control streams | Use synthetic demo inputs only in the prepared demo folders. |
| `simulation/sim_model_I_BW.mod` through `sim_model_III_AS.mod` | Three NONMEM dose-simulation controls | Model IV is represented by the structurally equivalent Model II simulation scenario. |
| `data/synthetic_demo.dat` | Artificial dose histories and concentrations for fictional IDs | Demonstrates the estimation input format; does not reproduce study data or estimates. |
| `data/simulation/` | Artificial dose-simulation input records | Demonstration grid only; `DV=0` is an input placeholder. |
| `data/simulation_demo_combined/` | Harmonized synthetic grid across four covariate scenarios | For inspection/filtering, not a direct input to one control stream. |
| `scripts/` | Synthetic-input generators, checks, and project-path helpers | Generators do not read clinical files. |
| `analysis/diagnostics/` | GOF, pc-VPC, and NPDE R scripts | Require original model/PsN result files, not included. |
| `analysis/simulation/` | C₂h boxplot, dose/PTA, and Ka sensitivity R scripts | Require NONMEM simulation outputs or PTA summary tables, not included. |
| `analysis/supplementary/` | Figure S2 ETA(CL) genotype comparison R script | Requires the patient-derived `.fit` input, not included. |
| `docs/` | Data definitions, verification notes, and this inventory | Describes scope and limitations of the demo package. |

## Portable input layout

By default, scripts read analysis inputs from `data/analysis_inputs/` and write results to `results/`. The input directory is Git-ignored. Use only locally authorized source outputs; do not copy patient-level files into a public repository.

| Analysis | Default input layout under `data/analysis_inputs/` | Environment overrides |
| --- | --- | --- |
| GOF | `GOF/{BASE/5-M1.fit,BW/7.fit,RS1065852/15.fit,AS/12-M1.fit}` | `ATX_GOF_INPUT_DIR`, `ATX_GOF_OUTPUT_DIR` |
| pc-VPC | `pcVPC/BW/pcVPC`, `pcVPC/RS1065852/pcVPC/pcvpc-overall`, `pcVPC/AS/pcVPC`, `pcVPC/RS1065852/pcVPC/pcvpc-strat` | `ATX_PCVPC_INPUT_DIR`, `ATX_PCVPC_OUTPUT_DIR` |
| NPDE | `NPDE/{BW,RS1065852,AS}/{vpc_original.npctab.dta,vpc_simulation.1.npctab.dta}` | `ATX_NPDE_INPUT_DIR`, `ATX_NPDE_OUTPUT_DIR` |
| Dose/PTA | `simulation/{BW/attainment_ModelI_BW_3strata.csv,RS1/attainment_ModelII_3strata.csv,AS/attainment_ModelAS_3strata.csv}` | `ATX_SIMULATION_ROOT`, `ATX_DOSE_PTA_OUTPUT_DIR` |
| C₂h boxplots | `simulation/{BW,RS1,AS}/...sim*.tab` | `ATX_SIMULATION_ROOT`, `ATX_BOXPLOT_OUTPUT_DIR` |
| Ka sensitivity | `ka_sensitivity/{0.5_sim.tab,0.5_sim1.3.tab,2_sim.tab,2_sim1.3.tab}` plus `simulation/RS1/{sim_RS1.tab,1.3/sim_RS1.3.tab}` | `ATX_KA_INPUT_DIR`, `ATX_RS1_INPUT_DIR`, `ATX_KA_OUTPUT_DIR` |
| Figure S2 | `supplementary/7.fit` | `ATX_FIGURE_S2_INPUT_FILE`, `ATX_FIGURE_S2_OUTPUT_DIR` |

The environment-variable path may be absolute or project-relative. The NPDE script also accepts an optional output-directory argument when called via `Rscript`.

## Deliberately excluded

- Raw individual-level clinical, demographic, dosing, concentration, and genetic data.
- Patient-derived NONMEM `.fit`, `.tab`, `.lst`, `.ext`, and NPDE tables.
- PsN VPC `vpc_results.csv` and `vpctab*` outputs, including the automatically generated VPC simulation control stream.
- Empirical PTA/attainment tables and original model run outputs.
- NONMEM/PsN executables, credentials, and machine-specific configuration.

Although empirical results are excluded, several scripts process individual-level input files when run locally. Those outputs remain protected data and must not be uploaded merely because their filenames contain no direct identifiers.

## Verification and limitations

The scripts have been changed to resolve the project root from their own location and to use project-relative default input/output folders. This removes dependence on the original author-computer directory layout; it does not create missing empirical inputs. The synthetic demo datasets are not substitutes for the original analysis results, and the diagnostic figures cannot be regenerated from the packaged synthetic inputs alone.

Before public deposit, confirm author approval for sharing the model/code files, select a license with author approval, run the static checks, and inspect the staged file list for patient-derived files, credentials, or personal information. NONMEM estimation and end-to-end figure reproduction are not claimed unless separately executed and documented.
