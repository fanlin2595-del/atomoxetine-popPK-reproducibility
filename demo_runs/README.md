# Synthetic concentration dataset with four estimation models

Each subfolder contains one NONMEM estimation control stream (`model.mod`) and a matching copy of `synthetic_demo.dat`. The data are synthetic, with invented concentrations for 96 fictional IDs; they are not patient data and do not reproduce the study estimates.

| Folder | Model | Input in the same folder |
| --- | --- | --- |
| `model_I_BW/` | Model I, body weight | `synthetic_demo.dat` |
| `model_II_rs1065852/` | Model II, rs1065852 genotype group | `synthetic_demo.dat` |
| `model_III_AS/` | Model III, activity score | `synthetic_demo.dat` |
| `model_IV_CYP2D6/` | Model IV, CYP2D6 predicted phenotype | `synthetic_demo.dat` |

The `$DATA` statement in every control stream is `synthetic_demo.dat`, so each model can be run independently from its own folder using a locally installed/licensed NONMEM launcher. Example, if `nmfe75` is installed:

```sh
cd model_I_BW
nmfe75 model.mod model.lst
```

Repeat in the other model folders. Do not run all four in one shared directory because their output filenames and working files are model-specific. Check the listing for estimation and covariance completion. NONMEM was not run during package preparation; the synthetic data and model starting values are for demonstration only.
