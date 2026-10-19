# Synthetic NONMEM simulation pairs

Each subfolder is a standalone run directory. Keep the control stream and its `.dat` file together, and run NONMEM from inside that subfolder. These are synthetic demonstration inputs, not the study's original simulation grid or clinical dosing recommendations.

| Run folder | Control stream | Matching input dataset | Expected simulation table |
| --- | --- | --- | --- |
| `sim_model_I_BW/` | `model.mod` (Model I, body weight) | `synthetic_sim_BW.dat` | `sim_model_I_BW.tab` |
| `sim_model_II_rs1065852/` | `model.mod` (Model II, rs1065852 group) | `synthetic_sim_RS.dat` | `sim_model_II_rs1065852.tab` |
| `sim_model_III_AS/` | `model.mod` (Model III, activity score) | `synthetic_sim_AS.dat` | `sim_model_III_AS.tab` |

The `$DATA` statement in each `model.mod` names the matching `.dat` file in its same folder. Each input has a header beginning with `@`, which is skipped by `IGNORE=@`. The three controls use `ONLYSIM NSUB=500`; they are simulation controls, not estimation controls. `DV=0` in the input is a placeholder, not an observed concentration; NONMEM generates predictions in the simulation output.

Run the locally licensed NONMEM installation separately in each folder. For example, if `nmfe75` is available:

```sh
cd sim_model_I_BW
nmfe75 model.mod model.lst
```

Repeat with the other folders. Check the listing for normal completion and verify that the expected `.tab` file is created. NONMEM was not run as part of preparing this package.
