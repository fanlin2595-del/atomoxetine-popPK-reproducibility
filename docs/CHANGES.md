# Stage 1, version 2: source mapping and changes

## Portability pass (2026-10-19)

- Prepared a separate `atomoxetine_reproducibility/` working copy; earlier candidate folders were left unchanged.
- Added `scripts/project_paths.R` and updated analysis scripts to locate the project from the script file, use project-relative default inputs/outputs, and accept environment-variable overrides.
- Removed hard-coded author-machine paths and the Ka script's `setwd()` behavior.
- Added an ignored local staging folder for authorized empirical inputs and directed generated files to `results/`.
- No model equations, statistical settings, figure styling, or numerical analysis choices were intentionally changed in this portability pass.
- Script parsing and end-to-end checks are pending for this copy; do not interpret path portability as a successful rerun of the original analyses.

Prepared 2026-10-18 from the latest four author-supplied files. Original source
files and the version-1 package were not modified. This remains a partial code
release candidate; the remaining analysis scripts and NONMEM tests are pending.

## Source mapping

| Supplied source | Prepared control stream | Model |
| --- | --- | --- |
| `Model I_BW.modmod` | `models/model_I_BW.mod` | I: body weight |
| `15.modmod` | `models/model_II_rs1065852.mod` | II: rs1065852 grouping |
| `12-M.modmod` | `models/model_III_AS.mod` | III: activity score |
| `11-M.modmod` | `models/model_IV_CYP2D6.mod` | IV: predicted phenotype |

## Changes already present in the latest author sources

These are differences from the previous supplied files, not new statistical
decisions made during packaging. All were retained as supplied.

- **Model I:** no change to active mathematical or estimation statements.
- **Model II:** SIGMA starting variance changed from 0.0711 to 0.0712.
- **Model III:** THETA(4) changed from `(0,0.585)` to `(-3,0.585,3)`.
  The initial estimate is unchanged, but the allowable search interval now
  includes negative exponents and has an explicit upper bound. This is an
  estimation-bound change, not merely formatting.
- **Model IV:** IWRES now uses `IRES/(IPRED+DEL)`. Its THETA, OMEGA and SIGMA
  starting values now match the latest Model II, including SIGMA=0.0712.
  The phenotype selector `CYP2D6.EQ.2` is retained.

Matching starts and mathematical structure do not establish the historical
final estimates or convergence. The corresponding final `.lst` files still
need to be checked.

## Packaging changes applied to the latest sources

1. Assigned descriptive `.mod` filenames and model-specific problem titles.
2. Replaced ambiguous inherited header comments with source/version notes.
3. Replaced the clinical input filename with `synthetic_demo.dat`.
4. Assigned distinct descriptive `.fit` output filenames to avoid collisions.
5. Normalized line endings to LF; input/table columns and all numeric settings
   remain unchanged.
6. Corrected **one full-width semicolon** in `11-M.modmod`, on the OMEGA(V) line:
   `0.329` followed by U+FF1B and `IIV IN V` is written as
   `0.329; IIV IN V`. This restores an ASCII comment delimiter, without changing
   the numeric value. The original file itself was not edited.

There are no demo-only changes to any model equation in v2. The obsolete
version-1 Model IV IWRES substitution has been removed from `prepare_demo.R`.
All run copies have the same active code as the corresponding prepared models.

## Synthetic demonstration update

The generator now reads the latest Model II SIGMA starting variance (0.0712).
Its synthetic observations were regenerated rather than copied from version 1.
The seed and fictional design are unchanged. No clinical data were read.

The demonstration script, input and static checks are newly prepared supporting
materials, not substitutes for the original study's analysis scripts. The input
cannot reproduce the manuscript's empirical results.
