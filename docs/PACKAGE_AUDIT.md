# Package audit

Audit date: 2026-10-19.

## Included

- Four final estimation control streams and three dose-simulation control streams.
- Synthetic demonstration data and generators, with explicit manifests distinguishing synthetic inputs from empirical study records.
- GOF, pc-VPC, NPDE, Ka-sensitivity, C₂h boxplot, dose/PTA, and Figure S2 R scripts.
- File inventory, project-relative input conventions, and Git ignore rules for local analysis inputs and generated results.

## Portability changes in this pass

- Added a shared project-path helper. Analysis scripts locate the package from their own source-file path, not the current working directory.
- Replaced author-computer absolute paths with project-relative defaults and optional `ATX_*` environment-variable overrides.
- Routed generated plots and tables to project-relative `results/` subfolders by default.
- Removed the Ka sensitivity script's working-directory change; its inputs and outputs are now explicit paths.
- Added a local staging guide and Git ignore rules so empirical input files and generated results are not included by default.

## Not claimed

- The actual patient-derived fit/VPC/NPDE/simulation inputs are not bundled, and the diagnostic figures cannot be recreated from synthetic demo data alone.
- The PsN-generated VPC simulation stream was not included.
- Portability edits do not verify model estimates, convergence, VPC/NPDE calculations, or published figure values. Only checks listed in `VERIFICATION.md` are claimed.
- The package has not been uploaded. A public license still requires author approval.
