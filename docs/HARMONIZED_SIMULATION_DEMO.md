# Synthetic dose-simulation input grids

Run `Rscript --vanilla scripts/prepare_simulation_demo.R` from any directory.
The script writes one comma-delimited combined dataset with an `@` header under
`data/simulation_demo_combined/`. Existing output is never overwritten.

These are fictional design grids, not patient records and not simulated
concentration results. The script reads no clinical data and assigns no
empirical concentration values. `DV=0` is a placeholder for records where
the simulation control stream is expected to generate an observation.

The compact design uses whole-kilogram weights from 17 to 60 kg, candidate
once-daily doses of 0.5, 0.8, and 1.0–1.4 mg/kg/day, and TAD observations at
2, 12, and 24 hours after the tenth daily dose. It demonstrates how a
simulation grid can be encoded and inspected; it does not reproduce the
original finer weight grid, random-effects design, 500 simulation replicates,
source-data distributions, or manuscript results.

The combined file stacks the BW, rs1065852, activity-score (AS), and CYP2D6
scenarios and identifies them with `SCENARIO`, `COVARIATE`, and
`COVARIATE_VALUE`. Model IV's genotype/phenotype equivalence is represented
only by arbitrary category codes for demonstration. The harmonized dataset is
for inspection and filtering; filter a single scenario and map its covariate
field to the matching control stream before use. This package does not yet
include runnable simulation control streams or execute NONMEM from this
generator.
