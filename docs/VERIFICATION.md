# Version-2 verification

Prepared and checked on 2026-10-18 using R 4.4.0. NONMEM was not invoked.

Passed:

- Compared all four prepared models with the latest supplied sources. After
  allowing the documented title/comment/path/format changes and the single
  full-width semicolon correction, all active statements and numeric settings
  matched, including input/table column order and Model III parameter bounds.
- Confirmed matching Model II/IV THETA, OMEGA, SIGMA, ERROR, ESTIMATION and COV
  blocks. Their PK blocks match after identifying the respective group selector.
- Confirmed that all four models use the source-supplied IWRES zero safeguard;
  generated demo copies introduce no active-code change.
- Regenerated the artificial input using Model II SIGMA=0.0712: 96 fictional
  subjects, 192 explicit dose records and 841 retained artificial observations.
  There were 23 artificial observations excluded below 5 ng/mL.
- Checked the 18-column schema, numeric values, chronology, ADDL/II dose
  histories, TAD, and conditional Model II/IV structural equivalence.
- Extracted the ZIP to a new temporary directory and reran the generator and
  checker. The regenerated input matched byte-for-byte.
- Synthetic input MD5: `811518bbce3bee7965c5a7921356f85e`.

Not verified:

- NM-TRAN compilation, NONMEM estimation or covariance completion.
- Final estimates or convergence of the original clinical-data runs.
- The actual clinical data coding and historical preprocessing.
- End-to-end execution of the statistical, diagnostic, simulation and plotting
  workflows. Their R source files are now collected under `analysis/`, but the
  scripts require empirical result files that are intentionally not included.

The demonstration seed and artificial design are not claims about the historical
study. No real patient records were needed or read during this preparation.

## Portability pass (2026-10-19)

All 13 R files parsed successfully after the portability edits. The model-package and simulation-module static checks both passed. A separate path-helper check confirmed project-root discovery and project-relative input/output resolution. The figure-generation scripts have not been executed end-to-end after these edits; empirical input folders remain empty by design, so successful figure regeneration requires authorized source outputs that are not part of the public package.
