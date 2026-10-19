# Input specification and boundaries

The supplied `$INPUT` lists **18 ordered fields**. The demonstration uses a
space-delimited NONMEM data file, not an export of clinical records. Its header
begins with `@` and is skipped by `IGNORE=@`.

| Position | Field | Role and demonstration convention |
| --- | --- | --- |
| 1 | ID | Fictional subject identifier, generated starting at 900001. |
| 2 | DAT2 | Dropped by the supplied `$INPUT`; zero placeholder, never a real date. |
| 3 | TIME | Elapsed hours from an arbitrary common origin; chronological within ID. |
| 4 | AMT | Oral dose in mg on dosing records; zero on observation records. |
| 5 | ADDL | Number of additional doses after the explicit dose; zero on observation records. |
| 6 | II | Additional-dose interval in hours; zero on observation records. |
| 7 | DV | Artificial concentration in ng/mL; zero placeholder when MDV=1. |
| 8 | MDV | 1 for dosing records without an observation; 0 for observed concentrations. |
| 9 | WT | Artificial body weight in kg; used by every model. |
| 10 | SEX | Arbitrary demo category 1 or 2; not used by these PK equations. |
| 11 | AGE | Artificial age in years; not used by these PK equations. |
| 12 | rs1065852 | Demo labels 0=CC, 1=CT, 2=TT; original study coding needs confirmation. |
| 13 | rs_grp | Demo 0=CC+CT, 1=TT. Model II applies its effect at 1. |
| 14 | rs5030865 | Zero placeholder only; no biological allele call is asserted. Unused by these equations. |
| 15 | rs16947 | Zero placeholder only; no biological allele call is asserted. Unused by these equations. |
| 16 | AS | Positive artificial activity score (2, 1.25, 0.5 by the demo's group). Model III uses AS/1.25. |
| 17 | CYP2D6 | Demo 1=reference category and 2=affected category; Model IV applies its effect at 2. Clinical labels need author confirmation. |
| 18 | TAD | Hours since the most recent dose, including implied ADDL doses; zero at dosing records. |

The synthetic mapping intentionally aligns the two affected-group indicators to
demonstrate the restricted Model II/IV equivalence. It does not imply that one SNP
determines activity score or phenotype in general, or establish original clinical
coding. A documented clinical translation table is still needed.

## Event history

The demo contains two dose epochs for each artificial subject. Additional doses
from the first epoch stop before the second epoch starts. TAD is computed from
the expanded artificial dose schedule, rather than from the preceding table row.
TIME is numeric elapsed time, not an HH:MM clock string. Because DAT2 is dropped,
this example does not rely on DAT2 to reconstruct multi-day chronology.

No EVID, CMT or SS field is added to the supplied `$INPUT`. ADVAN2's default oral
dose/central-observation compartments and the provided AMT/MDV encoding are
retained. Actual NONMEM parsing and execution must be verified on the intended
licensed installation.

With dose in mg, apparent volume in L and `S2=V/1000`, the intended model output
is ng/mL. The source reference weight is 30.7 kg; it is not replaced by the
artificial dataset's median weight. KA stays fixed at 1 h^-1.

## Missing documentation from the original workflow

The four control streams alone cannot establish the original mappings for SEX,
the three SNP columns, rs_grp, AS, and CYP2D6. Nor do they establish how elapsed
time was derived from clinical dates, how dose changes were encoded, or where
BLQ exclusion occurred. Those rules must be documented from the actual study
workflow, without releasing patient-level inputs.

Do not apply the Model III power relationship to AS=0 as though it were a
validated PM model. The synthetic demonstration intentionally avoids that case.
