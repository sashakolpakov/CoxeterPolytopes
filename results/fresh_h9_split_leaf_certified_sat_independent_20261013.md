# Independent replay of all 1,761 necessary escape constraints

The independent exact replay passed all **1,761 archived constraints**
against the **3,040 archived cap entries**. It evaluated **5,353,440
pairing signs** and found **814,787 positive cutter incidences**. Every
true cutter is present in its archived clause, and there are **no extra
allowed cutters** in this batch.

This audit checks necessary conditions only. The archived SAT run ended
without an UNSAT certificate or a certified finite-volume candidate.
A subsequent, separate
[exact endpoint-cover argument](fresh_h9_split_leaf_endpoint_independent_20261013.md)
excludes finite-volume completion from the fixed bank with the seed
facets retained; that conclusion does not come from these clauses alone.

The [new audit result](../data/enumerations/fresh_h9_split_leaf_certified_sat_independent_20261013.json)
is separate from the preserved
[earlier 1,041-constraint audit](../data/enumerations/fresh_h9_split_leaf_escape_constraints_independent_20261013.json).
The [independent checker](../scripts/fresh_h9_split_leaf_escape_constraints_independent_20261013.py)
now accepts `--source` and `--output`; an alternate source requires an
explicit output path to prevent overwriting the original audit.

## What the replay proves

The checker uses only the Python standard library and imports no project
producer or search module. It reconstructs the exact seed Gram matrix,
checks inertia `(9,1,6)` in its sixteen-generator presentation, and
verifies the strict timelike center of square `-64-128*phi`. It checks
every cap's literal normal against its covector, norm two, and strict
center pairing.

For each escape it requires exactly sixteen built-in integer
coefficients, a positive integer denominator, positive Lorentz norm,
correct time orientation, and all sixteen strict seed inequalities.
This gives 28,176 strict seed checks. It recomputes the archived norm
and center-pairing numerators and denominators exactly.

The cap-pairing signs are decided in `Q(phi,sqrt(D))` by exact comparisons
of squares in `Q(phi)`, rather than by the producer's interval signs.
Dependent or zero radicals are handled explicitly. Source files and
their hashes are checked before and after the audit.

If a selected set of caps contains no positive cutter of a certified
spacelike escape, that escape satisfies every selected cap inequality.
Perturbing toward the common strict timelike center makes those
inequalities strict while preserving positive norm for sufficiently
small perturbations. Thus finite volume requires at least one cutter
from each clause. This implication does not establish that all clauses
together are inconsistent.

The audited search archive has SHA-256
`9ff2e6029741634d71211a3c0b5b31d938bc38cc86ac72ac85e81936b93c7daf`.

## Input validation fixes

The [escape-constraint producer](../scripts/fresh_h9_split_leaf_escape_constraints_20261013.py)
now checks vector length and coefficient types before any truncated
`zip` calculation can occur. It rejects lengths 15 and 17, booleans,
integral `Fraction` objects, and floating-point coefficients. Built-in
integers of arbitrary size remain supported.

Its radical enclosure routine now requires the radicand to be exactly
nonnegative under the algebraic `Q(phi)` sign test. Checking only the
upper endpoint of a rounded interval was insufficient. The regression
uses `D=-1/(10^24)^2`: its interval upper endpoint is zero, but its exact
sign is negative, so the new test rejects it. Zero and unit radicands
remain accepted and correctly enclosed.

Producer replay passed both the original 1,041 records and the later
1,761 records, including every malformed-vector control, the negative
radicand control, and the zero-escape control. No archived escape or
cutting clause was changed by this hardening.

## Reproduction

```sh
python scripts/fresh_h9_split_leaf_escape_constraints_20261013.py --verify
python scripts/fresh_h9_split_leaf_escape_constraints_20261013.py --verify \
  --source data/enumerations/fresh_h9_split_leaf_certified_sat_20261013.json

python scripts/fresh_h9_split_leaf_escape_constraints_independent_20261013.py \
  --verify \
  --source data/enumerations/fresh_h9_split_leaf_certified_sat_20261013.json \
  --output data/enumerations/fresh_h9_split_leaf_certified_sat_independent_20261013.json

python scripts/fresh_h9_split_leaf_escape_constraints_independent_20261013.py --verify
```

The last command replays the original audit using its unchanged default
source and output. To regenerate the new audit, omit `--verify` from
the command that supplies both paths. These operations do not rerun the
SAT search and do not certify a polytope.
