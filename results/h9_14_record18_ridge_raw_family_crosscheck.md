# Record 18: literal ridge-family cross-check

Status: two independently implemented exact searches now agree on every
literal compact-incidence/ridge target in the high-addition branches.  This
cross-check does not produce a hyperbolic Coxeter example and does not test
finite-volume noncompact chambers.

## What was checked

The production search uses a collision-free byte-matrix symmetry quotient,
then expands every accepted quotient target under the complete stabilizer of
the prescribed missing four-faces.  Every expanded family is checked
literally for:

- no transversal of order at most four;
- exactly the scheduled number of order-five transversals;
- every six-set having ridge degree zero or two;
- the prescribed local incidence bounds; and
- the exact minimal-blocker inventory.

The normalized expanded families were then compared as literal sets with the
families obtained by expanding the standalone independent audit.  Equality
here means equality of every missing-five-face family, not merely equality of
counts or orbit sizes.

## Exact results

For `a=1/add6`, representatives `23` and `15` have no target.  Representative
`53` has eight quotient targets and exactly 1,080 literal targets.  The two
implementations return exactly the same 1,080 families.

For `a=2/add4`, all twelve representatives were replayed.  Ten have no target.
Representative `(53,198)` has six literal targets and `(53,202)` has two, for
eight in total.  Again, the two implementations return exactly the same
literal families.

The independent downstream audit already proves that all eight quotient
classes in the `a=1` branch and all three subgroup classes representing the
eight `a=2/add4` literal targets admit no compact Coxeter labeling.  Therefore
these incidence survivors do not yield compact hyperbolic Coxeter polytopes.
No statement about the finite-volume noncompact branch follows from this
compact ridge audit.

## Artifacts

- `scripts/explore_h9_14_record18_ridge_aware_additions.cpp`
- `scripts/certify_h9_14_record18_a1_ridge_raw_families.py`
- `scripts/certify_h9_14_record18_a2_add4_ridge_raw_families.py`
- `data/enumerations/h9_14_record18_a1_ridge_raw_families.json`
- `data/enumerations/h9_14_record18_a2_add4_ridge_raw_families.json`
- `tests/test_h9_14_record18_ridge_raw_families.py`
