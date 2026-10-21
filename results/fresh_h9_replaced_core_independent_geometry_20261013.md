# Independent geometry check of the final 45-wall chamber

Case 0 of the [round-four selection](../data/enumerations/fresh_h9_replaced_core_20261013_round4_selection.json)
passes a second exact geometry check. The displayed chamber has infinite
volume. Its 45 genuine facets and absence of affine subdiagrams do not make
it compact: an exact strict future spacelike escape survives.

The [checker](../scripts/fresh_h9_replaced_core_independent_geometry_20261013.py)
uses `Golden` arithmetic and a separate recursive quadratic-radical
implementation. It checks the literal rank-ten Gram matrix, signature `(9,1)`,
inverse, and all 45 normals' Gram pairings and positive norms. It then checks
all 45 strict inequalities for the common timelike center and the first saved
spacelike escape. The radical arithmetic handles dependent and zero auxiliary
radicands explicitly; it assumes neither that every tower is a field nor that
the auxiliary square roots are independent.

Writing `I` for the seed Gram inverse, the archived positive covectors `h,y`
represent the actual center and escape as `-I h,-I y`. Their exact values are:

```text
center norm = 7/2 - 65 phi < 0
escape norm = (684182955592023711 - 325126963008839837 phi)
              / 2500000000000000000 > 0
center–escape pairing = -1598251911/1000000000 - 2077989 phi/8000000 < 0
```

All 45 wall pairings at both points are strictly negative. A segment from the
strict timelike point to the strict spacelike point crosses the future light
cone while remaining strict in every wall. An open patch of ideal directions
therefore lies in the chamber, proving infinite volume.

The exact input extraction reuses the construction's source profile loader;
its frozen [input bundle](fresh_h9_replaced_core_independent_geometry_20261013_input.json)
is then replayed without importing any construction module. Source artifact
hashes are checked. The second audit does **not** repeat all 990 compositum
pair classifications or the complete affine census. Those checks are
separately covered by the construction agent's
[producer-free verifier](../scripts/fresh_h9_replaced_core_20261013_verify.py).

The facet argument is sound under those verified acute pairings. For a strict
timelike point `v`, projection to wall `r_i` gives
`v_i = v - (v,r_i) r_i/(r_i,r_i)`. Its norm remains negative. For every other
wall `r_j`, the hypotheses `(v,r_j)<0` and `(r_i,r_j)<=0` give
`(v_i,r_j)<0`. Thus every retained wall has a relative-interior point and is a
genuine facet; the reported facet count is not inflated by redundancy.

The [audit result](fresh_h9_replaced_core_independent_geometry_20261013.json)
also records rejection of the zero escape and regression checks for zero and
dependent radicals. This result concerns the displayed chamber only. Further
caps are not excluded.

Replay from the repository root:

```sh
PYTHONPATH=src python scripts/fresh_h9_replaced_core_independent_geometry_20261013.py --verify
```

Regenerate the frozen input from the source profiles if intentionally needed:

```sh
PYTHONPATH=src:. python scripts/fresh_h9_replaced_core_independent_geometry_20261013.py --extract
```
