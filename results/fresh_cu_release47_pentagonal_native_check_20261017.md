# Exact rejection of the six provisional pentagonal-cap improvements

**FINAL/FROZEN after generation, fresh saved replay and independent review.**
All six apparent nine-open-end candidates in the
[pentagonal-cap scout](fresh_cu_release47_pentagonal_cap_20261017.md) have
an inadmissible native K/R contact. The exact63-wall/199-finite-vertex/
ten-open-end infinite-volume control is unchanged.

The [source](../scripts/fresh_cu_release47_pentagonal_native_check_20261017.py)
and [saved certificate](fresh_cu_release47_pentagonal_native_check_20261017.json)
reconstruct only the necessary local geometry, using the frozen physical
interval engine. This is not an independent implementation of that engine.

## Exact local recipes

In all six records P3 equals old47. Starting from the exact control,
project the old normals off the positive S6 Gram. The following successive
three-frame norm recipes define the new normals:

| Normal | Frame | Native finite labels | A2 profile |
|---|---|---|---|
| L |30,52,old47 |3,2,2 |1 or4 |
| K |30,L,old47 |4,2,4 |0 |
| A |7,30,K |3,2,2 |0 |
| R |7,A,old47 |2,2,2 |0 |

At every step the frame Gram is inverted with outward rational intervals.
A strictly positive complementary line proves that the nondegenerate
three-frame is Lorentzian. The positive remaining norm gives exactly two
algebraic norm branches. The saved sign explicitly defines the chosen
branch; comparison with numerical products only identifies that sign and
checks that the other sign is separated. Native A2 offsets are subtracted
before the frame solve and restored for products. The norm-two property
and prescribed contact equalities follow from these algebraic recipes,
not from an interval merely containing a desired value.

The A frame is fixed even though the original scout used a moving C(u)
construction. Each of its three A families prescribes A7=3,A30=2,AK=2.
The independently reconstructed fixed frame selects the same branch in
each numerical record. Thus C(u) and the different A/C labels are not
needed for this necessary K/R gate. The exact check is repeated for each
of the six saved candidates, including both L profiles.

## Strict native angle gap

Every reconstructed branch gives the following outward enclosure:

\[
-1.999991758569355 < \langle K,R\rangle
< -1.999991758567004.
\]

The saved rational bounds establish strictly

\[
-2\cos(\pi/1095)<\langle K,R\rangle<-2\cos(\pi/1094).
\]

Since the finite Coxeter products are strictly ordered, no finite label
can equal this product. It is also strictly greater than−2, so it is
neither parallel nor dotted. The six candidate chambers therefore fail a
native Coxeter condition regardless of their numerical target-face cone
counts. This is an exact rejection of the six literal branches, not an
all-label or all-parameter exclusion of pentagonal cap configurations.

The floating classifier's absolute tolerance had labelled these products
1094, although the spacing between consecutive high-label products is
much smaller than that tolerance. Future discovery must treat such
unprescribed contacts as provisional and compare the tolerance with the
local finite-label spacing; exact promotion still requires an exact
contact construction or exclusion check. No saved predecessor is changed
to conceal these diagnosed false positives.

## Verification and limits

Generation93838 and fresh saved replay98435 both terminated with
exit0/PASS. A separate no-producer saved-data walk checked the six strict
rational angle gaps, branch signs, positive kernel/remaining squares and
bindings, alongside the scout's literal disk and saved incidence records.
It did not reconstruct the physical intervals or run a numerical cone.

Root read the complete source and proof: fixed-frame branch definitions,
interval inverses, positive complementary lines/remaining norms, profile
offsets and strict threshold separation. That review found no issue in
the stated rejection, and did not independently reconstruct the inherited
physical engine. Direct source bindings, Python AST and local report links
were checked. No new cone, ambient finite-vertex census or volume proof is
asserted by this local check. No old file or exact control was modified.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_cu_release47_pentagonal_native_check_20261017.py
OPENBLAS_NUM_THREADS=1 python scripts/fresh_cu_release47_pentagonal_native_check_20261017.py --verify
```
