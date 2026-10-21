# A four-duty cavity after releasing wall47

No improved chamber is promoted. The exact deletion audit identifies a cheaper
cavity than the previous release52 experiment, and a new coordinated four-cap
search tests its whole boundary. All 8,113 coherent numerical collections have
15 open ends, equal to the control. Their six extra intrinsic finite vertices
do not count as progress toward compactness.

The control is the [64-wall chamber with198 finite vertices](fresh_nonarithmetic_replace47_close12_20261017.md)
and its [exact S6 boundary](fresh_nonarithmetic_replace47_close12_boundary_20261017.md).
The [source](../scripts/fresh_nonarithmetic_release47_cycle_patch_20261017.py)
and [result](fresh_nonarithmetic_release47_cycle_patch_20261017.json) keep exact
integer/profile statements separate from the floating-point discovery stage.

## Exact old-geometry cost and conditional topology

Deleting47 loses precisely the two finite supports

`S6+[7,12,47]` and `S6+[12,47,52]`,

where `S6=[0,1,2,3,5,6]`. Repeating the frozen complete finite-incidence
closed-edge audit shows that the only old compact face of dimension at least2
that loses its certificate is the quadrilateral `S6+[12]`. In particular,
the8-heptagon and all higher-dimensional compact faces have no vertex on47.
This is an exact cost of deleting the old wall; it does not certify any new
floating-point geometry.

The numerical full boundary after deletion has23 finite and7 exterior rays,
with15 open ends. Its B endpoint duties form the four-cycle

`a=7 — b=12 — c=52 — d=30 — a=7`.

This differs from the seven-duty cavity and18 ends caused by deleting52.
A proposed four-cap disk has boundary cap order `L,P,R,Q`, where L terminates
both `ab` and `bc`, P terminates `cd`, Q terminates `da`, and R has one old
neighbor d. Its dual triangles are

`abL, bcL, cdP, daQ, cLP, dPR, dRQ, aQL, LPR, LRQ`.

The script checks all edge incidences, all vertex links, and Euler
characteristic1. The old–new/new–new edge counts are8 and5; the triangle
counts with one/two/three new vertices are4,4,2. Under the explicit retained
incidence, simple/no-ideal, genuine disk hypotheses, this gives

\[
 \Delta(2E-3V_f)=2(8+5)-3(4+4+2)=-4.
\]

Thus a geometrically realized filling of this disk would reduce15 ends to11.
The integer disk is a template, not a realization theorem. It avoids the
previous failed P–N/P–Q joins: when P and Q are dotted, R is solved to join
both instead of demanding that P and Q meet.

## Complete finite profile restriction on the left seed

The already proved [seven actual-S6 A2 profiles](fresh_nonarithmetic_h3_a2_profiles_20261016.md)
are retained. A cap L terminating both `7/12` and `12/52` must be right to12
and52. Indeed their normalized Gram is

\[
 \begin{pmatrix}2&-\sqrt3\\-\sqrt3&2\end{pmatrix},
\]

whose inverse has diagonal2 and nonnegative entries. Any nonzero Coxeter
contact has magnitude at least1 and therefore Schur cost at least2,
incompatible with a positive-definite finite corner.

The remaining finite test is extension of the A3 chain `5—6—7`. Its inverse is

\[
 \frac14\begin{pmatrix}3&2&1\\2&4&2\\1&2&3\end{pmatrix}.
\]

Let `(u,v,w)` be the nonnegative contact magnitudes to these three nodes.
Their inverse quadratic cost must be strictly below2. An endpoint contact
with label at least6 already costs at least9/4. Two nonzero contacts cost
at least2, with equality only for the two endpoints at label3. A center
contact of label at least4 costs at least2. Consequently the complete list is:

* no contact to5/6, with the7-label2,3,4,5;
* contact to5 of label3,4,5, with7 right;
* contact to6 of label3, with7 right.

These eight profiles yield at most16 normal branches. The exact code checks
all strict costs and exclusion inequalities using rational arithmetic and
the elementary squared comparison for the label5 cost
`(9+3sqrt5)/8 < 2`. This is a complete necessary finite-profile statement,
including all high labels. Checking the resulting normals against the other
old walls remains numerical in this package.

## Bounded coordinated search and full boundary result

The changed cavity is tested in all four rotations. Every solve includes
all seven A2 profiles, both real norm branches, the required positive-definite
projected endpoint corners, every retained native old-wall product including
white offsets, and the common old strict center. Assigned labels are every
integer2–12 together with16,24,32. No all-label claim is made for the other
frames or new–new contacts.

Only three multiple-endpoint seeds survive numerically. All terminate
`7/12` and `12/52`; their7-labels are3,4,5 and their12/52 contacts are right.
The label3 seed is exactly the old47 in the numerical coordinates, which is
explicitly checked rather than counted as a new wall. The label4/5 seeds have
44/41 P candidates but no compatible Q candidate in the assigned domain.

For the label3 seed,49 P and3 Q candidates pass their individual gates.
Their147 products give62 dotted pairs and85 non-Coxeter products; none is
finite. Each of the62 dotted pairs is then used as a new frame `(old30,P,Q)`
for R, with both required finite endpoint corners. This produces17,944 real
branches, of which8,473 pass the63 retained native-wall gates and center.
The remaining L/R relation is dotted in8,113 cases and non-Coxeter in360;
it is never finite in the saved domain.

The failure of the prescribed disk is not used to discard alternative net
improvements. **Every one of the8,113 coherent four-cap collections is scored
by the complete numerical S6 triple-intersection ledger.** All give exactly
31 finite rays,7 exterior rays,41 compact edges,11 one-open edges,2 two-open
edges, and3 entirely exterior edges. Hence all have15 open ends. The B
boundary has moved to new facets without reducing its number of open ends.
No new physical normals, ambient finite census, or finite-volume claim is
promoted from these scouts.

The current12–52 pi/6 bond is a concrete source of rigidity: it forces the
left cap right to both sides before other old-wall gates are considered.
Changing that bond through a coordinated replacement of12 is a different
remaining family. This report neither tests nor excludes it.

## Status and reproduction

**FINAL/FROZEN.** Final generation and a fresh saved replay both passed,
together with source hashes, syntax and report links.
The old-deletion and abstract-disk computations are exact integer
incidence checks; the eight-profile lemma is exact elementary algebra. All
normal construction, Coxeter recognition, and projective-boundary searches
are numerical diagnostics, not geometry certificates or nonexistence proofs.

```sh
python scripts/fresh_nonarithmetic_release47_cycle_patch_20261017.py
python scripts/fresh_nonarithmetic_release47_cycle_patch_20261017.py --verify
```
