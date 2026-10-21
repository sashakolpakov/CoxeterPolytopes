# Three K3-carried single doublings and an infinite family in `H12`

## Result

The certified 172-facet Jacobian--Kummer chamber has a much larger positive
doubling stock than its cusp count alone suggests.  Exactly 124 of its 172
facets meet every incident finite-angle facet with even Coxeter order.  Three
single doublings give pairwise nonisomorphic finite-volume noncompact Coxeter
polytopes

\[
  JK\text{-}H12\text{-}D2\text{-}W31: H^{12}/296,
  \qquad
  JK\text{-}H12\text{-}D2\text{-}W61: H^{12}/317,
  \qquad
  JK\text{-}H12\text{-}D2\text{-}W77: H^{12}/326.
\]

They are distinct without a graph-isomorphism assumption because their facet
counts differ.  Each has exactly twice the seed volume and stays on the same
rank-13 full Picard lattice

\[
 U(2)\oplus D_5(-1)^2\oplus\langle-4\rangle.
\]

Thus the same seven-dimensional very-general K3 family carries all three
descendants as exact integral divisor-reflection chambers.  This is a weak
full-Picard carrier statement, not a full-nef or automorphism theorem.

More strongly, seed walls 31 and 36 are parallel doubling walls.  Their
primitive roots satisfy

\[
 a^2=b^2=-4,\qquad (a,b)=4,\qquad (a+b)^2=0.
\]

They generate `D_infinity`, so every finite path of `N` seed chambers is a
finite-volume noncompact Coxeter polytope in `H12` with volume exactly `N`
times the seed volume and exactly

\[
  f(N)=124N+48
\]

facets.  Distinct `N` give distinct volumes and facet counts, hence infinitely
many pairwise nonisometric polytopes, all on the same K3 carrier.  The formula
is an all-parameter theorem, not an interpolation: a separate hostile audit
solves every possible wall-orbit collision for all integral translation
exponents.

The production certificate is
[`certify_jacobian_kummer_h12_single_doublings.py`](../scripts/certify_jacobian_kummer_h12_single_doublings.py),
its archive is
[`jacobian_kummer_h12_single_doublings.json`](../data/enumerations/jacobian_kummer_h12_single_doublings.json),
and the regression gate is
[`test_jacobian_kummer_h12_single_doublings.py`](../tests/test_jacobian_kummer_h12_single_doublings.py).

## Complete doubling-wall census

For a wall `w`, the script reconstructs every normalized relation from the
oriented divisor vectors.  It retains `w` exactly when no finite incident
relation has odd order.  The only finite relations in the seed are `m2`,
`m3`, and `m4`, so this is an exhaustive exact test, not a graph heuristic.

| selected-wall class | number of seed walls | seed `m2` count | seed `m4` count | parallel | ultraparallel | facets after doubling |
|---|---:|---:|---:|---:|---:|---:|
| represented by wall 31 | 12 | 46 | 20 | 45 | 60 | 296 |
| represented by wall 61 | 16 | 25 | 4 | 59 | 83 | 317 |
| represented by wall 77 | 96 | 16 | 4 | 13 | 138 | 326 |

All 124 walls are valid single-doubling inputs.  The certificate proves the
three displayed facet-count classes for the entire stock, but it only freezes
one exact representative in each class.  It does not yet assert that all
walls in one row give isomorphic labelled diagrams; a finer census may add
more than three types.

The count has a transparent check.  A double starts with two copies of the
171 nonselected seed facets.  A facet orthogonal to `w` is fixed by the
doubling reflection and the two occurrences coincide.  For the three
representatives the exact boundary-root reconstruction finds no further
identification, hence

\[
 2\cdot171-46=296,
 \qquad 2\cdot171-25=317,
 \qquad 2\cdot171-16=326.
\]

## Irredundancy and Coxeter geometry

Let `h` be the frozen strict seed interior and let `s_w` be reflection in the
selected wall.  The descendant interior is

\[
 h_D=h+s_w(h).
\]

For every nonselected seed facet `r`, the exact timelike seed witness

\[
 q_r=h-\frac{(h,r)}{r^2}r
\]

lies on `r` and strictly inside every other seed inequality.  The two-copy
boundary uses `q_r` and `s_w(q_r)`.  Direct rational replay proves that every
transported witness is timelike, lies on its target only, and has strictly
positive pairing with all other descendant roots.  The minimum off-target
slack is one in all three cases.  Hence every one of the 296, 317, or 326
displayed inequalities is a genuine facet; no boundary count relies on
genericity.

Every boundary root is primitive and its reflection is integral on the full
Picard lattice.  Every pair has nonnegative oriented product and is exactly a
Coxeter angle, parallel, or ultraparallel.  The coarse relation inventories
are:

| descendant | `m2` | `m3` | `m4` | parallel | ultraparallel |
|---|---:|---:|---:|---:|---:|
| `H12/296` | 5,051 | 264 | 1,472 | 3,397 | 33,476 |
| `H12/317` | 5,369 | 260 | 1,554 | 3,444 | 39,459 |
| `H12/326` | 5,511 | 276 | 1,585 | 3,487 | 42,116 |

The archive also records the exact normalized square of every ultraparallel
pair and every primitive oriented boundary vector.

## Volume certificate and Vinberg boundary

The seed `JK-H12/172` already has an independent GMP cone enumeration and a
literal Vinberg audit: every elliptic rank-11 edge has two finite/ideal
endpoints, every maximal affine component has its complete parabolic
neighbourhood, and there are no hyperideal rays.  The even-angle doubling
lemma now makes each displayed union a convex Coxeter polytope consisting of
exactly two seed chambers with disjoint interiors.  Consequently its volume
is exactly twice the finite seed volume.  Since the descendant reflection
group has index two in the nonuniform seed group, it remains nonuniform and
the polytope remains noncompact.

This is an equivalent exact finite-volume proof; the certificate does not
pretend to have enumerated the enormous descendant vertex sets or rerun the
literal Vinberg edge list on 296--326 facets.

## Novelty boundary and next step

The three diagrams are new to the project ledger.  Worldwide priority is
open pending a separate literature and corpus audit.  The mechanism is the
classical Coxeter doubling lemma; the positive addition is its exact use on
the Jacobian--Kummer divisor chamber, the complete 124-wall input census, the
three large boundary systems, and their common full-Picard K3 carrier.

The immediate continuation is to classify all 124 single doubles up to exact
labelled-graph isomorphism and exploit the much larger disjoint doubling-wall
graph for higher-product families.

## Independent all-parameter audit

The independent
[`parallel-redoubling hostile audit`](jacobian_kummer_h12_parallel_redoubling_hostile_audit.md)
reconstructs the rank-13 carrier and all 172 walls without importing this
production builder.  It reproduces the 124-wall doubling census and the exact
ordered W31/W61/W77 boundary lists, witnesses, Coxeter types, integral
reflections, and full-rank Smith data.

For W31/W36 it proves the `D_infinity` presentation and the unipotent identity
`(T-I)^3=0`.  An exhaustive symbolic solve over every integer translation
exponent finds 220 finite positive self-collisions and 28 identity cases, but
no negative, cross-orbit, endpoint--side, or other adverse collision.  The
surviving endpoint and side-wall orbits therefore give `124N+48` for every
`N>=1`; direct reconstruction for `1<=N<=12` agrees.  This audit also proves
that the non-endpoint roots already span the full Picard lattice, so every
descendant retains the same full-Picard weak carrier.

## Replay

```bash
PYTHONPATH=.:src python3 -u scripts/certify_jacobian_kummer_h12_single_doublings.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_jacobian_kummer_h12_single_doublings.py
PYTHONPATH=.:src python3 -u scripts/audit_jacobian_kummer_h12_parallel_redoubling_hostile.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_jacobian_kummer_h12_parallel_redoubling_hostile_audit.py
```
