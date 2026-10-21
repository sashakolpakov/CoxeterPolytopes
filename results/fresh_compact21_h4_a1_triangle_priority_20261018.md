# H4 five-edge triangles and an actual funded charge transfer

In the compact arithmetic dimension-21 golden-field setting, a vertex of type
A1+H4+2E8 admits a new charge rule. For any fixed rational 0 <= t <= 1, remove
t from the extra H4 payment on its label-five edge. Restore t at one specified
angle of each triangular face that this deletion would otherwise underpay.
Only the four source forests listed below can receive these restorations.
The resulting upper budgets are

| Forest | New budget in terms of its old total W |
|---|---:|
| A1 + H4 + 2E8 | W - t |
| A1 + A4 + 2E8 | W + t |
| A5 + 2E8 | W + 2t |
| A1 + B4 + 2E8 | W + t |
| B5 + 2E8 | W + t |

All other forest budgets are unchanged. This is a face-payment and vertex-budget
lemma. It is not yet a solution of the resulting weight model, and it gives
no new dimension ceiling by itself. It requires no E8 triangle activation.

## Assumptions and normalization

The old compact arithmetic dimension-21 field reduction gives k=Q(sqrt(5)),
with integer ring Z[phi], phi=(1+sqrt(5))/2. We inherit compact acute simplicity,
the finite Coxeter classification, and arithmetic cyclic integrality and
conjugate positivity. The pinned old weight model uses ordinary unit payment
at connected marked pairs of graph distance at most two, plus an H4 unit for
the five-edge and the endpoint pair. Added old bonus weights are nonnegative.
No fixed D5 half-unit is used.

Number H4 as 0--1--2--3, with labels 5,3,3, and let the isolated A1 root be 4.
The marked pair is (0,1). A triangular face at this pair has a third wall c,
numbered 5. Both other finite corners contain each full untouched E8. A full
finite E8 has no connected finite extension, so c is orthogonal to both E8
factors. Its contacts with H4 and A1 are considered together from the start.

Normalize each root to square two. Write c=p+sqrt(H)e, where e is the old
future unit timelike vertex vector. Its contact vector determines the spatial
projection p and H=p.p-2. A third wall containing the old vertex strictly in
its halfspace has H>0. Conjugate positivity requires H-conjugate <= 0. These
necessary tests do not assert that every surviving wall extends to a compact
polytope.

## High bonds are safe or impossible, without a sampled cutoff

In any finite Coxeter diagram a label at least six lies in an isolated I2(m).
If c has such a bond to white H4 root 2 or 3, an old neighbor is retained in
at least one corner, contradicting that isolation. A high bond to marked root
1 is likewise impossible because white root 2 remains. A high bond to the whole
A1 forces c orthogonal to every H4 root in the union of the two corners; then
H=(2 cos(pi/m))^2/2-2<0, impossible.

A high bond to marked leaf 0 is possible in principle. At the corner omitting
1 it gives the ordinary unit at the adjacent pair (0,c). The old H4 corner
retains its ordinary five-edge unit after the deletion. These two units already
pay the triangle. Removing any extra H4 five-edge payments at other corners
cannot remove either ordinary unit. Thus every high-bond case is safe for the
proposed transfer. All cases that could need a restoration lie in the complete
bounded domain {2,3,4,5}^5.

## The complete bounded domain and all simultaneous deletions

The producer checks all 1,024 contact vectors. It tests both finite corners,
positive height, field membership of every simple cycle inside the old H4,
conjugate height, and connected affine principal obstructions containing c.
Cycle products include old path factors in Q(sqrt(5)), so the displayed
contact-product field test is equivalent for this purpose. Components without
c were already positive definite. A connected affine Coxeter diagram has no
label-five edge, which justifies omitting those subsets in the affine check.

There are 46 finite-eligible vectors: four have nonpositive height, eight fail
cycle-field membership, and 34 survive the necessary arithmetic gates. The
JSON records every eligible vector and every surviving root. The test computes
all three triangle payments after deleting **every** H4 five-edge extra at all
three corners, including corners other than the original one. Exactly four
patterns then have payment below two:

| c contact with H4 | c contact with A1 | H | Forest at the corner omitting 0 |
|---|---|---|---|
| label 3 at root 3 | orthogonal | 1+sqrt(5) | A1 + A4 + 2E8 |
| label 3 at root 3 | label 3 | 3/2+sqrt(5) | A5 + 2E8 |
| label 3 at root 3 | label 4 | 2+sqrt(5) | B5 + 2E8 |
| label 4 at root 3 | orthogonal | 4+2sqrt(5) | A1 + B4 + 2E8 |

In each pattern the old baseline corner payments are (2,0,0); the only removed
H4 five-edge extra is at the original vertex. Restoring t at the corner
omitting 0 changes the sum from 2-t to exactly two. Every other surviving
pattern already has payment at least two after all possible deletions; partial
deletion for t<=1 is therefore safe. This also covers triangles with several
H4 corners and rules out double counting their losses.

The first pattern is the genuine Lanner chain 5--3--3--3. It cannot be excluded
as a spurious arithmetic root. All four patterns have integral length-two
products with both E8 triangle-wall types, so that simple gate cannot fund the
deletion. This packet does not claim that all longer cyclic products of every
pattern were independently certified, or that the full chamber is compact.
The transfer proof needs only the necessary complete list above.

## The restoration rule and source angle counts

At a triangle of one of these four types, choose its corner q obtained from the
H4 vertex p by omitting old root 0. The restored angle at q is marked by old
root 1 and c. The two E8 factors are untouched whole spectators at q, so the
stated source forest is exact. No assertion that E8 triangle activity persists
at q is needed. The restored angle has ordinary baseline zero in every case.

The finite diagrams identify its possible unordered marked pairs intrinsically:

* In A4, it is the unique pair of endpoints. The two possible incoming H4 edges
  use this same angle. An angle at a simple vertex belongs to one two-face, so
  its restoration is counted once. The capacity is one, rather than the edge
  indegree bound two.
* In A5, number the path 0--1--2--3--4. The possibilities are (0,3) and (1,4),
  hence capacity two.
* In B4, it is the endpoint pair, hence capacity one.
* In B5, number the label-four edge (0,1), followed by the three label-three
  edges. The pair is (1,4), hence capacity one.

These counts bound actual restored angles; they do not presume that every
allowed reverse extension occurs. Define the restoration only once at any
angle whose triangular face has the specified deficient H4 corner. No source
angle receives multiple additions. This gives exactly the five budget rows
at the start of the report.

## Small-face inequalities and compatibility with the previous successor

For a k-gon the required sum is at least 5-k. A modified triangle is covered
above. A quadrilateral containing a deleted H4 angle keeps the ordinary unit
at that angle, which alone pays its requirement one. Quadrilaterals without a
deletion retain the old valid payment. Faces of length at least five need only
nonnegative payments; deletion of t<=1 from the old two-unit angle preserves
nonnegativity. Restorations are nonnegative everywhere.

The earlier active/inactive D5+2E8 successor composes with this rule. Its deleted
D5 fork angle is never triangular, and the inactive E8 deletion concerns pairs
inside E8, whereas the untouched whole E8 factors are white throughout each
triangle considered here. Thus no earlier deletion affects these triangular
payments. For quadrilaterals, this H4 rule preserves every ordinary basic unit,
so it preserves the retained-unit arguments of the prior successor. The H4
endpoint extra, all I2(5)/H3 triangle routes and all old nonnegative bonuses
remain available.

The old unpaid quadrilateral bonus rows can be retained literally: any face
where the new H4 deletion acts is already paid by its retained ordinary unit.
The new conditional restorations are positive and require no additional face
constraint. Thus a support successor may add one variable t with 0<=t<=1,
decrease the H4 target forest coefficient by t, and increase the four source
forest bounds with capacities 1,2,1,1.

## Verification and the current numerical target

The standard-library producer passed in terminal chunk `54fa96` in 0.008
seconds, with 29,196,288 bytes peak RSS, under the shared lock and 59-second
alarm. It uses exact Fraction/ordered-field arithmetic and imports no optimizer.
The complete bounded-domain certificate is
`fresh_compact21_h4_a1_triangle_priority_20261018.json`.

The separate streaming budget probe evaluates the saved support primal at all
five forests and checks its reconstructed H4 and A4 rows against the original
literal rows. Its source-angle bonus values are all zero. It finds old totals
40.00447918 at H4 and A4 sources, 39.50447918 at A5/B5, and 38.00447918 at B4.

The exact repaired dual has new t-column load
35161798003/1000000000247>0, so it does not survive this added variable. Its
simple repair using the t<=1 row lowers its value to
39969317248112/1000000000247<40. This defeats that particular obstruction only.
No new optimization or full-model primal verification is claimed here.

```sh
python3 -B scripts/fresh_compact21_h4_a1_triangle_priority_20261018.py --verify
python3 -B scripts/fresh_compact21_h4_transition_budget_probe_20261018.py --verify
```

Verification status: internal exact enumeration plus written geometric and
charging proof; independent review pending. No Lean verification or dimension
improvement is claimed.
