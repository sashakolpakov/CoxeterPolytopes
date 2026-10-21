# D5 triangle exclusion beside A1, E7 and E8 in compact arithmetic dimension21

18 October 2026. **New complete local triangle-domain calculation and written
geometric proof.** At a compact arithmetic Coxeter H21 vertex with diagram
A1+D5+E7+E8, at most one of the D5 faces marked (0,2), (1,3), (1,4) is a
triangle. All other nonadjacent D5 marked faces are nontriangular.

The D5 numbering is 0-1-2 with leaves3,4 attached to2. This extends the
[previous D5+2E8 lemma](review_compact21_d5_two_e8_triangle_coupling_20261018.md)
to the other repeatedly supported forest in the exact dimension21 weight
dual. It does not lower its vertex allowance or the compact arithmetic
dimension ceiling. No old charging rule, unit payment, or certificate has
been modified.

The [producer](../scripts/fresh_compact21_d5_a1_triangle_coupling_20261018.py)
and [exact output](fresh_compact21_d5_a1_triangle_coupling_20261018.json) use
standard-library rational/ordered-field arithmetic and a complete six-root
contact domain. The only actual ground field at compact dimension21 is
Q(sqrt5), by the separate completed internal
[field reduction](fresh_compact21_single_field_20261016.md), which is pinned.
Finite Coxeter classification, compact acute simplicity, and Vinberg's
cyclic-integrality/conjugate-positivity criterion remain theoretical inputs.

## The untouched E7 and E8 really are orthogonal

Suppose a triangular D5 face has third wall c. Each other corner replaces
one of its two marked old D5 roots with c. Both full21-root corners are
positive definite. The old whole E8 occurs in both. Finite E8 has no
connected finite extension, so c is orthogonal to it.

If c meets the untouched whole E7, that connected finite E7+c diagram must
be E8. Maximality of finite E8 then forces c orthogonal to every other old
root in each new corner. The union of the two corners contains every old
root, so c has no other contacts. But the full old21 roots together with
c would have positive-definite Gram

  (E7+c) + A1 + D5 + E8,

of rank22, impossible in Lorentz signature (21,1). Thus c is orthogonal
to the E7 as well. This argument does not assume the remaining factors
are exceptional. It uses an entire untouched E7, and both actual corners.

Only the D5 roots, the one whole A1 and the timelike vertex direction
remain. The whole A1 is included from the start; its contact was not
assumed zero.

## Why the finite label domain is complete, including A1

Every contact is seen in a finite corner, so it is a finite Coxeter bond.
At a D5 root, a label at least six would require an isolated I2(m)
component in every corner containing it. For a white root this forces
its old neighbors into both distinct omitted singletons, impossible.
For a marked root, its old neighbors cannot all be the other mark because
the pair is nonadjacent. Thus all D5 contacts have labels at most five.

The whole A1 is isolated in the old forest, so this previous argument
alone would not bound its label. If c-A1 has label m at least six, the
finite I2(m) classification instead forces c orthogonal to every old D5
root retained in either new corner. Their union is the whole D5. After
the untouched-component result above, c then projects only to the A1
direction. With root square two its squared height is

  H=(2 cos(pi/m))^2/2 - 2 < 0,

contradicting the strictly positive height of a new wall whose halfspace
contains the old finite vertex. This excludes every high A1 label, with
no sampled cutoff or cyclotomic-degree assumption. The complete domain
is therefore {2,3,4,5}^6 for each of the six nonadjacent D5 marked pairs.

## Exact domain and arithmetic filters

The producer checks both finite corners for all6*4^6=24,576 vectors.
Its finite classifier recognizes disconnected finite diagrams through
trees, degrees, branch-arm lengths and the A/B/D/E/F/H/dihedral lists;
no numerical eigenvalue is used. The compact Gram, inverse, heights and
principal tests use exact Fraction arithmetic in Q(sqrt2,sqrt5).

There are 139 finite-eligible occurrences. The filters have the following
complete counts, in their actual execution order:

|Reason|Occurrences|
|---|---:|
|Nonpositive height|78|
|Cycle outside Q(sqrt5)|10|
|Conjugate Gram not positive semidefinite|45|
|Remaining necessary triangle roots|6|

Cycles within the connected D5 are exactly the new root plus the unique
old path between two contacts; the old path contributes only rational
label-three factors. The isolated A1 supplies no additional simple cycle.
After these cycle checks the height H belongs to Q(sqrt5), and the
conjugate full Gram requires H-conjugate at most zero by its Schur
complement. Every candidate remaining after these necessary arithmetic
gates is also checked against all connected affine principal obstructions
containing c. Subsets without c were positive definite already. A
connected affine Coxeter diagram cannot contain a label-five bond, which
justifies that subset pruning. None of the six remaining candidates has
such an obstruction.

Write s=sqrt5, phi=(1+s)/2. The complete remaining list is

|Name|D5 contact magnitudes|A1 contact|Squared height H|Marked face|
|---|---|---|---|---|
|C0|phi at0|0|(s-1)/2|(0,2)|
|L3|phi at3|0|(-1+5s)/8|(1,3)|
|U3|1 at1, phi at3|0|(23+13s)/8|(1,3)|
|L4|phi at4|0|(-1+5s)/8|(1,4)|
|U4|1 at1, phi at4|0|(23+13s)/8|(1,4)|
|C*|1 at0|phi|(s-1)/4|(0,2)|

The new C* is a necessary local arithmetic normal. Its occurrence in
this list is not an assertion that a compact H21 polytope exists or that
every local normal extends to a chamber.

## All six candidate walls are mutually incompatible

Let e be the common finite vertex's future unit vector. Each candidate
has the unique form c_i=p_i+sqrt(H_i)e, with projection p_i fixed by its
contacts. If two such actual walls occur, their product must be

  c_i dot c_j = p_i dot p_j - sqrt(H_i H_j).

The exact six-root inverse is used to compute p_i dot p_j. For every one
of the15 unordered candidate pairs it is positive and its square
exceeds H_i H_j. Thus the actual product would be strictly positive,
forbidden by acuteness.

The ten pairs among the old five agree with the earlier independent
D5 checker. For the new C*, the projection products with C0, an L root
and a U root are phi, (1+s)/4 and (5+s)/4, respectively. Their differences
of squares from the respective height products are

  3(1+s)/4,  (-7+5s)/16,  (9+5s)/16,

all strictly positive. C* therefore cannot rescue simultaneous D5
triangles. Because the three potentially triangular marked faces use
disjoint candidate sets, at most one is triangular.

This is a restriction on actual face incidence at this stated vertex
type. It does not classify its quadrilaterals or all adjacent marked
pairs. The known global charge inequality still pays quadrilaterals;
nontriangular does not imply uncharged. No change to the old exact
weight dual or its face-payment obligations has been made.

## Reproduction, resource limits and verification status

```sh
python3 -B scripts/fresh_compact21_d5_a1_triangle_coupling_20261018.py --verify
```

Generation passed in terminal chunk `5c2a9e`, exit0, in0.183 seconds.
The process's measured peak RSS was28,508,160 bytes (about27.2MiB).
It held the global proof flock, used a59-second wall alarm and a periodic
one-GiB RSS guard, and imported no NumPy, optimizer, or Lean process.
macOS rejected lowering Python's DATA/AS/RSS rlimits in this managed
environment; that first launch stopped before enumeration, and no OS
hard-memory-cap claim is made for the replacement periodic guard. Its
state is intrinsically bounded by the tiny six-root domain. No temporary
proof files or bytecode cache were produced.

Separate-process full replay passed in chunk `40c6b4`, exit0, in0.198
seconds, with peak RSS28,672,000 bytes. Independent mathematical review
is pending; generation and replay are not labelled independent review.
The full global field reduction and published theoretical inputs are
inherited, not replayed or proved in Lean here. No complete Lean or
external human verification is claimed.

Producer SHA-256:
`db7d38e57bcb34caff75a5b165cadf9ab0b78e720b3c05cb96f71ff11c5804c8`.
Certificate SHA-256:
`fad25113e98c81706a5305d3d9c8c9fa546b3035ed82ca54a5d00e43ec59202f`.
