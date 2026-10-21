# Connected white extensions of arithmetic triangle stars

This package enumerates the connected graph domain relevant to a proposed
compact arithmetic dimension21 weight rule over Q(sqrt5). It does not
validate a weight vector or prove a dimension bound. Common orthogonal
spectators and the marked containing facet must still be evaluated for
every retained state, including states that have further children.

The proposed pruning is conditional on exactly zero distance-two discount,
an adjacent3-edge discount between0 and2/3, and nonnegative compensation
and cross payments. These are hypotheses on the eventual exact rule;
they are not inferred from a numerical optimizer.

The [source](../scripts/fresh_sqrt5_triangle_extension_complete_20261016.py)
binds the frozen [structural reduction](fresh_sqrt5_discounted_face_domain_reduction_20261016.md),
[triangle baseline theorem](fresh_sqrt5_triangle_weights_20261015.md), and
[exact Q(phi) principal filters](fresh_sqrt5_quad_filters_20261015.md).
The earlier rank9 pilot remains separate bounded evidence.

## Why these seeds and this rank bound suffice

For a triangular two-face in a compact simple dimension21 Coxeter chamber,
write W for the nineteen containing white facets and b0,b1,b2 for its
three boundary facets. Each W together with two black roots is spherical.
The full star is indefinite. A minimal indefinite subset must therefore
contain all three black roots; compactness rules out parabolic subsets.
It is a Lanner diagram, of order3,4 or5. The frozen Q(sqrt5) catalogue
gives exactly the76 marked seeds used here.

The full star has twenty-two roots. This is a dimension-specific bound,
not an arbitrary cutoff on classical components. A connected extension
can be ordered so that each added white root touches an earlier root,
by following a spanning tree from the seed. Components orthogonal to
the entire grown graph are common spherical spectators. Their rank and
facet-mark choices are reserved for the later payment calculation.

Every current graph still contains its indefinite seed and has a
positive-definite principal submatrix obtained by omitting any black root.
Consequently its physical inertia is exactly `(nodes-1,1,0)`; in
particular it is nonsingular. The Q(phi) determinant remains nonzero at
the conjugate embedding, so conjugate positive semidefiniteness here is
equivalent to positive definiteness. The source checks these statements
with exact one-dimensional Schur remainders.

## Exhaustive extension choices and safe pruning

The structural reduction proves that every white vertex has degree at
most3. Degree4 would force an induced affine D4, whose positive null
combination is a feasible ideal point of the compact chamber. The finite
angle alphabet and spherical corner constraints give these choices for
an added vertex:

* one neighbor: label3,4,5,6 or10;
* two neighbors: labels33,34,43,35 or53;
* three neighbors including a white root: all labels3;
* the three black neighbors: all labels3, or one label4 or5 and two3s.

Every profile is then tested on all three full corner graphs. Their
connected components must be literal finite Dynkin diagrams. Remaining
profiles must pass the cyclic squareclass condition, the exact physical
Lorentzian test, and the conjugate positive-definite test. These are
necessary filters; passing a local graph does not assert its realization
inside a complete chamber.

Let B be the sum of distance-one/two baseline payments plus the two
entire-H4 corrections, and k the immutable number of black3-edges.
The source prunes only B>=4, B=3 with k<2, or B=2 with k=0. The inherited
baseline is extension-monotone. Three black3-edges would be affine A2,
excluded by compactness. All other states remain in the search even when
their present corner components have rank below6: extensions can make
them eligible for a discount.

## Quotients, persistence and obligations

All three black roots have the same color; all white roots have a different
color. A coarse degree/label bucket only selects possible matches. Every
identification then requires an exact full labeled graph isomorphism,
including right contacts, preserving the two colors. The literal map is
saved. Thus black permutations are allowed and no unproved graph
fingerprint is used as an equality test.

Each completed core and all newly registered children are committed in
one SQLite transaction. The export retains every state, the76 seed maps,
complete trial-domain accounting and accepted-child maps. A fresh saved
reconstruction checks each profile, each rejection, each literal map and
reachability of every registry state from a seed.

For every unpruned state of n nodes, the later evaluator must test all
common spectator forests of total rank `22-n`, with the same forest at
all three corners. It must consider each literal white facet mark and
every allowed external facet-mark type, using the actual corner component
counts in all denominators. **Expanded states are also zero-additional-
contact obligations.** Checking only maximal or rank22 states would be
insufficient. The safe prunes cover arbitrary spectators only under the
explicit discount/nonnegativity hypotheses above.

## Terminal enumeration

The [full export](fresh_sqrt5_triangle_extension_complete_20261016_rank22.json)
has996 colored graph classes and no pending state:

| Terminal role | Count |
|---|---:|
| Expanded, also a zero-additional-contact obligation | 737 |
| Hereditary-safe under the stated coefficient assumptions | 223 |
| Unpruned full rank22 stop | 36 |

All649,625 generated contact profiles are accounted for:648,062 fail a
spherical corner,2 fail the cyclic field condition,253 have a negative
conjugate direction, and1,308 lead to accepted children. The76 seed maps
and1,308 child maps are literal full-label/color-preserving maps.
Rejected nonspherical profiles are saved by their exact indices in each
core's exhaustive ordered profile list; all other trial records retain
their explicit contact profile and arithmetic decision.

There are773 unpruned payment obligations before common spectators and
facet marks are added. Of these,501 currently have `(B,d)=(2,1)` and94
have `(B,d)=(3,2)`; the others are retained because their black3 edges can
become eligible under connected growth. Their zero-addition payments must
still be checked rather than inferred from descendant payments.

The [replay source](../scripts/fresh_sqrt5_triangle_extension_complete_20261016_replay.py)
reconstructs all profile/filter decisions, all literal maps and reachability,
then checks the closed SQLite database against the complete JSON export.
Generation, the complete fresh reconstruction, and a fresh saved replay
of the [replay certificate](fresh_sqrt5_triangle_extension_complete_20261016_replay.json)
all passed. The latter checks all649,625 profiles,1,384 literal maps and
the terminal database/export equality. This is a same-producer
verification, not independent arithmetic or a completed payment proof.

The bounded pilot already exposed new shortfalls in a proposed numerical
weight vector. Those are failures to cover this necessary local domain,
not counterexamples to the earlier frozen compact bound. They do not
assert global realizability of every retained local star. No completed
payment or dimension21-exclusion claim follows from this catalogue.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_sqrt5_triangle_extension_complete_20261016_replay.py --verify
```
