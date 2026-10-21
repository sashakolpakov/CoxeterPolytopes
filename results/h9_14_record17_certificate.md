# H9/14 record 17: exact eight-orbit whole-shadow residual

## Result

For the second-layer missing-edge record

```text
index 17
graph6 M???????????o?g?_
edges  0-12, 1-12, 0-13, 2-13, 12-13
```

the complete blocker-and-ridge reduction is finite and exactly enumerable.
After every currently proved compact-H9/14 nerve condition is imposed, there
are **21 canonical-`X` labelled representations and exactly eight incidence
complexes modulo the automorphism group of the missing-edge graph**:

| missing four-faces `a` | labelled representations | incidence orbits |
|---:|---:|---:|
| 0 | 1 | 1 |
| 1 | 17 | 4 |
| 2 | 3 | 3 |
| 3 | 0 | 0 |
| 4 | 0 | 0 |
| 5 | excluded by strict link `g4` | 0 |
| 6 | excluded by strict link `g4` | 0 |

All eight residual complexes have every ridge in exactly two facets, connected
dual graph, the scheduled face vector, no minimal nonface of order above five,
and every face link has the mod-two homology of the appropriate sphere.

This is a sharp incidence theorem, not yet a Coxeter example.  A mod-two
homology sphere need not be PL; this certificate supplies neither convex
coordinates nor Coxeter labels nor a Lorentzian Gram matrix.

## 1. The graph, support, and schedules

Put

```text
W = {0,...,11},       I = {3,...,11},       S = {12,13}.
```

The three inclusion-minimal vertex covers of the missing-edge graph are

\[
C_0=\{12,13\},\qquad C_p=\{0,1,13\},\qquad
C_q=\{0,2,12\}.
\]

The independently reconstructed low-skeleton counts are 310 independent
triples, 735 independent four-sets, and 1,212 independent five-sets.  The
seven Dehn--Sommerville schedules are:

| `a` | all missing independent five-sets | facets | summed vertex-link `g4` |
|---:|---:|---:|---:|
| 0 | 8  | 110 | 42 |
| 1 | 14 | 108 | 36 |
| 2 | 20 | 106 | 30 |
| 3 | 26 | 104 | 24 |
| 4 | 32 | 102 | 18 |
| 5 | 38 | 100 | 12 |
| 6 | 44 | 98  | 6  |

Strict `g4>=1` in each of the fourteen vertex links excludes `a=5,6`.
The exact local-link interval at each of vertices 12 and 13 is `[0,0]`.
Consequently every missing four-face lies on `W` and must meet both
`{0,1}` and `{0,2}`; equivalently, it contains 0 or contains both 1 and 2.
The same inherited special-link bound places every additional minimal
five-face on `W`.  There are 495 four-sets on `W`, 201 satisfying the stated
support rule, and 792 five-sets on `W`.

## 2. All facet-complement types

Every five-element graph cover is of exactly one of the following forms:

\[
\begin{aligned}
Y_T&=\{12,13\}\cup T,&T&\in\binom W3 &&(220),\\
P_A&=\{0,1,13\}\cup A,&A&\in\binom{\{2,\ldots,11\}}2 &&(45),\\
Q_B&=\{0,2,12\}\cup B,&B&\in\binom{\{1,3,\ldots,11\}}2 &&(45).
\end{aligned}
\]

Thus the archived 310-element ridge-support core is not merely an upper
bound here: it is exactly the set of five-element graph covers.

Let `M` be the clutter of missing four- and additional minimal five-faces.
For `C_p`, let `Z_p` consist of the rows of `M` avoiding `{0,1}`.  If `T` is
a minimal transversal of `Z_p`, then `C_p union T` is a minimal transversal
of the full missing-face clutter.  The three elements of `C_p` remain
essential, witnessed respectively by the private graph edges `0-12`,
`1-12`, and `2-13`; elements of `T` are essential by minimality inside
`Z_p`.  Purity therefore forces every member of `b(Z_p)` to have order two.
The identical argument for `C_q` uses private edges `0-13`, `2-13`, and
`1-12`.

The one-special ridge equations say that every triple contains zero or two
selected `p` edges, and likewise for `q`.  Hence both selected graphs are
complete bipartite cuts.  Blocker involution identifies the two sides of
each cut with the rows of `Z_p` or `Z_q`.  Every such row has order at most
five, while the two side orders sum to ten.  Therefore both cuts are
necessarily `K5,5`.  This is an exact theorem, not a guessed balanced ansatz.

## 3. The complete ridge system

The replay enumerates all `binom(14,6)=3003` six-sets.  Exactly 735 groups
are nonempty:

* 120 `p` triangle equations;
* 120 `q` triangle equations;
* 495 equations containing both special vertices.

There are exactly 1,980 `y`, 405 `p`, and 405 `q` incidences.  For
`R in binom(W,4)`, the both-special equation is

\[
\sum_{T\in\binom R3}y_T
+\mathbf 1_{\{0,1\}\subset R}\,p_{R\setminus\{0,1\}}
+\mathbf 1_{\{0,2\}\subset R}\,q_{R\setminus\{0,2\}}
\in\{0,2\}.
\tag{R}
\]

Normalize the `p` cut by colouring vertex 2 zero and writing
`P subset I`, `|P|=5`, for the opposite-colour isolates.  Normalize `q`
similarly at vertex 1 with a five-set `Q subset I`.  Define a graph `X` on
the nine isolates by

\[
x_{ij}=y_{0ij}.
\]

Equation (R) is equivalent to the reconstruction formulas

\[
\begin{aligned}
y_{1ij}&=x_{ij}\mathbin\oplus p_{ij},&
y_{2ij}&=x_{ij}\mathbin\oplus q_{ij},\\
y_{12i}&=\mathbf1_{i\in P\triangle Q},&
y_{ijk}&=x_{ij}\mathbin\oplus x_{ik}\mathbin\oplus x_{jk},
\end{aligned}
\]

with `y_01i=y_02i=y_012=0`, together with these exact restrictions:

1. `X` is triangle-free and induced-`2K2`-free;
2. switching `X` by the cut `P` is triangle-free;
3. switching `X` by the cut `Q` is triangle-free;
4. on `D=P triangle Q`, `x_ij` equals the `P`-cut edge indicator.

The production and hostile replays reconstruct the selected five-complements
and check all ridge degrees directly, so the formulas are not used as an
unchecked shortcut.

## 4. Exact census

`nauty geng -q -t 9` supplies one representative of each of the 1,897
unlabelled triangle-free graphs of order nine.  Independent graph6 decoding
and property replay retain 156 induced-`2K2`-free graphs.  Every labelled
solution can be relabelled to one archived `X`, carrying `P,Q` with it, and
all `binom(9,5)^2` cut pairs are then tested.  Exactly 11,545 pairs satisfy
the coupled ridge characterization.

The complete attrition table is:

| gate | `a=0` | `a=1` | `a=2` | `a=3` | `a=4` |
|---|---:|---:|---:|---:|---:|
| target `|y|=60-2a` | 1434 | 1376 | 1455 | 1869 | 521 |
| every independent triple is a face | 654 | 584 | 159 | 136 | 21 |
| exactly `a` minimal four-nonfaces | 36 | 82 | 76 | 53 | 19 |
| exact missing-five count and support | 36 | 82 | 76 | 53 | 19 |
| complete blocker has no row above five | 1 | 17 | 3 | 0 | 0 |
| exact scheduled full face vector | 1 | 17 | 3 | 0 | 0 |

For each state at the blocker gate, all `2^14` subsets are classified as
faces or nonfaces, and the complete inclusion-minimal blocker is recovered.
This is what excludes `a=3,4`; no bounded solver result is involved.

### Correction of the exploratory uniqueness claim

An earlier exploratory filter incorrectly treated the schedule entry
`missing_five_count_offset = 8+6a` as the number of *inclusion-minimal* A5
rows.  It is the number of all missing independent five-sets.  When `a>0`,
that includes five-sets in the upper shadow of an A4.  The correct direct
face-table replay retains the 17 `a=1` and three `a=2` representations shown
above.  All of them have the scheduled full face vector and pass the local
and link audits.  The frozen certificate therefore supersedes the erroneous
“unique `a=0`” wording.

## 5. The eight automorphism orbits

The missing-edge graph has automorphism group `S9 x C2`: `S9` relabels the
isolates, and `C2` exchanges `(1,12)` with `(2,13)` while fixing vertex 0.
The graph `X`, defined by `x_ij=y_0ij`, is invariant under the `C2`.  Hence
distinct unlabelled `X` words below cannot become equivalent.  Within a
fixed word, an exact coloured-graph backtrack quotients `(P,Q)` by `Aut(X)`
and branch interchange.

| `a` | `X` graph6 | representative local `P=Q` | labelled multiplicity | facets | blocker histogram `(order:count)` | facet hash prefix |
|---:|---|---|---:|---:|---|---|
| 1 | `H???FB_` | 0,1,2,3,4 | 6 | 108 | 2:5, 4:1, 5:6 | `719e1b41557b` |
| 1 | `H???FB_` | 3,4,5,7,8 | 4 | 108 | 2:5, 4:1, 5:6 | `72b5451251f7` |
| 2 | `H???F?^` | 3,4,5,6,7 | 1 | 106 | 2:5, 4:2, 5:4 | `13e150f81621` |
| 2 | `H???Fbo` | 4,5,6,7,8 | 1 | 106 | 2:5, 4:2, 5:4 | `bb8789352097` |
| 0 | `H??EEBo` | 4,5,6,7,8 | 1 | 110 | 2:5, 5:8 | `3580d27d6ea9` |
| 1 | `H??EEBw` | 0,1,2,3,4 | 1 | 108 | 2:5, 4:1, 5:6 | `7f077b63f896` |
| 1 | `H??FFB_` | 0,1,2,3,4 | 6 | 108 | 2:5, 4:1, 5:6 | `267972a560ae` |
| 2 | `H??FF?^` | 3,4,5,6,7 | 1 | 106 | 2:5, 4:2, 5:4 | `d764cd7e0333` |

The local isolate labels 0 through 8 correspond to global vertices 3
through 11.  For the `a=0` representative the eight A5 masks are

```text
121, 122, 124, 481, 1921, 2078, 3590, 3968
```

or, as global vertex sets,

```text
{0,3,4,5,6}  {1,3,4,5,6}  {2,3,4,5,6}  {0,5,6,7,8}
{0,7,8,9,10} {1,2,3,4,11} {1,2,9,10,11} {7,8,9,10,11}.
```

## 6. Topology audit and exact scope

For every orbit representative the verifier independently checks:

* all selected complements lie in the archived 310-element support core;
* every active ridge has degree two;
* the dual graph is connected;
* the complete face vector equals its Dehn--Sommerville schedule;
* the complete minimal-nonface blocker has orders only 2, 4, and 5;
* every strict local `g` and local Macaulay inequality used by the global
  selection model;
* every link Euler characteristic; and
* every reduced Betti vector of every face link over `GF(2)`.

All eight pass.  This makes them genuine, sharply specified candidates for
the next primal step.  It does **not** prove any of the following:

1. PL-sphericity;
2. convex polytopality;
3. a compact Coxeter labelling or Lorentzian Gram realization;
4. finite volume in the noncompact sense.

The compact metric gate is now closed for all eight orbits by the
[all-orbit Coxeter certificate](h9_14_record17_all_coxeter_label_gate.md):
complete order-four/five Lanner joins and every actual rank-nine vertex leave
zero states on each orbit.  A finite-volume noncompact conclusion would still
require Vinberg's parabolic face/cusp criterion; the ridge, homology, and
compact ellipticity tests do not substitute for that criterion.

## 7. Companion primal and Coxeter results

Three frozen certificates continue the incidence theorem.

First, [the ordinary-polytopality certificate](h9_14_record17_polytopality.md)
constructs its inverse-multiwedge base from `C(8,5)` by two exact rational
stellar subdivisions and then restores the four clone pairs.  Thus the
selected `a=0` nerve is the boundary of a rational simplicial convex
9-polytope.  Its polar gives a rational simple polytope and hence a projective
`Q`-factorial toric carrier.  This is a genuine algebro-geometric incidence
object, but not a Coxeter example.

Second, the original [selected compact-Coxeter label certificate](h9_14_record17_coxeter_label_gate.md)
proves that this same selected nerve cannot be a compact hyperbolic Coxeter
nerve.  The three A5 supports with common core `{3,4,5,6}` have complete
Lanner join census `420 -> 1332 -> 4620`; exact actual-vertex restrictions
leave 80 states at rank five and zero at rank six.  Since the decisive union
contains no dotted pair, Perron--Frobenius monotonicity makes the obstruction
unrestricted in all remaining finite labels.  That certificate is
selected-nerve-only and compact-only.

Third, the
[all-orbit compact-Coxeter certificate](h9_14_record17_all_coxeter_label_gate.md)
reconstructs every one of the eight residual nerves and independently joins
the complete 108/420 order-four/order-five Lanner tables.  Five orbits die
after three supports and the remaining three after four.  This upgrades the
metric conclusion to the whole graph shadow, while remaining compact-only.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/enumerate_h9_14_record17_x_graphs.py --verify-only
PYTHONPATH=src:. python3 scripts/independent_h9_14_record17_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record17.py --verify-only
PYTHONPATH=src:. pytest -q tests/test_h9_14_record17_certificate.py
```

Frozen canonical digests:

```text
X archive       6264998333719c38b3406d3d12b0abef911d100c07be5ec01f4a97f19e084879
hostile summary d30dfc3d25b5e4bb0c5067f96fabfa863f843167d02774596b1efffd26adb06d
certificate     b98e58b390cbf8d9b40e08da16da29f78a2f04805bb0b13d7f2efb623b9afb51
```
