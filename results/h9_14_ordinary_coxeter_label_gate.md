# Exact Coxeter-label gate for the five earlier ordinary-polytopal H9/14 nerves

> **Scope update.** Record 35 is now a sixth ordinary-polytopal nerve.  It is
> not part of the five frozen families in this report and has a separate
> unrestricted compact-Coxeter obstruction in
> `results/h9_14_record35_coxeter_label_gate.md`.

## Result

The five ordinary-polytopal H9/14 nerves known at this proof stage admit Coxeter
labels satisfying every local finite/Lanner condition.  Thus none of the five
is eliminated at the discrete-label layer.

The most natural sparse families can be enumerated completely up to the full
nerve symmetry group.  They are nevertheless all globally impossible: their
right-angled cross-factor pairs split the Gram matrix into orthogonal pieces
containing at least two negative directions.  More generally, every genuine
Gram candidate needs at least `4,2,1,1,1` non-right edges across the
minimal-nonface components of records `30,42,43,44,47`, respectively.

This is a useful primal--dual movement of the boundary, but it is not a new
hyperbolic Coxeter polytope.  The unrestricted label CSP, dotted entries,
rank ten, and inertia `(9,1,4)` remain open.

**One-edge upgrade.**  The next cross-factor layer has now also been
exhausted exactly.  Although 6,342 one-edge symmetry orbits pass the local
finite/Lanner gate, every one has a principal Gram submatrix of negative
index at least two.  Hence at least two new cross-factor non-right edges are
necessary relative to each frozen family.  See
`results/h9_14_one_cross_edge_gate.md`.

## 1. The complete finite-label domain

The search uses labels

\[
                 m\in\{2,3,4,5\}.
\]

It is not restricted to the earlier golden alphabet `{2,3,5}`.  The exclusion
of larger finite labels is theorem-driven: Burcroff's Lemma 10.10 says that a
compact Coxeter 9-polytope with 14 facets has neither a missing triple nor a
multi-multiple edge.  A multi-multiple edge is precisely a finite Coxeter edge
with `m >= 6`.  See A. Burcroff,
[Near Classification of Compact Hyperbolic Coxeter d-Polytopes with d+4
Facets and Related Dimension Bounds](https://arxiv.org/abs/2201.03437),
Lemma 10.10.

The verifier recognizes every finite connected Coxeter diagram in the
classification `A,B,D,E,F,H,I2`, component by component, using integer labels
rather than floating-point eigenvalues.

## 2. Complete Lanner generation

The order-four and order-five Lanner tables were regenerated exhaustively,
not merely copied from a list.  For every complete pair-label word over
`{2,3,4,5}` the computation:

1. checks every one-node deletion with the complete finite Coxeter
   classifier;
2. rejects a finite full diagram;
3. computes the exact Gram determinant in
   `Q(sqrt(2),sqrt(5))`; and
4. separates negative-determinant Lanner words from zero-determinant affine
   words.

One fixed finite deletion is used as an exhaustive pruning step.  Interlacing
shows that a word with all one-node deletions positive definite can have at
most one non-positive direction, so the exact determinant dichotomy is
complete.

| order | all words | finite fixed-deletion prefixes | minimally nonfinite | affine | Lanner placements | isomorphism types |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 4,096 | 25 | 135 | 27 | 108 | 9 |
| 5 | 1,048,576 | 182 | 677 | 257 | 420 | 5 |

The generated placement sets agree exactly with the permutation closures of
the classical nine and five Lanner representatives.

## 3. Exact censuses inside the sparse ansätze

The phrase "symmetry orbit" below always means an orbit under the full
automorphism group of the complete minimal-nonface hypergraph, equivalently
of the completed nerve.  It does not mean that all labels in a pair orbit were
artificially equated.

| record | vertices checked | higher minimal nonfaces | sparse labelled solutions | full-nerve symmetry orbits | global fate of the whole sparse family |
|---:|---:|---:|---:|---:|---|
| 30 | 128 | `2 x A4` | 11,664 | 45 | inertia exactly `(9,5,0)` |
| 42 | 100 | `A4 + A5` | 46,448,640 | 6,120 | negative index at least 2 |
| 43 | 90 | `4 x A5` | 12,600 | 15 | negative index at least 2 |
| 44 | 90 | `4 x A5` | 60,480 | 25 | negative index at least 2 |
| 47 | 90 | `4 x A5` | 7,560 | 5 | negative index at least 2 |

Here `A4` and `A5` denote supports of order four and five, not finite
Dynkin types.

### Records 30 and 42: product ansätze

Record 30 is the nerve of
`Delta_1^3 x Delta_3^2`.  Put an arbitrary Lanner diagram on each
four-facet simplex block and make every cross-block meeting pair right.  The
minimal-nonface automorphism group is

\[
 (S_2\wr S_3)\times(S_4\wr S_2),\qquad |G|=55,296.
\]

Internal `S4` permutations reduce each block to one of nine Lanner types, and
the outer `S2` makes the pair unordered.  The complete orbit count within the
ansatz is therefore

\[
 \binom{9+2-1}{2}=45.
\]

Record 42 is the nerve of `C_5 x Delta_3 x Delta_4`.  Choose one of nine
order-four types, one of five order-five types, and label the five meeting
edges of the polygon factor arbitrarily by `{2,3,4,5}`.  The full group is
`D5 x S4 x S5`, of order 28,800.  Burnside gives 136 four-color bracelets,

\[
 \frac{4^5+4\cdot4+5\cdot4^3}{10}=136,
\]

so the complete ansatz has `9*5*136 = 6,120` orbits.

In both records every polytope vertex diagram is an orthogonal union of
proper deletions of the chosen Lanner blocks, plus rank-one or rank-two finite
factors.  This proves the local gate for every member; the archived verifier
also checks all 128 and 100 vertices directly.

### Records 43, 44, and 47: multiwedge ansätze

Each multiwedge nerve has four missing five-faces.  One is the isolated
`boundary(Delta_4)` join factor; the other three overlap.  The exact search
puts one of all 420 order-five Lanner placements on every support, joins the
three overlapping tables by equality on their common pairs, labels every
remaining meeting pair by two, and tests all residual vertex diagrams.

The complete 90-vertex list is independently recovered as five deletions of
the isolated block times 18 residual vertices.  Hence the table join loses no
constraints and multiplying by the 420 independent isolated placements is
exact.

| record | successive overlapping-table joins | passing cluster states | all labelled solutions | automorphism group | orbit-size histogram | orbits |
|---:|---|---:|---:|---|---|---:|
| 43 | `420 -> 1,332 -> 4,212` | 30 | 12,600 | `S5 x S3 x C2` | `360^3, 720^8, 1440^4` | 15 |
| 44 | `420 -> 1,332 -> 4,620` | 144 | 60,480 | `S5 x S4 x C2` | `1440^12, 2880^11, 5760^2` | 25 |
| 47 | `420 -> 1,332 -> 222` | 18 | 7,560 | `S5 x S3 x S3` | `1080^3, 2160^2` | 5 |

Every displayed group generator is replayed on both the complete minimal
nonface clutter and all 90 vertices, and every orbit is traversed explicitly.
The exported canonical representatives are therefore ready for independent
Gram-gate regression.

## 4. The global no-go theorem for all five sparse families

The local success is not a metric success.

For record 30 the normalized Gram matrix is an orthogonal direct sum of two
order-four Lanner blocks, each of inertia `(3,1,0)`, and three dotted
two-facet blocks, each of inertia `(1,1,0)`.  For every legitimate choice of
the three dotted values the full inertia is therefore

\[
                  (9,5,0),
\]

whereas an H9 Coxeter Gram matrix on 14 facets must have inertia `(9,1,4)`.
This recovers and generalizes the earlier golden block-product
countercertificate to all nine order-four Lanner types and all dotted values
inside this ansatz.

For record 42 the mutually orthogonal order-four and order-five Lanner blocks
already supply two negative directions, independently of every polygon
dotted value.

For records 43, 44, and 47 the isolated `boundary(Delta_4)` Lanner block is
orthogonal to the residual block, which contains another Lanner principal
submatrix.  Interlacing again gives negative index at least two.

There is also an ansatz-independent strengthening.  Make a hypergraph whose
hyperedges are *all* minimal nonfaces.  Every hypergraph component is already
connected in the nonzero Gram graph: dotted pairs are connected, Lanner
diagrams are connected, and overlapping supports glue connected diagrams.
Moreover each such component contains an indefinite principal block.  If
there are `c` components, fewer than `c-1` non-right meeting edges cannot
connect them, leaving at least two orthogonal indefinite blocks and therefore
negative index at least two.  The component counts and forced lower bounds
are:

| record | minimal-nonface components | necessary cross-component non-right edges |
|---:|---:|---:|
| 30 | 5 | at least 4 |
| 42 | 3 | at least 2 |
| 43 | 2 | at least 1 |
| 44 | 2 | at least 1 |
| 47 | 2 | at least 1 |

Thus:

> In every one of these five earlier ordinary-polytopal H9/14 nerves, the tested
> factor-orthogonal label ansatz is locally Coxeter-admissible but globally
> impossible.  Any genuine solution must satisfy the unrestricted
> cross-component lower bounds `4,2,1,1,1` above.

These are exact forced-edge cardinality clauses for the next unrestricted
CSP, not merely observations about the enumerated sparse families.

## 5. Algebro-geometric interpretation

The algebro-geometric track has not been abandoned.  The toric
anticanonical-boundary constructions already produce these five divisor
incidence nerves automatically.  The present calculation says precisely what
the divisor/root lattice must add beyond the Stanley--Reisner incidence
algebra:

- the local divisor clusters can carry every required elliptic/Lanner
  decoration;
- a decomposable intersection pairing is impossible because it has too many
  negative directions; and
- the desired variety or surface must create nonzero cross-intersections
  between Lanner-bearing divisor clusters while preserving every rank-nine
  finite restriction.

So the next inverse-design variable is not merely the variety's incidence
complex.  It is an indecomposable Neron--Severi/root coupling across the
toric-looking factorization.  That is a concrete algebraic design constraint,
not a fishing expedition.

## 6. Reproducibility and scope

The exact local verifier and export APIs are in
`src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py` and
`src/coxeter_polytopes/combinatorics/h9_14_sparse_coxeter_labels.py`.

- `iter_product_orthogonal_label_representatives(...)` exports the 45 and
  6,120 product-family orbit representatives as complete meeting-pair label
  maps.
- `audit_sparse_h9_14_labels(i).orbit_representative_label_maps()` exports
  the 15, 25, or 5 multiwedge representatives.

Frozen certificates:

- `data/enumerations/h9_14_ordinary_coxeter_label_gate.json` (unified
  five-record theorem);
- `data/enumerations/h9_14_product_coxeter_labels.json`;
- `data/enumerations/h9_14_sparse_coxeter_labels.json`.

Replay scripts:

- `scripts/certify_h9_14_ordinary_coxeter_label_gate.py`;
- `scripts/certify_h9_14_product_coxeter_labels.py`;
- `scripts/certify_h9_14_sparse_coxeter_labels.py`.

The certificates distinguish three logically separate statuses:

1. the displayed sparse discrete families are `SAT`;
2. every member of those sparse families is globally Gram `UNSAT`; and
3. unrestricted labels with cross-factor non-right edges, followed by dotted
   variables and the full Lorentzian rank/inertia equations, remain open.

No claim of a new compact hyperbolic Coxeter polytope is made.
