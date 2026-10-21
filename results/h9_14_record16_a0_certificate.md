# Record 16, branch `a=0`: analytic ridge exclusion

## Result and exact scope

The `a=0` Dehn--Sommerville branch of H9/14 second-layer survivor 16,

```text
M???????????o?o?_
```

is impossible.  The proof is analytic and solver-independent: blocker
duality reduces the facets to two Boolean families, the complete ridge
system reduces those families to two ordinary graphs, and a
Cauchy--Schwarz bound contradicts the required 110 facets.

This is a branch theorem, not yet an exclusion of record 16.  Strict
vertex-link `g4` already excludes `a=5,6`; the branches `a=1,2,3,4` remain
open.  Consequently the global graph-shadow partition remains

\[
121\text{ independently impossible}
+7\text{ explicit ordinary-polytopal}
+14\text{ undecided}.
\]

## 1. The graph shadow and the schedules

Write

\[
U=\{0,1\},\qquad I=\{2,\ldots,11\},\qquad P=\{12,13\}.
\]

The five missing edges are

\[
\{0,12\},\{1,12\},\{0,13\},\{1,13\},\{12,13\};
\]

thus the active graph is `K4-e` on `U union P`, with ten isolates.  Its
three inclusion-minimal vertex covers are

\[
C_P=\{12,13\},\quad
C_{12}=\{0,1,12\},\quad
C_{13}=\{0,1,13\}.
\tag{1.1}
\]

The exact schedule data are:

| `a=|A4|` | missing-five offset | facets | sum of vertex-link `g4` |
|---:|---:|---:|---:|
| 0 | 8  | 110 | 42 |
| 1 | 14 | 108 | 36 |
| 2 | 20 | 106 | 30 |
| 3 | 26 | 104 | 24 |
| 4 | 32 | 102 | 18 |
| 5 | 38 | 100 | 12 |
| 6 | 44 | 98  | 6  |

Every one of the fourteen vertex links has strict `g4>=1`, by the inherited
Lanner-order bound and the equality case of the generalized lower-bound
theorem.  Hence `a=5,6` are impossible.

Before adding missing four- or five-faces, the vertex-link values are

\[
g_3=(4,4,6^{10},1,1),\qquad
g_4=(5,5,7^{10},1,1).
\tag{1.2}
\]

In the branch `a=0`, every added minimal five-face through vertex 12 or 13
would reduce its baseline `g4=1` to zero.  All eight missing five-faces
therefore lie on

\[
W=U\sqcup I=\{0,\ldots,11\}.
\]

Call their clutter `H`, and let `Z` be the subclutter consisting of the
members of `H` that avoid both vertices of `U`.

## 2. Blocker duality gives two facet families

In a valid pure eight-sphere, a facet complement is an inclusion-minimal
transversal of the complete minimal-nonface clutter and has order five.
In particular, no transversal has order at most four.  Using (1.1), this
implies

\[
\tau(H)\geq3,\qquad \tau(Z)\geq2.
\tag{2.1}
\]

Every five-transversal belongs to exactly one of the following families.
Indeed, it contains one of the three minimal graph covers in (1.1).  If it
contains `C_P`, its other three vertices form a three-set `T` in `W`.  If it
does not contain `C_P`, it contains exactly one of `C_12,C_13`; the other
special vertex is absent, and the two remaining vertices must be isolates.
This exhausts all 310 five-element covers of the missing-edge graph.

1. For a three-set `T` in `W` hitting every member of `H`, there is the
   complement
   \[
   P\cup T.
   \]
   Write its indicator as `y_T`.

2. For a two-set `S` in `I` hitting every member of `Z`, there are the two
   complements
   \[
   \{0,1,12\}\cup S,\qquad \{0,1,13\}\cup S.
   \]
   Write their common indicator as `z_S`.  The indicators for 12 and 13 are
   equal, rather than merely equinumerous: `U` already hits every member of
   `H\Z`, while the same pair `S` must hit `Z`, independently of which one of
   12 or 13 occurs.  Thus either both displayed five-sets are transversals or
   neither is.

The condition `tau(Z)>=2` rules out `y_{01i}` for every isolate `i`: the set
`P union {0,1,i}` hits `H\Z` automatically but hits `Z` exactly when the
singleton `{i}` does, which (2.1) forbids.  This argument is independent of
the later `z`-indicator graph and therefore still applies in its `r=0` and
`r=1` cut cases.  It also prevents overlap between the two displayed
families.  Since no transversal has order at most four, every
five-transversal is inclusion-minimal.  Conversely every facet complement is
such a minimal transversal, so there is no third facet type.  Therefore the
required facet count is exactly

\[
|y|+2|z|=110.
\tag{2.2}
\]

## 3. The complete ridge system

The complement of a ridge is a six-set.  Its degree is the number of selected
five-complements contained in that six-set, and hence must be zero or two.
There are only two nonempty kinds of symbolic ridge equation.

For every four-set `Q` in `W`, the six-set `P union Q` gives

\[
\sum_{T\in\binom Q3}y_T
+2\,\mathbf 1_{Q=\{0,1,i,j\}}z_{ij}\in\{0,2\}.
\tag{R1}
\]

For either `p` in `P` and every isolate triple `i,j,k`, the six-set
`{p,0,1,i,j,k}` gives

\[
z_{ij}+z_{ik}+z_{jk}\in\{0,2\}.
\tag{R2}
\]

The independent replay first classifies all 310 graph-cover five-sets: 220
are `y`-type, including the ten forbidden `y_{01i}`, while 90 are the two
copies of 45 paired `z` variables.  It then constructs all 3,003 six-sets
directly.  Exactly 735
have nonempty symbolic equations: 495 of type (R1), including 45 with the
doubled `z` term, and 240 of type (R2).  It checks 1,980 `y`-facet/six-set
incidences and 810 `z`-facet/six-set incidences.  Thus (R1)--(R2) are the
complete ridge system, not a selected subsystem.

## 4. The `z` graph is a complete bipartite cut

Regard `z` as a graph on the ten isolates.  Equation (R2) says that every
three vertices span zero or two `z`-edges.

Nonadjacency is consequently transitive: if `ij` and `jk` are nonedges but
`ik` is an edge, the triple would span exactly one edge.  Hence nonadjacency
is an equivalence relation and `z` is complete multipartite.  Three distinct
parts would produce a `z`-triangle, also forbidden by (R2).  Therefore

\[
z=K_{r,10-r}
\tag{4.1}
\]

for some `0<=r<=5`, up to exchanging the two sides.  The case `r=0` is the
empty indicator graph; it does not mean that the clutter `Z` is empty.

The replay is exhaustive here.  The nine adjacencies from one fixed isolate
can be chosen in 512 ways.  For every other pair `i,j`, (R2) on the triple
containing the fixed isolate uniquely forces its indicator to be the XOR of
those two seed adjacencies.  All 512 reconstructed graphs satisfy (R2), are
distinct, and are complete bipartite.  Their edge-count histogram is

\[
0:1,\quad9:10,\quad16:45,\quad21:120,\quad24:210,\quad25:126.
\]

## 5. The `y` family is controlled by one graph `X`

Equation (R1) is even modulo two for every `Q`.  Define a graph `X` on the
isolates by

\[
x_{ij}=y_{0ij}.
\]

Apply parity first to `Q={0,1,i,j}` and then to `Q={0,i,j,k}`.  Since
`y_{01i}=0`, one obtains

\[
y_{1ij}=x_{ij},\qquad
y_{ijk}=x_{ij}\mathbin\oplus x_{ik}\mathbin\oplus x_{jk}.
\tag{5.1}
\]

The exact, rather than merely parity, part of (R1) has three consequences.

- If `z_ij=1`, then (R1) for `{0,1,i,j}` forces `x_ij=0`.
- A triangle of `X` makes the ridge equation on `{1,i,j,k}` equal four, so
  `X` is triangle-free.
- An induced `2K2` in `X` makes all four `y`-triples of its four-set equal
  one, again giving ridge degree four.  Thus `X` is induced-`2K2`-free.

Every cross-pair of the cut (4.1) has `z_ij=1`, so the first bullet forces all
cross-edges of `X` to vanish.  If `X` had an edge `ij` inside one cut side and
an edge `kl` inside the other, the four cross-pairs would all be absent;
hence the induced graph on `{i,j,k,l}` would be exactly `2K2`, contradicting
the third bullet.  Consequently all edges of `X` lie in one cut side, on at
most

\[
s=\max(r,10-r)=10-r
\tag{5.2}
\]

vertices.

## 6. The degree bound contradicts all six cuts

Let `e` be the number of edges of `X`, with degrees `d_i`.  Formula (5.1) and
triangle-freeness show that the isolate triples with `y=1` are precisely the
triples containing exactly one edge of `X`.  Therefore

\[
|y|=2e+\sum_{ij\in E(X)}(10-d_i-d_j)
=12e-\sum_i d_i^2.
\tag{6.1}
\]

Only `s` vertices can have nonzero degree.  Cauchy--Schwarz gives

\[
s|y|\leq12es-4e^2
=9s^2-(2e-3s)^2,
\]

and hence

\[
|y|\leq9s.
\tag{6.2}
\]

Combining (2.2), (4.1), and (6.2) gives the complete case table.

| `r` | `s=10-r` | `|z|=r(10-r)` | required `|y|=110-2|z|` | upper bound `9s` | gap |
|---:|---:|---:|---:|---:|---:|
| 0 | 10 | 0  | 110 | 90 | 20 |
| 1 | 9  | 9  | 92  | 81 | 11 |
| 2 | 8  | 16 | 78  | 72 | 6  |
| 3 | 7  | 21 | 68  | 63 | 5  |
| 4 | 6  | 24 | 62  | 54 | 8  |
| 5 | 5  | 25 | 60  | 45 | 15 |

Every gap is strictly positive.  This proves that the `a=0` branch cannot
satisfy the ridge equations and the required facet count simultaneously.

## Trust and reproducibility

The frozen manifest is
`data/enumerations/h9_14_record16_a0_certificate.json`.

- manifest file SHA-256:
  `54e4e81dea471849bbd7af090e5c2aa1a01f681f5fb9bf1d3ae49d082004c58f`;
- canonical internal certificate SHA-256:
  `883e2bff979d1c6681950c4e2027b1c5a5fdb5856b09542af56a166a63975965`;
- proof-module SHA-256:
  `fd2d02c223d683d9e72ad2e6e72edd9d928c6b0a4fe2525dc254b3b5bab8e564`.

Replay with

```text
PYTHONPATH=src:. python3 scripts/certify_h9_14_record16_a0.py --verify-only
PYTHONPATH=src:. pytest -q tests/test_h9_14_record16_a0_certificate.py
```

The archived bounded Z3 timeout is retained only as untrusted corroboration.
No missing-five census, Coxeter labeling, or Gram-matrix feasibility claim is
used in this branch proof.
