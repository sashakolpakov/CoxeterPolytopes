# Exhaustive golden Lanner census: twelve complete compact K3 chambers

**Status (2026-08-09): complete theorem for the stated class.**  Consider
symmetric norm-two Coxeter Gram matrices over `O=Z[phi]`, `phi^2=phi+1`,
whose off-diagonal entries are exactly `0,-1,-phi`.  Require the diagram to
be connected, the identity Gram to have inertia `(k-1,1)`, the conjugate Gram
to be positive definite, and every proper principal subgram to be positive
definite.

There are exactly

\[
  \boxed{4+5+3=12}
\]

unlabelled diagrams: four of rank three, five of rank four, and three of
rank five.  There are none of rank at least six.  Every one of the twelve is
already the complete chamber of its full norm-two `O`-reflection
arrangement, is compact by the literal finite-vertex criterion, and returns
unchanged as a K3 eigen-nef section.

The earlier eight handpicked seeds therefore missed three compact `H3`
chambers and one compact `H4` chamber.  All four missing rows are positive;
none is a hidden-wall failure.

## 1. The theorem

### Theorem 1.1 -- exhaustive equal-root golden Lanner class

Up to simultaneous permutation of rows and columns, every matrix in the
stated class is one of the following:

| hyperbolic dimension | complete list |
|---:|---|
| 2 | `(2,5,5)`, `(3,3,5)`, `(3,5,5)`, `(5,5,5)` |
| 3 | `[3,5,3]`, `[5,3,5]`, `D4(5,3,3)`, `C4(one 5)`, `C4(opposite 5s)` |
| 4 | `[5,3,3,3]`, `[5,3,3,5]`, `D5(terminal 5)` |

Here brackets denote a path, `D4(5,3,3)` is the three-legged star with edge
labels `5,3,3`, the two `C4` rows are four-cycles with respectively one
`5`-edge and two opposite `5`-edges, and `D5(terminal 5)` is the fork

```text
5--3--o
      |\
      3 3
```

with the `5` on the long terminal branch.

Every row has the following additional properties.

1. Its seed simplex is the complete chamber of all norm-two roots in its
   `O`-lattice: no extra mirror cuts the interior.
2. Every simplex vertex is finite; hence the chamber is compact and has
   finite volume.
3. Its even trace lattice embeds primitively in the K3 lattice by Nikulin's
   stable criterion.
4. After the even-golden visible-shell collapse and K3 Weyl framing, the
   projectivized eigen-nef section is the same closed Coxeter simplex.

### Corollary 1.2 -- exact dimensional closure

There is no member of this class in rank `k>=6`, equivalently no such compact
simplex in `H5` or above.

Indeed the hypotheses say precisely that the diagram is connected and
minimal non-elliptic/non-parabolic: it is a Lanner diagram.  Lanner's
classification has maximum order five.  This is a theorem boundary for
compact simplices, not a boundary for general compact Coxeter polytopes with
more than `n+1` facets.

## 2. Exhaustive finite diagram enumeration

For rank `k`, the certifier begins with every one of the
`2^(k choose 2)` support graphs.  A safe exact prefilter forgets the `3/5`
decoration, replacing every nonzero off-diagonal entry by `-1`.  If a
decorated proper deletion were positive definite, the spectral radius of
its nonnegative weighted adjacency matrix would be below two.  Entrywise
decreasing every `phi` weight to one can only decrease that spectral radius.
Thus every surviving support deletion must be simply-laced elliptic.  This
Perron--Frobenius test discards no possible answer.

Every surviving edge is then independently decorated by `3` or `5`.  The
script applies exact identity and conjugate inertia tests and canonicalizes
the full edge-label tuple under all `k!` node permutations.  The literal
search ledger is:

| rank | all support graphs | supports after safe prefilter | decorated labelled Grams tested | labelled survivors | unlabelled orbits |
|---:|---:|---:|---:|---:|---:|
| 3 | 8 | 4 | 20 | 10 | 4 |
| 4 | 64 | 19 | 176 | 54 | 5 |
| 5 | 1024 | 137 | 2384 | 240 | 3 |

For every orbit representative, the certifier separately checks all
`2^k-2` nonempty proper principal subgrams at both real places.  The twelve
named representatives have exactly the twelve canonical permutation keys
found by the unlabeled search, so the naming step neither loses nor adds an
orbit.

Ranks at least six are closed by Lanner's theorem rather than an arbitrary
computational cutoff.  Compact Coxeter simplices are exactly Lanner
diagrams; the modern literature continues to use Lanner's complete list in
this form.  See [Alexandrov's discussion and references](https://arxiv.org/abs/2203.07248).

## 3. Full-root completion

For every surviving Gram `G`, use the exact incenter

\[
  x=-G^{-1}(1,\ldots,1)^t,
  \qquad (x,e_i)=-1,
  \qquad x^2<0.
\]

If `r=sum c_i e_i` is a norm-two root, then `r^perp` meets the closed seed
simplex exactly when the real identity-place coefficients contain both weak
signs.  For a root with both strict signs, split its support into
`I={c_i>0}` and `J={c_i<0}`.  Since every off-diagonal Gram entry is
nonpositive,

\[
 h(c_I,c_I)+h(c_J,c_J)\le h(r,r)=2.
\]

Both supports are proper and elliptic.  Their inverse diagonals bound every
identity coefficient, while positive definiteness of `G^sigma` and
`h^sigma(r,r)=2` bound every conjugate coefficient.  A weakly one-sign
boundary contact has proper support and the same two-place bound.  Finally,
for `c=a+b phi`,

\[
 \operatorname{Tr}(c^2)=2a^2+2ab+3b^2
\]

turns those bounds into exhaustive finite integer boxes.  This is the same
support/sign compiler used for the original eight rows; there is no
heuristic height.

The complete result is:

| diagram | facets | norm-two mirrors meeting closed seed | interior cuts | complete chamber |
|---|---:|---:|---:|---|
| `(2,5,5)` | 3 | 9 | 0 | original seed |
| `(3,3,5)` | 3 | 8 | 0 | original seed |
| `(3,5,5)` | 3 | 10 | 0 | original seed |
| `(5,5,5)` | 3 | 12 | 0 | original seed |
| `[3,5,3]` | 4 | 25 | 0 | original seed |
| `[5,3,5]` | 4 | 27 | 0 | original seed |
| `D4(5,3,3)` | 4 | 26 | 0 | original seed |
| `C4(one 5)` | 4 | 28 | 0 | original seed |
| `C4(opposite 5s)` | 4 | 44 | 0 | original seed |
| `[5,3,3,3]` | 5 | 64 | 0 | original seed |
| `[5,3,3,5]` | 5 | 114 | 0 | original seed |
| `D5(terminal 5)` | 5 | 108 | 0 | original seed |

The contact mirrors beyond the simple walls touch lower-dimensional boundary
faces and do not subdivide the open chamber.  The exhaustive mixed-sign root
list is empty in every row.  Therefore every remaining root can be oriented
with weakly nonnegative seed coefficients, making its inequality a
nonnegative combination of the simple-wall inequalities.  This proves full
root-completeness.

## 4. Literal compactness, not inferred finite volume

The JSON records all homogeneous vertices `v_i=-G^{-1}e_i`, their exact
negative squares, and all wall pairings.  It also records the exact
determinant and two-place inertia of each deleted-wall Gram.  In every row:

- `G` has identity inertia `(n,1,0)` and conjugate inertia `(n+1,0,0)`;
- every rank-`n` deleted-wall diagram is elliptic at both places;
- all `n+1` vertices are finite;
- there are zero ideal and zero hyperideal vertices.

Thus Vinberg's simplex criterion proves compactness directly.  Finite volume
is a consequence, not an inference from a wall prefix.

## 5. The four newly recovered positive rows

The rows absent from the handpicked list have the exact invariants:

| new row | dimension | contacts/cuts | trace determinant | discriminant invariants | period dimension |
|---|---:|---|---:|---|---:|
| `D4(5,3,3)` | 3 | `26/0` | -16 | `(2,2,2,2)` | 12 |
| `C4(one 5)` | 3 | `28/0` | -31 | `(31)` | 12 |
| `C4(opposite 5s)` | 3 | `44/0` | -55 | `(55)` | 12 |
| `D5(terminal 5)` | 4 | `108/0` | -16 | `(4,4)` | 10 |

All four trace lattices are even, have signatures `(1,7)` or `(1,9)`, and
satisfy

\[
 \operatorname{rk}N^\perp\ge\ell(A_N)+2.
\]

Hence their primitive K3 embeddings are automatic in the stable range.
Together with the original eight, the full trace census is:

| diagram | `det N` | discriminant invariants |
|---|---:|---|
| `(2,5,5)` | -16 | `(4,4)` |
| `(3,3,5)` | -20 | `(2,10)` |
| `(3,5,5)` | -36 | `(6,6)` |
| `(5,5,5)` | -100 | `(10,10)` |
| `[3,5,3]` | -11 | `(11)` |
| `[5,3,5]` | -19 | `(19)` |
| `D4(5,3,3)` | -16 | `(2,2,2,2)` |
| `C4(one 5)` | -31 | `(31)` |
| `C4(opposite 5s)` | -55 | `(55)` |
| `[5,3,3,3]` | -4 | `(2,2)` |
| `[5,3,3,5]` | -20 | `(2,10)` |
| `D5(terminal 5)` | -16 | `(4,4)` |

## 6. Exact K3 return

For each field lattice `M`, put

\[
 N(M)=\left(M,
 \operatorname{Tr}_{F/\mathbf Q}
 \left(-\frac{\phi}{\sqrt5}h\right)\right).
\]

The certificate constructs its full integral Gram, Smith invariants, the
self-adjoint integral golden action `tau` with `tau^2=tau+1`, and the paired
K3 facet-root Gram.  Field evenness leaves the two visible K3 root shells
`2` and `4-2phi`; multiplication by the unit `phi` takes the second to norm
two without changing its mirror.  Thus the visible K3 arrangement is exactly
the full norm-two arrangement already completed above.

A generic lattice-polarized K3 period has exact Picard lattice `N(M)`, and a
K3 Weyl change of marking gives

\[
 \boxed{
 \mathbf P\bigl(\operatorname{Nef}(X)\cap E_\phi^+\bigr)
 =\overline P
 }
\]

for the same simplex `P`.  Each facet is supported by the orthogonal K3 root
pair `r_i,(phi-1)r_i`.  The period-family dimensions are 14 in `H2`, 12 in
`H3`, and 10 in `H4`.  The golden action is a divisor-lattice endomorphism;
no surface automorphism is claimed.

## 7. Normal-scaling closure

The frozen 12-row census is now the input to
`golden_normal_scaling_universality.json`.  Every diagram has a connected
`3/5` graph and hence exactly one normal-scaling class modulo oriented vertex
units and common form homothety.  Thus the census is not restricted to an
arbitrarily favored equal-root normalization: the symbolic scaling theorem
returns that normalization exhaustively for this class.

The scaling theorem remains scoped.  It does not assert equal norms across a
dotted nonsquare bridge, and it does not make K3-carrier success universal
for general compact arithmetic polytopes.  Here it closes the normalization
quantifier for precisely the twelve classified diagrams.

## 8. Lineage

The abstract diagrams are classical compact Lanner simplices; no
worldwide-new abstract-polytope claim is made.  The positive project result
is sharper:

- the earlier list of eight was not exhaustive even in its own stated
  equal-root golden class;
- the exhaustive list contains twelve, adding four explicit compact
  `P -> K3 -> same P` returns;
- all twelve original simplices, rather than unspecified subdivisions, are
  exact full-root chambers; and
- Lanner's theorem closes every higher rank in this simplex class.

## 9. Reproducibility

- Census certificate:
  `data/enumerations/golden_lanner_exhaustive_census.json`
- Census certifier:
  `scripts/certify_golden_lanner_exhaustive_census.py`
- Census tests:
  `tests/test_golden_lanner_exhaustive_census.py`
- Updated scaling certificate:
  `data/enumerations/golden_normal_scaling_universality.json`

Replay with

```text
PYTHONPATH=src:. python3 scripts/certify_golden_lanner_exhaustive_census.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_golden_lanner_exhaustive_census.py
PYTHONPATH=src:. python3 scripts/certify_golden_normal_scaling_universality.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_golden_normal_scaling_universality.py
```

References:

- F. Lanner, *On complexes with transitive groups of automorphisms*,
  Meddelanden fran Lunds Universitets Matematiska Seminarium 11 (1950),
  1--71.
- S. Alexandrov, *Lanner diagrams and combinatorial properties of compact
  hyperbolic Coxeter polytopes*,
  [arXiv:2203.07248](https://arxiv.org/abs/2203.07248).
- V. Nikulin, *Integral symmetric bilinear forms and some of their
  applications*,
  [Math. USSR Izv. 14 (1980)](https://doi.org/10.1070/IM1980v014n01ABEH001060).
