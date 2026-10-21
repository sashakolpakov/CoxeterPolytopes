# Nikulin as a hostile regression target: conceptual dependency audit

**Date:** 2026-08-10  
**Scope:** arithmetic hyperbolic reflection finiteness, the narrow-face/edge-
graph field funnel, and the relation between field bounds and dimension bounds.  
**Policy:** none of the results below is imported as a dependency of the
project's AG/cohomological-packet route.

## Executive verdict

There is one explicit false implication in the published narrow-edge reduction:
in the `u<2` branch of Proposition 3 of *On ground fields of arithmetic
hyperbolic reflection groups*, a minimal hyperbolic subgraph need not be
Lanner.  The Coxeter orthoscheme `[3,3,6]`, with the middle two nodes marked as
the endpoint walls, is a four-node V-arithmetic finite edge graph with `u=1`,
and is minimal indefinite but quasi-Lanner.  Its affine proper subdiagram is an
exact counterexample to the printed sentence.

The proposition's **field-set conclusion is nevertheless repairable**.  A
minimal indefinite diagram is either Lanner or quasi-Lanner.  In the latter
case it is a finite-volume noncompact V-arithmetic simplex, hence an arithmetic
noncompact simplex, and its ground field is `Q`.  Since `Q` is already in the
small field sets used by Nikulin, later arguments use only the repaired
field-set conclusion and do not acquire a new numerical branch.

The other major conceptual boundary is not an error in Nikulin but is an easy
error for this project to make:

```text
finite field list or bounded [k:Q]  -/->  bounded hyperbolic dimension.
```

Nikulin's final finiteness chain does **not** make that implication.  It uses a
separate dimension cutoff due to Vinberg.  The narrow-face, four-node, and
Fekete arguments are a field funnel, not a global dimension theorem.

The finite five-shape enumeration, all parameter lists, and most decimal
extremum computations in the transition paper remain unverified by this audit.
They are finite and project-certifiable, but cannot yet be treated as native
certificates.

## 1. Frozen primary artifacts

The following local artifacts were read directly.

| Work | Local artifact | SHA-256 |
|---|---|---|
| Nikulin, 1981/82 classification paper | `/private/tmp/nikulin_1981_eng.pdf` | `696652b6a67c6ed420de6c1cb61d9422e3619e28e081a78131a5323ce2edee9d` |
| Nikulin, ground fields I, arXiv v2 | `/private/tmp/nikulin1.pdf` | `bf5e66703c20245ab229b199aa464b3955f3645716d7ae1b967566cadce082d0` |
| Source of ground fields I, arXiv v2 | `/private/tmp/nikulin_ground_fields_I_src/arrefield2v2.tex` | `486855e2702e2ad0bc1e7a8fc828d9bc38ce948ea7dc5ff0e4d98ea0dcdc54b7` |
| Nikulin, 2007 finiteness paper | `/private/tmp/nikulin_finiteness.pdf` | `9a820fef0b88158b283057f7d97acc89130610a9a4aacc6c53a8215e23ce874a` |
| Transition paper, arXiv v4 | `/private/tmp/nikulin_transition.pdf` | `e5391ab018e6e3f92ba9e2af8c35aff79332f973f3021f9c085148f52336f54b` |
| Transition paper, final journal | `/private/tmp/nikulin_transition_journal.pdf` | `716e6e7333efbfaa22060b60dc97e8303363a4914409ccdd56ef7a90c923c7af` |
| Belolipetsky survey | `/private/tmp/belolipetsky_survey.pdf` | `28f6a68799028258e023964c88f5d70cf01df79d45a96ab526c8e9e01c530200` |
| ABSW | `/private/tmp/absw.pdf` | `974ac7c5141ca1c45b515f145215ad115aba867e42551a9b103233ea52930ee4` |

Primary public records:

- [1980 arithmetic/Galois-lattice paper](https://www.mathnet.ru/eng/im1742),
  Math. USSR-Izv. **16** (1981), 573--601;
- [1981/82 classification paper](https://www.mathnet.ru/eng/im1550),
  Math. USSR-Izv. **18** (1982), 99--123;
- [2007 finiteness paper](https://doi.org/10.1070/IM2007v071n01ABEH002349),
  Izv. Math. **71** (2007), 53--56;
- [ground fields I, arXiv:0708.3991v2](https://arxiv.org/abs/0708.3991),
  published in CRM Proc. Lecture Notes **47** (2009), 299--326;
- [final transition paper](https://doi.org/10.1070/IM2011v075n05ABEH002561),
  Izv. Math. **75** (2011), 971--1005.

## 2. Verdict ledger

| ID | Exact locus | Classification | Consequence |
|---|---|---|---|
| C1 | arXiv:0708.3991v2, Proposition 3, printed pp. 8--9, `u<2` branch | **False literal implication** | `[3,3,6]` is minimal indefinite and quasi-Lanner, not Lanner. |
| C2 | Same proposition, asserted field-set conclusion | **Compressed but repairable** | Add the quasi-Lanner branch; it has ground field `Q`, so the stated field union remains valid. |
| C3 | 1982 Theorem 2.1.3(2), p. 104 | **Sound and independently reconstructible** | A hyperbolic principal subgraph carries the full Vinberg ground field. |
| C4 | 1982 Fundamental Theorem 4.1.1, pp. 113--117; 2009 Theorem 5, arXiv pp. 10--13 | **Sound only after the explicit compactness split** | Compact chambers are simple; finite-volume noncompact chambers cannot be put into the same h-vector proof. |
| C5 | 2009 Section 3.3, arXiv pp. 7--8; 2011 Section 3.3, p. 976 | **Unverified finite enumeration** | “Easy to draw” is not a completeness certificate for the five four-node shapes or their label tuples. |
| C6 | 2011 Theorems 4.6--4.9, pp. 989--998 | **Unverified computation** | Most interval extrema and parameter branches still need exact symbolic replay. |
| C7 | 2009 Lemma 11, arXiv p. 14 | **False literal notation, harmless repair** | `P(M)={delta_1,...,delta_4}` must mean the wall set of the quadrilateral `M_4`, not that of the ambient chamber. |
| C8 | 1982 Theorem 4.1.1 versus Theorem 2.3.1 | **Threshold mismatch, harmless repair** | The theorem says degree `>=N(14)` while the invoked edge theorem assumes degree `>N(14)`; enlarge or replace the threshold. |
| C9 | 2011 Remark 5.1, p. 1000 | **Acknowledged historical error** | The 2009 claim `N(14)<=56` had only established `<=120`; later `<=25` supersedes it. |
| C10 | 2011 Theorem 3.2, p. 977 | **Unsupported equality wording** | The proof gives `N(14)<=25`; it does not exhibit a degree-25 field proving equality. |
| C11 | 2010 arXiv v4 equation (24), p. 15 | **False numerical formula, journal-repaired** | Exponent `2` gives the wrong norm; journal (4.19), p. 986, uses exponent `3`. |
| C12 | 2007 Theorem 2.2, p. 56 | **Logically separated and sound in form** | Global finiteness explicitly appends Vinberg's independent `n<30`; no field-to-dimension implication is asserted. |

## 3. Exact counterexample to the Lanner sentence

### 3.1 The printed implication

Proposition 3 of arXiv:0708.3991v2, pp. 8--9, considers the unique
hyperbolic connected component of a finite edge graph with distinguished
endpoint walls `e_1,e_2` and `u=e_1.e_2`.  In the branch `u<2` it says:

```text
all one- and two-vertex subgraphs are negative definite
    => a minimal hyperbolic subgraph L is Lanner.
```

The implication omits the possibility that a proper subgraph is parabolic.

### 3.2 The four-node matrix

Take the Coxeter chain `[3,3,6]`, in Nikulin's square-minus-two
normalization:

```text
        1          1        sqrt(3)
  o -------- o -------- o ------------ o
  1          2          3              4
```

Its Gram matrix is

```text
G = [ -2       1       0        0      ]
    [  1      -2       1        0      ]
    [  0       1      -2      sqrt(3)  ]
    [  0       0     sqrt(3)   -2      ].
```

The leading principal determinants are

```text
-2,  3,  -4,  -1.
```

Hence the exact `LDL^T` pivots are

```text
-2,  -3/2,  -4/3,  1/4,
```

so `G` has signature `(1,3)` and determinant `-1`.

The vertex-deletion data are:

| deleted node | determinant | spectrum/type |
|---:|---:|---|
| 1 | `0` | `{-4,-2,0}`; affine chain `[3,6]` |
| 2 | `-2` | `{-2,-2-sqrt(3),-2+sqrt(3)}`; negative definite |
| 3 | `-6` | `{-3,-2,-1}`; negative definite |
| 4 | `-4` | `{-2,-2-sqrt(2),-2+sqrt(2)}`; negative definite |

All two-node subgraphs are negative definite because their off-diagonal
entries are `<2`.  The whole graph is therefore minimal indefinite, but the
deletion of node 1 is parabolic.  It is quasi-Lanner, not Lanner.

Mark nodes 2 and 3 as the endpoint walls.  Deleting node 2 or node 3 is
negative definite, so the corresponding two vertices of the simplex are
finite and the simplex is a valid finite edge chamber around their connecting
edge.  But

```text
u = e_2.e_3 = 1 < 2.
```

The affine deletion shows that the simplex has one ideal vertex and finite
volume.  Its ground field is `Q`: the graph is a tree, so the nontrivial cyclic
products are the edge squares `1,1,3`.  Its coefficient extension is
`Q(sqrt(3))`.  Thus it is also a V-arithmetic example, not merely an abstract
signature matrix.

This is an exact counterexample to the sentence in the proof, within the
proof's own finite-edge hypotheses.

The independent machine-readable replay is frozen in:

- `scripts/certify_nikulin_336_quasilanner_counterexample.py`;
- `data/enumerations/nikulin_336_quasilanner_counterexample.json`;
- `tests/test_nikulin_336_quasilanner_counterexample.py`.

It certifies the full inertia, all deletion inertias, the marked value `u=1`,
the cyclic-product ground field, and the primitive affine null vector
`(1,2,sqrt(3))` using exact arithmetic only.

### 3.3 Exact repair

Choose a minimal indefinite principal subgraph `L`.  Since all proper
principal subgraphs of `L` are negative semidefinite, precisely two cases are
possible.

1. Every proper subgraph is negative definite.  Then `L` is Lanner and the
   compact Lanner/triangle field lists apply.
2. Some proper subgraph is parabolic.  Then `L` is quasi-Lanner and defines a
   finite-volume noncompact Coxeter simplex.

By 1982 Theorem 2.1.3(2), `L` is V-arithmetic with the same ground field as the
original hyperbolic edge graph.  By Theorem 2.1.3(1), finite volume upgrades it
to an arithmetic reflection simplex.  A noncompact arithmetic quadratic-form
reflection lattice has ground field `Q`: a cusp supplies a nonzero isotropic
vector over the ground field, while every nongeometric real conjugate of the
form is definite; a nonzero isotropic vector cannot survive such a conjugate.

Therefore the quasi-Lanner branch contributes only `Q`, already contained in
both `FL^4` and `FT`.  Proposition 3's stated **field union** survives, although
its printed proof does not.

### 3.4 Downstream dependency audit

No later dependency located here needs the false assertion that the `u<2`
diagram is compact Lanner.

- Lemma 6 of arXiv:0708.3991v2 uses only the contrapositive field statement:
  if the ground field is outside `FL^4 union FT`, then `u>2`.
- Theorems 5 and 9 use Lemma 6 only to enter the marked four-node `u>2`
  families.
- Journal Proposition 3.1 (2011, p. 976) cites the 2009 proof, while Theorems
  3.2 and 5.2 use only the partition into low-degree `FL^4,FT` fields and the
  five `u>2` families.
- The transition bound is unchanged because the repaired branch has degree
  one.

Thus C1 is a real local proof error, but not a counterexample to the later
field-degree bounds once C2 is inserted explicitly.  Any census or AG carrier
that records actual diagram type, rather than only the ground field, **must**
retain the affine/quasi-Lanner branch.

## 4. Ground-field inheritance is a sound small-probe principle

Nikulin's 1982 Theorem 2.1.3(2), p. 104, says that a hyperbolic principal
subgraph of a V-arithmetic graph is V-arithmetic with the same ground field.
This can be reconstructed without treating it as an oracle.

Let `K_0` be the field generated by cyclic products of the hyperbolic
subgraph and `K` the ground field of the whole graph.  If `K_0` were a proper
subfield of `K`, total reality and separability give a nonidentity embedding
of `K` that fixes the geometric embedding of `K_0`.  Extend it to the totally
real coefficient field.  On the connected subgraph, all two-cycles fix the
absolute values of Gram entries and all longer cycles make the resulting sign
cocycle a vertex-sign coboundary.  The conjugate submatrix is therefore
diagonally congruent to the original hyperbolic submatrix and has the same
signature.  But V-arithmeticity requires the full conjugate Gram form to be
negative definite.  A negative-definite matrix cannot have a hyperbolic
principal submatrix, contradiction.

This is one of the sound insights to translate into the project:

```text
small hyperbolic Gram packet + full cyclic-product algebra
    => the packet remembers the global ground field.
```

The qualifier “full cyclic-product algebra” is essential.  A surface or Belyi
carrier that remembers only a few labels, but permits the remaining field in a
moving period or cross-ratio, does not inherit this theorem.

## 5. Compactness and simplicity: where the topological argument stops

### 5.1 The published compact branches make the split

The 1982 proof of Fundamental Theorem 4.1.1, stage 4, p. 114, says that the
degree is larger than one, hence the arithmetic chamber is compact, and only
then invokes simplicial corners.  The 2009 proof of Theorem 5 likewise begins:

```text
noncompact chamber => ground field Q;
otherwise assume compact => the chamber and narrow facet are simple.
```

Within that scope, the h-vector/face-average argument is legitimate.  It gives
for the narrow `(n-1)`-facet the strict average bounds

```text
a < 4+4/(n-2)   for even n,
a < 4+4/(n-3)   for odd n.
```

Together with the separately classified small Coxeter components, this forces
a small local field diagram in high compact dimension.

### 5.2 The noncompact export is invalid

A finite-volume noncompact Coxeter chamber need not be simple at an ideal
vertex.  Its horospherical link may be a reducible affine Coxeter chamber, a
product of Euclidean simplices.  In dimension three, an ideal right-angled
octahedral vertex has four incident walls, not three.  Thus the compact
simple-polytope h-vector and incidence identities cannot be applied blindly.

Strict terminology also matters: for a finite-volume hyperbolic polytope an
entire positive-dimensional face does not live at infinity.  “Ideal face” in
this context means an ideal vertex; ordinary positive-dimensional faces may be
noncompact and have ideal vertices.

The implication “noncompact arithmetic => `K=Q`” is only a field statement.
It does not restore simplicity, bound the number of affine components, or give
a dimension theorem.

## 6. Finite graph and numerical parts not yet certified

### 6.1 Five four-node shapes

Ground fields I, Section 3.3, arXiv pp. 7--8, and the transition paper,
Section 3.3, p. 976, say that it is “easy to draw” all connected four-node
finite edge graphs with `u>2`, obtaining five shapes, and that the integer
parameters are easily listed from the two spherical endpoint deletions.

That is a finite claim, but the article does not provide a completeness proof
that meets the project's standard.  A native certificate must generate all
rank-three spherical endpoint-deletion diagrams, glue them along the common
two-node subdiagram, quotient by marked graph isomorphism, and recover exactly
the five shapes and every allowed label tuple.  Until that is done, the shape
and tuple lists are **unverified computation**, not false.

### 6.2 Interval extrema

Theorems 4.6--4.9 of the 2011 paper, pp. 989--998, contain the finite
parameter branches and numerous decimal interval maxima.  The exponent-two
polynomial printed in arXiv v4 equation (24), p. 15, is exactly wrong; the
final journal equation (4.19), p. 986, repairs it with exponent three.  The
project's exact audit is frozen separately in
`results/nikulin_transition_polynomial_hostile_audit.md`.

The final journal branch checked there is sound.  The remaining branches still
need exact determinant derivations, algebraic interval isolation, and Sturm or
resultant certificates.  Their publication is evidence, not a native proof.

### 6.3 Correction history

The correction trail is material:

- The 1982 Section 1 Fekete/hyperbolic-number formulas contained arithmetic
  inaccuracies; the qualitative existence argument survives corrected
  constants.
- The 2009 paper claimed `N(14)<=56` but, as the author records in 2011 Remark
  5.1, p. 1000, its calculation established only `N(14)<=120`.
- The same incorrect condition entered the 2008/09 low-dimensional papers;
  the author says their much looser terminal bounds `138` and `909` were not
  affected.  This audit has not replayed that assertion.
- The 2011 journal proof gives the current published upper bound `N(14)<=25`.
  The wording `N(14)=25 (or <=25)` in Theorem 3.2 should not be read as an
  equality theorem: no degree-25 witness is supplied.

## 7. Field finiteness is not dimension finiteness

The exact logical directions are:

```text
n >= 10
  => a compact narrow facet has a forced small local Gram packet
  => the ground field lies in finite small-diagram families
  => [K:Q] is bounded.
```

This is not reversible.  In particular,

```text
[K:Q] bounded
  -/-> n bounded,

finite list of K
  -/-> n bounded.
```

The 1982 theorem is a high-dimension-to-field statement (equivalently, very
large field degree forces low dimension).  It is not an absolute dimension
cutoff.  The 2007 all-dimensional finiteness paper keeps the inputs separate:

1. the 1980 theorem gives finiteness for fixed dimension and fixed field-degree
   data (its core fixed-rank/fixed-extension statement is Theorem 3.1.1,
   p. 582 of the English translation);
2. the 1981 narrow-face method and the 2007 induction bound field degree in
   each remaining dimension;
3. Vinberg's separate theorem gives `n<30`;
4. only then is the union over dimensions finite.

The 2007 induction itself is short but still depends on the old edge theorem's
exceptional-diagram classification.  In Theorem 2.1, pp. 54--56, the facet and
codimension-two branches inherit the same field by Theorem 2.1.3(2); the final
branch chooses an edge containing a connected negative-definite triple and
invokes Theorem 2.3.1 to exclude its three exceptional diagrams.  The
inductive architecture is coherent, but that finite exceptional-diagram step
has not been independently replayed here.

Maximality enters the final classification theorem, not the local edge-field
bound.  Nikulin supplies no implication of the form

```text
reflection maximality => bounded labels, rank, Picard number, or dimension.
```

An AG translation of maximality may still be valuable, but a boundedness
theorem must be proved for the translated condition; maximality cannot be used
as a name for the missing bound.

## 8. AG/cohomological-packet translation

The useful pieces of the regression target translate as follows.

| Coxeter mechanism | Project-native AG/packet analogue |
|---|---|
| Narrow facet of minimality 14 | A geometric construction giving an absolute conjugate/height envelope for marked integral coordinates |
| Hyperbolic subgraph inherits the full field | A functorial bounded integral action whose full trace algebra is already visible on a small cohomological summand |
| Fekete/norm argument | An integral marked coordinate plus strict bounds at all nongeometric embeddings gives an effective degree bound |
| h-vector face average | A global geometric/combinatorial theorem forcing one bounded local packet when carrier rank/dimension is high |
| Four-node families | Finite marked arithmetic dessins, provided shape and label completeness is independently certified |

The `[3,3,6]` counterexample adds a mandatory cusp clause.  A compact Lanner
diagram corresponds, in the surface analogy, to a negative-definite curve
configuration that can be contracted.  A quasi-Lanner diagram contains an
affine null direction.  It behaves instead like a fibre/cusp or log-canonical
boundary cluster and is **not** interchangeable with a compact
negative-definite exceptional configuration.  A carrier grammar must therefore
record at least:

```text
elliptic / parabolic / hyperbolic type,
the primitive affine null vector,
and the attachment of the affine cluster to the marked hyperbolic packet.
```

Ground field `Q` at a cusp does not permit this data to be discarded.  It
removes a field extension, not the affine geometry.

Finally, a Belyi surface alone does not recreate the Fekete mechanism.  The
carrier must prove integrality, one global finite return grammar, and a strict
conjugate envelope.  Nor does a finite dessin field list bound dimension.  A
dimension theorem still requires a uniform cohomological-rank bound, an
intrinsic full trace-module injection, or a separate global theorem forcing
and excluding local packets.

## 9. Belolipetsky attribution

No public primary source was found in which M. Belolipetsky claims a deeper
failure of Nikulin's narrow-face or transition proof.  The public record located
here instead uses the results:

- [Belolipetsky's survey](https://arxiv.org/abs/1506.03111), Section 3,
  states the minimality-14 theorem and the published `N(14)<=25` bound;
- [Belolipetsky--Linowitz](https://arxiv.org/abs/1211.4443), Corollary 5.1,
  combines their three-dimensional result with Nikulin's higher-dimensional
  result;
- [ABSW](https://arxiv.org/abs/math/0612132), Introduction, describes its
  method as independent of Nikulin's contemporaneous proof, not as a
  refutation.

Therefore the explicit `[3,3,6]` defect should be reported on its own merits.
No private criticism should be attributed to Belolipetsky without the exact
email, comment, referee report, or link.

## 10. Final classification

### False

- The `u<2 => minimal hyperbolic subgraph is Lanner` sentence in 2009
  Proposition 3.
- The literal ambient-wall assertion `P(M)={delta_1,...,delta_4}` in 2009
  Lemma 11.
- The obsolete numerical claims `N(14)<=56` and arXiv-v4 equation (24), both
  acknowledged or repaired in the 2011 journal paper.

### Compressed but repairable

- Replace Lanner by “Lanner or quasi-Lanner”; the quasi-Lanner branch has
  ground field `Q`.
- Make compactness explicit before every simple-polytope incidence argument.
- Replace equality language for `N(14)` by the proved upper bound.
- Fix the strict/non-strict threshold mismatch by enlarging `N(14)`.

### Unverified computation/classification

- Completeness of the five marked four-node shapes and every label tuple.
- Exact replay of all determinant intervals and maxima in Theorems 4.6--4.9.
- The old Theorem 2.3.1 exceptional-diagram branch used by the 2007 induction.
- The claim that the incorrect low-dimensional Fekete condition did not affect
  the terminal `138` and `909` bounds.

### Sound insight worth translating

- The narrow-face conjugate envelope.
- Full-field inheritance by a hyperbolic principal Gram packet.
- Integral norm/Fekete conversion of an interval packet into an effective field
  bound.
- Compact h-vector averaging as a way to force a bounded local packet.
- The strict separation between a field funnel and a global dimension theorem.
- The necessity of carrying affine cusp/null-root data separately from compact
  negative-definite AG configurations.

**Bottom line:** Nikulin supplies a valuable regression architecture, not the
project's terminal theorem.  The concrete Lanner/quasi-Lanner error is real and
locally repairable.  The dimension problem remains genuinely ours.
