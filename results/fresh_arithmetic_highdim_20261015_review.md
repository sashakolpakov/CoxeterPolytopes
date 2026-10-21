# Independent review of the compact arithmetic dimension-28 candidate

**PASS: independent exact generation, saved replay and source/proof audit.** The [standalone checker](../scripts/fresh_arithmetic_highdim_20261015_review.py) and [certificate](fresh_arithmetic_highdim_20261015_review.json) bind the frozen [primary source](../scripts/fresh_arithmetic_highdim_20261015.py) and [primary certificate](fresh_arithmetic_highdim_20261015.json). No producer arithmetic or combinatorial code is imported. The [primary report](fresh_arithmetic_highdim_20261015.md) is read but not hash-bound, so its historical pending-review sentence is not a statement about this completed review.

The reviewed conclusion is: a **compact** hyperbolic Coxeter polytope whose Lannér subdiagrams all have order at most three has dimension at most 28. The arithmetic corollary applies when its ground field is not contained in \(\mathbb Q(\sqrt2,\sqrt5)\), in particular to the cubic-seven field. This is not an unrestricted compact dimension-28 theorem and makes no claim about noncompact chambers or priority.

## Independent calculations

The checker constructs 160 connected underlying graphs through rank 80, using a different vertex numbering for the forked diagrams. Breadth-first distances agree with the tree identity

\[
q_3(T)=|E(T)|+\sum_v\binom{\deg v}{2}
 +\sum_{uv\in E(T)}(\deg u-1)(\deg v-1).
\]

The three terms count paths of lengths one, two and three. This gives the all-rank formulas \(3r-6\) for paths with \(r\ge3\), \(3r-5\) for \(D_r\) with \(r\ge5\), and \(3r-4\) for \(E_6,E_7,E_8\). The boundary cases are checked separately. Every one of the producer's 100 graph records is also reconstructed and compared. Arbitrary finite \(I_2(m)\) is covered by its two-node graph; this is not a label cutoff.

For a forest of rank at least two, either a component of rank at least two contributes deficit at least four, or there are at least two isolated nodes and the total deficit is at least six. Thus \(q_3\le3n-4\). An additional integer partition dynamic program gives the sharper value 82 at rank 29, but the proof needs only the weaker bound 83.

The local graph replay independently enumerates the four connected black triples and all 212 minimal quadrilateral connection configurations. Every triangle has at least two direct black edges. Each quadrilateral configuration has a cross-pair path of length at most three with only white internal vertices. Extra edges can only shorten the path.

The rank-29 symmetric h-vector coefficients are reconstructed independently: their minimum is 196 and the endpoint coefficient is larger. Consequently the total required face weight is strictly greater than \(84f_0\), contradicting the vertex bound \(83f_0\). The script also verifies all higher dimensions directly. At odd \(n=29+2t\) and even \(n=30+2t\), the strict h-vector lower bound minus \(3n-4\) is respectively

\[
1+7t+\tfrac12t^2,
\qquad 4+\tfrac{15}{2}t+\tfrac12t^2,
\]

both positive for \(t\ge0\). Hence this review does not rely on first importing the general compact upper bound 29.

## Geometric hypotheses checked

The global compactness assumption is essential. For an acute chamber, a positive-definite wall subset defines an actual face: projecting a strict interior point onto its common orthogonal complement remains future timelike and strict at every other wall, because the inverse of its positive-definite Gram Z-matrix is nonnegative. A positive affine dependence instead supplies a feasible future null point. Compactness excludes the latter, so every minimal nonspherical subset is Lannér.

For a triangular two-face, adjoining any two black sides to its white support gives an actual finite vertex. A minimal indefinite subset of its full star must therefore contain all three black nodes; the order hypothesis leaves exactly their connected triple.

For a quadrilateral, adjoining either opposite black pair to the white support is indefinite. Positive definiteness would give an actual intersection of opposite sides, while semidefiniteness would supply a forbidden ideal point. A minimal indefinite subset contains its two opposite black nodes and at most one white node. Two such subsets cannot span orthogonal indefinite subspaces in a Lorentz space of negative index one. Their shared white node or a cross edge therefore yields the required path of length at most three. Its black endpoints are adjacent sides of the quadrilateral, so the path belongs to an actual spherical vertex diagram. The graph census checks this topology but does not replace these geometric premises.

Thus triangles receive at least two units, quadrilaterals at least one, and larger polygons receive at least \(5-k\) by nonnegativity. In particular, the earlier invalid localization from a merely compact facet of a noncompact chamber is not used here.

## Arithmetic scope

I checked Proposition 17 directly on printed p.91 of [Vinberg's original paper](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt). It states that a Lannér determinant generates the arithmetic ground field and becomes positive under its other embeddings. The preceding arithmeticity criterion and its short proof justify the claimed general field restriction, not only the cubic example.

Independently, a connected Lannér diagram on at least four nodes cannot contain an edge of label at least six: some adjacent nonright edge produces a proper triple whose determinant is already zero at the relaxed path \((6,3)\), and becomes nonpositive when any edge increases. The checker expands the full three-edge determinant and verifies this monotonic obstruction exactly. Parallel or dotted edges already violate a proper two-node positive-definite minor. Hence all labels are at most five and the determinant belongs to \(B=\mathbb Q(\sqrt2,\sqrt5)\).

The determinant also lies in the cyclic-product ground field \(k\), since its permutation expansion consists of products of Gram cycles. If \(k\ne\mathbb Q(\det L)\), a nonidentity real embedding of \(k\) fixing the negative determinant extends to the totally real Gram field, contradicting conjugate positive semidefiniteness. This reconstructs the consequence of Proposition 17 used here: such a diagram forces \(k\subseteq B\). For the cubic-seven case the separate irreducible polynomial \(x^3-x^2-2x+1\), for the alternative generator \(2\cos(\pi/7)\), has discriminant 49 and three real roots. Its degree cannot divide four.

The statement that the only remaining compact arithmetic dimension-29 ground fields are \(\mathbb Q(\sqrt2)\) and \(\mathbb Q(\sqrt5)\) additionally imports the high-dimensional field theorem cited by the primary report. It is not a new field-classification census in this review.

## Verification limits

The optional seven-triangle arithmetic catalogue was read for consistency but was **not independently re-enumerated** here. It is not used by the dimension argument. The checker's finite computations are accompanied above by the all-rank formulas and face proofs; no finite graph cutoff is promoted to an all-rank theorem without those arguments. Literature priority and external peer review remain outside scope.

```sh
python3 scripts/fresh_arithmetic_highdim_20261015_review.py
python3 scripts/fresh_arithmetic_highdim_20261015_review.py --verify
python3 scripts/fresh_arithmetic_highdim_20261015.py --verify
```

All three commands pass. Source hashes, syntax and local report links also pass. The review files are final and frozen.
