# A compact dimension-28 candidate for the cubic-seven arithmetic field

The exact [checker](../scripts/fresh_arithmetic_highdim_20261015.py) and [certificate](fresh_arithmetic_highdim_20261015.json) support the following research proof candidate:

> A compact hyperbolic Coxeter polytope whose Lannér subdiagrams all have order at most three has dimension at most 28. Consequently, a compact arithmetic reflection polytope with ground field \(\mathbb Q(\cos(2\pi/7))\) has dimension at most 28.

The arithmetic conclusion applies more generally whenever the ground field is not a subfield of \(\mathbb Q(\sqrt2,\sqrt5)\), in particular to every nontrivial odd-degree ground field. Compactness of the **whole polytope** is required. No noncompact bound, construction, priority claim, or exhaustive literature-status claim is made. Independent review is being requested; the primary exact generation and saved replay pass.

The practical high-dimensional reduction is that a compact arithmetic chamber of dimension 29 can only have ground field \(\mathbb Q(\sqrt2)\) or \(\mathbb Q(\sqrt5)\), conditional on acceptance of the proof below. This excludes one of the three previously retained high-dimensional fields; it does not exclude compact arithmetic dimension 29 altogether.

## Published inputs and their scope

Vinberg's original paper supplies the finite spherical diagram classification, the compact face/Lannér correspondence, and the following two propositions. Proposition 2, printed pp.70–72, turns angle weights with vertex total at most \(cn\) and a \(k\)-gon total at least \(5-k\) into \(n<8c+6\). Proposition 17, printed p.91, says that every Lannér determinant generates the arithmetic ground field. These are theorem inputs, not new results here. [Vinberg, original paper](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).

For \(n\ge22\), the three compact arithmetic fields are \(\mathbb Q(\sqrt2),\mathbb Q(\sqrt5),\mathbb Q(\cos(2\pi/7))\), as stated in Vinberg's Theorem 2, printed p.69. Nikulin's corrected version explicitly lists these together with \(\mathbb Q\) in the general arithmetic statement. The rational case does not give a compact chamber in this dimension range. [Vinberg, Theorem 2](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt), [Nikulin, PDF p.2](https://arxiv.org/pdf/0708.3991).

The proof below gives a different weighting for the stated subclass. The previous repository [compact29 reduction](fresh_compact_dimension_bound_20261015.md) and its [independent review](fresh_compact_dimension_bound_20261015_review.md) remain unchanged. Their unresolved general compact29 case is not being declared solved.

## The arithmetic restriction on minimal indefinite diagrams

Use normals of squared length two. A connected Lannér diagram on at least four nodes has only labels \(2,3,4,5\). Indeed, a parallel or dotted edge would already violate positive definiteness of a proper two-node subdiagram. If a finite edge had label \(m\ge6\), connectedness would give a third node adjacent nontrivially to one of its endpoints. The proper three-node Gram would be nonpositive at the relaxed path \((6,3)\): its determinant after unit normalization is

\[
1-\cos^2(\pi/6)-\cos^2(\pi/3)=0.
\]

Increasing either edge or adding the third edge cannot restore positive definiteness, by the nonnegative Perron eigenvector of the corresponding adjacency matrix. This contradicts the definition of a Lannér diagram.

Thus every entry of this Gram matrix belongs to \(B=\mathbb Q(\sqrt2,\sqrt5)\): the possible nonzero off-diagonal magnitudes are \(1,\sqrt2,(1+\sqrt5)/2\). Its determinant belongs to \(B\). Multiplying a unit-normalized determinant by the rational factor \(2^{|L|}\) does not change its generated field. Proposition 17 therefore forces the entire arithmetic ground field to be a subfield of \(B\) if such a diagram exists. In particular a cubic field cannot occur. This argument does not need a numerical upper bound on labels or the explicit rank-four/five Lannér list.

For a compact chamber there are no parabolic subdiagrams, so every minimal nonspherical subdiagram is Lannér. Under the resulting hypothesis, every minimal nonface has size two or three. This is weaker than the flag condition that all minimal nonfaces have size two; the flag dimension bounds are not being applied.

## Local face weights

Let \(P\) be a compact Coxeter \(n\)-polytope satisfying the Lannér-order hypothesis. Its finite vertices are simple. A planar angle consists of a vertex \(v\) and a two-face \(F\) through it. In the rank-\(n\) spherical diagram at \(v\), mark the two nodes released to obtain \(F\). Give this angle weight one if their graph distance is at most three, and zero otherwise. Distances ignore edge labels; disconnected nodes have infinite distance.

The following elementary incidence facts specify why the argument applies to actual faces. If a wall subset has positive definite Gram \(G\), orthogonally project a strict interior future point \(H\) to its common wall intersection. Since \(G^{-1}\) is entrywise nonnegative and all wall products are nonpositive, the projection is future timelike and satisfies every remaining wall inequality strictly. Hence the intersection is an actual face. A positive affine kernel similarly gives a nonzero future null vector satisfying all wall inequalities, which is impossible for compact \(P\). Thus a set of facets with empty intersection contains a minimal indefinite diagram, and every proper subset of such a minimal diagram is spherical.

For a triangular two-face, let \(W\) be its \(n-2\) white facet nodes and let \(a,b,c\) be its three bounding black nodes. Each \(W\) together with two black nodes is a finite vertex diagram. The whole set is nonspherical. A minimal indefinite subdiagram must contain all three black nodes; the order hypothesis forces it to be exactly \(\{a,b,c\}\). It is connected, so at least two of the three black pairs are edges. The triangle therefore receives total weight at least two.

For a quadrilateral, write its successive black nodes \(a,b,c,d\). The sets \(W\cup\{a,c\}\) and \(W\cup\{b,d\}\) are nonspherical: a positive definite Gram would give a forbidden actual intersection of opposite sides, while a semidefinite affine block would give an ideal point of the compact chamber. Choose minimal indefinite diagrams \(L,M\) in these two sets. They contain their respective opposite black pairs and at most one white node each.

Two indefinite physical subspaces cannot be orthogonal in a Lorentz space with negative index one. Thus \(L\) and \(M\) either meet at a white node or have a cross edge. In either case, some black node of \(L\) and some black node of \(M\) are joined by a path of length at most three whose internal nodes are white. These two black nodes are adjacent sides of the quadrilateral, and the path lies in their actual spherical corner diagram. Hence at least one quadrilateral angle has weight one. The checker enumerates all 212 connected graph possibilities for this argument, but the preceding reasoning supplies its geometric justification.

For \(k\ge5\), nonnegative angle weights automatically sum to at least \(5-k\). We have therefore established the needed face inequality for every two-face. No assertion about a merely compact facet of a noncompact chamber is used.

## The spherical-forest bound and the dimension inequality

Write \(q_3(S)\) for the number of unordered pairs of nodes at distance at most three in a spherical diagram \(S\). For its connected components the counts are:

| Component | \(q_3\) | \(3r-q_3\) |
|---|---:|---:|
| \(A_1\) | 0 | 3 |
| \(I_2(m)\), every finite \(m\) | 1 | 5 |
| Any path on \(r\ge3\) nodes | \(3r-6\) | 6 |
| \(D_4\) | 6 | 6 |
| \(D_r\), \(r\ge5\) | \(3r-5\) | 5 |
| \(E_6,E_7,E_8\) | \(3r-4\) | 4 |

This includes \(B_r,F_4,H_3,H_4\) in the path row. For a path the counts at distances one, two, three are \(r-1,r-2,r-3\). For \(D_r\), \(r\ge5\), they are \(r-1,r-1,r-3\). Extending the long arm introduces exactly one new pair at each of the three distances, proving these formulas at all ranks. The exceptional counts are \(14,17,20\). The program reconstructs 100 literal graphs as checks of the formulas and boundaries; the all-rank conclusion is analytic.

A connected component of rank at least two has deficit at least four. A forest of total rank at least two consisting only of isolated nodes has deficit at least six. Additivity consequently gives

\[
q_3(S)\le3n-4.
\]

The vertex total of our angle weights is exactly \(q_3(S_v)\). Vinberg's combinatorial Proposition 2, applied with \(c=3-4/n\), yields

\[
n<30-\frac{32}{n},\qquad n^2-30n+32<0.
\]

At \(n=29+u\), \(u\ge0\), the left side is \(u^2+28u+3>0\). Thus \(n\le28\). This proves the subclass bound directly, without importing the general compact ceiling 29. As a separate check at dimension 29, the strict odd-dimensional face count gives a total greater than \(84f_0\), whereas the displayed vertex bound is \(83f_0\); the certificate reproduces the full symmetric h-vector identity.

## Auxiliary complete cubic-field triangle catalogue

This section is not needed for dimension 28. It supplies explicit finite-angle information for subsequent arithmetic refinements.

Let \(t=2\cos(2\pi/7)\), so \(t^3+t^2-2t-1=0\), and put \(s=t^2+t-1=2\cos(\pi/7)\). A finite label \(m\) in a diagram with cyclic-product field \(K=\mathbb Q(t)\) must satisfy \(\mathbb Q(\cos(2\pi/m))\subseteq K\). The real cyclotomic degree formula gives \(\varphi(m)\in\{2,6\}\) for \(m\ge3\). This reduces analytically to \(m\in\{3,4,6,7,9,14,18\}\): every prime divisor has \(p-1\mid6\), and its exponent is bounded by \(p^{e-1}\mid6\). The cubic fields at 9 and 18 differ from \(K\). Their defining cubic has discriminant 81, whereas the cubic for \(K\) has discriminant 49; a common nontrivial field discriminant cannot divide both. Thus the complete finite-label alphabet is \(2,3,4,6,7,14\). [Nikulin, real cyclotomic degree formula (23)](https://arxiv.org/pdf/0708.3991).

For these labels the squared norm-two edge magnitudes are \(0,1,2,3,t+2,t^2+t+1\). The square classes of \(2,3,2+s\) are independent: the norms of all seven nonempty products are nonsquare rational integers. For a triangle without a right edge, its three-cycle belongs to \(K\) exactly when each of the labels 4, 6, 14 occurs an even number of times. For a right triangle the three-cycle is zero. Each proper two-node diagram is positive definite at all three embeddings.

The checker examines all 56 unordered label triples and tests hyperbolicity, the cyclic-product field, and the determinant sign \((-,+,+)\), using only exact fractions and isolating intervals for the three roots of the cubic. Exactly seven triples survive:

| Labels | Determinant of the norm-two Gram |
|---|---|
| \((2,3,7)\) | \(2-2t\) |
| \((2,3,14)\) | \(4-2t-2t^2\) |
| \((2,4,7)\) | \(-2t\) |
| \((2,7,7)\) | \(-4t\) |
| \((2,7,14)\) | \(2-4t-2t^2\) |
| \((3,3,7)\) | \(2-4t-2t^2\) |
| \((7,7,7)\) | \(-2-12t-4t^2\) |

These are the complete necessary cyclic-field/signature possibilities for a Lannér triangle in the stated arithmetic field; no finite-volume chamber is constructed by this list. Each contains a 7 or 14 edge, which must be an isolated rank-two component in any spherical vertex diagram containing it. No finite label cutoff or sampled embedding sign is used.

## Verification and literature status

```sh
python3 scripts/fresh_arithmetic_highdim_20261015.py
python3 scripts/fresh_arithmetic_highdim_20261015.py --verify
```

The certificate records 100 finite graph checks, all 212 minimal quadrilateral graph connections, the complete 56-triple arithmetic decision table, and the scalar and h-vector identities. It binds the previous compact29 primary and independent arithmetic review. Geometric incidence, the all-rank graph induction, and the published arithmetic inputs are explicit proof obligations in this report, rather than being silently identified with a finite computation.

A bounded primary-literature check on 2026-10-15 covered the original Vinberg paper, Nikulin's corrected field paper, and searches directed at Lannér-order restrictions and dimension 28. It did not establish priority or the absence of an earlier equivalent result. The result is presented as a repository proof candidate with auditable dependencies, not as a claimed new published best bound.
