# A coarse dimension bound from a compact facet

**Theorem.** Let P be a finite-sided, finite-volume hyperbolic Coxeter polytope in H^n. If P has a compact facet, then **n≤143**. The facet need not itself be a Coxeter polytope. No narrowness assumption is required for this lemma.

The proof uses a coarser version of Vinberg's graph-distance weighting. It allows parabolic subdiagrams in the ambient face stars. This report establishes the compact-facet branch separately; any combination with the [conditional narrow-facet product-cusp bound](fresh_nikulin_narrow_20261015_products.md) requires review of both proofs and their common hypotheses.

The [small exact checker](../scripts/fresh_nikulin_narrow_20261015_compact_facet.py) and [certificate](fresh_nikulin_narrow_20261015_compact_facet.json) verify the final scalar inequalities, the h-vector coefficient minima, classical distance-count examples, and an explicit counterexample to the stronger localization premise discussed below. The geometric argument and the stated classification inputs are not replaced by a numerical test.

## Published inputs and the localization issue

Prokhorov defines quasi-Lannér diagrams by minimal indefiniteness and records both the maximal order ten and extraction from any hyperbolic diagram in §2, Propositions 1–2, printed pp.407–408. These statements supply the classification input here. [Primary paper](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&paperid=1484&what=fullteng).

Vinberg's §4 Proposition 10, printed p.77, bounds the number of unordered node pairs at graph distance at most c in a spherical diagram of order r by cr. His §1, pp.70–71, provides the symmetric nonnegative h-vector estimate for a simple polytope. [Primary paper](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).

The finer compact-chamber weights cannot simply be imported: Vinberg's §3.2 and Proposition 9 use absence of parabolic subdiagrams, and §6.4 invokes that fact. Compactness of one marked face does not imply this absence. The Coxeter tetrahedron with path labels 3–3–6 has a compact triangular facet at node 0 but an affine G2 subdiagram on nodes 1,2,3. The checker verifies determinant −1/16, dual-ray squares (0,−4,−12,−8), a strict timelike center, and all three spherical corners of facet 0. This diagnoses a gap in the finer transplant; it does not contradict Vinberg's compact theorem or the coarser proof below.

## Minimal indefinite diagrams

All Gram matrices below have diagonal one and nonpositive off-diagonal entries. A minimal indefinite principal subdiagram L is connected: otherwise its indefinite component is a smaller such subdiagram. Every proper principal submatrix is positive semidefinite.

If every proper subdiagram is positive definite, L is Lannér. Otherwise choose a connected affine component A of a proper positive-semidefinite subdiagram. Its null vector has strictly positive coefficients. Connectedness of L supplies an outside node v adjacent to A. The Gram of A+v is already indefinite: evaluating on the positive null vector plus a sufficiently small positive multiple of v gives a negative value. Minimality therefore forces L=A+v. In particular the affine component has codimension one; there can be no additional disconnected affine or spherical part. Its nonzero coupling to the affine null vector makes L nonsingular with inertia (|L|−1,1).

Consequently the general minimal-indefinite situation has exactly the Lannér/quasi-Lannér scope of the cited classification. This also follows directly from Prokhorov's Proposition 2. It gives |L|≤10; the Lannér case has the stronger bound five, which is not needed. A strict dotted edge, if present in a minimal indefinite subdiagram, already constitutes a two-node indefinite diagram and causes no exception to the size bound.

## Compact two-face stars

Let Q be a compact two-dimensional face of P. Let W be the n−2 facet normals whose facets contain Q, and mark the other facets that cut its polygonal sides black. All corners of Q are finite. Therefore W together with each pair of consecutive black nodes is a positive-definite diagram of order n.

The whole star has a negative direction. One way to see this is to project its black normals off the positive subspace spanned by W: their hyperplanes bound the compact polygon Q in a Lorentzian three-dimensional space. The polygon side normals span that space. Every star subdiagram consequently has negative index at most one, as does the ambient Gram matrix.

Two elementary incidence facts will be used. A positive-definite set of wall normals of an acute chamber has an actual common face: project a strict interior timelike point orthogonally off their span. The inverse of a positive-definite Z-matrix is entrywise nonnegative, so all other wall inequalities remain strict and the projected vector stays timelike. If instead a principal block is positive semidefinite and singular, one connected affine component supplies a positive null combination c. This combination is nonzero since its pairing with a strict interior point is negative. It is future null, orthogonal to every normal in that positive-semidefinite block, and nonpositive against every other chamber normal. It is therefore an ideal point of the chamber. These facts do not require the entire chamber to be compact.

For a triangle, any indefinite subdiagram of the star must contain all three black nodes: omitting a black node puts it inside a positive-definite corner diagram. Choose a minimal indefinite subdiagram L. It has at most ten nodes. In a connected graph with three marked nodes, at least two marked nodes can each be deleted while the other two remain joined by a path. For example, use the leaves of the minimal subtree joining the marks in a spanning tree. The resulting two paths lie in L with one black node removed, hence have length at most eight. They lie in two actual corner diagrams of Q.

For a quadrilateral, split the black nodes into its two pairs U,V of opposite sides. The full principal block W+U is indefinite. If it were positive definite, the actual-face projection argument would make the opposite sides intersect in Q. If it were singular positive semidefinite, the null-vector argument would produce an ideal point on every white wall, hence in Q, contradicting compactness. The same reasoning applies to W+V. This is the required replacement for the global no-parabolic assertion.

Choose minimal indefinite subdiagrams L⊆W+U and M⊆W+V. Each must contain its two black nodes, since W plus either single black node is positive definite. Both have at most ten nodes and at most eight white nodes. They cannot be orthogonal: their two negative directions would give ambient negative index at least two. Thus they overlap in a white node or are joined by an edge. A shortest path between the black sets U and V in this connected union has only white internal nodes. It has length at most 8+1+8=17. Its two black endpoints are consecutive sides of Q, so the path is contained in their actual corner diagram.

Assign a planar angle of a compact two-face weight one if its two marked nodes have graph distance at most 17 in the full ambient finite-vertex diagram, and zero otherwise. The triangle argument gives total weight at least two; the quadrilateral argument gives at least one. Polygons with k≥5 need only nonnegative total weight. In all cases the total is at least 5−k.

## Counting on a compact facet

Let F be a compact facet, d=n−1. It is a compact simple d-polytope because all of its vertices are finite vertices of the acute chamber. Its two-faces are also two-faces of P, so the preceding weights apply. At each vertex of F these weights count a subset of the unordered pairs in the full ambient spherical diagram, whose order is n. The published spherical distance bound therefore gives a vertex subtotal at most 17n=17(d+1). No reflection structure on F is being asserted.

Write f_i for the face numbers of F. Summing the two-face lower bounds gives

    total weight ≥ 5 f_2 − binom(d,2) f_0.

For a simple polytope its nonnegative h-vector is symmetric with h_0=h_d=1. Pairing the coefficients binom(k,2) and binom(d−k,2) gives the strict bounds

    f_2/f_0 > d(d−2)/8          if d is even,
    f_2/f_0 > (d−1)^2/8        if d is odd.

Strictness follows from the nonzero endpoint coefficients, which exceed the middle coefficient minimum. It follows that the average vertex subtotal is strictly greater than d(d−6)/8 in even dimension, and (d−1)(d−5)/8 in odd dimension. Comparing with 17(d+1) gives a contradiction whenever d≥143. Explicitly, writing the odd dimensions as d=143+2t and the even ones as d=144+2t, eight times the difference of lower and upper bounds is respectively

    4t²+288t+12 > 0,
    4t²+292t+152 > 0.

Thus d≤142 and n≤143, as claimed.

The proof is deliberately coarse and makes no current-best or novelty claim. Its useful feature is that it permits the affine subdiagrams that obstruct the direct use of the finer compact-chamber weighting. The checker imports no earlier producer arithmetic. Generation and saved replay are run with:

```sh
python scripts/fresh_nikulin_narrow_20261015_compact_facet.py
python scripts/fresh_nikulin_narrow_20261015_compact_facet.py --verify
```
