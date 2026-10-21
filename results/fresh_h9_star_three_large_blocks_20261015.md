# Three large root blocks force a forbidden spherical triangle

In the saved star incidence, **at most two of the four internal labels can be at least6** in a compact Coxeter realization. No restriction is placed on the eight leaves' own labels or on other labels and dotted distances. This excludes the newly saved119,877-visit necessary assignment, whose four internal labels are all in the unbounded category≥7.

The [exact checker](../scripts/fresh_h9_star_three_large_blocks_20261015.py) and [certificate](fresh_h9_star_three_large_blocks_20261015.json) check all four placements directly against the378 actual vertices and complete minimal-nonface list. They do not use a global incidence automorphism or a previous parameter bank.

| Block | Roots | Cut leaf |
|---|---|---|
|0|3,4|14|
|1|6,7|15|
|2|9,10|16|
|3|12,13|17|

Choose any three blocks with internal labels at least6. Every needed cross-root pair, and every cut-leaf contact with a foreign block, is forced right by an actual three-node face. If such a contact were nonright, its normalized Gram would be entrywise no larger than the path with labels6–3,

\[
\begin{pmatrix}
1&-\sqrt3/2&0\\
-\sqrt3/2&1&-1/2\\
0&-1/2&1
\end{pmatrix},
\]

whose strictly positive null vector is \((\sqrt3,2,1)\). Thus the actual face Gram would have a nonpositive quadratic form, contradicting positive definiteness. The source verifies **96 literal path occurrences** across the four placements. Equality at internal label6 is excluded as well.

The six roots consequently span a positive-definite orthogonal sum of three finite dihedral diagrams. Each cut leaf is right to the two foreign blocks. Its own block together with that leaf is a literal size-three minimal nonface. In a compact Coxeter nerve this is a Lannér triple, so its determinant is strictly negative. With normalized own block cosine \(c<1\) and leaf contact cosines \(a,b\ge0\), the residual square after projecting off the own block is

\[
\rho=\frac{1-c^2-a^2-b^2-2cab}{1-c^2}<0.
\]

The exact determinant expansion is checked. This uses the full Lannér condition, rather than an assumed lower bound on an own angle; it covers right own contacts, arbitrarily large finite own labels, and dense own-contact patterns.

The three projected cut-leaf normals are therefore timelike. The projections onto their distinct owner blocks are mutually orthogonal, so a right ambient pair of leaves would also be a right residual pair. Two orthogonal timelike vectors give two negative directions, impossible in a Lorentzian space. Hence all three leaf pairs are nonright.

But these three cut leaves form an actual spherical face in every placement. Its three Coxeter labels are now all at least3. For the all-ones vector its normalized Gram has quadratic form

\[
3-2\bigl(\cos(\pi/m_1)+\cos(\pi/m_2)+\cos(\pi/m_3)\bigr)\le0,
\]

again contradicting positive definiteness.

This is a dimension-independent Lorentz obstruction under the stated literal face and Lannér premises. It is **not an unconditional dimension bound**: no occurrence theorem is asserted for arbitrary high-dimensional polytopes. On compact star0 it gives four early predicates, each requiring only three internal lower bounds≥6. The search file is unchanged by this proof, and no continuation is performed here.

The [independent review](fresh_h9_star_three_large_blocks_20261015_review.md) reconstructs all96 path guards,12 missing triples and four actual leaf triangles, and independently checks the positive kernel and Schur determinant. Generation and saved replay pass. It uses no producer arithmetic or previous parameter bank.

Generation and saved replay use the source above, with the optional --verify argument.
