# Coupled cusp restrictions for a stated three-type ideal RA5 facet family

18 October 2026. **New necessary compatibility result; full independent internal review completed.** Assume an entirely ideal finite-volume right-angled hyperbolic five-polytope has every four-facet isometric to one of: the regular 24-cell C, its reflected facet double D, or the specific sixteen-cell positive carrier T constructed in the frozen [doubling packet](ideal_ra4_positive_kappa_doubling_20261018.md). Then it must have at least 30 T facets and 66 D facets. A particular cusp of every T forces an eight-T ambient cusp and one of two possible widths is excluded. These statements are conditional on this explicit facet spectrum. They do not exclude the spectrum, construct an ideal RA5 polytope, or restrict all unlisted four-facet types.

## Complete input cusp stars

The [producer](../scripts/fresh_ideal_ra5_known_carrier_cusp_compatibility_20261018.py) reconstructs all24+42+258 cusp stars directly from the frozen exact wall/vertex sets. Each has three opposite pairs of unit normals. If the normalized ideal vertex is v and an opposite pair has sum $\lambda v$, its horospherical side length is proportional to $\lambda$. This follows in the upper-half-space model with v at infinity: unit normals have horizontal Euclidean norm one, and their sum coefficient is the separation of their parallel boundary planes, up to the common horosphere scale.

The three normalized cusp shapes are denoted $C_0=(1,1,1)$, $A=(1,1,2)$ and $B=(1,2,2)$. Their exact multiplicities are:

|Four-facet type|Number of $C_0$ cusps|Number of A cusps|Number of B cusps|Kappa|
|---|---:|---:|---:|---:|
|C|24|0|0|0|
|D|36|6|0|−2|
|T|195|49|14|6|

For each of the324 cusps, the full saved record retains the three widths and the number of vertices in each of its six incident three-facets, with the actual facet indices. These are complete stars of the three specified polytopes, not a catalogue of all ideal RA4 polytopes. In particular T has128 octahedral,48 nine-vertex,18 fourteen-vertex and4 twenty-four-vertex three-facets.

## Global count and orthogonality constraints

Write $x,y,z$ for the numbers of C,D,T facets. The frozen universal ideal5 carrier identity gives

\[
 -2y+6z=48,\qquad y=3z-24. \tag{1}
\]

Every ambient cusp is a rectangular four-box. Since all its three-box facets have one of the three listed shapes, its side lengths, up to scale and permutation, must be $1111,1112,1122$ or $1222$. Indeed every pair of lengths belongs to a triple, so their ratio is at most2; three distinct lengths would form a forbidden triple. Let $a,b,c,d$ be the numbers of these four ambient shapes. The eight box facets have respectively the following three-box type counts:

\[
\begin{array}{c|ccc}
 &C_0&A&B\\\hline
1111&8&0&0\\
1112&2&6&0\\
1122&0&4&4\\
1222&2&0&6.
\end{array}
\]

Counting actual cusp occurrences gives

\[
8a+2b+2d=24x+36y+195z,\quad
6b+4c=6y+49z,\quad 4c+6d=14z. \tag{2}
\]

Subtracting the last two equations gives $6(b-d)=6y+35z$, so $6\mid z$. This is a count of actual vertices of facets; native orientations are not additional occurrences.

The22 large three-facets of each T must meet22 distinct T neighbors: neither C nor D has a three-facet with14 or24 vertices. Distinctness follows because the intersection of two convex facets is a single face. Form the simple graph on the z T facets, joining any two sharing a ridge. Its minimum degree is at least22, so it has at least11z edges. It has no clique on six vertices: six mutually adjacent right-angled facets would give six mutually orthogonal unit spacelike normals in Lorentz space of positive index five.

For completeness, the elementary weighted clique bound gives $e\le2z^2/5$ for any graph without a six-clique. Maximize $\sum_{ij\in E}t_it_j$ over nonnegative weights summing to one, and choose a maximizer with smallest support. If two supported vertices are nonadjacent, concentrating their combined weight on the one with the larger weighted neighbor sum cannot reduce the objective, contradicting minimal support. Thus the support is a clique of size at most five. Its objective is at most $(1-1/5)/2=2/5$. Applying this to uniform weights proves the edge bound. Consequently $z\ge55/2$, and with $6\mid z$:

\[
 \boxed{z\ge30,\qquad y\ge66.} \tag{3}
\]

This uses only necessary orthogonality, without assuming any abstract graph has a Lorentz realization.

## A forced local eight-T cusp

T has exactly one cusp, saved vertex85, whose incident three-facet sizes are $(24,24)$, $(24,24)$ and $(14,14)$ on three opposite pairs. Their widths are exactly $(4,4,8)$ in the saved normalization. Thus, if T is a facet of the hypothetical P, its six neighboring facets at this cusp must all be T. The ambient four-box has widths $(t,4,4,8)$ for some positive t; the original T is perpendicular to the t direction.

Consider a side facet perpendicular to one of the width4 directions. Its three-box has widths $(t,4,8)$. To have one of the allowed shapes, t must be4 or8. The opposite of the original T has shape A and hence is D or T; it cannot be C. This proves completeness of the four tests below.

For each box facet, the exact domain consists of all its actual cusp records, with all permutations and sign flips of coordinate axes that preserve the three widths up to a common scale. On any two adjacent box facets, the recorded vertex count of their common three-face must agree. Allowing all such cusp symmetries is a necessary overdomain: no assertion is made that every local symmetry extends to a global symmetry of the four-polytope.

|Fourth width t|Opposite four-facet|Necessary decorated-box matching|
|---:|---|---|
|4|D|Excluded by complete finite search|
|4|T|An explicit matching survives|
|8|D|Excluded by complete finite search|
|8|T|Excluded by complete finite search|

The t=8 obstruction can also be seen directly in the saved cusp table: at every B cusp of T, the four sides corresponding to the two long coordinate directions are nine-vertex three-facets. A required neighbor perpendicular to an old width4 direction has shape B, but its side facing the original T lies in the long t direction and is a twenty-four-vertex three-face. These labels cannot agree. The t=4 opposite-D case requires the coupled adjacent-facet checks; it is not inferred just from one three-box shape.

The surviving t=4,T assignment is only a consistent width/ridge-size labeling. It does not give a complete boundary complex, isometries of the full three-faces, mutually consistent facet placements, a Lorentz Gram realization, or a chamber. The exclusions are valid because all these additional conditions can only shrink the necessary domain.

## What these constraints leave open

The arithmetic data retain the nonnegative count assignment

\[
 (x,y,z)=(361,66,30),\qquad(a,b,c,d)=(2051,241,105,0).
\]

It satisfies (1), (2), the elementary facet lower bound $f_4\ge1+198+258=457$ imposed by any T facet, and the ideal5 Euler/facet-incidence identities, with457 facets and2397 vertices. This is not a geometric example and does not prove simultaneous feasibility of all local decorated-box assignments. It shows why the currently checked aggregate equations alone do not exclude the three-type spectrum. Arbitrary unlisted facet types remain entirely outside this conditional theorem.

The [standalone checker](../scripts/check_ideal_ra5_carrier_cusps_saved_20261018.py) imports no producer. It reconstructs the 324 stars from frozen coordinates, checks every saved domain, independently decides the four matching problems by arc-consistency plus complete branching, and verifies the count equations and all no-claim flags. Full original four-polytope completeness and the universal carrier identity remain inherited. The universal orthogonality/clique and actual-incidence arguments are written proofs, not inferred from the finite searches. Generation and separate-process replay completed in session `46215`, terminal chunk `f8ef0a`, both exit 0/PASS. The final independent checker, including an explicit guard for the displayed B-cusp observation, completed in session `5302`/chunk `ffa62e`, exit 0/PASS.

The root agent completed a full mathematical read of the producer, report and standalone checker. That review checked the width/normal-sum relation, complete four-box shapes, cusp count and divisibility, distinct-neighbor and positive-index clique arguments, all four necessary decorated-box domains, and the limited status of the surviving count assignment. It did not execute the certificate or reprove the inherited four-polytopes or carrier theorem. No additional review, external human review or full formal-kernel verification is claimed.
