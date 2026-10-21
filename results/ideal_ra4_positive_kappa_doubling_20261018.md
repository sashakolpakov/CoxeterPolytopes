# A positive-kappa ideal right-angled four-polytope

18 October 2026. **Final exact construction, with two independent full mathematical reviews.** This packet disproves the proposed universal inequality $3f_0-4f_3+24\le0$ for ideal right-angled hyperbolic four-polytopes. It does not construct an ideal right-angled five-polytope or improve a dimension bound.

The construction is an iterated reflection double of the regular ideal 24-cell. Its final face vector and volume are

\[
 (f_0,f_1,f_2,f_3)=(258,1032,972,198),\qquad
 V=\frac{64\pi^2}{3},\qquad \kappa=3f_0-4f_3+24=6.
\]

Thus the positive carriers required by the [five-face incidence identity](right_angled_five_face_carrier_20261018.md) really exist as four-polytopes. Whether they can occur together as the facets of an ideal right-angled five-polytope remains open here. In particular, their existence does not realize any incidence arrangement in dimension five.

## The actual doubling operation

Use the Lorentz metric $\operatorname{diag}(-1,1,1,1,1)$ and outward unit facet normals. Let $n$ be a facet normal of an ideal right-angled polytope $P$, and let $s_n$ denote reflection in its wall. The polytope

\[
 D=P\cup s_nP
\]

is exactly the intersection of the halfspaces with normals $m$ and $s_nm$ for all old facet normals $m\ne n$.

Indeed, for $x\in P$, acuteness gives $\langle m,n\rangle\le0$ and $\langle x,n\rangle\le0$, and hence

\[
 \langle x,s_nm\rangle
 =\langle x,m\rangle-2\langle m,n\rangle\langle x,n\rangle\le0.
\]

The same holds on the reflected copy. Conversely, a point of the asserted intersection with $\langle x,n\rangle\le0$ is in $P$; a point on the other side reflects into $P$. This proves the halfspace-union identity and convexity. The two interiors are disjoint, so volume doubles. In the Klein ball the vertex set is the union of the two ideal vertex sets: the union is the convex hull of those sets, and every distinct point on the strictly convex sphere remains extreme. No finite vertex is created.

Facets adjacent to the gluing facet are orthogonal to its wall and merge with their reflected copies; all other facets occur in two copies. In this packet the saved normals and vertex incidences explicitly verify every merging and the absence of other identifications or redundant walls. Each remaining wall contains four linearly independent null vertex rays and is therefore a genuine facet. Every cusp has six incident walls with three disjoint asymptotic pairs and all other products zero. These checks also directly verify the right-angled local structure after each operation.

If the gluing facet has $b$ ideal vertices, it is an ideal right-angled three-polyhedron with $b+2$ polygonal faces. Consequently

\[
 v'=2v-b,\qquad f'=2f-b-4,\qquad
 \kappa'=2\kappa+b-8. \tag{1}
\]

For the bounded level comparison below, these formulas hold for an arbitrary facet double. Orthogonal neighbors merge; a nonorthogonal facet has a finite interior point strictly on its original side of the gluing wall, so its original and reflected facets remain distinct. There are no further cross-copy wall identifications: if $m=s_n\ell$ for two old unit outward normals, then $\langle m,n\rangle=-\langle\ell,n\rangle$; both sides are nonpositive by acuteness, so both vanish and $m=\ell$. Thus only orthogonal neighbors merge. No other old facets coincide. An adjacent facet with $a$ vertices meeting the gluing facet in a $k$-gon becomes a facet with $2a-k$ vertices; a nonadjacent facet retains $a$ vertices in each copy.

## An explicit three-step path

Start from the frozen [exact double-24-cell certificate](fresh_ideal_ra4_double24_volume_20261018.json), whose actual halfspace and complete original vertex census are inherited. It has eight nine-vertex double-octahedral facets. Choose the length-two path with unit normals

\[
 F=(1,0,1,-1,0),\quad G=(1,0,1,0,-1),\quad
 H=(1,0,1,1,0).
\]

The products are $\langle F,G\rangle=\langle G,H\rangle=0$ and $\langle F,H\rangle=-1$. Both consecutive intersections are quadrilaterals. Double first across $F$, then the retained wall $G$, then the retained wall $H$.

After the first double, $G$ has $2\cdot9-4=14$ vertices, while the chosen copy of $H$ still has nine: $H$ and $F$ are asymptotic, not ridge-adjacent. After doubling across $G$, the chosen $H$ wall has fourteen vertices. The complete exact sequence is:

|24-cell copies|Gluing-facet vertices just used|Vertices|Facets|$\kappa$|
|---:|---:|---:|---:|---:|
|2|Initial frozen double|42|38|−2|
|4|9|75|63|−3|
|8|14|136|108|0|
|16|14|258|198|6|

The final volume is sixteen times $4\pi^2/3$. Independently, the universal ideal RA4 volume formula $(v-f+4)\pi^2/3$ gives the same result. This formula and the original regular 24-cell are the published inputs from [Kolpakov, 2012](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf); the new reflection and incidence data are checked here.

## The first positive level within whole-polytope doubling

There is a small complete level comparison, without a claim of minimality among all ideal RA4 polytopes. Start with one regular 24-cell and repeatedly double the entire polytope across one facet. Facet transitivity makes the two-cell result unique. Its facets have six or nine vertices.

If the next double is across a six-vertex facet, (1) gives $\kappa=-6$. Every resulting facet has at most $2\cdot9-3=15$ vertices, since a polygonal gluing face has at least three vertices. Therefore every ensuing eight-cell double has $\kappa\le2(-6)+15-8=-5$.

If the next double is across a nine-vertex facet, (1) gives $\kappa=-3$. Two adjacent nine-vertex facets of the frozen double meet in a quadrilateral, as verified for the complete eight-facet domain. Thus a merged nine-vertex facet has fourteen vertices; a merged octahedron has at most nine, and unmerged facets have at most nine. Every ensuing eight-cell double has $\kappa\le2(-3)+14-8=0$.

The initial 24-cell has $\kappa=0$ and its double has $\kappa=-2$. Hence sixteen cells is the first whole-polytope doubling level at which positive kappa can occur, and the explicit path realizes it. This says nothing about arbitrary gluings, other reflection chambers, the least possible positive kappa, or the minimum volume of an arbitrary positive-kappa ideal RA4 polytope.

## Positive kappa is unbounded

The example also rules out a universal constant upper bound on kappa. Every ideal right-angled three-dimensional facet has at least six vertices. Therefore any further facet double obeys

\[
 \kappa_{t+1}=2\kappa_t+b_t-8\ge2\kappa_t-2.
\]

Starting from $\kappa_0=6$, induction gives $\kappa_t\ge4\cdot2^t+2$. Every operation is again an actual ideal right-angled polytope by the same halfspace-union proof. No new numerical computation is required for this consequence. The five-face identity consequently needs compatibility among its facets or other incidence restrictions; a universal constant bound on an individual positive carrier is unavailable.

## Verification and literature scope

The [producer](../scripts/fresh_ideal_ra4_positive_kappa_doubling_20261018.py) reconstructs the three reflection steps from the frozen double, checks the complete resulting normal and ideal-vertex sets, facet rank, strict center, all known pair products, cusp links, exact merging counts and formulas (1). The [standalone saved checker](../scripts/check_ideal_ra4_positive_kappa_saved_20261018.py) imports no producer code and checks the same saved transformations with independently written arithmetic and reverse-column rank elimination. It inherits the original 24-cell census and the written universal doubling argument; it does not enumerate all four-subsets of the final 198 facets.

The literature search on 18 October 2026 did not locate a claim of the proposed inequality. That search is not a novelty proof. The recent [Whitehead preprint](https://arxiv.org/abs/2609.05389) concerns three-dimensional ideal right-angled polyhedra and does not supply this four-dimensional inequality. The present construction directly refutes the inequality regardless of its literature status.

Generation `6ba641`, separate-process fresh replay `1d5a60`, and standalone verification `5738ef` all completed with exit0/PASS. Root and construction independently read the complete producer, written proof and standalone helper: both mathematical/accounting reviews PASS. Their scope includes the halfspace-union and vertex-completeness proof, cross-copy identification argument, genuine facets and cusp links, the three-step path, the bounded first-positive doubling level and the unbounded-kappa corollary. Neither reviewer executed the scripts or reproved the frozen original 24-cell census. The universal doubling argument is a written proof supported by exact finite data, not a formal kernel verification. No optimizer, historical classification census or dimension-bound computation is used.
