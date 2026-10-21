# A corrected small-face weighting for the cubic-seven arithmetic field

**Independently reviewed result:** a compact arithmetic hyperbolic Coxeter polytope with ground field **Q(cos(2π/7)) has dimension at most20**. The argument below uses a simultaneous quadrilateral-star rank constraint, rather than excluding individual marked spherical angles. It does not improve the unconditional compact bound29, and it makes no novelty or priority claim.

The [exact checker](../scripts/fresh_compact_global_bound_20261015.py) and [certificate](fresh_compact_global_bound_20261015.json) verify43,836 relaxed local graph cases, both symbolic determinant identities, their cubic-field remainders,101 finite-diagram graph records, the exceptional four corners, and the final h-vector arithmetic. The geometric reduction and all-rank arguments are written here; the finite computations do not replace them.

A [separate independent review](fresh_compact_global_bound_20261015_review.md) passes generation and saved replay. Its standalone Fraction polynomial engine reconstructs the determinant identities and remainders, an independent35,904-graph enumeration leaves the same nine supports, and122 finite graphs and all four corner matrices are checked. The hostile proof review separately checks arbitrary extra white contacts, rational branch elimination, the full-star rank condition and the assertion that the exceptional D₅ components are entire. No producer arithmetic is imported.

The separate [arithmetic field package](fresh_arithmetic_highdim_20261015.md) supplies the complete cubic-seven finite-label/Lannér restrictions and the weaker dimension28 fallback. Its source/data are unchanged. The present package does not alter any frozen construction or search cursor.

## Inputs and normalization

Normalize facet normals to square1, so a finite bond labelled m has product −cos(π/m). We use the standard compact Coxeter facts: vertex diagrams are spherical forests; elliptic subsets correspond to intersecting facets; affine subdiagrams are absent; a minimal indefinite subdiagram is Lannér. We also use arithmetic cyclic products and the symmetric nonnegative h-vector. These are the published inputs, not new computational conclusions. [Vinberg, §§1–3 and7, especially Proposition17](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).

For this ground field the finite labels are2,3,4,6,7,14. Lannér diagrams have order at most3, and every Lannér triangle has a7 or14 bond. The exact field package separately checks the seven allowed triangle types. In a finite diagram any bond of label at least6 is an isolated rank-two component: adding a nonright adjacent bond gives a three-node path determinant at most1−3/4−1/4=0.

Set c=cos(π/7), t=1−c². Both c and c² have degree3 overQ, with irreducible polynomials

\[
 P(c)=8c^3-4c^2-4c+1=0,
 \qquad 64(c^2)^3-80(c^2)^2+24c^2-1=0.
\]

We have c>9/10. Indeed c>cos(π/6)>4/5, P is strictly increasing on[4/5,1], and P(9/10)=−1/125. No upper bound on an unexamined integer label is being imposed: the complete label list is an arithmetic field consequence.

## A quadrilateral with no short marked corner

Let W be the spherical white support of an actual compact quadrilateral two-face, and let a,b,c₀,d be its four black sides in cyclic order. The symbol c₀ denotes a node; c continues to denote cos(π/7). Each W plus an adjacent black pair is positive definite. The opposite extensions W+a+c₀ and W+b+d are indefinite. They cannot be positive definite by the face correspondence, and cannot be affine in a compact chamber. Their minimal indefinite subsets L,M include both respective black nodes; hence they have order2 or3.

Assume no adjacent black pair is connected by a path of length at most2 in its finite corner diagram. The following reduction is exhaustive.

1. All adjacent black pairs are right. A Lannér2 opposite pair would be orthogonal to the other minimal indefinite diagram: any cross black bond gives a length-one marked path, while a cross black-to-white bond gives a length-two marked path. Two orthogonal indefinite diagrams violate Lorentz index one. Thus both L and M are triangles.
2. Their white nodes u,v are distinct. Their only possible connecting bond is u–v. A high bond involving a white node would be isolated in each finite corner and would force u–v right, again making L and M orthogonal. Therefore the high bonds are the two opposite black pairs a–c₀ and b–d, and u–v is nonright.

The checker enumerates every binary support of the two connected diagrams, all cross-bond subsets, shared or separate white nodes, and all choices of distinguished high bond. It only imposes the necessary finite-corner forest/isolation conditions and absence of a short marked path. The nine remaining supports are precisely the nonempty attachment patterns on either side of the two high opposite bonds, joined by u–v. This is a relaxed complete finite verification, not a numerical Gram test.

### Eliminate label14 and dense attachments

The two opposite I₂ blocks have a positive-definite direct sum. Eliminate them by a Schur complement. A singleπ/3 contact to an I₂(7) block costs more than

\[
 \frac{1}{4(1-(9/10)^2)}=\frac{25}{19}.
\]

A singleπ/3 contact to I₂(14), or twoπ/3 contacts to both roots of I₂(7), costs more than

\[
 \frac1{2(1-c)}>5.
\]

Larger contacts only increase these costs, since the inverse I₂ Gram has positive entries. If either opposite bond is14, or either white node touches both endpoints of its opposite7 bond, the two residual white diagonals are less than−4 and−6/19. Their off-diagonal remains their old product: the two white nodes have no cross contacts to the other black block. Its square is at most1. Thus their determinant exceeds

\[
 4\frac6{19}-1=\frac5{19}>0.
\]

This supplies a negative-definite two-plane, impossible in a Lorentz Gram. Consequently both opposite bonds are7, and every white neighbor touches exactly one endpoint of its associated opposite pair.

Let U consist of all white neighbors of a or c₀, and V of all white neighbors of b or d. They are nonempty and disjoint. Each triple consisting of an opposite7 pair and one of its white neighbors is Lannér. Applying the preceding nonorthogonality argument to every such pair of triples forces every U–V bond nonright. Since W is a forest, min(|U|,|V|)=1. Interchange the two pairs if necessary and write U={u}.

Each of b,d has at most one white neighbor: two neighbors v,w would give the cycle b–v–u–w–b, or its d analogue, inside the finite diagram W+b. Thus |V| is1 or2. These statements cover all white attachments, not merely the initially selected L,M.

Every bond touching the resulting white core has label3 or4. A label6 or higher next to another nonright bond violates spherical isolation, and label5 is excluded by the ground field. Any remaining white component attaches to this core at at most one node, because W is a forest. Eliminate all such branches and discard orthogonal white components. Leaf elimination has rational diagonal corrections: every squared branch bond is rational, and the branch is a finite tree. The surviving white diagonals p,q,r are rational, positive and at most1. Any nonempty attached branch decreases its attachment diagonal strictly.

The full star has rank rank(W)+3: its black projections span the Lorentz3 normal space of the compact quadrilateral. Thus the reduced six- or seven-node Gram is singular. Merely checking a local Lorentz Gram without this rank requirement would be insufficient.

## The rank condition leaves exactly one exception

### One white neighbor on each side

Write the two surviving white diagonals p,q, their bond magnitudeδ, and their contacts to black endpoints α,β. Put A=α²,B=β²,d₀=δ². The full six-node determinant is

\[
 (pq-d_0)t^2-(pB+qA)t+AB.
\]

All coefficients are rational and pq−d₀>0 because the reduced white Gram is positive definite. This is a nonzero polynomial of degree2 in c², so cannot vanish. This excludes, in particular, the apparent six-node path7–3–3–3–7 as an actual quadrilateral star. That path can have Lorentz signature; its excessive rank prevents it from being the required two-face star.

### Two white neighbors on one side

Name the white nodes v,u,w along a path, with c₀–u of magnitudeα, b–v of magnitudeβ, and d–w of magnitudeγ. Write δ,ε for u–v,u–w. After branch elimination the white Gram is

\[
 H=\begin{pmatrix}p&-\delta&-\varepsilon\\-\delta&q&0\\-\varepsilon&0&r\end{pmatrix}>0.
\]

Use A=α²,B=β²,C=γ²,d₀=δ²,e₀=ε² and L=δεβγ. The full seven-node determinant is

\[
 N(c)=(pt-A)(qrt-qC-rB+BC)
       -d_0(rt-C)t-e_0(qt-B)t-2Lct.
\]

This identity is checked directly from the literal7×7 matrix. All its coefficients are rational. The only point needing arithmetic beyond squared bonds is L: the cycle b–v–u–w–d–b has cyclic product a nonzero rational multiple of Lc, which belongs to the ground field. Since c belongs to that field and L belongs toQ(√2), the cubic field's absence of a quadratic subfield forces L rational.

The coefficient of c in the remainder of N modulo P is

\[
 \frac{pqr-d_0r-e_0q-8L}{8}=\frac{\det H-8L}{8}.
\]

Singularity therefore implies detH=8L. On the other hand, positive-definite determinant monotonicity under increasing diagonal entries gives

\[
 \det H\le 1-d_0-e_0\le\frac12\le8L.
\]

The last inequality uses δ,ε,β,γ≥1/2. Equality throughout forces p=q=r=1 and δ=ε=β=γ=1/2. The equality of white diagonals removes every off-core white branch. Substitution into the remaining remainder gives

\[
 \frac{(4A-1)(16c^2-9)}{64}=0,
\]

so α=1/2 as well. The entire affected diagram is now fixed: white A₃, both opposite bonds7, and all five other bonds3. Its determinant factors as

\[
 N(c)=\frac{(4c+3)(8c^3-4c^2-4c+1)}{64}=0.
\]

This is an exception to the plain distance-two rule, not a contradiction. Two corners have type A₄+A₁ and the other two have type D₅. At each D₅ corner the marked nodes are the far leaf of the length-two arm and one of the two short leaves, at distance3. Because the white diagonals stayed1 and all other black–white contacts were exhausted, these D₅ diagrams are **entire connected components** of the corresponding full vertex forests. Arbitrary other white components are orthogonal to them.

## Corrected weights and the global contradiction

Assign a planar-angle weight from its finite vertex diagram as follows:

* weight1 if its two marked nodes have graph distance at most2;
* one additional unit if they form an isolated I₂(7) or I₂(14) component;
* weight1/2 if they are the far leaf of the length-two arm and a short leaf in an entire D₅ component;
* weight0 otherwise.

The last rule concerns exactly two pairs in each D₅ component; they have distance3 and hence receive no base unit.

A triangular two-face has its three black nodes forming a Lannér triangle: every proper black subset with W is spherical, and the minimal indefinite subset has order at most3. It therefore has a high7/14 bond. The corresponding corner has that bond as an isolated I₂ component and receives weight2. A quadrilateral either has a short marked corner of weight1, or is the classified exception and receives1/2+1/2 from its D₅ corners. A face with k≥5 sides has nonnegative weight, which is at least5−k. Thus every two-face has total weight at least5−k.

For a connected finite diagram of rank r, the distance-at-most-two pair count is2r−3 for a path of rank at least2, and2r−2 for D/E trees. Rank1 contributes0. The extra high-I₂ unit makes its total2=2r−2. D₅ alone has total8+1=9=2r−1; every other component contributes at most2r−2. These formulas follow directly by counting edges and pairs of neighbors; the checker verifies101 representatives and their exceptional marked pairs.

For a vertex of total rank n≥6, either its diagram is one component other than D₅, or it has at least two components, each with deficit at least1 from2r. Hence every vertex total is at most2n−2.

Double-counting planar angles and using the symmetric nonnegative h-vector gives a strict average lower bound

\[
 \frac{n(n-6)}8\quad(n\text{ even}),\qquad
 \frac{(n-1)(n-5)}8\quad(n\text{ odd}).
\]

At n=21 the average weight is strictly greater than40, while the vertex upper bound is40. Explicitly,

\[
 f_2-50f_0=\sum_{i=0}^{10}(10-i)(11-i)h_i>0,
 \qquad 5\cdot50-\binom{21}{2}=40.
\]

For n=21+2u the non-strict lower-minus-upper difference is u²/2+5u; for n=22+2u it is u²/2+11u/2+2. The strict face-count inequality therefore excludes every n≥21. This proves the stated bound using the published inputs and the independently reviewed geometric argument above.

The consequence for the previously established high-dimensional arithmetic field list is specific: the cubic-seven field disappears in dimensions at least22, leaving the two quadratic candidates. This does not eliminate either quadratic field or settle unconditional compact dimension29.

## Verification scope

Generation and saved replay pass. The replay command is:

```sh
python scripts/fresh_compact_global_bound_20261015.py --verify
```

The checker imports SymPy for symbolic polynomial determinants. It imports no earlier producer arithmetic. It binds the frozen cubic-field package and compact h-vector package. The separate arithmetic/proof review now passes; the primary certificate's historical status string still says `PROOF_CANDIDATE` because its frozen source/data were preserved. No finite sample of stars is being substituted for the all-white-attachment reduction, and no external priority search completeness is claimed.
