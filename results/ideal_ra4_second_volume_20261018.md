# Second volume for entirely ideal right-angled hyperbolic 4-polytopes

Date: 18 October 2026. **Status: internally checked mathematical proof of the second distinct volume, with two independent full agent reviews and exact construction checks.** Equality classification remains open. The universal argument below is separate from the finite construction computation.

The conclusion is about the **second distinct volume value**:

\[
 V_1=\frac{4\pi^2}{3},\qquad V_2=\frac{8\pi^2}{3}.
\]

The first value has the unique ideal right-angled 24-cell as its minimizer, by Kolpakov. The second is attained by reflecting and adjoining another 24-cell across an octahedral facet. **This note does not classify every polytope attaining the second value, and does not prove that the double is its unique minimizer.** All polytopes here have dimension four, finite volume, every vertex ideal, and every dihedral angle right.

## Literature audit and precise inputs

Targeted searches for second-smallest/second-minimal volume, the doubled 24-cell, and the four-dimensional ideal right-angled volume spectrum found no primary source asserting the second value. This is a search result, not a claim that no such result exists. Three-dimensional second-volume results and four-manifold volume spectra concern different classes.

The published inputs used here are particularly small:

1. Kolpakov, *On the optimality of the ideal right-angled 24-cell*, Algebraic & Geometric Topology **12** (2012), 1941–1960, published 27 October 2012. Lemmas 1–2 give the incidence and volume identities; Theorem 3 gives the unique minimum. Proposition 5 states that an ideal right-angled 3-polyhedron with a k-gonal face has at least 2k vertices, with equality precisely the k-antiprism. Section 4.3, Case D, pages 1951–1952 and Figure 6, explicitly gives the complete admissible octahedrite list through ten vertices; in particular its nine-vertex member is unique. [Published paper](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf).
2. The nine-vertex classification is also a consequence of the complete small octahedrite catalogue of Deza and Shtogrin cited there. We use only the nine-vertex entry, not an unbounded census or an assumed global vertex limit. The complete-list assertion in Kolpakov's Case D is the direct published input being invoked.
3. The regular ideal octahedron is the unique ideal right-angled 3-polyhedron with six vertices. The eight-vertex case is the square antiprism and there is no seven-vertex case; these also follow from Proposition 5 and the defect identities below. A later published discussion explicitly records these small cases: Alexandrov–Bogachev–Vesnin–Egorov, *On volumes of hyperbolic right-angled polyhedra*, Sbornik: Mathematics **214** (2023), 148–165, page 150. [Paper](https://www.mathnet.ru/php/getFT.phtml?jrnid=sm&option_lang=eng&paperid=9740&what=fullteng).

No external human review, publication, or Lean formalization of the new argument is claimed.

## Exact discrete volume reduction

Write \(v,e,r,f\) for the numbers of vertices, edges, ridges and facets of a four-dimensional polytope \(P\). A horospherical vertex link is a Euclidean rectangular box. It has eight vertices, twelve edges and six faces. Consequently

\[
 e=4v,\qquad r=3v+f,\qquad \sum_{R\in\Omega_2(P)} k(R)=12v.
\]

Here the second identity also uses Euler's equation for the boundary 3-sphere, and \(k(R)\) is the number of vertices of a ridge. Define

\[
 W(P)=\sum_R(k(R)-3)=3(v-f)\in3\mathbb Z_{\ge0}.
\]

Gauss–Bonnet, or the reflection-orbifold Euler characteristic, gives

\[
 \chi=1-\frac f2+\frac r4-\frac e8=\frac{v-f+4}{4},\qquad
 \operatorname{Vol}(P)=\frac{\pi^2}{3}(v-f+4)
 =\frac{4\pi^2}{3}+\frac{\pi^2}{9}W(P).
\]

Thus a volume strictly between the 24-cell and its double would require exactly \(W=3,6,9\), corresponding to \(5\pi^2/3,2\pi^2,7\pi^2/3\). No stronger divisibility is assumed. The remainder excludes these three values for **every** polytope in the class.

## Facet defect and the opposite-cusp inequality

A facet \(F\) is itself an ideal right-angled hyperbolic 3-polyhedron. Its vertices have valency four. If it has \(v_F\) vertices, then it has \(2v_F\) edges and \(v_F+2\) polygonal faces, so

\[
 d(F):=v_F-6=\sum_{R\subset F}(k(R)-3)\ge0.
\]

Every ridge lies in two facets; therefore

\[
 \sum_F d(F)=2W(P). \tag{1}
\]

If \(d(F)=0\), all its faces are triangles and it is the regular ideal octahedron. If \(d(F)>0\), it has a face of size \(k\ge4\); the antiprism bound \(v_F\ge2k\) implies \(d(F)\ge2\). In particular there is no defect-one facet.

Suppose \(P\) is not the 24-cell. A positive-defect facet exists: otherwise (1) gives \(W=0\), the volume formula gives the minimum, and Kolpakov's Theorem 3 identifies \(P\) with the 24-cell. Choose a facet of **minimum positive defect** \(m\). Let \(q(F)\) count its nontriangular faces and let \(s(F)\) count its ideal vertices whose intrinsic horospherical rectangle is not a square. Square/nonsquare is independent of the horosphere's scale.

Each nontriangular face of \(F\) has a distinct neighboring facet. It is non-octahedral because it contains that same nontriangle, so its defect is at least \(m\).

At each cusp \(x\) of \(F\), the face of the ambient rectangular box opposite its \(F\)-face belongs to another facet \(F_x\). Those two Euclidean rectangles are congruent. If the rectangle in \(F\) is nonsquare, then so is the one in \(F_x\); hence \(F_x\) is not a regular octahedron and \(d(F_x)\ge m\).

These opposite facets are distinct for distinct cusps and none is a ridge-neighbor of \(F\). More geometrically, the two supporting hyperplanes for opposite box faces are asymptotic, meeting only at that ideal point. The same hyperplane therefore cannot be opposite at a different cusp or intersect \(F\)'s hyperplane along a ridge. Equivalently, the local box links plus convexity of the facet intersection give the same distinctness. This is precisely the facet-versus-cusp distinction used in Kolpakov's Section 4.3; here it is applied only to the nonsquare cusps.

Counting these \(1+q(F)+s(F)\) distinct non-octahedral facets in (1) yields

\[
 \boxed{\quad 2W(P)\ge m\,[1+q(F)+s(F)].\quad} \tag{2}
\]

This bound does not bound the total number of facets or vertices. All other facets simply contribute nonnegative defect.

## Two exact cusp-shape facts

**Antiprisms.** For the ideal right-angled k-antiprism \(A_k\), every cusp rectangle has aspect ratio

\[
 2\cos(\pi/k).
\]

Here is a direct normal calculation, so the nonsquare assertion does not depend on a picture. Put \(c=\cos(\pi/k)\), \(h^2=(1-c)/(1+c)\), \(r^2=2c/(1+c)\), and \(b^2=2/[c(1+c)]\), taking positive square roots. Use Lorentz coordinates \((t,x,y,z)\) with metric \((-1,1,1,1)\). Place the two regular k-gons at heights \(\pm h\), radius \(r\), with angular offset \(\pi/k\). At an upper vertex \(u=(1,r\cos\theta,r\sin\theta,h)\), \(\theta=\pi/k\), the four incident wall normals are

\[
 N=(h,0,0,1),\quad M_0=(1,b,0,h),\quad
 M_2=(1,b\cos2\theta,b\sin2\theta,h),\quad
 L=(1,b\cos\theta,b\sin\theta,-h).
\]

They have squares \(r^2,\mu^2,\mu^2,\mu^2\), where \(\mu^2=2(1-c)/c\). Adjacent products are zero, the two opposite pairs are parallel, and \(bc=r\), \(\mu=2h/r\). After unit normalization,

\[
 M_0/\mu+M_2/\mu=(2/\mu)u,\qquad
 N/r+L/\mu=(h/r+1/\mu)u.
\]

For two opposite unit normals whose sum is \(\kappa u\), their Euclidean separation on a fixed horosphere is proportional to \(\kappa\). This follows by taking \(u\) to infinity in upper-half-space: the walls become parallel Euclidean planes with unit horizontal normals. Therefore the ratio of the two separations is

\[
 \frac{2/\mu}{h/r+1/\mu}
 =\frac{2r^2}{1+h^2}=2c.
\]

For \(k=3\) this is one. For every \(k\ge4\) it is strictly greater than one. The antiprism symmetries act transitively on its \(2k\) vertices, so \(s(A_k)=2k\).

**Nine vertices.** If \(v_F=9\), the antiprism bound excludes pentagons or larger faces. Thus \(F\) is an octahedrite; the complete published small list cited above has exactly one nine-vertex type. It is two regular ideal octahedra glued across a triangular face. Right-angled rigidity makes that combinatorial identification geometric. Its three gluing cusps have rectangles obtained by adjoining two squares along an edge, with aspect ratio two. Its other six cusps remain squares. Hence

\[
 d(F)=3,\qquad q(F)=3,\qquad s(F)=3. \tag{3}
\]

The rational doubled-24-cell certificate independently realizes eight such facets and checks all their cusp ratios: exactly three ratios two and six ratios one per facet. This verifies the claimed geometry of the model; the **completeness/uniqueness** of the nine-vertex case is supplied by the published classification, not by that example computation.

## Universal gap proof

First, any non-octahedral facet has at least two nontriangular faces. Indeed if there were just one, of size \(k\), its defect would be \(k-3\), contradicting the antiprism inequality \(d\ge2(k-3)\). If there are precisely two, of sizes \(k\ge l\), then

\[
 (k-3)+(l-3)=d\ge2(k-3),
\]

forcing \(k=l\) and equality. Such a facet is exactly \(A_k\). Thus a non-antiprism non-octahedral facet has \(q\ge3\).

Apply (2) to a facet of minimum positive defect \(m\ge2\).

* If it is an antiprism, then \(m=2k-6\), \(q=2\), \(s=2k=m+6\). Consequently \(2W\ge m(m+9)\ge22\).
* Otherwise \(m=2\) is impossible: the bound \(k\le v_F/2=4\) forces exactly two quadrilaterals, hence an antiprism.
* If \(m=3\), equation (3) gives \(2W\ge3(1+3+3)=21\).
* If \(m=4\), then \(v_F=10\) and \(k\le5\). A pentagon would give equality in the antiprism bound and hence \(A_5\), already treated. Otherwise its four units of defect are four quadrilaterals, so \(q=4\) and \(2W\ge4(1+4)=20\).
* If \(m\ge5\) and it is not an antiprism, then \(q\ge3\), giving \(2W\ge4m\ge20\).

Every case yields \(W\ge10\). Since \(W\) is a multiple of three, in fact

\[
 P\not\cong C_{24}\quad\Longrightarrow\quad W(P)\ge12
 \quad\Longrightarrow\quad\operatorname{Vol}(P)\ge\frac{8\pi^2}{3}. \tag{4}
\]

This is a universal local-incidence argument, not a census of four-polytopes below an arbitrary size. The only finite classification input is the complete nine-vertex **three-dimensional facet** case, forced by the minimum-defect argument.

## Exact double and attainment

The producer [fresh_ideal_ra4_double24_volume_20261018.py](../scripts/fresh_ideal_ra4_double24_volume_20261018.py) uses the Lorentz form \((-1,1,1,1,1)\). The original 24 unit outward normals are \((1,a)\), where \(a\) runs through the D4 roots with two entries \(\pm1\). In the Klein time-one section its vertices are the eight \(\pm e_i\) and the sixteen vectors \((\pm1/2,\pm1/2,\pm1/2,\pm1/2)\). Every vertex is null. All 10,626 active-four-wall subsets are checked, establishing that these are all vertices of the bounded halfspace intersection.

Reflect across \(n=(1,1,1,0,0)\) by \(x\mapsto x-2\langle x,n\rangle n\), omit the gluing wall on both sides, and identify repeated walls. This gives 38 distinct unit walls and 42 distinct ideal vertices. The exact halfspace products are nonpositive for every vertex and every retained wall. The resulting intersection equals the union of the two cells: on either side of the gluing hyperplane the corresponding original 24 inequalities hold, while all retained inequalities hold on both cells. Thus there are no unexamined extra vertices or escape directions.

All faces are reconstructed as nonempty intersections of facet vertex sets, with exact affine rank. Every ridge belongs to exactly two orthogonal walls, every edge to three facets, and every cusp to six walls forming three parallel pairs and orthogonal cross-pairs. The certificate gives

| Quantity | 24-cell | Double |
| --- | ---: | ---: |
| Vertices/cusps | 24 | 42 |
| Edges | 96 | 168 |
| Ridges | 96 | 164 |
| Facets | 24 | 38 |
| Triangular ridges | 96 | 152 |
| Quadrilateral ridges | 0 | 12 |
| Octahedral facets | 24 | 30 |
| Doubled-octahedron facets | 0 | 8 |
| Ridge defect W | 0 | 12 |
| Reflection orbifold Euler characteristic | 1 | 2 |
| Volume | \(4\pi^2/3\) | \(8\pi^2/3\) |

The double has six cusps on the gluing facet and 36 off it. Every cusp has eight incident edges, twelve incident ridges and six incident facets. Each doubled-octahedron facet has nine vertices, eighteen edges, eight triangular faces and three quadrilateral faces. The facet–vertex incidence check is \(30\cdot6+8\cdot9=6\cdot42=252\). Volume also follows directly by addition of the two disjoint interiors.

## Verification status and remaining equality problem

Final rational construction generation `f42b71` and separate replay `a81264` exited zero and passed, including complete original vertex enumeration, reflected walls, all face incidences and cusp Gram/shape tests. The standalone [saved checker](../scripts/check_ideal_ra4_second_volume_saved_20261018.py), run `6fe807`, independently repeated these checks with reverse-column elimination, reconstructed the full face intersection lattice, and verified seven rational-function identities in the antiprism calculation. After adding an explicit gluing-wall membership guard, its final run `cffdc9` exited zero and passed. It imports no producer. The universal geometric counting argument and the published facet-classification input are not replaced by this arithmetic check.

Root's independent full written-proof, construction-producer and standalone-helper mathematical read passed. It independently checked the precise published nine-vertex completeness and antiprism equality inputs, the normal/horosphere calculation, all minimum-defect cases, and the complete construction/checker accounting by reading. Construction's second full written-proof/producer/helper mathematical read also passed, with an independent reading of those same published premises. Neither reviewer executed the computation, reconstructed the historical classification census, or claimed external human review. The exact rational generation/replay/helper runs are the author's separate computational checks.

Equation (4) plus the double establishes the second **distinct value**, using the stated published inputs. It does not establish that every W=12 polytope is the double. That requires a new equality analysis; the coarse inequalities above are not all equalities at the double and do not themselves determine its full face lattice. The construction JSON's negative proof-status flags refer to what that construction computation establishes by itself; the universal theorem is the written argument above.
