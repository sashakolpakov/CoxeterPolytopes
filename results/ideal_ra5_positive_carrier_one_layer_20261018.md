# The positive carrier survives every first-layer facet and cusp condition

18 October 2026. **Universal attachment lemma and exact example; internally reviewed.** The positive-kappa sixteen-cell ideal RA4 polytope admits a simultaneous right-angled attachment in hyperbolic five-space with 457 genuine facets and all 258 original cusps closed by rectangular four-box links. The resulting chamber has infinite volume. Thus its first-layer Gram and cusp conditions do not obstruct its occurrence as a facet; a finite-volume completion is neither constructed nor excluded.

This is not merely a positive-semidefinite local diagram assignment. Every proposed wall has an explicit future timelike relative-interior facet witness, all pair products are checked, the original four-polytope is the exact intersection with its wall, and an explicit open ideal escape proves the volume statement. The input is the frozen [positive-kappa construction](ideal_ra4_positive_kappa_doubling_20261018.md).

## A universal one-layer lemma

Let Q be any finite-volume entirely ideal right-angled hyperbolic four-polytope. Choose Lorentz coordinates in its five-dimensional ambient vector space so that $p=(1,0,0,0,0)$ is a strict interior point, and normalize every ideal vertex to $v=(1,x)$. Let its outward unit facet normals be $m_i$; then $m_{i0}>0$. Embed this space as the orthogonal complement of a new unit spacelike vector n. Define

\[
 h_v=L v-n
\]

for every ideal vertex v, with a common real $L>1$. The old normals are extended by a zero final coordinate. The proposed five-dimensional chamber has walls n, all $m_i$, and all $h_v$.

These are exactly the possible forms of the cusp-opposite walls, up to an individual positive parameter. In a genuine rectangular four-box cusp containing Q, the opposite of n is asymptotic to n and orthogonal to the six old ridge-neighbor walls. Those six old normals span a rank-four subspace with null radical $\mathbb Rv$; their orthogonal complement in the old Lorentz space is $\mathbb Rv$. Thus the new normal lies in $\mathbb Rv+\mathbb Rn$. Its product with n is −1, so it is $\lambda_vv-n$, and strict negativity on the interior of Q forces $\lambda_v>0$.

There are finitely many vertices and facets. The quantities

\[
 a=\min_{\langle v,m_i\rangle<0}(-\langle v,m_i\rangle)>0,
 \qquad b=\min_{v\ne w}(-\langle v,w\rangle)>0
\]

are positive. Take L such that $La>1$ and $L^2b>2$. Then

\[
 h_v^2=1,\quad \langle n,h_v\rangle=-1,\quad
 \langle h_v,m_i\rangle=L\langle v,m_i\rangle\in\{0\}\cup(-\infty,-1),
\]

and

\[
 \langle h_v,h_w\rangle=1+L^2\langle v,w\rangle<-1\quad(v\ne w).
\]

Together with the old right-angled wall products, every pair is orthogonal, asymptotic or ultraparallel. Hence every finite dihedral angle of the actual intersection is right. The point $p-n/2$ is future timelike and strictly satisfies every wall inequality. No algebraic-field restriction is required by this lemma.

## Every wall is a genuine facet

The intersection with $n^\perp$ is exactly Q: its old inequalities define Q, and $\langle x,h_v\rangle=L\langle x,v\rangle<0$ for every future timelike $x\in Q$. The point p is a relative-interior witness for the n facet.

For each old ridge-neighbor wall $m_i$, choose a finite relative-interior point q of its facet in Q. Then $q-\epsilon n$, for sufficiently small positive $\epsilon$, is timelike, has zero product with $m_i$, and has strictly negative products with every other wall. In the exact example q is the average of that facet's normalized ideal vertices; the certificate chooses and checks a rational epsilon separately for each wall. This proves that all old facets extend to genuine four-dimensional facets of the new chamber.

For each new wall $h_v$, use the explicit vector

\[
 x_v=v+\frac{p-Ln}{L^2}.
\]

It is future timelike, with $x_v^2=-L^{-2}-L^{-4}<0$, and $\langle x_v,h_v\rangle=0$. Its product with n is $-1/L$. Every old-wall product is
$\langle v,m_i\rangle-m_{i0}/L^2<0$. For every other vertex w,
$\langle x_v,h_w\rangle=L\langle v,w\rangle<0$. Thus $x_v$ is a strict relative-interior witness for $h_v$, proving genuineness without a global vertex census.

## Original cusps and the global escape

At each original ideal vertex v, exactly eight walls meet: n, $h_v$, and the six old incident walls. Their Gram diagram is the product of four affine $A_1$ pairs. The new pair has positive sum $n+h_v=Lv$; the three old pairs have positive multiples of the same v. All remaining walls are strictly negative at v. Consequently a sufficiently small horospherical section is a bounded rectangular four-box. Each original ideal vertex is therefore an actual isolated cusp of the new chamber, with the required bounded local cusp link.

Nevertheless,

\[
 u=p-n
\]

is future null and has strictly negative product with **every** wall:
$\langle u,n\rangle=-1$, $\langle u,m_i\rangle=-m_{i0}<0$, and $\langle u,h_v\rangle=1-L<0$. By continuity, an open neighborhood of its ideal direction also strictly satisfies every inequality. This proves infinite volume: in the Klein ball there is an open spherical patch with an inward collar inside the chamber, and the hyperbolic volume density diverges as the radial coordinate tends to one. Equivalently, small outward spacelike perturbations of u still satisfy all wall inequalities, displaying a genuine open escape rather than an unclosed single ideal vertex. This argument does not assert that every possible finite-volume extension has the same escape; further walls could cut it.

The lemma therefore realizes the entire required one-layer facet/cusp star of every ideal RA4 Q, but does not solve its finite-volume extension problem.

## Exact positive-carrier instance

For the frozen sixteen-cell Q, $f_0=258$, $f_3=198$ and $\kappa=6$. Its saved rational coordinates have

\[
 a=1/25,\qquad b=1/1150.
\]

Choose $L=48$. Then $La=48/25>1$ and $L^2b=1152/575>2$. The new wall set has exactly

\[
 1+198+258=457
\]

members. The [producer](../scripts/fresh_ideal_ra5_positive_carrier_one_layer_20261018.py) checks all 104,196 pair products, all 457 genuine-facet witnesses, all 258 complete eight-wall cusp links, and the strict timelike center and ideal escape. Its [independent saved checker](../scripts/check_ideal_ra5_one_layer_saved_20261018.py) imports no producer and uses reverse-column rank elimination. Both retain the full inherited Q geometry and merely extend it to the new coordinate; no historical census or final high-dimensional cone enumeration is claimed.

The construction has rational coordinates in standard Lorentz space. This statement is not an arithmetic-lattice or finite-covolume claim. All local and actual-facet conditions checked here survive; completion of the opposite/new facets away from the original cusps remains the obstruction to investigate.

Generation (session `57676`, terminal chunk `160d1b`) and the separate fresh replay plus standalone helper (session `17764`, terminal chunk `d89074`) completed with exit 0/PASS. The root agent completed a full mathematical read of this report, the producer and the helper: the universal opposite-wall form, norm/products, genuine-facet witnesses, exact Q section, all eight-wall cusp boxes and open ideal escape were checked by source and proof review. That reviewer did not execute the certificate or reprove the frozen input Q. A further peer read has not been performed and is not claimed. The universal geometric argument is a written proof, not a formal-kernel proof. No new right-angled dimension bound, finite-volume RA5 example, or impossibility of completion is claimed.
