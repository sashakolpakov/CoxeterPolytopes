# A new H20 projection and three bounded cusp links

The construction retains the preceding 43 walls, adds a cap that cuts their saved escape and closes another finite edge, then adds three mutually orthogonal null-source caps. The final chamber has **47 actual facets and three explicitly bounded ideal-vertex links**. Four distinct finite vertices and two finite edges from the two stages remain certified. A renewed strict spacelike escape proves that the chamber still has **infinite volume**.

This is a geometric, noninducing extension of the complete H18 source with an A1² support. It is outside the previously closed direct inducing-wall families. The retained four-wall conditional nonquasi-arithmeticity witness persists, but the current infinite-volume chamber is not asserted to be a lattice. No arithmetic reflection search in dimension20 is used.

The three exact producers and certificates passed generation and saved replay:

- [New projection and finite edge](fresh_h20_new_projection_20261015_exact.json), [producer](../scripts/fresh_h20_new_projection_20261015_exact.py).
- [Three bounded cusp links](fresh_h20_new_projection_20261015_cusps.json), [producer](../scripts/fresh_h20_new_projection_20261015_cusps.py).
- [All-profile gate on31 parallel-derived null rays](fresh_h20_new_projection_20261015_null_gate.json), [producer](../scripts/fresh_h20_new_projection_20261015_null_gate.py).

These are producer checks; independent reviews are separate artifacts. All earlier checkpoint files remain unchanged. The exact constructions have no numerical-scout dependency.

## The escape-cutting cap

Retain the [frozen43-wall certificate](fresh_h20_joint_route_20261015_exact.json), with source normals0–39, A1 support normals40–41 and the preceding cap42. The physical form is $Q_0\oplus2I_2$. Solve for $Y\in\mathbb Q(\sqrt2)^{19}$ with zero products at

\[
1,2,3,4,5,7,8,9,10,11,13,14,15,16,18,20,21
\]

and products $-1$ at source walls12 and37. The exact solution is future timelike with $Y^2=-15/2$. Add

\[
R_{43}=(Y,-\sqrt{19}/2,0),\qquad R_{43}^2=2.
\]

Its43 pair relations are18 right, two π/3, one parallel (to source wall6), and22 dotted. In particular it is right to support41 and dotted to support40 and cap42. The saved43-wall escape has new product

\[
-\frac{3001081}{287500}+\frac{292431}{57500}\sqrt2
+\frac34\sqrt{19}>0.
\]

All946 pairs and44 actual facet witnesses are checked exactly. The field basis is
$1,\sqrt2,\sqrt3,\sqrt6,\sqrt{19},\sqrt{38},\sqrt{57},\sqrt{114}$,
stored as eight rational coefficients. The complete finite-neighbor census on the new cap checks20 possible20-root supports and finds exactly one:

\[
\{1,2,3,4,5,7,8,9,10,11,13,14,15,16,18,20,21,37,41,43\}.
\]

Replace43 by40 to obtain an old finite vertex. The two physical vertex vectors are future timelike, zero exactly on their respective supports and strict at all other walls. Their common19-root support is positive definite; each exchanged wall is strict at the opposite endpoint. Hence they bound an actual finite edge. The two previously stored vertices are also rechecked, and all four active sets are distinct. This is not a complete ambient vertex census.

## Three null-source caps at the first admissible scale

Purely timelike source projections are strict at every source cusp. They cannot change its incident link: the old bounded source link still has an unbounded A1² transverse quadrant. This explains why completing source cusps requires a different projection type.

The three source normals $N_{17},N_{23},N_{24}$ have norm1 and mutual product$-1$. Thus

\[
c_0=N_{17}+N_{23},\quad c_1=N_{17}+N_{24},\quad
c_2=N_{23}+N_{24}
\]

are future null source vectors with $(c_i,c_j)=-2$ for $i\ne j$. Set

\[
W=(0^{19},-\sqrt2/2,-\sqrt2/2),\qquad U_i=c_i+W.
\]

Each $U_i$ has norm2 and meets both support walls at π/4. The three $U_i$ are mutually orthogonal because $W^2=2$. They become walls44–46.

All135 new pair relations are exactly63 right,12 π/3, six π/4 and54 dotted. Each new cap has four π/3 source contacts:

| Cap | Null source pair | π/3 source contacts |
|---|---|---|
|44|(17,23)|0,8,16,20|
|45|(17,24)|4,12,16,20|
|46|(23,24)|0,4,8,12|

For a cap $t c_i+W$ with $t>0$, any one of these norm-two source walls has product$-t$. The first possible nonright Coxeter cosine is$1/2$, so $t\ge1$. The chosen scale1 is therefore the least admissible positive scale, rather than a large-scale dotted-only repair.

At each $c_i$, the producer checks all47 inequalities and reconstructs the complete23-wall zero set. Its Gram has rank19 and its physical normal span has rank20. Its four connected components have ranks8,8,1,2; each has a strictly positive kernel and a positive-definite principal deletion. Each positive kernel's physical sum is a positive multiple of $c_i$. The last component is exactly the affine C2 triangle on walls40,41 and$44+i$. These facts certify a bounded full19-dimensional Euclidean cusp link and an actual proper ideal vertex.

The common center is strict and timelike. All47 facets have explicit timelike projection witnesses and exact strict signs at the other46 walls. All four stored finite vertices, and therefore both previously certified finite edges, survive strictly against the new caps. No complete finite-vertex census on the three new cusp caps is claimed. These are three specified source cusps, not a complete source cusp census.

## Renewed escape and a uniform31-ray obstruction

Keep the source part of the preceding43-wall escape and change its transverse coordinates to $(-7/10,-1/10)$. The resulting $E$ is strict at all47 walls, points to the same future sheet as the center, and has square

\[
E^2=\frac{19945826604799}{661250000000}
-\frac{87506613029}{4132812500}\sqrt2
\approx0.2197964>0.
\]

This is an exact infinite-volume certificate, not an unsuccessful numerical volume test.

There is also a uniform obstruction to filling this escape with further caps whose source projections lie on a source parallel-pair ray. The40 retained source walls have exactly31 parallel pairs. The separate gate reconstructs their normalized sums $c$, checks that all31 are future null and source-feasible, and verifies all465 distinct-ray products are negative. Thus the31 directions are distinct.

For an arbitrary proposed cap normalized to square2, write

\[
r=t c+w,\qquad w\in\mathbb R^2,\quad w^2=2.
\]

Acuteness excludes negative$t$. At$t=0$, acuteness to the two support normals puts$w$ in their negative orthant, which prevents a strict common interior. Hence$t>0$. Every nonincident source wall$N_j$ gives the necessary all-label inequality

\[
t\ge
\frac{\sqrt{N_j^2/2}}{-(c,N_j)},\qquad
L(c)=\max_j\frac{\sqrt{N_j^2/2}}{-(c,N_j)}.
\]

The transverse part of$E$ has square1, so Euclidean Cauchy–Schwarz yields

\[
(r,E)\le t(c,E_{\rm src})+\sqrt2
\le L(c)(c,E_{\rm src})+\sqrt2.
\]

The certificate proves $(c,E_{\rm src})<0$ and
$L(c)^2(c,E_{\rm src})^2>2$ for every one of the31 rays. The smallest exact squared margin is

\[
\frac{38065153644889}{330625000000}
-\frac{134418676}{1653125}\sqrt2
\approx0.1385579>0.
\]

Thus **every real transverse profile**, with any finite Coxeter labels or disjoint-wall relations, retains$E$ on these31 source rays. Any finite coherent extension using only those projection rays remains infinite volume. The three explicit cusp repairs fall within this theorem.

The31 rays exhaust the literal parallel pairs, **not a proved complete cusp list**. A complete source-cusp census would also have to rule out maximal parabolic diagrams without an affine A1 component. General mixed source residuals, other source-cusp directions, or changed retained source walls are not excluded.

## Bounded numerical provenance

The [swap diagnostic](fresh_h20_new_projection_20261015_scout.json), [producer](../scripts/fresh_h20_new_projection_20261015_scout.py), is separate from every exact proof. It swaps the two transverse contacts of the preceding1,856 saved numerical recipes and checks all43 old relations with finite labels2–60 and tolerance$10^{-7}$. Of1,833 acute proposals,1,694 pass that bounded numerical angle test; six cut the preceding escape. They form two triples with positive products within each triple; across triples there are three parallel and six dotted pairs. The one simplest cutter above is exactified from its literal linear equations, not from a floating vector.

The renewed$E$ is numerically strict for all1,694 accepted swapped proposals, with maximum pairing about$-0.148849$. This is only a diagnostic for the saved bank. It is not an exact bank theorem or an all-profile exclusion.

All four new scripts support saved replay with the --verify option; use OPENBLAS_NUM_THREADS=1 for the numerical scout. No larger cap bank or further source census was run in this package.
