# Discriminant59: an explicit infinite-order chamber period

20 October 2026. Paper derivation; no scientific program, worker, census,
optimizer, Lean or build ran in this lane. The global conclusion uses
the **complete first-wall premise** in Section4. Its independent census
and exact execution are separate interfaces; a compatible prefix alone
does not discharge that premise.

## 1. Literal host and the two corner frames

Put O=Z[phi], phi^2=phi+1, b=7+5phi and delta=5+12phi=phi*b. Use

    L=E8(O) orthogonal-sum [[-2phi,-1],[-1,4+2phi]].

Write p=(1,0;0), c=(1,1;0), and let e_i denote the native integral E8
basis, with outward chamber roots alpha_i=-e_i. The vectors p,c,e_i
are an integral basis of the FULL lattice. In the equivalent binary
basis p,c the Gram is

    [[-2phi,-phi^3],[-phi^3,2]].

Let

    w=omega0=(4,5,7,10,8,6,4,2),  w^2=4,
    (w,alpha0)=-1,  (w,alpha_i)=0 for i=1,...,7,
    h=2p+phi^3*c,   h^2=-2delta,   (h,c)=0.

The retained alpha_1,...,alpha_7 form D7. Its orthogonal complement in
the ten-dimensional space is U=span(p,c,w), of signature(2,1).

The actual long Bugaenko facet has normal

    r=phi^2*p-2c,       r^2=2b.

The archived Bugaenko side wall with index0 is

    a=phi*p+phi*w,      a^2=2,      (a,r)=0.

Its vertex with native D7 is the literal archived facet vertex
[0,2,3,4,5,6,7,9]. Here the retained native E8 roots1,...,7 correspond
to archived wall indices9,2,3,4,5,6,7. Thus D7+r+a is an ACTUAL full
finite vertex germ of the ambient chamber C, by the already established
full Bugaenko r-facet interface. A future-pointing apex on its ray is

    P=2p+phi*w=h-phi^3*c+phi*w,       P^2=-4/phi.

The inherited data are
`data/known_examples/bugaenko_compact_h8.json`, SHA256
`c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c`.
The actual-facet construction is in
`compact_h9_golden_discriminant59_edge_route_20261020.md`; determinant
similarity alone is not used to recognize the new facet.

Put

    t=38+61phi, z=133+213phi, A=10+18phi,
    v=t*h+z*w,          H=A*h+t*w.

The exact identities

    z=delta*(5+9phi),   A=2z/delta,
    delta*t^2-2z^2=-b

give

    v^2=2b, (v,c)=0, (v,H)=0, H^2=-4/phi.

Both v and H are perpendicular to the retained D7. The apex H is
future-pointing since A>0 and (h,H)=-2delta*A<0. The two orthogonal
frames (r,a,P) and (v,c,H) in U therefore have the SAME diagonal Gram
(2b,2,-4/phi), with the same time orientation.

## 2. An integral isometry, with no genus assumption

Define G to fix every alpha_i for i=1,...,7 and set

    Gp=(56+92phi)p+(119+193phi)c+(99+160phi)w,
    Gc=(36+59phi)p+(76+123phi)c+(63+103phi)w,
    G(alpha0)=alpha0+(14+23phi)p+(30+48phi)c+(25+40phi)w.
                                                               (1)

These are integral columns in a full integral basis. Since
w=-(4alpha0+5alpha1+7alpha2+10alpha3+8alpha4+6alpha5+4alpha6+2alpha7),
they imply

    Gw=-(56+92phi)p-(120+192phi)c-(99+160phi)w.

Substitution, using phi^2=phi+1, now gives

    Gr=v,       Ga=c,       GP=H.                              (2)

For example G(p+w)=(phi-1)c, so Ga=c; the other two identities follow
by r=phi^2*p-2c and P=2p+phi*w. Thus G preserves the orthogonal-frame
Gram in U. It fixes D7 pointwise, so it preserves the full lattice Gram
over K. All entries in(1) are integral and det(G)^2=1, so its determinant
is a unit and its inverse is integral. It is an actual O-isometry of L,
not merely a K-similarity or an isometry of an index sublattice. Equation
(2) and the future choice of H show it preserves time orientation.

For direct checking in the ORIGINAL native basis (x,y;e_i), the same
columns are

    Gp=(175+285phi,119+193phi;(99+160phi)w),
    G(e_y)=(-63-103phi,-43-70phi;-(36+57phi)w),
    G(e0)=(-44-71phi,-30-48phi;e0-(25+40phi)w),
    G(e_i)=e_i, i=1,...,7.

They follow from e_y=c-p and e0=-alpha0, and independently expose the
integrality on the original host rather than suppressing the binary
half-glue. The automatic short-facet glue lemma is useful motivation,
but is not needed as an unverified integrality shortcut here.

## 3. Infinite order

On U, in the basis p,c,w, the three diagonal entries are
56+92phi,76+123phi,-99-160phi. On D7 the map is the identity. Hence

    tr(G)=40+55phi>10.

A finite-order real ten-dimensional matrix has roots of unity as
eigenvalues and absolute trace at most10. Consequently G has infinite
order. No asymptotic height argument, numerical eigenvalue or a product
of the known order-two Bugaenko symmetry is used.

## 4. The actual first-wall interface and chamber preservation

**Required premise.** The complete outgoing D7+c edge from the actual
E8+c vertex has v as its first simple wall. In particular its endpoint
has the FULL finite germ D7+c+v, with the displayed outward normals.
An all-height exclusion up to a proven upper endpoint, including every
finite profile, is needed. Post-selection compatibility with known caps
alone is insufficient. The independent complete finite-profile census
is to be bound separately; no first-prefix inference is made here.

Under this premise, (2) maps the full old actual germ D7+r+a onto the
full new actual germ D7+v+c, with future apices P and H. Every integral
isometry of L permutes its ENTIRE integral-mirror arrangement: it
conjugates any integral reflection to another integral reflection.
Therefore G takes the full chamber C to a chamber having precisely C's
inward germ at the new vertex. Chambers with that same germ coincide.
Thus G(C)=C.

This step is what a bare explicit Lorentz isometry or a pair of finite
prefix vertices would fail to prove. The actual old vertex comes from
the full calibrated r facet; the actual new vertex must come from the
complete edge premise. No classification or genus uniqueness for
v-perp is necessary.

## 5. Global conclusion once the premise is discharged

Let W be the subgroup of O^+(L) generated by ALL integral reflections.
It acts freely on its chambers: W intersect Stab(C) is trivial. Since
every G^k stabilizes C and G has infinite order, the cosets G^k W are
pairwise distinct. Indeed G^(k-l) in W would both stabilize C and lie
in W, forcing G^(k-l)=1 and k=l.

Hence W has infinite index in O^+(L): the L59 host is nonreflective.
This is a direct full-group obstruction, not the failure of a finite
wall prefix. It does not require an assertion that a particular partial
polytope is compact, nor an orthogonal-lattice nonreflectivity transfer
that omits oblique ambient mirrors.

Its scope is the LITERAL displayed integral discriminant59 host. It is
not an exclusion of all golden binary extensions, all quadratic forms
over Q(sqrt5), or their entire commensurability classes. It supplies no
new dimension ceiling or compact H9 example. Until the first-wall
interface has its independent proof and bound source, Section5 is the
conditional consequence of that explicitly stated premise.
