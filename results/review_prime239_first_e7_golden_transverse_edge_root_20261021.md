# Root review: the first norm239 transverse edge

21 October 2026. This is a complete manual mathematical review, not a
checker run, a Lean compilation or an external referee report.

Reviewed source:
`results/compact_h9_prime239_first_e7_golden_transverse_edge_candidate_20261021.md`

SHA-256:
`2a207fb82414203e34625d323238abc4c297541b6e34f23c4d28b16b7c3efa15`.

The root is not a contributing author of this candidate. I read all six
sections and independently checked the profile reduction, original-module
conditions and every finite row below. The inherited first binary edge,
primitive norm alphabet and finite spherical classification remain the
source's explicitly named inputs. Acceptance of this written theorem
also requires the separate outside-author review; this document does
not substitute for it.

## 1. The actual edge and its integral upper

The stated E8 dual columns have native covectors at nodes7 and6,
respectively. Their scalar products are2,6,3. Thus releasing node7 at
the inherited actual E8+A1 vertex gives the ray
`q+tau Omega7`, with square `-2D+2tau^2`; its products with the retained
walls are zero and with the released wall are `-tau`. An initial segment
belongs to the original chamber by the actual starting germ.

Expanding `phi(q-m)/2` in the ORIGINAL binary basis gives
`(13+22phi,-10-17phi)`. In particular the displayed half does not
introduce a lattice extension. Adding `phi^3 Omega7` gives a norm2
integral root R. Its products with q, m and Omega7 are
`-D phi`, `-phi` and `2phi^3`. Hence its positive collision is
`tau*=D/(2phi^2)`. The identity

    4phi^4-D=7-4phi>0

puts that collision strictly inside hyperbolic space. A norm2 vector
in this even integral module defines an integral reflection. It is
primitive: division by a nonunit would make its squared norm divided
by2 a nonunit square in O, whereas that value is1. This supplies an
actual finite mirror upper without assuming that R is first.

## 2. Why the list concerns every possible first chamber wall

Local finiteness on the compact segment from the starting vertex to
the upper ensures a first chamber boundary by the upper time. The
retained E7+m system has rank8. At its finite first endpoint the
spherical local reflection chamber has a rank9 simple system. The
finite classifications leave E7+I2(k), or an E8 completion attached
by a simple edge to node6 and perpendicular to m. No connected
rank9 finite diagram contains E7. This classifies first SIMPLE walls;
it does not wrongly classify every arrangement mirror meeting the ray.

For a norm2 first root, any nonzero finite contact with m is an
algebraic integer kappa with `1<=kappa<2`, `|kappa'|<2`. Writing
`kappa=A+Bphi` gives B0 or1 and precisely kappa1 orphi. A norm2b
root has a contact p with a retained norm2 root divisible by b, by
reflection integrality. If p is nonzero, `p^2/b=b eta^2` has positive
integer rational norm at least239. The physical finite angle and the
positive definite conjugate give rational norm strictly below16.
Thus every such contact is zero. These arguments justify all five
listed profiles, including the absence of a norm2b E8 completion.

All simple chamber walls are pairwise nonpositive. This supplies
`h=-r.(-e7)>=0` even though -e7 has been released at this edge.
No compactness or finite-volume assumption about the full chamber is
used here.

## 3. Original dyadic coordinates and the norm equations

The representation

    r=(Tq-kappa m)/2+h Omega7+epsilon Omega6

follows from the full retained covector conditions. E8 unimodularity
gives h in O. For completeness, the binary coordinates of q generate
the unit ideal: `det(q,m)=-2` puts2 in that ideal; q modulo2 is the
nonzero vector `(1+phi)(1,1)`, and O/(2) is the field with four
elements. Thus an original-integral binary vector forces T in O.
Reduction modulo2 then gives exactly `T=kappa mod2O`, since q and m
have the same displayed nonzero reduction. The converse integrality
also follows from that congruence. This does not treat (q,m) as an
O-basis; their additive index is4.

The norm2 equation and positive collision are consequently

    D T^2=4h^2+12epsilon h+12epsilon^2+kappa^2-4,
    tau=D T/(2h+3epsilon), T>0.

I checked both using q squared -2D, m squared2 and the dual-column
products2,6,3. In particular kappa0 forces T=2U with U in O.

## 4. Complete finite rows for epsilon0 and norm2

The physical upper together with the negative conjugate D' gives

    h^2 <= (4-kappa^2)phi^4/(7-4phi),
    |h'|^2 <= (4-kappa'^2)/4.

These inequalities include every wall at or before the upper. Using
`8/5<phi<81/50`, the physical bounds are h<7.3 and h<6.3 for
kappa0 and1. For kappaphi the bound is exactly h<=phi^3, since
`(7-4phi)phi^2=3-phi`.

The integer B in `h=A+Bphi` is `(h-h')/sqrt5`. For kappa0 the
strip gives B0..3, and its seven positive rows are precisely

    1, phi, 1+phi, 1+2phi, 2+2phi, 1+3phi, 2+3phi.

The rational norms of h squared minus1 are respectively

    0, -1, -5, -16, -11, -9, -45.

The zero row forces U0 and has no positive collision. Every other
row fails the necessary divisibility by239 in `D U^2=h^2-1`.

For kappa1 the strip `|h'|<.87`, h<6.3 gives exactly

    phi, 1+phi, 1+2phi, 2+2phi, 1+3phi.

I expanded `4h^2-3`; the pairs of golden coefficients are
`(1,4),(5,12),(17,32),(29,48),(37,60)`. Their rational norms,
using `N(A+Bphi)=A^2+AB-B^2`, are

    -11, -59, -191, -71, -11.

None is divisible by239. For kappaphi the smaller strip gives
exactly phi,1+phi,1+2phi. The coefficients of `4h^2+phi-3`
are `(1,5),(5,13),(17,33)`, with norms -19,-79,-239. Only the
last survives; it equals `D phi^2`, hence T=phi at the physical
embedding. It satisfies the actual parity constraint and gives
exactly R, at tau*. These are complete interval classifications,
not samples of a search.

## 5. The other two profiles

For the E8-completion profile the equation reduces to
`D U^2=(h+1)(h+2)` and
`tau^2=D(1-1/(2h+3)^2)`. The known upper yields h<2.15;
the conjugate equation gives strictly -2<h'<-1. Therefore B=1
and A=-1, so the only possible h is phi-1. The right side is
phi cubed, of rational norm -1, contradicting divisibility by239.

For a norm2b wall, integrality against the native E8 roots forces
h=bj with j in O. The binary parity gives T=2U. Its norm equation
is `D U^2=h^2-b`. The physical and conjugate upper bounds become

    0<j<1.1, |j'|<.45.

Here `b(7-4phi)=44-13phi` and `b'=33-17phi>5` justify the
stated rational estimates. These two strips force the golden
coefficient B=0. Then j is a positive rational integer of absolute
value below .45, a contradiction. No norm2b cut remains.

## 6. Full endpoint and exact scope

Every possible first simple wall at or before the integral upper is
now R itself. Thus the endpoint is the future ray

    Q=2phi^2 q+D Omega7,
    Q^2=-2D(7-4phi)<0.

It is feasible against ALL original walls because it is reached along
the actual chamber edge before its first collision. The eight
retained walls and R have the positive definite diagram E7+I2(5).
The finite local chamber therefore has exactly this simple facet
germ; additional mirrors in its finite local reflection group are
not additional chamber facets. The conclusion uses no binary-only
exclusion on the transverse ray and transports no norm79 identity.

**Root verdict:** no mathematical blocker in the complete written
first-transverse-edge theorem, at its explicitly inherited inputs.
It proves a second actual finite edge for this literal norm239 host.
It does not prove finite volume, chamber completion, reflectivity,
nonreflectivity, a new polytope or any new dimension bound. There is
no runtime or Lean credit from this manual review. A source-bound
outside-author full review is still required before the new acceptance
receipt can record this theorem.
