# Quadratic13: exact short-facet inheritance and a compact C3 prefix

20 October 2026. **Paper proof; no scientific execution.** This continues
the literal family L_n=<-eta>+I_n, eta=(3+sqrt13)/2, over O=Z[omega],
omega=(1+sqrt13)/2. The complete norm alphabet1,2 and the first two
edge calculations are proved in
`compact_h9_sqrt13_b9_two_edge_route_20261020.md`, SHA256
`b922647902d32beccd98fbe5a2781ddcbb9cafaeb10a12902eb3ebbce809547f`.
That proof has independent written review. Here C_n means the chamber
of the ENTIRE integral reflection arrangement with the standard B_n
germ at p=(1;0^n).

The new conclusions are: the actual coordinate-short facet of C_n is
exactly C_(n-1), and L_3 admits an explicitly compact nine-wall Coxeter
polyhedron. Thus low rank3 does not exclude the dimension9 host. There
is no compact4 or compact9 claim, and no novelty claim for the abstract
three-polyhedron.

## 1. Equality of the restricted and intrinsic arrangements

Let f=e_n, of norm1, and H_f=f-perp. Its intrinsic lattice is literally
L_(n-1). Normalize any integral reflection normal r to norm q=1 or2,
using the complete norm alphabet and unit rescaling. Put t=(r,f) in O.
If H_r meets H_f inside hyperbolic space, then |t|<sqrt(q). At the
other real embedding the entire form is positive definite, so Cauchy
gives |t'|<=sqrt(q). Writing t=a+b*omega therefore gives

    |b|=|t-t'|/sqrt13 <=2sqrt2/sqrt13 <1.

Thus t is an ordinary integer. The crossing possibilities are exactly

    q=1: t=0;             q=2: t=0,+1,-1.

The projection r-t*f is integral in H_f and has norm q-t^2, hence1 or2.
It is an intrinsic integral reflection normal. In the only oblique case,
q=2 and |t|=1, the intrinsic normal has norm1. This argument concerns
ALL intersecting mirrors and does not assume fundamental dihedral
angles for an arbitrary mirror pair.

Conversely every intrinsic norm1 or2 normal extends by a zero final
coordinate to an ambient integral reflector orthogonal to f. Therefore
the restriction of the complete ambient mirror arrangement to H_f is
EXACTLY the complete intrinsic arrangement, with repetitions ignored.
No omitted oblique restrictions or index-submodule identification is
being used. A noncrossing mirror has no zero on H_f and cannot subdivide
an intrinsic chamber.

To pass from arrangements to actual facets, choose a point z in the
relative interior of the f-facet of C_n, avoiding every other mirror.
Let D be its cell in the restricted arrangement. Every ambient wall
inequality has constant sign on D, unless it is identically zero on H_f;
its sign agrees with z. Consequently D is exactly the open facet cell,
not a proper subset or a union of intrinsic chambers. Taking closures
gives the actual facet. The B_n germ projects to the standard B_(n-1)
germ: the adjacent long root -e_(n-1)+e_n projects to the short root
-e_(n-1), and the earlier chain is unchanged. Thus that facet is C_(n-1).

In particular compactness of C_n implies compactness of C_(n-1).
If one literal C_m is noncompact, no higher C_n can be compact. This
is a theorem for the displayed integral family, not its entire
commensurability class or arbitrary Coxeter polytopes. No assumption
that every metric facet of an arbitrary Coxeter polytope is Coxeter
enters the proof.

## 2. Explicit C3 normals and the two new finite uppers

Work in coordinates (x;z1,z2,z3) with form -eta*x^2+z1^2+z2^2+z3^2.
Use the five-wall seed

    a1=(0;-1,1,0),      a2=(0;0,-1,1),     a3=(0;0,0,-1),
    v =(1;omega,0,0),   d =(eta;eta+1,eta+1,1).

The last short coordinate facet is the compact C2 quadrilateral. The
first two finite edge arguments specialize to this rank: the only
finite B2 extension is B3, and the six one-root extensions of the
retained short A1 plus the long v are precisely the previous six rows.
F4 cannot intervene: two norm1 roots could not have an integral product
-1/2, so a pair of adjacent short simple roots is unavailable.

At p release a3, retaining A2 on a1,a2. Its edge is
p+tau*(0;1,1,1). The primitive norm1 root

    r=(2eta;2eta+1,2eta,2eta)

has r^2=1 and products (r,a1)=-1, (r,a2)=0, (r,a3)=-2eta and
(r,v)=(r,d)=-eta-1. It is a genuine finite upper, since

    (6eta+1)^2-3eta*(2eta)^2=1,
    tau_r=2eta^2/(6eta+1)<sqrt(eta/3).

There are only six finite A2-extension profiles: isolated norms1,2;
an A3 extension at either end, of norm2; and a B3 extension at either
end, of norm1. Put Z equal to the sum of the three spatial coordinates.
Their equations are Z^2-3eta*x^2=R, with R=3,6,4,1 respectively.
The proved upper gives x^2<=4eta^2R and x'^2<=Reta/3. All rows have
0<x<17, |x'|<13/5. If x=a+b*omega, this gives 0<=b<=5 and six
possible integer a around floor((sqrt13-1)b/2). Thus a complete future
first-wall calculation is finite; it is unnecessary for the compact
prefix proof below, and r is not assumed first in the full arrangement.

At q=(omega;eta,0,0), release a3 and retain the orthogonal norm2 pair
a2,v. The edge is q+tau*(0;0,1,1). A second norm2 upper is

    s=(eta+2;2eta+1,eta+1,eta+1).

It is orthogonal to a2,v and cuts at
tau_s=eta*omega/(eta+1)<sqrt(eta). Its remaining products are
(s,a1)=-eta, (s,a3)=-eta-1 and (s,d)=1-eta. All are native/ultraparallel.

For completeness this edge also has a tiny bounded complete problem.
Write a root as

    (X/2)q-(ell/2)v+(Z/2)(e2+e3)+(m/2)(e2-e3).

The seven finite profiles have (norm,ell,m)
(1,0,0),(2,0,0),(1,1,0),(1,0,1),(2,1,0),(2,0,1),(2,1,1).
The last is A3; two mixed double bonds would be affine, so the omitted
(1,1,1) is not finite. In every row

    Z^2-eta*X^2=R=2norm-ell^2-m^2,
    X^2<=R*omega^2,       X'^2<=R*eta,
    X congruent eta*ell modulo2,    Z congruent m modulo2.

The coarse bounds 0<X<5 and |X'|<4 give b=-1,0,1,2 and a=-5,...,6
for X=a+b*omega, a complete48-pair superset. Again no unbounded prefix,
unit exponent or numerical angle tolerance is involved.

## 3. A canonical nine-wall Coxeter polyhedron

Use the proven involution G from the first two germs, fixing e3. Put

    r'=Gr=(3eta+1;4eta+2,3eta+1,2eta),
    s'=Gs=(2eta-1;2eta,2eta,eta+1).

The norms of (a1,a2,a3,v,d,r,r',s,s') are
(2,2,1,2,2,1,1,2,2). Here is the complete strict upper triangle of
their Gram matrix; zero entries are included.

| | a2 | a3 | v | d | r | r' | s | s' |
|---|---|---|---|---|---|---|---|---|
| a1 | -1 | 0 | -omega | 0 | -1 | -eta-1 | -eta | 0 |
| a2 | | -1 | 0 | -eta | 0 | -eta-1 | 0 | 1-eta |
| a3 | | | 0 | -1 | -2eta | -2eta | -eta-1 | -eta-1 |
| v | | | | -1 | -eta-1 | -1 | 0 | -eta |
| d | | | | | -eta-1 | 0 | 1-eta | 0 |
| r | | | | | | -2eta | -1 | 0 |
| r' | | | | | | | 0 | -1 |
| s | | | | | | | | -1 |

All negative entries of magnitude1 give label3 for equal norm2 and
label4 for mixed norms1,2. Every remaining nonzero pair is strictly
ultraparallel: omega>2, eta>3, eta-1>2. Thus the nine halfspaces are
an acute native Coxeter system, not a numerical near-Coxeter proposal.
Every root is integral of norm1 or2, so every reflection is integral.

There is a common strict future point p+(0;3,2,1)/100. Substitution
using 3<eta<7/2 gives all nine inequalities strictly and its square is
-eta+14/10000<0. The future affine section x=1 is bounded, because
the first four walls imply

    0<=z3<=z2<=z1<=eta/omega.

Moreover the homogeneous cone has no nonzero point with x=0: these
inequalities then force all three spatial coordinates0. A negative x
is impossible as well. Thus this affine section sees every nonzero
feasible ray; no direction at infinity is hidden by the chart.

## 4. Fourteen finite vertices force compactness

There is a useful general elementary fact. Suppose a finite acute set
of normals has a strict future interior point x, and a selected subset
has positive definite Gram A. Since A has nonpositive off-diagonal
entries, A^-1 is entrywise nonnegative. For example write A=D-B,
with D positive diagonal and B entrywise nonnegative. Positive
definiteness makes the spectral radius of D^(-1/2)BD^(-1/2) less than1;
the convergent geometric series proves the inverse claim.

Project x orthogonally off the selected normals. Their product vector
with x is strictly negative, so the projection ADDS a nonnegative
linear combination of those normals to x. All unselected inequalities
remain STRICT, by acuteness. The square becomes
x^2-y^T A^-1 y<0, and its product with x is this same negative square,
so it stays future. Therefore a rank-three positive definite subset
gives an actual finite vertex with exactly its three active walls.
Different triples give distinct vertices. This is an actual feasibility
proof, not an assertion that any spherical subgraph of an arbitrary
nonacute configuration is realized.

The following fourteen triples all have positive definite Gram:

| Triple | Finite type |
|---|---|
| a1,a2,a3 | B3 |
| a2,a3,v | B2+A1 |
| a3,v,d | B3 |
| a1,a3,d | B2+A1 |
| a1,a2,r | B3 |
| a2,v,s | 3A1 |
| v,d,r' | B3 |
| a1,d,s' | 3A1 |
| a2,r,s | B2+A1 |
| a1,r,s' | B2+A1 |
| v,r',s | B2+A1 |
| d,r',s' | B2+A1 |
| r,s,s' | B3 |
| r',s,s' | B3 |

Each wall appears in one of these actual vertices, so all nine are
actual facets of the bounded three-dimensional affine polyhedron.
Every bounded3-polyhedron with F facets has at most2F-4 vertices:
Euler gives V-E+F=2 and minimum vertex degree3 gives2E>=3V.
Here F=9, so V<=14. The fourteen constructed finite vertices therefore
exhaust ALL affine vertices; in particular no additional spacelike or
null vertex has been silently ignored.

The entire bounded affine polyhedron is the convex hull of these
vertices. Its squared Euclidean spatial norm is convex, and every
vertex has norm strictly below eta. Hence the whole polyhedron lies
strictly inside the hyperbolic ball. It is a compact Coxeter
three-polyhedron with exactly9 facets and14 finite vertices.

The full integral chamber C3 is a subchamber of this compact integral
Coxeter polyhedron. Consequently it too is compact, and its full
reflection group has finite index in O^+(L3). Alternatively the
nine-wall reflection group already has finite covolume and is contained
in the full integral group; the uniform arithmetic quotient gives the
same conclusion. This proves reflectivity of the literal L3 but does
NOT identify the exact fundamental chamber of the full integral group.

## 5. Research consequence and verification boundary

The short-facet inheritance argument is exact, but rank3 is alive:
it cannot provide a hereditary obstruction to L9. The compact prefix
gives explicit caps and a finite local geometry that may guide rank4;
it is not a claim that the same nine walls close a higher-dimensional
chamber. Rank4 and rank9 remain unresolved.

All calculations here are paper expansions in eta^2=3eta+1. A small
fixed algebra checker could check the nine norms, all36 pairs, strict
center and14 positive Gram triples; its purpose would be to corroborate
the universal projection-and-Euler proof, not to assume a historical
census or run an unbounded Vinberg recursion. No such checker has been
prepared or executed in this lane. Independent source review is still
required before treating this successor as settled evidence.
