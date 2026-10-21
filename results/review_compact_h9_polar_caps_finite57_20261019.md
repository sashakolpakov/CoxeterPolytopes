# Independent review of the discriminant19 polar caps

19 October 2026. **Written exact algebra and geometric review: PASS within
the archival scope below.** No matrix program, checker, scientific search,
optimizer, Lean process or build was executed. The small archived coordinate
and incidence files were read; their completion verifier was not rerun.
All reviewed source bytes were preserved.

## Reviewed inputs and scope

* New report `results/compact_h9_discriminant19_polar_caps_research_20261019.md`:
  `792634e935866ed78d95d0a93c04fdd7a8484f3605b5c50a66c1d376c9599f1b`.
* Predecessor including its21st-wall addendum,
  `results/compact_h9_discriminant19_mixed_norm_seed_20261019.md`:
  `70fa010732d837d5c044701b4c18ed5e68cb36c18cc0954dc7d7f369250a5c87`.
* Native coordinate archive `data/candidates/h19_vinberg_prefix.json`:
  `91c16b581e5763d6a283b6818f0e087e941fbb987b6b05a7e283febf1e32d0f1`.
* Complete H8 example `data/known_examples/bugaenko_compact_h8.json`:
  `c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c`.

Here H19 is the rank-ten discriminant19 form, with hyperbolic dimension9.
The result is a25-wall arithmetic Coxeter seed containing two complete
compact H8 facets and at least83 actual finite vertices. The24-wall
predecessor is infinite-volume. Finite volume of the25-wall seed is
unresolved. This is neither a compact H9 example nor a complete mixed-norm
mirror census.

## 1. Native covectors, norms and reflection integrality

I checked the native E8 index convention, Ctheta=e6 and Comega0=e0.
The sparse covectors of k,w,a,d,z follow from the displayed coordinates
and the stated binary block. For instance the E8 coefficient vector of w
divided by b is

    phi^3 theta-phi e6.

Its C-image is phi e5+e6, exactly as reported. For a the C-image is
phi e4+e5+phi^2 e6; for d it is
e0+(3+4phi)e5+phi^5 e6+2phi e7. The listed binary covectors follow because
all four new vectors have phi*x+2y=0. For z, Comega0=e0 gives its one-entry
E8 covector immediately.

Pairing each vector with its sparse covector gives

    w^2=z^2=2b,       a^2=d^2=k^2=2.

In particular the new short roots really have norm2, independently of
their asserted transport from the H8 example. The long covectors are
divisible by b, so reflections in w,z are integral. Reflections in the
norm-two roots are integral because both root and Gram entries are
integral. Their displayed y coordinates are units, proving primitivity
for the four newly added vectors w,a,d,z. The earlier k has norm2 and is
the transported primitive H8 root.

The polar relations and heights also check:

    v^2=(b/2)phi^(-8),   w=2phi^4 v,
    q=phi^(-4)z,        q^2=2b phi^(-8)=2(115-71phi)>0,
    height(w)=height(z)=phi^9.

The equality b phi^(-2)=7-3phi verifies the displayed E8 coordinates of q.

## 2. Complete product and strict-center ledger

For any archived vector u=(x,y;z0,...,z7), the four product columns reduce
to these linear functionals:

    (u,w)=b[-phi^5 x+phi z5+z6],
    (u,a)=-b phi^3 x+phi z4+z5+phi^2 z6,
    (u,d)=-b phi^6 x+z0+(3+4phi)z5+phi^5 z6+2phi z7,
    (u,z)=b[-phi^5 x+phi^2 z0].

Substituting every retained old coordinate row, plus r and k, agrees with
the ENTIRE compatibility table. The remaining new/new entries also check,
including a.d=0 and

    z.w=-4b phi^4,    z.a=-2b phi^3,    z.d=-2b phi^6.

For two norm-two roots, every nonzero finite product is -1 or-phi, giving
pi/3 or pi/5. All other nonzero short/short entries are strictly below-2.
All long/short nonzero entries are -b times a positive value at least1,
so normalized magnitude exceeds sqrt(b)/2>1. The listed long/long entries
are strictly ultraparallel after normalization as well. No parallel or
forbidden finite-angle product was silently substituted for a dotted one.

The strict-center values use an important normalization: its E8 part is
C^(-1)1/10. Thus the spatial correction in pairing with a root is the
SUM OF THAT ROOT'S E8 COORDINATES divided by10, rather than the sum of its
covector entries. Using this gives precisely all four displayed negative
new-wall products and the earlier (c,k)=-(83+139phi)/10. For example
sum(a_E8)=86+113phi and sum(d_E8)=304+509phi; these give the stated a,d
center products after adding their binary terms.

Acute projection of this SAME strict timelike center onto any one wall
is timelike and strict on every other wall. Therefore every inequality
is an actual facet, including newly appended inequalities. The mirrors
really cut the predecessor's hyperbolic interior, rather than merely
removing a spacelike direction outside hyperbolic space. The statement
that such mirrors are outside the predecessor reflection subgroup also
uses the standard Coxeter-polyhedron fundamental-chamber theorem; it is
not an orbit enumeration or a claim about an arbitrary generating set.

## 3. The H8 basis transport includes the integral coordinates

The nine-wall correspondence in Section3 matches the archived Gram:
the eight images form E8 with fork at old3, and old9 has its pi/5 bond
at the long-arm endpoint old0. These are independent Lorentz bases,
because the archived eight E8 simple roots together with positive root8
are independent. Root8 has control coordinate1, so they form an integral
basis of the archived H8 lattice.

I checked the two missing-root combinations directly against the archived
coordinates. Root0 minus phi*root8 has the eight simple-root coefficients

    2-2phi, 4-3phi, 6-4phi, 5-3phi,
    4-2phi, 3-phi, 2, 3-2phi

in order1,2,3,4,5,6,7,9. Root10 minus phi^4*root8 has coefficients
phi on root4,2phi on root5,1+3phi on root6,2+4phi on root7, and-phi on
root9. These reproduce exactly the report's native expressions for a,d.
Negative coefficients cause no integrality problem: this is a basis
isometry, not a positive root combination or a word in the old group.

Consequently the eleven normals on w-perpendicular space are the full
archived Bugaenko H8 system, with its complete compact chamber. The strict
center selects the correct halfspaces and time component. Every other
wall is strictly ultraparallel to w. With the common strict center, the
whole w hyperplane lies in the retained side of each such wall, so none
can cut down that compact facet. The same reasoning preserves the full
original r-facet. The k-coordinate/covector and original spatial-index
permutation also match the archive, so the21st-wall interface is genuinely
included in this review rather than inferred from the earlier20-wall PASS.

## 4. The83 vertices and the volume distinction

The complete H8 archive contains41 finite vertices. Completeness and their
certified feasibility remain inherited from that archived exact example;
this review did not rerun its Vinberg/vertex verifier. Each eight-wall
intrinsic elliptic vertex extends with its orthogonal r or w normal to an
actual ambient nine-wall elliptic vertex. Other walls are strict there
because they are ultraparallel to the corresponding facet hyperplane.

The r and w hyperplanes are ultraparallel, so these two41-element finite
vertex sets are disjoint. The initial p is preserved, remains finite, and
is strict against r,w and every new wall. Thus there are at least83
actual finite vertices. Appending z preserves both entire facets and p,
so the same count is valid at25 walls. It is not a complete ambient
vertex count.

For24 walls, q is a positive multiple of z and has nonpositive products
with every retained wall by the complete ledger. Its norm is strictly
positive. Perturbing toward the strict center gives a STRICT spacelike
feasible point. The segment from the future timelike center crosses the
future light cone strictly inside all finitely many wall halfspaces,
exposing an open ideal patch. This proves infinite volume at24 walls,
despite the83 finite vertices.

The new z wall cuts q by a strictly positive product. Nothing proves that
all other escapes have been removed. The25-wall chamber's final cone,
vertices/cusps and volume remain unresolved. No complete root search,
finite-volume conclusion, K3 nef marking or dimension record should be
inferred from this scoped construction.
