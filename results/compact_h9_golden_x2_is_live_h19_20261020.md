# The x=2 golden binary extension is the existing H19 host

20 October 2026. Exact paper identification; no scientific computation,
worker, census, Lean or build was run. This is a duplicate-host stop rule,
not a nonreflectivity theorem for the still-live H19 host.

Consider the proposed family

    a=1+x+phi*x^2,   B_x=[[-2phi,-1],[-1,2a]],
    phi^2=phi+1.

At x=2 this gives B_2=[[-2phi,-1],[-1,6+8phi]]. The live discriminant19
host has binary Gram B_19=[[-2phi,phi],[phi,2]]. The explicit matrix

    U=[[1,-1],[0,-phi^2]]

satisfies

    U^t B_19 U=B_2,       det(U)=-phi^2 in O-units.

Indeed its first column has square-2phi, its product with the second is
2phi-phi^3=-1, and the second square is6+8phi. Thus this is an INTEGRAL
lattice isometry with integral inverse, not merely a K-similarity. The
E8 block is matched by the already calibrated native permutation
(0,1,2,3,4,5,6,7)->(0,7,1,2,3,4,5,6); alternatively use the same E8
Cartan ordering on both sides and the identity there.

The exact relation of the long classes is

    b_2=11+17phi=phi^4*(4+phi).

The isometry preserves the actual chamber sector as well. It fixes the
future point p=(1,0;0), and maps

    r_family=(phi-1,-2;0) --> phi^2*(1,2;0).

The target is a POSITIVE unit multiple of the known H19 long wall r.
The native E8 simple walls map with their calibrated outward signs.
Consequently the full inward E8+A1(long) germ at p is exactly the
already-used live H19 germ. Equality of integral mirror arrangements
then identifies the whole chamber, without a finite-prefix inference.

The family's suggested short vector c_2=(2,1;0) maps to
(1,-phi^2;0). It is not an additional actual wall of this chamber:
its product with the existing H19 wall old0=(0,-1;0) is2+phi>0.
The earlier binary short root is already

    d=phi^-2*(1,1;0),       d^2=2,
    U*d=(0,-1;0)=old0.

Thus the suggested new binary endpoint is blocked by a known actual
wall. The unit basis change identifies both the host and its selected
germ, so there is no reason to repeat the H19 edge search under the x=2
family name. The H19 construction status itself is unchanged.

The literal target form, native permutation and actual old0/r germ are
recorded in `compact_h9_discriminant19_mixed_norm_seed_20261019.md` and
`src/coxeter_polytopes/search/h19_vinberg.py`; no earlier H19 full-host
exclusion is being assumed.
