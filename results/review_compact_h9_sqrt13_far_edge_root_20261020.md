# Root review: quadratic13 far edge and chamber period

20 October 2026. **PASS, complete paper argument.** The new fixed
checker has not yet been executed at this review stage. The reviewed
source is `results/compact_h9_sqrt13_a9_far_end_and_period_20261020.md`,
SHA256 `86ffab56f4696298fac4ae34b8bc8772dd1e7a6daf2e95e48d56291e810bbcf2`.
The earlier B9 route b9226479 and A8 interval e3e62a6 are the explicit
inputs; both earlier workers now have generation and separate replay PASS.

I independently derived the orthogonal-plane quadratic form before the
producer supplied its diagonal basis. In coordinates x,z its entries are
`7eta+5`, `-(15eta-3)`, `13eta+34`, with determinant `-9/eta`.
I checked that P0,R0 are perpendicular with norms -eta,9 and that their
coordinate determinant is the unit eta. The isolated-profile equation
therefore has integral coefficients A,B, so the two separate split-prime3
exclusions apply at every height.

For each attached placement I independently expanded the native linear
forms

    A=2omega*x-(eta+4)z-(k+eta-5),
    Z=-9eta*x+9(eta+2)z+(eta+3)k-9.

They give the claimed mod9 congruence. The displayed reconstruction via
the integral P0,R0 basis proves sufficiency as well as necessity; there
is no missing denominator or assumption that the old A8 lattice glue is
unchanged. The projection norm k(9-k)/9 is the standard A8 fundamental
weight norm, giving exactly R=1,4,10 in the six attached placements.

I derived the finite upper independently from the integral point
inversion at P, rather than proposing it from a height search:
`w0 alpha8=a+alpha1+...+alpha8`, `b=J_P w0 alpha8`.
The unit norm of P makes J_P integral with integral inverse. Its future
orientation follows from fixing the future timelike P. The two endpoints
are `(eta+3)P0 +/- eta R0`, both finite with square -10eta; the plus
endpoint is forward from Q on the released edge. Thus b supplies an
ordinary finite upper before being declared first.

Strict feasibility at Q and contact by Q' give both sides of the
inequality `Z>=(eta+3)|A|`. Together with the norm equation this yields
`10A²<=R`. The conjugate bound and `R<=10` force the omega coefficient
of A to vanish. Thus the exhaustive domains are {0} for R1/4 and
{-1,0,1} for R10. Mod9 excludes the R1/4 rows; sqrt10 is not in K;
R10 permits k1 only. The A=-1 root is through Q and fails the STRICT
outgoing condition, leaving b alone. This exclusion is essential: it
cannot be silently removed by a cap compatibility filter.

The first-contact argument retains all eight actual walls and identifies
the full A9 endpoint germ. Consequently H maps actual inward germs of
one arrangement cell, proving H(C)=C. The argument would not be valid
for two partial corner lists. The inherited G has already been proved
to preserve the same full chamber and integral arrangement.

I separately expanded all four `(gamma_j,G gamma_j)` values and
`(P,GP)`. They agree with the source and give `trace(HG)=268eta+80>10`.
The four orthogonal gamma reflections do give the longest A8 element.
The strict trace excess proves infinite order, despite H and G each
being involutions. Distinct powers stabilize C, have trivial intersection
with the chamber's reflection group, and give distinct cosets. The
admissible anisotropic integral orthogonal group is an arithmetic
lattice, so the full reflection subgroup has infinite covolume.

The complete written conclusion is nonreflectivity of the literal
`<-eta>+I9` host. It excludes a compact reflection chamber whose walls
preserve this lattice. It does not exclude every commensurable lattice,
all quadratic13 forms, or compact H9 Coxeter polytopes generally; no
uniform dimension ceiling or new example follows. The pending checker
will corroborate finite algebra only. Full geometric Lean verification
is not claimed; the existing shared ChamberPeriod module checks the
final group/coset implication with explicit geometric and trace inputs.
