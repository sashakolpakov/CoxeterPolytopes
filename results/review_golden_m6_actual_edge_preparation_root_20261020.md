# Root pre-execution review: one complete M6 interval

20 October 2026. **PASS for the finite reduction and code interface.**
This entry precedes execution; it does not claim the finite census has
passed. No scientific job, build or formalization was performed for
this review.

Read in full:

* `results/golden_dyadic_m6_actual_short_edge_preparation_20261020.md`,
  SHA256 `6a23b543e11d30d49810779d166c7cd4294ce4cc3a75172398805089bcd9dccc`;
* `scripts/check_golden_dyadic_m6_actual_edge_20261020.py`,
  SHA256 `bd3ec12320cf6e5e141cf87883b3ada2a9d4d5d4e5ad1a76c3c2a2891463f01c`;
* the pinned helper's field arithmetic, integrality, inverse and guarded
  main interface. Importing it does not execute its separate model check.

## Actual geometry and complete finite domain

With spatial coordinate e1=-c, the old actual a1 is precisely
`u=(phi-1;phi,0^5)`. The fourth-vertex full seal supplies the entire
literal M6 arrangement, so no norm4 transverse root or H3 contact can
be an unlisted continuation. The unimodular M6 primitive pool is 1,2.

After releasing beta6, every first wall meets a finite retained A5
corner before the explicit finite upper h. The only whole-component
finite extensions are isolated A1, endpoint A6, D6 at nodes2/4,
E6 at node3 and a norm1 endpoint B6. The length squareclasses and
field exclude other native products. Multiple contacts create a
nonfinite cycle. These are exactly the nine literal profiles.

For contact index k the norm identity is
`Z^2-6phi*a^2=6q-k(6-k)`; k=0 has all coordinates equal.
The upper has `ZH^2-6phi*aH^2=1`. Comparing collision parameters
therefore gives `a^2<=R*aH^2`, and the other embedding gives
`a'^2<=phi*R/6`. These bounds have no arbitrary height parameter.

The coefficient box is complete: `|a'|<9/5`; if B>=14 then
`a=B*sqrt(5)+a'>14*sqrt(5)-9/5>59/2`, whereas the exact upper
is below 59/2. Negative B is incompatible with a>0. The conjugate
bound then gives -1<=A<=9 for 0<=B<=13. There are 154 pairs
per R, six distinct R values, hence exactly 924 bound tests.

## Code and reconstruction

The quadratic helper uses exact Fraction arithmetic. Its sign and
O_phi-integrality tests are correct. The square-root routine is
complete: for nonrational x, a root's rational part has square
`(x.a +/- sqrt(N(x)))/2`; both signs are tried. Its two returned
roots give all roots in the field. The rational branch correctly
allows rational or sqrt(5)-rational roots.

The worker tests the exact inequalities inside the conservative box,
both square roots, positive Z, all six reconstructed coordinates,
the norm equation and all known actual beta/u inequalities. A true
first wall must pass these tests even if further walls are unknown.
It aborts on any additional accepted candidate. Norm1/2 vectors of
this unimodular lattice are automatically primitive; h also has a
unit spatial coordinate. No optimizer or restricted root prefix is
substituted for this finite domain.

I checked the endpoint formula and the six positive Gram minors:
the retained A5 plus h is B6, with final determinant1. The full
seven-dimensional map is computed from an actual O-basis, then
checked for integrality, square identity, all Gram products and
future orientation. Its graph interpretation is not assumed in
place of these checks. The child has a 45-second alarm and bounded
certificate/output; the parent must hold the shared lock and apply
the stronger process-tree/resource guard. Separate generation and
replay remain required.

## What a successful census would establish

After the finite result, h is an actual wall and the reversal maps
the actual full B6+H3 germ to the actual full germ at its opposite
endpoint. The integral orthogonal split extends this map by identity
on H3, so it preserves the full ambient arrangement and its chamber.
It makes w=A(u) actual as well. This is a finite-order involution,
not a nonreflectivity certificate.

The independently known compact P6 double lineage is an explicit
external compactness input for the nine-wall cone. Once all nine
walls are actual, acute projection places every finite vertex of
that compact cone in the full chamber. Convexity then proves equality
with the full M6 chamber, without assuming absence of hidden roots
from compactness alone. No compact dimension9 example follows.

## First execution: the expected candidate count was wrong

The initial bounded generation at worker pin `bd3ec123...` stopped at
the fail-closed assertion that exactly one integral acute candidate
survived. It took 0.046 seconds and wrote neither certificate nor
success ledger. A diagnostic-only assertion-message change, worker
SHA256 `f8698a35d2fa2fa940ea0da8ff95f7e3b3911b1a9d4785416ea92dc4cedd8b63`,
ran the same domain and stopped in 0.039 seconds, showing exactly two
candidates. No domain, glue or acuteness filter was narrowed.

The second candidate is the norm2 D6-profile vector

    g=(4phi^3;(4phi+3)^2,(4phi+2)^4)=2h+beta1.

Its time and spatial-sum coordinates are twice h's, so its collision
parameter is exactly the same. It is a positive nonsimple B6 root.
The full chamber lies in the h halfspace, since h is an integral
arrangement mirror oriented from the known old vertex. With the
actual beta1 halfspace this makes g's inequality redundant; equality
forces both h and beta1 to vanish, so g cannot define a facet.

This preserves the proposed unique SIMPLE h conclusion but falsifies
the original expected count of all surviving integral normals. The
source and checker must explicitly retain both candidates and verify
their tie and positive-root relation. Original pre-execution PASS
does not certify successful enumeration. Corrected generation and
separate replay are still required before promoting firstness or the
reversal.

## Corrected two-root interface

I read and rederived the complete correction at preparation SHA256
`2afd545c3c3166cd07381949b972fb012280b18ed673a090025b6aac30537dd2`
and worker SHA256
`9e48ab80a337f43925a1be7ee6761b02da5eaa0c60df6dfb649bc640ec5e6949`.
The nine profiles, coefficient box, square tests, reconstruction and
all acuteness filters are unchanged. The expected set is now exactly
the two specified profiles and vectors. The worker verifies g's
positive-root identity, explicit coordinates, norm and exact tie,
and retains both candidates in the certificate with h selected as
the sole simple root. The source proves the necessary codimension-two
redundancy without first assuming h actual. **PASS for this corrected
source/code interface; successful execution remains a separate gate.**
