# Root source review: cubic B16 logic kernel

22 October 2026. **Source-only mathematical/interface PASS; uncompiled.**
I read all137 Lean lines and the complete assumption map. This is not
elaboration, an actual axiom audit or a new geometric formalization.
The module's author is sharp_ra_foundations_audit. Root authored the
written cubic proofs and reviews this separate formal source at that
declared application boundary.

Source `formal/CoxeterProofs/CubicB16ShortCapKernel.lean`:
SHA-256 `ce8b53cafa9cd18a1520545618bb854578d0e740a1c2c23a811fa16ac0901d74`.
Map `docs/cubic_b16_short_cap_kernel_logic_20261022.md`:
SHA-256 `4f680190280224617a239cc9e3ab5f6e5c1f91c0121c9e55a48e03ad6a490181`.

The first two declarations use genuine Int coefficients. The trace
inequality4T<=S+15 with S<=28 bounds T by10; the explicit cases0..10
then combine with17<=S and4S<=T squared to leave exactly S25,T10.
The sign-dependent coefficient strips have the correct directions:
positive n uses the lower v and upper u bounds; negative n reverses
the multiplication signs. Their cleared constants2205,1845,2470
match .2205,.1845,.247. All integral n=-6..6 are covered and only
m0,n=-1 remains. Field integrality, the full basis, the divided
differences, |n|<=6 and the remaining l6 identification are outside,
as the map says. This does not assume the desired pair as an input.

The ordered-ring portion derives24eta<1 from positive B-24, eta
and B*eta1. With4h=2-eta this gives96h>47. Multiplying its positive
gap by the positive square(1-r) squared, then using the half-tangent
identity, gives `47r^2-478r+47<0`. The displayed exact completion
checks algebraically:

    (99-1000r)(478000-47(1000r+99))
      =1000000(47r^2-478r+47)-138647.

If1000r<=99 both left factors are nonnegative, while the right
side is negative. Thus the energy margin follows from primitive
height and half-tangent identities; no rational-function monotonicity
axiom or final energy margin is supplied.

The funded-reference inequality divides to
`eB>=base+.26(S-f)-.0338`. With base>=.89026, f<=3.311,
S>3.5 and eI>.099, its exact total is strictly greater than1.00460.
The final declaration derives the I margin from the preceding logic
before combining it with that reference. No output inequality is
assumed as its own premise. Geometry of the funded vector and its
positive coefficient, actual orthogonality and the link from field
coefficients to B/eta/J remain explicit written application inputs.

There is only a Std import and no custom axiom or proof placeholder.
The six `#print axioms` commands are requests, not six actual outputs.
I found no mathematical or assumption-map blocker in the prepared
source. Its tactic elaboration and actual dependency output still
require a resource-admitted single-module compilation followed by a
SEPARATE audit. Existing successful modules must not be rerun. No
worker, compiler, import execution, test, Lean or TeX process ran for
this review, and no runtime verification credit is granted.
