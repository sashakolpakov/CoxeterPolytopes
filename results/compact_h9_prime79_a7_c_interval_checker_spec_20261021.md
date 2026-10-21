# Norm79 A7+c interval: exact checker specification and correctness

21 October 2026. SOURCE PREPARATION ONLY. The new worker has not been
imported, executed or compiled; its fixtures, counters, minimum and germ
have no execution credit. Root alone may admit one guarded run after full
independent source review. This is one literal-lattice interval certificate,
not an H9 example, period, nonreflectivity theorem or dimension bound.

Worker: `scripts/check_compact_h9_prime79_a7_c_interval_20261021.py`,
SHA256 `59926019e0c228bb8364038660bed8e37b601bc771f4894fe4e024952390eda6`
(21,943 bytes, 514 lines). The author previously reviewed the mathematical
domain b4de9185 and contributed inherited exact-checker foundations; this
algorithm statement is not an independent review of its own code.

## 1. Mathematical inputs and original module

The exact domain proof is
`results/compact_h9_prime79_a7_c_finite_upper_and_complete_interval_contract_20261021.md`,
SHA `b4de91853f6fd625831d591554cb47f7a61400fde49c233017ef905e8b24c3fa`.
Its Sections 1--5 establish the actual starting germ, all-height exclusions,
finite upper and complete first-SIMPLE model/box theorem. Its inherited
geometric input is
`results/compact_h9_nonunit_coupling_prime79_full_base_and_first_transverse_20261021.md`,
SHA `33f6cb9dd20ba91e7e7974586075309f2f437564f2e8dd978a16594349e48da3`.
The complete primitive mirror alphabet, finite spherical classification,
arrangement-angle list, anisotropy, local finiteness and simplicity at a
finite corner retain precisely those written scopes. The worker does not
replace their proofs by a table assertion.

All ten-vector computations use the ORIGINAL basis X,Y,e0,...,e7 over
O=Z[phi], phi^2=phi+1. The binary Gram is [[-2phi,a],[a,2]], with
a=3+3phi and b=13+18phi; the orthogonal E8 block has diagonal 2 and
products -1 on 02,13,23,34,45,56,67. Put D=18+31phi,
p1=(2,-a;0), c=(0,-1;0), Q=omega1. Then p1^2=-2D, Q^2=8 and
the actual outgoing edge, releasing -e1, is p1+tau Q, tau>0.
The retained eight normals are -e0,-e2,-e3,-e4,-e5,-e6,-e7,c.

The worker solves the E8 Cartan system over Q(phi), rather than trusting
a copied inverse. It requires Q=(5,8,10,15,12,9,6,3) and checks

| j | 0 | 2 | 3 | 5 | 6 | 7 |
|---|---:|---:|---:|---:|---:|---:|
| I_j=Q.omega_j | 5 | 10 | 15 | 9 | 6 | 3 |
| omega_j^2 | 4 | 14 | 30 | 12 | 6 | 2 |
| I_j^2-8omega_j^2+16 | 9 | 4 | 1 | 1 | 4 | 9 |

The eleven inherited base normals are loaded from the small pinned H8
certificate, lifted from (X;e0,...,e7) to (X,0;e0,...,e7). All 121 saved
Gram entries are compared. The known actual-wall list adds

    r=(3phi,2;0), c=(0,-1;0),
    m=(phi-1,-phi;phi^2 omega7),
    n=(phi/2)p1-c/2+omega7+phi^3 omega0.

The worker checks original-coordinate integrality, norm 2 (r has norm 2b),
integral reflection covectors and orientation at p1 for all fifteen walls.
Their ACTUAL wall status is a geometric premise from 33f6, not inferred
from these algebraic checks.

The upper uses U*=1206+1953phi, Z*=19925+32240phi,
Y*=2490+4030phi and h*=(U*/2)p1-c/2+omega0+Y*Q. Its original
binary coordinates, norm, reflection, Z*^2-4DU*^2=5 and collision
tau*=DU*/Z* are explicitly checked, including D-4tau*^2>0. By b4de,
this integral mirror forces an actual finite next cut at or before tau*;
h* need not itself be a simple wall or satisfy known-wall acuteness.

## 2. Exactly nine models and 38,784 raw trials

The all-height exclusions of independent norm2/norm2b and c-only5 remain
mathematical premises from b4de, including its original-module congruence
argument. Every remaining first simple wall is one of these nine models:

| Models | Normal h | Square equation and reconstruction |
|---|---|---|
| c-only3 | (U/2)p1-c/2+YQ | z^2=DU^2+3; Y=z/4 |
| A8 at 0 or 7 | Xp1+omega_j+YQ | z^2=9+16DX^2; Y=(z-I_j)/8 |
| D8 at 2 or 6 | same | z^2=4+16DX^2; Y=(z-I_j)/8 |
| E8 at 3 or 5 | same | z^2=1+16DX^2; Y=(z-I_j)/8 |
| A9 bridge at 0 or 7 | (U/2)p1-c/2+omega_j+YQ | z^2=5+4DU^2; Y=(z-I_j)/8 |

Controls U,X and Y belong to O; the control is strictly positive.
Y>=0, and Y>0 in the c-only row. The c-only and bridge rows additionally
require (1-aU)/2 in O, equivalently U=phi mod 2O. Original coordinates,
norm exactly 2, integral reflection covector and all eight retained
products are then checked on the reconstructed vector, not just inferred
from a square equation. An integral norm2 vector in this even module is
primitive: a nonunit common divisor d would give 1=d^2 k in O.

Closed collision bounds are, respectively,

    U^2 <= (12/5)U*^2,    X^2 <= (R/20)U*^2,    U^2 <= U*^2.

The other embedding gives (-D')U'^2<=3,
16(-D')X'^2<=R, or 4(-D')U'^2<=5. These are exact non-strict
inequalities, so equality and all collision ties are retained.
The proven conservative strict strips and integer loops are:

| Models | physical cap | conjugate radius r | 0<=B<= | A offsets |
|---|---:|---:|---:|---|
| c-only3 | 6800 | 161/100 | 3050 | -1,0,1,2 |
| two A8 | 3000 | 7/10 | 1350 | 0,1 |
| two D8 | 2000 | 93/200 | 900 | 0,1 |
| two E8 | 1000 | 233/1000 | 450 | 0,1 |
| two A9 | 4400 | 26/25 | 1970 | -1,0,1,2 |

Here a control is A+Bphi and the offsets are relative to
n=floor(B(phi-1))=(isqrt(5B^2)-B)//2. This is an exact integer floor:
flooring sqrt(5B^2) before subtracting integer B and dividing by 2
preserves the final floor. Since r<sqrt5, B=(t-t')/sqrt5>-1 gives
B>=0. The bounds in b4de include every B and A; no saved survivor list,
approximate root order or unproved height cutoff is used. The loops visit

    4*3051 + 4*1351 + 4*901 + 4*451 + 8*1971 = 38784

raw trials, before any filter. The worker requires this exact total at
completion. All nine model counters are retained separately.

## 3. Complete O-square test and exact order

Every field element is stored as two Fractions in the basis (1,phi).
Multiplication uses phi^2=phi+1; conjugation sends (a,b) to (a+b,-b).
The physical sign of a+bphi is the sign of (2a+b)+b sqrt5, found by
comparing the rational squares (2a+b)^2 and 5b^2 when signs differ.
There are no floats, square-root approximations or symbolic imports.

For an integral radicand w, suppose v=x+yphi in O and v^2=w.
N(w)=N(v)^2 must be a nonnegative integer square. For EACH possible
n=N(v)=+sqrt(N(w)) or -sqrt(N(w)), the identities

    Tr(v)^2=Tr(w)+2n,    y^2=(Tr(w)-2n)/5

give necessary integer-square tests. The algorithm takes both signs of
each integer square root, requires Tr(v)-y even, reconstructs
x=(Tr(v)-y)/2, and requires exact resquaring and positive physical sign.
Conversely any returned element is in O, positive and squares to w.
Every possible integral root appears in these sign choices. Zero is
handled separately; trace-zero and rational roots are included. The
runtime fixtures include sqrt(5)=-1+2phi, phi^2, zero, and the nonsquare
3 whose rational norm is nevertheless a square. No fixture has yet run.

Only a positive root is needed: for a c-only row it is 4Y>0, and in
the other rows z=8Y+I_j>0. Divisibility and original binary parity are
checked after the square test. Exact resquaring prevents a rational-norm
false positive from entering the candidate list.

## 4. Known-wall filters, minimum and germ certificate

For each reconstructed mirror h the worker checks every one of the
fifteen known actual walls w. It rejects (h,w)>0. A vector equal to w
is recognized as the same mirror, not tested against itself. A zero
product is orthogonal. Otherwise set B=4(h,w)^2/(h^2 w^2). Both
reflection covectors are integral, so B is an algebraic integer; this is
also checked. B=4 is rejected by anisotropy, B>4 is dotted, and 0<B<4
must be exactly 1,2,phi^2,3,2+phi, corresponding to native labels
3,4,5,6,10. In particular the nonnative value 2-phi is NOT admitted
merely because its square root lies in K. This enforces the released
-e1 constraint as well as the retained and other known walls.

Each integral mirror is streamed with its original vector, reflection
covector, control, positive square root, Y, exact tau and complete
known-wall result (or first rejecting wall/reason). Rejections do not
disappear from the certificate. The collision is recomputed as
tau=-(p1,h)/(Q,h), required positive and <=tau*. Among admitted mirrors
the worker compares exact field differences, resetting only on a strict
improvement and retaining EVERY equal minimum. Within one model tau is
strictly increasing in its positive control, by the three rational
fractions in b4de Section 5. Thus at most one minimum occurs per model,
so retaining at most nine ties is proved, not an arbitrary truncation.

An actual first simple mirror is in the complete domain and passes all
known-wall tests. Any enumerated integral mirror oriented negative at p1
cannot cut the actual chamber edge before its actual endpoint. Therefore
the admitted minimum is the actual endpoint, conditional on the pinned
geometric premises. This does not assert that every admitted mirror is
a simple wall. A missing admitted minimum fails the checker, preserving
its partial evidence instead of silently declaring an empty interval.

To recognize the whole finite germ, for EACH minimum tie h the worker
forms the nine-root basis (eight retained roots,h). It checks strictly
positive exact Schur pivots, giving positive definite Gram, then solves
for EVERY other tied normal in that basis. It requires exact vector
reconstruction and records whether all nine coefficients are nonnegative.
No candidate is chosen because of its position in a loop.

The actual new simple normal s is among the ties. Every tied arrangement
root is a positive root of the actual finite tangent system (orientation
is negative just before the endpoint), hence is in cone(retained,s).
If all ties lie in cone(retained,h), then s also lies in that cone.
Substituting the two nonnegative representations forces the retained
coefficients to vanish and h to be a positive multiple of s. Their common
norm2 forces h=s. This proves the cone recognition rule. Only when the
accepted proposals are nonempty and all have the SAME exact normal does
the output claim a unique actual simple germ at the pinned geometric
boundary. Otherwise it reports a complete necessary interval with germ
recognition unresolved, retaining all proposals and ties. No full-chamber
metric isometry, native period or chamber stabilizer is constructed.

## 5. Pins, execution protocol, failure and scope

Eight fixed input hashes appear in the worker: the two mathematical
sources, the small H8 data and its transcript, plus the bounded runner,
RSS monitor, signal-safe runner and resource-threshold module. All eight
were streamed and matched during source preparation. No local Python
module is imported: the only decoded inherited JSON is the 7,648-byte H8
file. The other inputs are streamed for SHA binding; their mathematical
status is not replaced by that hash test. The worker's own SHA and this
specification's SHA are mandatory CLI inputs, checked with all fixed pins
before and after the finite calculation. Root must obtain the spec SHA
after freezing this file; it is intentionally not embedded in itself.

Prepared root-only entrypoint (NOT EXECUTED; SPEC_SHA is the frozen SHA):

```text
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock \
  --ledger results/prime79_a7_c_interval_execution_20261021.json -- \
  scripts/check_compact_h9_prime79_a7_c_interval_20261021.py --verify \
  --expected-source-sha 59926019e0c228bb8364038660bed8e37b601bc771f4894fe4e024952390eda6 \
  --spec-sha SPEC_SHA --out results/prime79_a7_c_interval_20261021.jsonl
```

The unchanged runner admits only above 40GiB disk and 6GiB available RAM,
aborts below 8/3, retains the shared exclusive lock, and monitors one
process tree for <59s, <1GiB RSS and one mathematical thread. The child
requires that actual wrapper parent/lock, -B, a 45-second alarm and
45-second CPU limit. Finite loop size is proved; runtime/RSS admission or
success is NOT predicted or certified by this source note.

The output must be a NEW JSONL directly under results, created exclusively.
It starts INCOMPLETE, flushes each integral-mirror record, and ends with
a complete result only after the final pin/parent checks. Failures append
an incomplete failure record when possible; a signal/watchdog interruption
leaves an incomplete prefix. No earlier output is overwritten or removed.
Stdout is only a compact summary with the certificate byte length and SHA;
root's separate terminal guarded ledger must bind the actual command,
source, resources and stdout. No certificate or PASS exists yet.

The finite domain is complete for this ONE literal A7+c edge. Even a
successful unique-germ result would not establish finite volume of the
entire chamber, an infinite-order chamber period, nonreflectivity, an H9
example, a family/form-class exclusion, a new dimension bound or Lean
verification. Those conclusions are absent from the checker.
