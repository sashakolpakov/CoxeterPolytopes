# Golden31 unit-bisector proof and fixed checker: outside review

20 October 2026. Reviewer: batch_guard. Complete paper/source review,
read-only apart from this record. No worker import/body, generation,
replay, scientific computation, Lean, or TeX was executed by the reviewer.

**Verdict: PASS. No gap found in the literal-host nonreflectivity proof
or in the fixed checker's stated verification interface.** Execution is
not claimed by this review and remains root-owned.

## Exact bindings

* Proof: `results/compact_h9_golden31_unit_bisector_preparation_20261020.md`,
  SHA256 `7259c05e994fd48d8bb5f70c172e98f9c0120bbfa13fc58f3a6b08174e5f9ea1`.
* Worker: `scripts/check_compact_h9_golden31_unit_bisector_20261020.py`,
  SHA256 `5b5af55734b6799fb3acce349bdad2e1fe08b64ebe7abcf2d8a3264f2f647245`.
* Fixed contract: `results/compact_h9_golden31_fixed_checker_preparation_20261020.md`,
  SHA256 `c199ff6aa72cb74ae21b1989d7c4258684b318ce707abeb3d699c68c96fd55fe`.
* Prior root review: `results/review_compact_h9_golden31_unit_bisector_root_20261020.md`,
  SHA256 `4cf76e054d1b8461a46d86fce9421a0b7139b741f6f45938761a6a15a71a80f9`.

The pinned golden binary audit, corrected L59 isolation argument and native
Bugaenko fixture/transcript were read at the interfaces used. This reviewer
previously checked those inherited interfaces, but did not author the new
golden31 proof. Its bisector germ, involution and trace were rederived here.

## Complete mathematical audit

The determinant and conjugate signatures are correct. Evenness,
primitivity, the golden PID and positive-unit square classification give
the complete necessary norm pool2,2b. Prime31 isolation uses the entire
finite-angle divisibility list, including arbitrary arrangement angles
2pi/5 and3pi/10. The long--long radical calculation also excludes
orthogonality. The nonexistence of K-null vectors excludes tangency.
Thus the base r facet is the FULL calibrated Bugaenko arrangement;
it is not inferred from an acute finite seed or a discriminant alone.

The binary basis(r,c) has unit determinant. I checked the physical
opposite-sign inequalities and the conjugate Schur bound for both norm
classes. In the long class b divides mu because b is odd. Both cases
give the same strict norm bound4/[31(4-b')]<1. Finite extension of the
retained E8 forbids a nonbinary first wall. This is an all-height
argument, so c and the full E8+A1(c) endpoint germ are actual.

The timelike vector t=c+g has square-2/phi, is future, and is feasible
by actual all-wall acuteness. A zero wall at t must be orthogonal to
both c,g. Solving those two equations integrally gives exactly the
displayed (2u,-phi^2*u;z) and z.omega7=b*u. For short roots, BOTH
embedding Cauchy bounds give N(u)^2<=16/341<1. For long roots,
unimodularity gives z=b*z0; then N(u)^2<=16/11<2 forces a unit,
while the conjugate bound requires k>0 and the physical bound excludes
every k>=1. The u0 long case is impossible by the nonunit b.
The remaining short roots are precisely ordinary E7: the separate
even-E8 coefficient argument excludes hidden O-valued norm2 roots.
Consequently the actual zero germ is EXACTLY E7, not merely a rank7
subset of an unknown larger germ.

W and the point inversion commute, are integral isometries, and fix t.
Their product H fixes each actual E7 normal with its inward sign.
Because t is an interior hyperbolic point with exactly that zero
arrangement, this identifies the full local chamber germ and hence
H(C)=C. No assumption that an arbitrary integral isometry stabilizes C
is used. The inherited S stabilizes the full base-facet germ separately.

In the native basis I independently recovered tr(S)=2,
tr(S P8)=-8-16phi and (omega7,S omega7)=-10-18phi. Expanding the
given St yields binary contribution-46-77phi and E8 contribution
26+41phi, so (t,St)=-20-36phi. The resulting trace is28+42phi>10.
The finite-order trace bound, free action of the full reflection group
on chambers and distinct-coset argument therefore give infinite index
in the integral arithmetic lattice. The explicitly imported arithmetic
lattice theorem then gives infinite chamber volume.

## Fixed checker audit and scope

I read the entire worker and contract. The Fraction-pair sign, norm,
conjugation and division routines are exact; the matrix conventions
use columns for images and the Gram adjoint for inverse verification.
The native E8 ordering, ten S columns, W/H rank-one formulas,
involution/Gram/integrality/future checks and trace contractions agree
with the independent paper derivation. The finite-angle list and both
Cauchy seals match the complete written proof. The program correctly
leaves the classical norm-pool, root-arrangement geometry, actual-germ
and covolume implications as written premises; fixed scalar checks do
not masquerade as enumeration of all mirrors.

There is no height, Pell or cone search. Pinned inputs and the worker
are streamed and rechecked before returning. Root-owner enforcement,
the45-second child alarm,128KiB certificate limit, exclusive creation
and fresh exact replay interface are consistent with the outer bounded
runner. No change is requested.

This excludes ONLY the literal golden31 host. It does not exclude an
entire similarity/commensurability class, the different silver31 host,
or the remaining infinite binary family. The Lucas/unit condition
correctly isolates exponent2 in the timelike range; it is not a
uniform-family extension. No new example or dimension ceiling follows.
