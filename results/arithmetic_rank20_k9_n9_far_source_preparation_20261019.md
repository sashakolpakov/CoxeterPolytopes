# Source preparation: one K9, N=9, far-endpoint3 arithmetic batch

19 October 2026. **Source prepared only. No import, AST validation,
calibration, enumeration or scientific worker has been executed.**
There is no result JSON and no arithmetic-absence claim.

The two new sources are:

* `scripts/exact_totally_real_field_20261019.py`, SHA256
  `fbd896a6160f232e0f6b4e8a715c990b05498d74be5c846e91aa699c95d598d2`;
* `scripts/fresh_rank20_k9_n9_endpoint3_far_20261019.py`, SHA256
  `2310f437a8210d34891a3e919965ffa35b5cde1a5102a405cd6a64826706ba5d`.

The batch pins the helper, written height theorem and enumeration
design, checks their streaming hashes before work, and rechecks all
sources before publication. Its own source hash is retained too and
must be pinned by the root's external scheduled-job record.

## Mathematical and implementation boundary

The pure helper has no I/O, driver, imported scientific dependency or
import-time computation. It supplies exact Fraction polynomial
arithmetic, positive-content Sturm chains, certified real-root
refinement, exact field inverses, Newton traces, rational trace duals,
and rational-rank generated-degree tests. It assumes the caller's
declared polynomial really is minimal. It does not claim to prove
the caller's integral basis.

The one batch is fixed to K9=Q(theta), theta=2cos(2pi/9), N=9,
r=11, far-endpoint3. It accepts no alternate field, rank or role.
The monic cubic X^3-3X+1 is irreducible by the rational-root test:
its only possible rational roots are +/-1 and neither is a root.
The three embeddings are certified by separate Sturm intervals,
ordered increasingly, with the physical root largest.

The external integral-ring input is explicitly

    O_K9=Z[theta].

This uses the maximal-real-cyclotomic integer-ring theorem. Real
root counts and cubic irreducibility alone do not supply that input.
All geometric implications, corank and cyclic-product integrality
remain the typed written inputs of the pinned height source.

For A=2x, the proved physical bound is 2<=A<24. Every nonphysical
conjugate satisfies |9A_i+1|<sqrt10. Since 10<361/36, the closed
rational interval [-25/54,13/54] contains each such A_i. The batch
constructs all integer trace values and complete trace-dual
coefficient boxes from these outer bounds. It recovers c0 from the
exact trace equation and traverses every integer tail in every box.
It does not choose a heuristic coefficient or height cutoff.

For every remaining vector it checks, exactly,

    Z=9A+1,
    C=(2Z^2-110)/9=18A^2+4A-12,
    delta=40/C in O_K9,
    z=4+delta,
    delta_phys>0, 0<z_i<8/11 for both i nonphysical,
    qc=4/9+Z^2/18, qc_phys>1, qc_i<1,
    (Z^2/180-11/36)*delta=1,
    Q(z)=K9.

The general helper is reusable, but this source does not contain
an eight-field search or a multi-role producer. In particular there
is no inaccurate requirement that A generate K in a future e4 row.

## Future scheduling and result interpretation

Root review must precede execution. One future invocation must run
under the shared serial lock and repository resource gates, with an
external 59-second/1GiB limit. The script additionally has a
55-second internal deadline, a complete-box visit cap and a bounded
512KiB certificate. It never spawns another process. A preexisting
output is never overwritten.

The single script argument is a new JSON output path directly under
`results`. The scientific payload source is the fixed batch above;
there is no proposed unrestricted shell or build command.

Complete empty and complete nonempty necessary-condition sets have
different explicit statuses. Every retained passing vector is marked
`UNRESOLVED_NECESSARY_CONDITION_HIT`, with its exact A,C,delta,z,qc,R
coordinates. Such a hit does not prove an actual face or a polytope.
Timeouts, excessive boxes, unresolved signs, source changes, internal
identities that fail, or output overflow give incomplete status and
block an absence claim. Process exit zero means only that this
necessary-condition traversal finished, even when unresolved hits
remain; it must never be interpreted as a universal geometric PASS.

Fixed exact Sturm/inverse/trace/sign/degree calibrations are included
for a future scheduled run. They have not been called.

An independent checker that reconstructs and traverses the complete
trace boxes is **not yet written or run**. Rechecking only saved hits
would be insufficient for absence. The producer records this gap in
its certificate. Any future sharing of the pure helper must be named
as a shared arithmetic dependency, rather than advertised as two
independent implementations.

Even a future independently replayed empty set would settle only
this field/N/role under the named inputs. It would not cover c
meeting the spectator, all ten remaining consumer families, their
triangle composition, or a simultaneous rank20 allowance of35.
