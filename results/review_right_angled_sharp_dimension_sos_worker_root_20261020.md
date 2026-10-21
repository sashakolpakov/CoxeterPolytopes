# Root source review: exact cubic tensor-square checker

20 October 2026. Full source review PASS; no execution is asserted here.

Reviewed all154 lines of
`scripts/check_right_angled_sharp_dimension_sos_20261020.py`, SHA256
`bd0ed1c4c1b7e18a219a2a8aca2b178dbe1d559d7e9ec90ad609150aedf8b2c3`.
The worker pins the fully reviewed dimension-eight proof9470fa1f...babdbe.

Its arithmetic is a literal polynomial ring over Fraction in twelve
independent variables u1,...,u6,v1,...,v6. Equalities compare complete
monomial/coefficient dictionaries. It constructs all36 entries of the
rank-two tensor and216 entries of the rank-three tensor for BOTH
variable sets, using the exact Kronecker-delta formula. This is not
a sampled substitution or a numerical positive-semidefinite test.

I independently checked the two generic contraction targets:

    A(u).A(v)=t²-(su+sv)/6+1/6,
    T(u).T(v)=t³-(3/8)t(su+sv)+(3/8)t,

and the diagonal/trace ideal remainders. In particular the rank-three
trace is u_a(su-1), and its squared norm minus5/8 is
(su-1)(su²+su/4+5/8). The tensor indexing36a+6b+c agrees with
the construction's nested loops.

The complete kernel minus the cubic polynomial is exactly

    (2-su-sv)((3/8)t+1/8).

The worker proves this as an unrestricted polynomial identity, then
also proves its explicit membership in the ideal(su-1,sv-1).
Thus unit-vector specialization is honest and universal. The same
method verifies the optional normalized harmonic-kernel formulas.

The factored cubic, all positive square coefficients, constant and
diagonal ratio27, and negative squared norm-2 for four equality-case
vectors are exact rational checks. The worker does not enumerate
spherical codes or prove all-wall projection, face inheritance, real
norm positivity, or Dufour's facet minimum. Those are written and
published mathematical inputs with their own reviews.

The CLI is exclusive write/verify/stdout, refuses certificate overwrite,
uses a45-second alarm and bounds the retained certificate by256KiB.
It streams before/after hashes of itself and the proof. The standard
root-owned bounded wrapper must retain the shared lock, one-thread
environment, resource gates,59-second watchdog and1GiB process-tree
cap. No imported optimizer, scientific runtime or detached worker is
present. A later actual generation and replay must have separate
execution ledgers.
