# Hostile audit: rational isotropic weak-log period gate

**Verdict (2026-08-10): PASS AFTER REPAIR, within the stated smooth complex
surface scope.**  This is a manual internal proof audit.  It is not represented
as an independent external referee report.

The audited theorem is
[`rational_isotropic_weak_log_period_gate.md`](rational_isotropic_weak_log_period_gate.md).
The audit tried to break the cohomological equivalence, the exact section-count
formula, the two nine-point blowup examples, and the claimed compiler
consequence.  It did not audit or certify any Coxeter chamber.

## 1. Hypothesis audit

The final theorem explicitly works over \(\mathbf C\).  Thus all uses of an
elliptic \(\operatorname {Pic}^0\), the countability of the torsion subgroup,
and the existence of a point outside a countable set are legitimate.  The
surface is smooth, projective, and rational; hence

\[
 q(S)=p_g(S)=0,\qquad \chi(\mathcal O_S)=1,
 \qquad \operatorname {Pic}(S)=\operatorname {NS}(S).
\]

The divisor \(F\) is explicitly effective, Cartier, smooth, connected,
genus one, and primitive.  It is also assumed nef with
\(F^2=K_S\cdot F=0\).  Primitivity and nefness are the compiler hypotheses;
the proof of the period theorem does not covertly use either one.

**Attack:** was degree zero of the normal character merely asserted?

**Result:** no.  The report records
\(\deg\mathcal O_F(F)=F^2=0\), so
\(N=\mathcal O_F(F)\) genuinely lies in \(\operatorname {Pic}^0(F)\).

## 2. Cohomology audit

For \(r\geq1\), the Cartier restriction sequence is exactly

\[
0\to\mathcal O_S((r-1)F)\to\mathcal O_S(rF)\to N^r\to0.
\]

Riemann--Roch gives \(\chi(\mathcal O_S(rF))=1\).  Serre duality reduces
\(H^2(rF)\) to \(H^0(K_S-rF)\).  An effective representative of
\(K_S-rF\), added to the effective divisor \(rF\), would be an effective
canonical divisor, contradicting \(p_g(S)=0\).  Thus the claimed \(H^2\)
vanishing is valid for every \(r\geq0\).

**Attack:** could a nontrivial degree-zero elliptic line bundle contribute a
section and corrupt the recursion?

**Result:** no.  An effective degree-zero divisor on a smooth connected
elliptic curve is zero, so a degree-zero line bundle has a nonzero section
exactly when it is trivial.

Suppose \(N\) has exact order \(m\).  Until \(r=m\), the restriction term has
no section, whence \(h^0(rF)=1\).  Since \(\chi=1\) and \(H^2=0\), this gives
\(H^1((m-1)F)=0\).  At \(m\), restriction is onto a one-dimensional target,
so \(h^0(mF)=2\).

**Attack:** do two sections really prove base-point-freeness?

**Result:** yes.  A lift \(u\) of \(1\in H^0(F,N^m)\) is nowhere zero on
\(F\), so its zero divisor \(G\sim mF\) is disjoint from \(F\).  The section
with zero divisor \(mF\) and \(u\) have disjoint zero sets.  They span
\(H^0(mF)\), so the complete system \(|mF|\) is a base-point-free pencil.

**Attack:** did the formula
\(h^0(rF)=\lfloor r/m\rfloor+1\) silently assume connected fibres or a
trivial Stein factorization?

**Result:** no.  At every \(r=qm\), the explicit section \(u^q\) restricts
nontrivially to \(F\), so the restriction map has rank one.  At every other
positive \(r\), its target has no section.  Hence the dimension increments
exactly at the multiples of \(m\).  This direct argument proves the floor
formula before any connected-fibre statement.  Riemann--Roch and \(H^2=0\)
then give \(h^1(rF)=\lfloor r/m\rfloor\).

If \(N\) is non-torsion, every positive target has zero sections, giving
\(h^0(rF)=1\), and then \(h^1(rF)=0\).  Conversely, at the least \(r\) for
which \(|rF|\) moves, the new section has nonzero restriction to \(F\), so
\(N^r\) is trivial.  If a multiple is base-point-free, its degree-zero
restriction is globally generated and hence trivial.  All directions of the
equivalence are therefore covered.

## 3. Blowup audit

For nine distinct ordinary points on a smooth cubic, the blowup has the
marked lattice \(I_{1,9}\), and the strict transform satisfies the equality

\[
[F]=3H-\sum_{i=1}^9E_i=[-K_S]
       \quad\text{in }\operatorname {Pic}(S)=\operatorname {NS}(S).
\]

This is a linear-equivalence statement, not a claim inferred only from
intersection numbers: it is an equality of divisor classes.  The vector is
primitive and has
\(F^2=K_S\cdot F=0\).

**Attack:** does the Bézout argument miss a curve on the blowup?

**Result:** no.  An irreducible curve is either exceptional, the strict
transform \(F\), or the strict transform of an irreducible plane curve
\(Q\ne C\).  These three cases give respectively \(F\cdot E_i=1\),
\(F^2=0\), and

\[
F\cdot\widetilde Q
=3\deg Q-\sum_i\operatorname {mult}_{p_i}Q\geq0
\]

by local intersection multiplicities with the smooth cubic.  Thus nefness is
proved for both choices of points.

The normal bundle formula

\[
\mathcal O_F(F)\simeq\mathcal O_C(3)(-p_1-\cdots-p_9)
\]

is the usual strict-transform normal bundle and has degree zero.  If the
points are a transverse complete intersection of \(C\) with another cubic,
one blowup at each point resolves the pencil and the character is trivial.
For the opposite example, after fixing eight points, the ninth-point map to
\(\operatorname {Pic}^0(C)\) is a translate of the Abel--Jacobi isomorphism.
The complex torsion subgroup is countable, so a distinct ninth point with
non-torsion character exists.

The examples therefore have the same marked \((\operatorname {NS},K,F)\)
and the same answer to the question “is \(F\) nef?”, but opposite answers to
semiampleness.  The report correctly does **not** assert that their complete
nef cones are marked-isomorphic.

## 4. Computation and scope audit

The deterministic producer checks the \(I_{1,9}\) Gram matrix, determinant,
signature, primitive vector, canonical and exceptional pairings, adjunction,
Riemann--Roch values, exact torsion recurrences, and the non-torsion
recurrence.  The archive records that this finite replay does not replace the
geometric proof.

The following possible upgrades were searched for and are all explicitly
denied by both report and archive:

* coverage of every isotropic weak-log class;
* equality of the two complete nef cones;
* realization of a prescribed Coxeter polytope;
* return of a complete nef chamber or all physical facets;
* finite hyperbolic volume;
* reflection saturation or maximality; and
* exhaustivity for raw arithmetic polytopes.

No reflective-lattice classification is imported.  The valid conclusion is
exactly the local compiler gate: in the square-zero smooth genus-one lane,
the Néron--Severi packet and nefness do not determine semiampleness; the
normal period, a base-point-free certificate, or an equivalent morphism
certificate must be supplied.

## 5. Frozen audited artifacts

| artifact | SHA-256 |
|---|---|
| theorem report | `10335b97010d4e3197e5f84060a3f38dfd3775ed849097eba00cb86ea86841e5` |
| exact certifier | `77acb4b484e2b8a1f6b39b414d15ba39616c4d3ef221122121f93369fc49a6e5` |
| deterministic archive | `518eccdf3798c197d031d512a064cc4ec47638788af1010be8787cfedc0ce254` |
| regression tests | `e660eb820b36cb16b6494f33fba2c848da0025a7e90eeabe3092fc24d2f88bcb` |

The archive replay and the dedicated regression suite must both pass before
these pins are used.
