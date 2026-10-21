# Independent hostile audit: the rational isotropic period gate

**Verdict (2026-08-10): PASS, in exactly the stated smooth complex
genus-one scope.**  This audit was reconstructed from the restriction
sequence and the two blowup constructions; it does not use the conclusions
of the earlier manual audit as premises.  It certifies no Coxeter chamber.

The audited package consists of
[`rational_isotropic_weak_log_period_gate.md`](rational_isotropic_weak_log_period_gate.md),
its exact producer, deterministic archive, and regression suite.  The target
claim is deliberately narrow: for a smooth connected genus-one divisor of
square zero on a rational surface, semiampleness contains a normal-bundle
period datum that the marked intersection lattice does not retain.

## 1. Reconstruction of the cohomology argument

Let \(F\subset S\) satisfy the hypotheses of the report and put
\(N=\mathcal O_F(F)\).  Since \(F\) is an effective Cartier divisor,
tensoring its ideal sequence by \(\mathcal O_S(rF)\) gives, for \(r\geq1\),

\[
 0\longrightarrow\mathcal O_S((r-1)F)
 \longrightarrow\mathcal O_S(rF)
 \longrightarrow\mathcal O_F(rF)=N^r\longrightarrow0.       \tag{1.1}
\]

There is no exponent shift in the last term.  Riemann--Roch gives

\[
 \chi(\mathcal O_S(rF))
 =1+\frac{r^2F^2-rK_S\cdot F}{2}=1.                         \tag{1.2}
\]

Also \(H^2(S,rF)=H^0(S,K_S-rF)^\vee=0\): an effective
\(K_S-rF\), together with the already effective \(rF\), would produce an
effective canonical divisor, whereas a rational surface has \(p_g=0\).
This also works at \(r=0\).

On a smooth connected genus-one curve, a degree-zero line bundle has a
section exactly when it is trivial.  Consequently the target in global
sections in (1.1) has dimension one precisely when \(N^r\simeq\mathcal O_F\),
and otherwise has dimension zero.

If \(N\) has exact order \(m\), then \(h^0(rF)=1\) for \(0\leq r<m\).
Equations (1.1)--(1.2) and \(H^2=0\) give
\(H^1((m-1)F)=0\), hence the restriction map at \(r=m\) is onto.  A lift
\(u\in H^0(S,mF)\) of a nonzero constant on \(F\) has a zero divisor disjoint
from \(F\).  Together with the canonical section vanishing on \(mF\), it
globally generates \(\mathcal O_S(mF)\).  Thus 

\[
                     |mF|\text{ is a base-point-free pencil}. \tag{1.3}
\]

For every \(q\geq1\), \(u^q\in H^0(S,qmF)\) restricts nontrivially to
\(F\).  This independently forces the restriction-map rank to be one at
every multiple of \(m\); between multiples the target has no sections.
Therefore

\[
 h^0(rF)=\left\lfloor\frac r m\right\rfloor+1,
 \qquad
 h^1(rF)=\left\lfloor\frac r m\right\rfloor.                \tag{1.4}
\]

This step does not assume connected fibres, a pre-existing fibration, or a
trivial Stein factorization.  If \(N\) is non-torsion, all positive targets
in (1.1) have no sections, so \(h^0(rF)=1\) and \(h^1(rF)=0\) for all
\(r\geq0\).  Conversely, if a positive-dimensional \(|rF|\) exists, a
section outside the one-dimensional image of \(H^0((r-1)F)\) restricts
nontrivially, forcing \(N^r\simeq\mathcal O_F\).  This recovers all three
equivalent conditions and closes the edge cases \(m=1\) and \(r=0\).

## 2. Reconstruction of the two nine-point blowups

For nine distinct points \(p_i\) on a smooth plane cubic \(C\), the blowup
has marked lattice \(I_{1,9}\),

\[
 K_S=-3H+\sum_iE_i,
 \qquad
 F=3H-\sum_iE_i=-K_S,
 \qquad F^2=K_S\cdot F=0.                                  \tag{2.1}
\]

The vector \(F\) is primitive.  Nefness is not being inferred from the
lattice: it is proved geometrically.  One has \(F\cdot E_i=1\), and for an
irreducible plane curve \(Q\neq C\), local Bézout on the smooth curve \(C\)
gives

\[
 3\deg Q=C\cdot Q
 \geq\sum_i I_{p_i}(C,Q)
 \geq\sum_i\operatorname {mult}_{p_i}Q,
\]

which is exactly \(F\cdot\widetilde Q\geq0\).  Exceptional curves,
\(F\) itself, and strict transforms of such \(Q\) exhaust the irreducible
curves on the blowup.

Restriction to the strict transform identifies

\[
 N\simeq\mathcal O_C(3)(-p_1-\cdots-p_9).                  \tag{2.2}
\]

If the nine points are the transverse intersection of \(C\) with a second
cubic, (2.2) is trivial and blowing up the reduced base locus resolves the
cubic pencil.  Thus \(|F|\) is base-point-free.  For the opposite example,
fix eight distinct points and vary the ninth.  The map

\[
 p_9\longmapsto\mathcal O_C(3)(-p_1-\cdots-p_8-p_9)
\]

is a translate of the negative Abel--Jacobi isomorphism
\(C\to\operatorname {Pic}^0(C)\).  The torsion subgroup of a complex
elliptic curve is countable, while \(C(\mathbf C)\) is uncountable; excluding
those preimages and the finitely many repeated points leaves a valid
non-torsion choice.  The theorem then makes every \(|rF|\) fixed.

Thus the examples really have the same marked 
\((\operatorname {NS},K,F)\) and both have \(F\) nef, while their
semiampleness answers differ.  The restriction map to the embedded cubic,
not the abstract intersection form, carries the missing period.  Nothing in
this argument proves that the two complete nef cones are marked-isomorphic,
and the audited report does not claim that.

## 3. Attacks on the certificate and on scope

The producer exactly checks the diagonal \(I_{1,9}\) form, the vectors in
(2.1), primitivity, exceptional pairings, adjunction, and the
Riemann--Roch polynomial.  Its torsion tables are a finite replay of the
restriction-map ranks proved above; they are not computational evidence for
the elliptic-line-bundle theorem or for the moduli choice of points.  The
archive says so explicitly.

The following possible scope upgrades were checked and are all expressly
excluded:

* arbitrary isotropic weak-log classes without a smooth genus-one member;
* equality of complete nef cones;
* realization of a prescribed Coxeter polytope;
* return of all physical walls or of the complete nef chamber;
* finite-volume, reflection-saturation, or reflection-maximality claims; and
* exhaustivity for raw maximal arithmetic polytopes.

Accordingly, the result supplies a necessary extra datum for the isotropic
compiler lane; it is not itself a polytope-to-surface existence theorem.
No Nikulin, Vinberg, or ABSW classification is imported.

## 4. Frozen artifacts

The independent tests recompute the lattice and recurrence without using the
producer for those core assertions, replay the archive byte-for-byte, enforce
the scope exclusions, and verify all pins below.

| artifact | SHA-256 |
|---|---|
| theorem report | `10335b97010d4e3197e5f84060a3f38dfd3775ed849097eba00cb86ea86841e5` |
| exact certifier | `77acb4b484e2b8a1f6b39b414d15ba39616c4d3ef221122121f93369fc49a6e5` |
| deterministic archive | `518eccdf3798c197d031d512a064cc4ec47638788af1010be8787cfedc0ce254` |
| original regression tests | `e660eb820b36cb16b6494f33fba2c848da0025a7e90eeabe3092fc24d2f88bcb` |
| independent hostile tests | `1f838231b9038a2afafe49d47a3bd60e827c0f7935b9907ef2877fa045c3aedb` |
