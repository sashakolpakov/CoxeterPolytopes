# Hostile audit of weak-incidence reflection-maximal finiteness

## Verdict

**PASS after repairing TeX damage and narrowing one scope word.**  The
carrier-relative theorem in
[rational_wrong_log_reflection_maximal_finiteness.md](rational_wrong_log_reflection_maximal_finiteness.md)
is correct on its declared input.  A literal complete finite-sided,
finite-volume rational full-nef chamber satisfying the complete
\((-1),(-2),(-4)\) Mori alphabet and the weak incidences has

\[
                         \rho(S)\leq26,\qquad n\leq25.
\]

If its facet group is maximal in any of the three stated reflection senses,
then it equals the full integral reflection group of the adapted host.  Each
host has at most one unmarked complete chamber orbit, only finitely many
hosts occur, and every resulting polytope has minimal Vinberg field
\(\mathbf Q\).

The audit found no theorem defect.  It did find two concrete presentation
defects in the owner report:

1. most inline mathematics after Section 2 had lost its TeX delimiters;
2. the status said “classification-free”, although the host count uses the
   elementary classification of rational Lorentzian Picard lattices.

The repaired status says “reflective-classification-free”, and the theorem
now names all three valid maximality modes literally.  No hypothesis or
conclusion was narrowed.

The audit does **not** import the producer's functions.  Its independent
certifier, archive, and tests are

* [audit script](../scripts/audit_rational_wrong_log_maximal_finiteness_hostile.py);
* [audit archive](../data/enumerations/rational_wrong_log_maximal_finiteness_hostile_audit.json);
* [hostile tests](../tests/test_rational_wrong_log_maximal_finiteness_hostile_audit.py).

## 1. Audited owner artifacts

| artifact | SHA-256 |
|---|---|
| results/rational_wrong_log_reflection_maximal_finiteness.md | a682712bf212aa79c15547aeba0f59c279612d07bd428df7c5fd9d5339cf2768 |
| scripts/certify_rational_wrong_log_maximal_finiteness.py | 22d551ca0ead5d362c6ca194a8d60aa3257958f806498da38bae6cfcac430f3d |
| data/enumerations/rational_wrong_log_maximal_finiteness.json | cf95f1f7ba6d3ed3f60a521419f817874ccb55f3b27d0a0ecde2fb5c0e8c16d7 |
| tests/test_rational_wrong_log_maximal_finiteness.py | 10ced3cc9b611724d8b172c2a3b5cec5f839dfc2cd20b70d72b10cdbe0efb20f |

## 2. The AG elimination, including the rank-26 endpoint

The two inputs are

\[
 D^2=10-\rho+b\geq0,\qquad
 \widehat c_2=\rho+2-\frac{7b}{4}\geq0.
\]

They are exactly

\[
                    b\geq\rho-10,\qquad 7b\leq4\rho+8.
\]

Eliminating \(b\), rather than searching a bounded rank horizon, gives

\[
          7(\rho-10)\leq4\rho+8,\qquad3\rho\leq78.
\]

Thus every rank at least \(27\) is impossible.  At the endpoint, both
inequalities force the unique integer value

\[
                         (\rho,b)=(26,16).                  \tag{2.1}
\]

This distinction matters.  The argument eliminates ranks **above** \(26\);
it does not prove that rank \(26\) itself is empty.  At (2.1),

\[
 D^2=0,\qquad\widehat c_2=0,\qquad\rho(Y)=\rho-b=10,
\]

and the imported semipositivity theorem permits equality.  Therefore the
audited conclusion is exactly \(n\leq25\), not \(n\leq24\).  The hostile
certifier derives the symbolic inequality and then uses ranks through
\(1000\) only as a poison check.

## 3. The adapted lattice and all three physical wall types

Let \(L=\operatorname{NS}(S)\), give it the doubled form, and set

\[
 M=\left\langle L,B_1/2,\ldots,B_b/2\right\rangle,\qquad
 q_M=2q_L.
\]

Disjointness and \(B_i^2=-4\) give

\[
 q_M(B_i/2,B_j/2)=-2\delta_{ij},\qquad
 q_M(B_i/2,x)=B_i.x\in\mathbf Z.
\]

So \(M\) is integral.  For a root \(r\), use

\[
                       s_r(x)=x-\frac{2q_M(x,r)}{q_M(r,r)}r. \tag{3.1}
\]

The preservation checks are then exact.

1. A boundary root \(r=B_i/2\) has norm \(-2\).  Equation (3.1) has
   integral multiplier \(q_M(x,r)\) on every generator of \(M\).
2. A physical \((-1)\)-curve \(E\) has norm \(-2\).  On \(L\),
   \(q_M(x,E)=2x.E\), and on a half-boundary generator,
   \(q_M(B_i/2,E)=B_i.E\).  Both are integral.  No divisibility beyond
   integrality is needed.
3. A physical \((-2)\)-curve \(R\) has norm \(-4\).  On \(L\),
   \(q_M(x,R)=2x.R\), so the multiplier in (3.1) is integral.  On
   \(B_i/2\), weak incidence gives
   \(q_M(B_i/2,R)=B_i.R=0\).

The independent script also retains the hostile control \(B_i.R=1\):
the last multiplier becomes \(1/2\), so the reflection need not preserve
the adapted host.  This confirms that (W2), rather than a silent parity
assumption, supplies the \((-2)\)-wall check.  The bound \(B.E\leq2\) is
used by the AG nef argument, not by host integrality; any integral
\(B_i.E\) is harmless for a norm-minus-two reflection.

## 4. Finite adapted-host interval and exact coarse count

For a rational surface, the standard minimal-model/blowup calculation gives

\[
 L\simeq U\text{ or }I_{1,1}\quad(\rho=2),\qquad
 L\simeq I_{1,\rho-1}\quad(\rho\geq3).
\]

This is elementary surface classification input, which is why the repaired
header promises only freedom from **reflective-lattice** classification.
Unimodularity gives the finite exponent-two interval

\[
                 L(2)\subseteq M\subseteq L(2)^*=\frac12L.  \tag{4.1}
\]

Every intermediate group is specified by a binary subspace
\(H\leq L/2L\), and

\[
 M_H=L+\frac12\widetilde H
\]

is integral precisely when

\[
                         x.y=0\pmod2\quad(x,y\in H).          \tag{4.2}
\]

The boundary span satisfies (4.2); counting every binary subspace only
overcounts.  The hostile replay reconstructs the Gaussian coefficients by
the recurrence

\[
 {n\brack h}_2={n-1\brack h}_2+
                      2^{n-h}{n-1\brack h-1}_2,
\]

not by the producer's product formula.  A literal enumeration of every
subspace through rank five independently gives

| rank | 1 | 2 | 3 | 4 | 5 |
|---:|---:|---:|---:|---:|---:|
| all subspaces | 2 | 5 | 16 | 67 | 374 |
| integral/totally isotropic controls | 1 | 2 | 4 | 11 | 31 |

Using two Picard types in rank two and one in ranks \(3,\ldots,26\), the
independent recurrence recovers

\[
2S_2(2)+\sum_{\rho=3}^{26}S_2(\rho)
=5{,}515{,}561{,}289{,}142{,}524{,}504{,}797{,}715{,}176{,}506{,}813{,}442{,}521{,}589{,}215{,}285{,}286.
\]

This counts nonintegral subspaces and isometric repetitions, exactly as the
owner states.  Those duplications enlarge a finite upper bound and cannot
invalidate finiteness.

## 5. Saturation under each maximality mode

Every facet reflection preserves \(M\), so

\[
                             W(P)\leq R(M)\leq O^+(M).
\]

If an integral reflection \(s\) is missing, then

\[
                       G=\langle W(P),s\rangle              \tag{5.1}
\]

is a proper discrete reflection overgroup: it is discrete because it lies
inside the discrete group \(O^+(M)\).  Therefore either of the following
immediately excludes (5.1):

1. absolute maximality among discrete hyperbolic reflection groups;
2. maximality among reflection subgroups of \(O^+(M)\).

The arithmetic-reflection version also works, with one necessary sentence
made explicit.  The finite-volume group \(W(P)\) and the arithmetic group
\(O^+(M)\) are lattices with \(W(P)\leq O^+(M)\), so \(W(P)\) has finite
index in \(O^+(M)\).  Every intermediate discrete group \(G\) in (5.1) is
therefore finite index and arithmetic.  Maximality among arithmetic
**reflection** groups excludes it.  In all three cases,

\[
                              \boxed{W(P)=R(M)}.             \tag{5.2}
\]

Maximality of the integral lattice, completeness of a known curve list, or
absence of another currently effective wall is not a substitute for any of
these three group-theoretic hypotheses.  The hostile archive records those
negative controls.

## 6. Exactly one complete chamber orbit

The mirrors of the discrete group \(R(M)\) are locally finite.  Indeed, if
infinitely many met a compact set \(K\), their reflections would satisfy
\(sK\cap K\ne\varnothing\), contradicting proper discontinuity.

Choose regular points in two components of the complement and connect them
by a generic compact path.  It crosses finitely many mirrors, one at a time.
Successively reflecting across them walks through adjacent chambers and
carries the first chamber to the second.  Thus the complete arrangement has
one \(R(M)\)-orbit.

After (5.2), the supplied Coxeter fundamental chamber is a chamber of this
complete arrangement: a further \(R(M)\)-mirror through its interior would
be a \(W(P)\)-mirror through the interior of a fundamental chamber.
Consequently one host contributes at most one unmarked chamber isometry
class.  It may contribute none if its full reflection group has no
finite-sided finite-volume chamber.

As a finite negative control, the hostile script realizes all six chambers
of the complete \(A_2\) arrangement as one orbit under its two simple
reflections.  Keeping only one mirror reaches two chambers, confirming that
the proof requires the **complete** mirror arrangement.

## 7. Vinberg field and exact scope

Every physical wall normal belongs to the rational quadratic space
\(M_{\mathbf Q}\).  Hence every Cartan entry

\[
                         \frac{2(r_i,r_j)}{(r_i,r_i)}
\]

is rational, as is every cyclic Cartan/Gram product.  These products generate
the minimal Vinberg field, so

\[
                                k(P)=\mathbf Q.
\]

This is a statement about the abstract polytope field.  It is not a claim
that the carrier surfaces have field of definition \(\mathbf Q\).

The complete audited theorem is carrier-relative.  It requires:

* a literal complete full-nef rational surface carrier;
* the complete \((-1),(-2),(-4)\) Mori alphabet and weak incidences;
* a finite-sided finite-volume Coxeter chamber;
* one of the three reflection-maximality hypotheses above.

It proves neither raw Coxeter-to-carrier essential surjectivity nor a rank
bound for hostile incidences.  It also proves qualitative finiteness, not a
terminating reflectivity decision for every candidate host and not a surface
isomorphism census.

## 8. Frozen replay

Run

    /Users/sasha/.venvs/global/bin/python -m pytest -q \
      tests/test_rational_wrong_log_maximal_finiteness.py \
      tests/test_rational_wrong_log_maximal_finiteness_hostile_audit.py

The hostile suite pins all owner artifacts and independently rebuilds the
AG endpoint, the three wall-reflection divisibility checks, binary host
counts, all maximality implications, complete chamber transitivity, the
rational field calculation, and the scope controls.

The terminal hostile artifacts are frozen as follows:

| artifact | SHA-256 |
|---|---|
| scripts/audit_rational_wrong_log_maximal_finiteness_hostile.py | 7f07c94d75a3c9124b1512a857615440f5bfe04a036288e6992f28d2aa7bdf65 |
| data/enumerations/rational_wrong_log_maximal_finiteness_hostile_audit.json | 43a3c1383e1bd408e3561d258d7118fba2a54b838ac4619c314457b539ab4eb2 |
| tests/test_rational_wrong_log_maximal_finiteness_hostile_audit.py | bb646eecf1fcdfa7bf080ef55400d2ac482b1ae68e046c849db78a680b1d2b97 |
