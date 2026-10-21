# Hostile audit of the minimal H2 adjoint core

**Status (2026-08-10): PASS after one presentation repair and one scoped
conditional addition.**  This is an independent audit of the H2 statement

\[
       B^2=-4,\qquad R^2=-2,\qquad B.R=2
\]

on a smooth projective rational surface.  It checks the adjoint argument,
the contraction iteration, the exact rank-12 identity, the transverse and
tangential singularity scopes, the free exceptional tail, and both
conditional closures.  It does not claim an original-rank theorem, a
reflection-maximal counterfamily, or a census result.

## 1. Adjoint effectivity survives the length-two intersection scheme

For \(D=B+R\), adjunction gives

\[
 K.B=2,quad K.R=0,quad D^2=-2,quad K.D=2,quad p_a(D)=1,
\]

and \((K+D).B=(K+D).R=0\).  The restrictions of \(\omega_D\) to the two
rational components therefore have degree zero and are trivial.  This does
not require the two contacts to be distinct.

The potentially delicate point is whether the unique dualizing section
could vanish on one component when the intersection scheme has length two.
Tensoring the union sequence by the invertible sheaf \(\omega_D\) shows that
the two component restrictions must agree on \(B\cap R\).  If one restriction
were zero, the other would be a section of \(\mathcal O_{\mathbf P^1}\)
vanishing on the nonempty length-two intersection scheme, hence would also
be zero.  Thus the nonzero section supplied by Gorenstein Serre duality is
nowhere zero and \(\omega_D\simeq\mathcal O_D\).

Since a rational surface has \(H^0(K)=H^1(K)=0\), the canonical sequence

\[
 0\longrightarrow\mathcal O_S(K)
  \longrightarrow\mathcal O_S(K+D)
  \longrightarrow\omega_D\longrightarrow0
\]

gives \(h^0(K+D)=1\).  Its unique divisor \(A\) is effective and disjoint
from \(D\).  The report's effectivity, uniqueness, and disjointness claims
are therefore valid for both two transverse contacts and one contact of
length two.

## 2. The contraction iteration and the number 12

If \(A\) is not nef and \(A.E<0\), effectivity makes \(E\) a component of
\(A\), while disjointness gives \(K.E=A.E<0\).  Writing \(E^2=-q\),
adjunction says

\[
                       K.E=2p_a(E)-2+q<0.
\]

Thus \(p_a(E)=0\) and \(q=1\): \(E\) is a smooth rational \((-1)\)-curve
disjoint from \(D\).  Its blowdown is an isomorphism near the pair, and the
same argument can be repeated.  Rank strictly decreases, so the iteration
terminates with \(K+D\) nef.

On that terminal core, Noether's identity and the unchanged pair
intersections give exactly

\[
 (K+D)^2=(10-\rho)+2K.D+D^2=10-\rho+4-2=12-\rho.
\]

Nefness forces \(\rho\leq12\).  The stated blowup of two transverse
\((1,1)\)-curves on \(\mathbf P^1\times\mathbf P^1\), at six private points
on one and four on the other, has \(K+B+R=0\) and \(\rho=12\).  It proves
sharpness for smooth rational pairs only; the owner correctly attaches no
complete-Mori or reflection-maximal claim to it.

## 3. Contraction and discrepancy scope

The pair matrix and canonical-degree vector are

\[
 Q=\begin{pmatrix}-4&2\\2&-2\end{pmatrix},\qquad
 \det Q=4,qquad
 Q\binom{a_B}{a_R}=\binom20.
\]

The unique solution is \((a_B,a_R)=(-1,-1)\).  With two distinct transverse
contacts, the exceptional divisor is an SNC two-component cusp, so the
contraction is lc and not klt.  This does not license the quotient-orbifold
\(c_2\) formula used in the klt lane.

The same matrix does not control a tangency.  After blowing up an ordinary
order-two tangency, the two strict transforms and the first exceptional
curve have discrepancy \(-1\) and meet at one triple point.  Blowing up that
point produces discrepancy \(1-1-1-1=-2\).  Hence the tangential contraction
is not lc.  The owner's separation of the SNC and tangential scopes is
necessary and correct.

## 4. The free tail and what it does not prove

For \(m\) successive free blowups away from \(D\), total-transform
exceptional coordinates give the physical chain

\[
 e_1-e_2,\ e_2-e_3,\ldots,e_{m-1}-e_m,\ e_m.
\]

In the standard host \(q_M=2q_{\operatorname {NS}}\), its Gram matrix has
diagonal \((-4,\ldots,-4,-2)\) and adjacent entries \(2\).  Reflection in
\(e_i-e_{i+1}\) swaps the two coordinates, while reflection in \(e_m\)
changes the last sign.  These generate every signed permutation, hence

\[
             W(B_m)=O(\langle-2\rangle^m),\qquad |W(B_m)|=2^m m!.
\]

Thus there is no omitted *pure-exceptional* reflection.  This local
saturation neither proves finite volume or reflection maximality of a full
carrier nor bounds mixed core/tail roots.  The original rank remains
\(\rho(S)=\rho(S_0)+m\), with \(m\) uncontrolled.

## 5. Conditional consequences are correctly scoped

If every component of the contracted tail is a physical facet and an
integral ample marking has an absolute total physical-facet degree

\[
                     T=\sum_C H.C,
\]

then each of the \(m\) components contributes a positive integer, so

\[
       m\leq T,\qquad \rho(S)\leq12+T,
       \qquad\dim\mathbf P(\operatorname {Nef}S)\leq11+T.
\]

Neither a per-facet degree bound nor a tail not known to consist of facets
implies this conclusion.

Separately, if a full trace lattice of rank \(d(n+1)\) genuinely descends
through this contraction, then it embeds in \(\operatorname {NS}(S_0)\), so
\(d(n+1)\leq12\), and \(n\leq11\) over \(\mathbf Q\).  Descent is not a
formal consequence of finite orbit or reflection maximality: in a literal
full-nef carrier the tail components can themselves lie in and span part of
the visible trace block.

## 6. Repairs and frozen scope

The hostile pass found pervasive lost inline-TeX delimiters in the report;
these were repaired without changing the proof.  It also added the
total-marked-degree consequence above, including both necessary premises
and the warning that an individual-facet degree bound is insufficient.

The independent regression test recomputes the pair determinant,
discrepancies, core-square constant, equality blowup arithmetic, and
\(B_m\) reflection actions without importing the owner's certifier.  It
also pins the archive and owner-file bytes, checks the archived scope fields,
and scans the report for control bytes and the known delimiter-loss
signatures.

The audited owner files have SHA-256 hashes

```text
1464b327d588f9a5e8c5b2a6809667194d27c09a69aad5805250df4053c9634b  results/rational_h2_minus4_minus2_genus_one_core.md
f75dc268fbaeea9efd23f39ebca8c7a29ac9903a0f2d20422c6399cad2f09adf  scripts/certify_rational_h2_minus4_minus2_genus_one_core.py
d240a14f49df3672ce604aae0f7ddb67cf3700eca966b96b3a83e59917a5b7e5  data/enumerations/rational_h2_minus4_minus2_genus_one_core.json
f0153a604090ac2a65e2179183a075ca6a4093ad68b2f5b370a02ac0621fb80d  tests/test_rational_h2_minus4_minus2_genus_one_core.py
```

Subject to those exact bytes, the hostile verdict is **PASS** for the
rank-12 adjoint-core theorem and its explicitly conditional consequences.
It is not a PASS for an absolute original-rank or dimension theorem.
