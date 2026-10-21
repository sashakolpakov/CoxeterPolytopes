# Hostile audit of the genus-three vertical-blowup invariant no-go

## Verdict

**PASS after one scope-precision repair.**  The
[vertical-blowup note](rational_genus3_vertical_blowup_invariant_no_go.md)
correctly proves that the generic curve, relative canonical algebra, Hodge
bundle, relative canonical model, relatively minimal model, and stabilized
model cannot detect vertical blowup depth.  The exact raw changes are
\((-1,+1)\) in \((K_f^2,e_f)\), so relative Noether cancels every such
blowup.

The local physical tail also has the claimed \(B_m\) signed-permutation
group, canonical colors, and full mixed-orbit span.  These are necessary
local and finite-volume attachment checks only.  The note does not claim a
complete-nef or reflection-maximal counterfamily.

The hostile audit repaired one terminology error: a smooth free centre
produces an unstable tail, whereas a nodal centre produces a bridge.  The
revised text now says explicitly that this description is made after the
usual semistable reduction and refers to the stabilized dual graph.

## 1. The relative invariants and blowup shifts are exact

For a genus-three fibration over \(\mathbf P^1\),

\[
 K_f=K_S+2F,\qquad K_S.F=4,\qquad F^2=0.
\]

On a rational surface, \(K_S^2=10-\rho\) and
\(e(S)=\rho+2\).  Since the smooth genus-three fibre has Euler number
\(-4\), direct substitution gives

\[
 K_f^2=26-\rho,\qquad e_f=\rho+10,
 \qquad K_f^2+e_f=36.
\]

Moreover

\[
 \chi_f=\chi(\mathcal O_S)+(3-1)=3,
\]

so the last equality is exactly \(12\chi_f\).  A point blowup adds one to
\(\rho,c_2,e_f\), subtracts one from \(K^2,K_f^2\), and leaves the generic
fibre fixed.  The independent tests verify these identities symbolically
over large unrelated ranks and blowup counts.

The warning about nefness is essential.  For the newest vertical
exceptional curve \(A\), one has \(A^2=-1\), \(F.A=0\), and
\(K_f.A=-1\).  Hence a relatively minimal slope or nefness inequality
cannot be applied to the nonminimal total space without first contracting
the very tail under investigation.

## 2. Canonical and stable models really erase the blowups

For a composition of point blowups between smooth surfaces,

\[
 K_{S/\mathbf P^1}=\pi^*K_{X/\mathbf P^1}+D_{\rm exc}
\]

with effective exceptional discrepancy.  For every \(n\geq0\), normality
of \(X\) gives

\[
                       \pi_*\mathcal O_S(nD_{\rm exc})=
                       \mathcal O_X.
\]

Indeed, a rational function with no divisorial pole downstairs is regular
away from finitely many codimension-two points and hence extends across
them.  Projection formula therefore identifies every graded piece and the
multiplication maps of the relative canonical algebras.  Relative Proj and
the Hodge summand \(n=1\) are unchanged.

The generic curves are identical, so uniqueness also identifies their
relative minimal models and, after any required finite base change, their
stable models.  What stays fixed is the stabilized fibre data.  Raw dual
graphs on deliberately blown-up smooth models of course acquire extra
vertices; stabilization contracts them.

## 3. The physical free tail and its colors are correct

Start from a smooth rational vertical curve \(C_0\) with \(C_0^2=-1\).
After the first blowup on \(C_0\), make every later blowup at a free point
of the newest exceptional curve, away from its intersections.  In
orthogonal total-transform coordinates the actual curves are

\[
 C=C_0-e_1,\qquad R_i=e_i-e_{i+1},\qquad A=e_m.
\]

They form a simple chain with squares \(-2,\ldots,-2,-1\).  Writing the
canonical discrepancy in these coordinates as \(e_1+\cdots+e_m\) gives

\[
 K.C=K.R_i=0,\qquad K.A=-1.
\]

After doubling the form, the root norms are \(-4\) and \(-2\), while the
canonical colors are \(0\) and \(-2\).  Because the construction occurs in
a fibre different from the distinguished affine-H3 fibre, its classes are
orthogonal to every displayed half-boundary component; the pure-tail
reflections therefore preserve the adapted host.

## 4. Signed permutations and mixed orbit do not bound depth

The reflections in

\[
 e_1-e_2,\ldots,e_{m-1}-e_m,e_m
\]

act as adjacent coordinate transpositions and a sign change of the last
coordinate.  Conjugating the sign change by the transpositions changes any
one sign, so these reflections generate all
\(2^m m!\) signed permutation matrices: the Weyl group \(W(B_m)\).

The exceptional projection of the mixed curve \(C_0-e_1\) is \(-e_1\).
Its orbit is exactly

\[
                         \{\pm e_1,\ldots,\pm e_m\},
\]

which has linear rank \(m\).  The independent hostile tests enumerate the
entire group in small ranks and verify the orbit and Gram/color formulas in
larger ranks without importing the producer certifier.

This proves only that canonical coloring, pure-tail saturation, and the
known orbit-span necessity do not control \(m\).  It does not prove that all
other Mori walls exist, that the resulting nef cone has finite volume, or
that its reflection group is maximal.  The construction is also
conditional on the stated local \((-1)\)-component; it is not advertised as
a universal fibre configuration.

## 5. Frozen producer artifacts

The repaired producer bytes pinned by this audit are:

    a2ee22bea4b5e9ccfafcc2ad2e99f12da9d2d26162c5e2e9c10ef58d7907adbb  results/rational_genus3_vertical_blowup_invariant_no_go.md
    6a8d5bafb99ed1f279f48716ee16bacf6bc5d2d68ed7231547a6153760651a06  scripts/certify_rational_genus3_vertical_blowup_invariant_no_go.py
    9f45165f365e24c08ffd44bb61cfe2f51c2fa169054a8fa8382f9454fa2ccc9b  tests/test_rational_genus3_vertical_blowup_invariant_no_go.py

The surviving gate is genuinely global: a total marked-degree theorem or
an equivariant contraction/descent theorem must use complete-nef and
maximality information which stable genus-three geometry forgets.
