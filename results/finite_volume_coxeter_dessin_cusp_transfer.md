# Finite-volume Coxeter dessins: rational cusp transfer

**Status (2026-08-09): proved general horospherical finiteness theorem
with an exact three-cusp \(H^6/9\) calibration.**  This extends the compact
hidden-Heegner-sheet gate to finite-volume noncompact arithmetic quadratic
Coxeter input.  It neither asserts that a requested Heegner divisor is
principal nor produces a new polytope by itself.

## Theorem

Let \(L\) have signature \((1,n)\), let \(P\subset\mathbb H(L)\) be a
finite-volume rational polyhedral Coxeter chamber, and fix finitely many
root norms and fixed-denominator discriminant cosets
\(\mu+L\subset L^\vee\).  Fix the determinant/overlattice extension, or
choose it from an independently proved finite list, so all denominators are
uniformly bounded.  Require every cusp of the finite-sided chamber to have
full parabolic rank \(n-1\), equivalently a compact Euclidean horosphere
cell under its cusp translation lattice.

At a rational cusp choose primitive isotropic
\(u\), fixed rational isotropic \(v\) with \((u,v)=1\), and let
\(K=(u,v)^\perp\), with positive norm \(\|z\|^2=-(z,z)\).
Then the projections of the chosen lattice cosets to \(K\), and the allowed
\(\alpha=(r,v)\), remain discrete.

Use horospherical coordinates

\[
x(h,z)=\frac{h^2+\|z\|^2}{2h}u+\frac1h v+\frac1h z,
\qquad (x,x)=1.
\]

Write a candidate negative root as

\[
r=\alpha u+\beta v+w,\qquad
\beta=(r,u),\quad\alpha=(r,v),\quad w\in K,
\]

and \(r^2=-\rho\), with \(\rho>0\).  The universal identity is

\[
\beta^2h^2+\|\beta z-w\|^2-\rho
=2h\beta\,(x(h,z),r).
\]

Consequently, when \(\beta\ne0\),

\[
(x(h,z),r)=0
\quad\Longleftrightarrow\quad
\boxed{\beta^2h^2+\|\beta z-w\|^2=\rho}.
\]

1. if \(\beta\ne0\), then
   \(h\le\sqrt{\rho}/|\beta|\).  On an integral or fixed-denominator coset,
   nonzero \(\beta\) is bounded away from zero, so all such mirrors lie in
   a computably truncated region;
2. if \(\beta=0\), the completed-square equality is only the root-norm
   identity and is not equivalent to the mirror equation.  Directly,
   \[
   (x(h,z),r)=\frac{\alpha-\langle z,w\rangle}{h}.
   \]
   Thus \(\|w\|^2=\rho\) and the mirror equation is independently
   \(\langle z,w\rangle=\alpha\).  The positive-definite shell for \(w\) is
   finite, and a compact cusp-translation cell bounds the discrete values
   of \(\alpha\).

There are finitely many cusps.  Hence the special-divisor sheets in any
fixed finite norm/discriminant support which meet \(P\) are exactly and
finitely enumerable.  Combining this enumeration with Vinberg's literal
maximal-parabolic-neighborhood test distinguishes listed vertical branches,
bounded thick-part sheets, harmless boundary contacts, and forbidden hidden
interior cuts.

## Proof

The Witt decomposition gives

\[
r^2=2\alpha\beta-\|w\|^2=-\rho.
\]

Multiplying \((x,r)=0\) by \(2h\) gives

\[
\beta h^2+\beta\|z\|^2+2\alpha-2\langle z,w\rangle=0.
\]

Multiply by \(\beta\), substitute
\(2\alpha\beta=\|w\|^2-\rho\), and complete the
positive-definite square.  Without assuming the mirror equation this gives
the universal identity above.  Division by \(\beta\) proves the equivalence
only in the nonvertical case.  Setting \(\beta=0\) in the original, unmultiplied
pairing gives the separate vertical equation.  Discreteness and bounded-shell
finiteness then prove the two cases.

## Field split

For an admissible arithmetic quadratic (simple-type) form over a totally real field
\(k\), every non-distinguished real form is definite.  If
\(k\ne\mathbb Q\), a \(k\)-isotropic vector would remain isotropic at a
definite real place, which is impossible.  Thus the quotient is compact.
Finite-volume noncompact input lies over \(\mathbb Q\), exactly where the
rational cusp theorem applies.  Meyer isotropy further leaves compact
rational quadratic input only in \(n\le3\).

A finite-covolume reflection subgroup has finite index in its ambient
arithmetic lattice, so the compact/isotropic alternative transfers to its
Coxeter chamber.  No claim is made here for a hypothetical non-quadratic
arithmetic carrier without this admissible reflection representation.

## Exact calibration

The certifier replays the known toric \(H^6/9\) chamber with 18 finite
vertices, three ideal vertices, and no hyperideal ray.  At each cusp it uses

\[
u=H-E_i,\qquad v=(H+E_i)/2
\]

and decomposes all nine walls.  Exactly seven have \(\beta=0\), and they are
literally the seven active facets of the independently certified maximal
parabolic neighborhood.  The other two have exact finite height caps.
Thus all three cusps have the predicted \(7+2\) split.

Finite volume is inherited only from the independent 18+3+0 ray and full
parabolic audit, never from finite sidedness.

## Remaining obligations

This theorem closes the missing finite-volume hidden-sheet finiteness step.
The general Coxeter-dessin programme must still compile each finite
special-divisor support, replay every face and parabolic gate, prove the
decorated-surface faithfulness theorem, and obtain nonzero classification
pushback on a non-singleton P-only class.

Artifacts:

- producer: scripts/certify_finite_volume_coxeter_dessin_cusp_transfer.py
- archive: data/enumerations/finite_volume_coxeter_dessin_cusp_transfer.json
- tests: tests/test_finite_volume_coxeter_dessin_cusp_transfer.py
