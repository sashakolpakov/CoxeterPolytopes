# Hostile audit: Vinberg integral points on quadric-arrangement moduli

**Audited owner bytes:**
`b2111131709ebce49dccfcb6b7009a994894a413c8e1432ff195a4c5ce2f63eb`

**Verdict:** **PASS, at the stated scoped strength.**  The report proves an
exact arithmetic-moduli translation, all-field finiteness on the Lambert
curve, and counterexamples to automatic global boundedness.  It does not
claim global arithmetic-reflection finiteness.

## Checks that matter

1. **Invariant coordinates.**  Under (r_i\mapsto\lambda_i r_i), the Cartan
   entry has weight (lambda_j/\lambda_i), so closed directed products are
   invariant.  Directed cycles generate the invariant monomials for this
   diagonal torus action.  The report correctly distinguishes Vinberg
   integrality of these cycles from integrality of individual Cartan entries
   and from good reduction of the marked arrangement.

2. **Lambert intervals.**  With (u=A-4), (v=B-4), one has (uv=4).
   Swapping the dotted pairs gives (0<u_0\leq2).  At every nonselected
   place, ambient definiteness makes the two relevant two-root blocks
   positive definite, hence (A^\sigma,B^\sigma\in(0,4)).  Therefore
   (u^\sigma,v^\sigma\in(-4,0)), and their product sharpens both intervals
   to ((-4,-1)).  The patched text no longer confuses the rank-three
   (4\times4) facet Gram, which is only positive semidefinite, with the
   positive-definite conjugate ambient form.

3. **Field primitivity.**  The selected and nonselected intervals are
   disjoint, so the selected embedding of \(\mathbf Q(u)\) cannot have two
   extensions to (k).  Equivalently the nonzero Lambert cycle graph is a
   tree and
   (k=\mathbf Q(A,B)=\mathbf Q(u)), since (B-4=4/(A-4)).  Thus the norm
   exponent really is ([k:\mathbf Q]-1), not merely
   ([\mathbf Q(u):\mathbf Q]-1).

4. **The cutoff (92).**  Writing
   (p=(T+1)(T+2)(T+3)(T+4)) gives
   (p([-4,-1])=[-1,9/16]), hence
   \(\|p(p+1)\|=225/256\).  On ([0,2]) the maximum is
   (360\cdot361=129960).  Since (F(u)\neq0),
   
   \[
   1\leq129960(225/256)^{d-1}.
   \]
   
   The exact comparisons at exponents (91) and (92) make (d=92) the
   last degree not excluded.  Bounded degree plus (|u^\sigma|\leq4) bounds
   every coefficient of the monic minimal polynomial, proving finiteness of
   (u), fields, and labelled Lambert Grams.

5. **Escape controls.**  The cyclotomic packet genuinely has bounded height
   and unbounded degree.  The Salem-power packet genuinely has fixed-degree
   trace field and unbounded distinguished height; it is explicitly and
   honestly confined to the (n=1) obstruction role.  For the Farey
   polygons, the displayed integral normals have norm one, consecutive
   pairings (-1), (f=N+2), area (N\pi), index (N) in the ideal
   triangle group, and an unbounded end-to-end two-cycle.  The determinant
   formulae exclude triple line concurrence.  At each cusp the conic and two
   half-weight lines give an lc triple point, and
   (K+\Delta=(N/2)H), so the pairs are stable.  The cross-ratio (N)
   forces a new bad prime along prime values of (N).  Their nonmaximality is
   stated, so they are controls on missing hypotheses, not counterexamples to
   maximal finiteness.

## Replay and scope

The owner certificate and test have hashes
`927967d97ea1b90fb399cb631d47b98e04eec72356ba27de598f30b60e5aa025`
and
`1c4f61fed2864af763b1b3c41746fe50a364cebfa801ec536964c06cfad7c5c6`.
The test passes.  It checks the symbolic identities; Vinberg's criterion,
the reflection-subgroup index, and the moduli interpretation remain
mathematical inputs proved in the text rather than things a symbolic test can
certify.

No substantive defect remains.  In particular, the report never promotes
the Lambert contraction to a universal coordinate and explicitly identifies
the missing implication from absolute arithmetic reflection maximality to
bounded complexity, degree, and height.
