# All lengths of a single attached white chain

The [checker](../scripts/fresh_sqrt5_single_chain_completion_20261015.py) and
[certificate](fresh_sqrt5_single_chain_completion_20261015.json) solve the
rank-completion equation for every possible length of one simply-laced
white chain attached to a nonsingular finite core. Only its first node
touches the core. There is no length cutoff.

Let H be the scaled core Gram matrix, v the first-node contact column,
d its positive squared norm, and c=vᵀH⁻¹v. An attached A_l chain has block
(d/2)Cartan(A_l), whose first inverse diagonal is2l/[d(l+1)]. The full Gram
can have the required corank1 only if

\[
(l+1)d-2lc=0,\qquad l=\frac{d}{2c-d}.
\]

The denominator-zero case has no solution because d>0. In all other cases
the displayed field element must be a positive rational integer. Thus each
contact profile has at most one possible length, and the two embeddings
satisfy the same exact rank equation.

Among all458 accepted first-node contact profiles,37 have zero denominator
and258 have no positive integer length. The remaining163 give lengths1
(147 profiles),2 (9),3 (5), and4 (2). All163 pass exact finite-corner,
physical-inertia and conjugate-inertia checks and remain unpaid by the
baseline short-path/H4 rule. They are principal stars requiring further
analysis, not new finite-volume polytopes. The longer chains add16 literal
stars beyond the previously studied first extensions.

Entire external components orthogonal to both the core and this chain do
not affect the rank condition. Multiple touching chains, branching,
two-ended attachments and dotted cores remain outside this lemma. The
calculation is a complete all-length result for its specified family,
not a complete quadrilateral classification. The
[independent review](fresh_sqrt5_single_chain_completion_20261015_review.md)
passes, separately checking all458 costs, all163 whole-star inertias and
all652 corner matrices in both embeddings.
