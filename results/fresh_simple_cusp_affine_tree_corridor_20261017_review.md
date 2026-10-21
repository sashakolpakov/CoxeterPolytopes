# Review of the long affine-tree corridor lemma

**PASS for the stated central-corridor branch.** Root read the
[proof](fresh_simple_cusp_affine_tree_corridor_20261017.md) and
[source](../scripts/fresh_simple_cusp_affine_tree_corridor_20261017.py).
This is an independent mathematical/source review, not an independently
implemented arithmetic replay. The package's own exact checks and saved
replay retain their separate scope.

Deleting either marked root leaves two whole finite B/D components of
rank at least nine. A finite extension can only lengthen one exposed long
end by a simple bond; the exceptional diagrams are too small. An affine
extension touching just one component would leave a disconnected corner,
incompatible with a simple ideal vertex. Joining both components yields
exactly the two-long-port profile of the original affine B/C/D diagram.
These observations justify the complete four-profile list in this scope.

The small-affine exclusions were checked with the end positions explicit.
A new simple leaf at a marked root, the outer exposed port, or the first
inner white root admits the E8 arms of lengths 1,2,5. At marked distance
four, the last embedding continues through the other marked root into
the opposite long endblock. The rank-nine hypotheses provide enough
simple bonds before a double bond or a fork; the old graph is a tree, so
there is no unrecorded closing edge. Every excluded affine order is below
the ambient n=r+s+P+1 >= 23.

Direct multiplication confirms the stated all-rank inverse port columns
for B and D, with diagonal entry one. White elimination then has
alpha=1/P, outer-single constant one and moment zero, inner-single constant
one and moment one, and both-single constant zero and moment one. The
only cross term comes from the intervening A path. Substitution into
the reviewed four-by-four determinant formula cancels the quadratic
cross term exactly. Each of the two remaining single-root contributions
is at most one on its height interval. Therefore the determinant is at
least 1-2/P >= 1/2, contradicting the actual star's zero determinant if
the opposite marked roots were orthogonal.

The finite rank-nine-through-100 checks confirm an explicit symbolic
identity; they are not the reason for an all-rank conclusion. Cases near
the decorated ends, small endblocks, and the full quadrilateral-payment
theorem remain outside this result. No new ambient dimension bound is
claimed.
