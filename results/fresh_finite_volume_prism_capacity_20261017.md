# Compact prism capacities in a finite-volume chamber

The [shared prism/simplex capacity theorem](fresh_prism_simplex_shared_capacity_20261017.md)
extends unchanged to compact faces of a finite-sided, finite-volume acute
Coxeter polytope P in hyperbolic n-space. Its other vertices may be ideal
and nonsimple. Select maximal **compact** prism faces
Q = Delta1 x Delta(k-1), k >= 2, and let S be a maximal compact simplex
face of dimension d. For k >= 3 count the two simplex bases of Q; for a
square count all four edge bases. Then

\[
 \sum_{k\ge2}\sum_Q k\,m(Q,S,k)\le(d+1)(n-d),
\]

where m counts the selected bases contained in S. Maximality here is among
compact prism faces, and all faces and base incidences are actual. No claim
is made for a general collection of overlapping diagram aggregates.

Work in the convex projective closure of P. Every finite vertex is simple.
The [compact-simplex transfer](fresh_finite_volume_simplex_truncation_20261017.md)
therefore gives a unique maximal compact simplex containing each compact
simplex base, including a compact edge that lies in no larger compact
simplex. The remaining steps of the compact proof use only the following
local facts.

1. In the triangle-and-two-squares completion, all six named vertices are
   finite because the three supplied faces are compact. The four vertices
   whose three neighbors are already fixed are simple in the generated
   three-face. If any further vertex existed, whether finite or ideal,
   deleting the other two named vertices would disconnect its graph.
   Balinski's theorem excludes it. The resulting triangular prism has only
   the six finite vertices and is compact.
2. In the higher-dimensional completion, every displayed vertex belongs to
   a supplied compact square or to its compact simplex base. Each is finite
   and already has the full degree of the generated face. Connectivity
   excludes every additional vertex, including ideal vertices. The facet
   reconstruction in the original proof then gives a compact prism.
3. The face interval at a finite vertex of S is Boolean. Hence S has exactly
   n-d immediate cofaces, even if some of those cofaces are noncompact.
   Two selected prisms in one such coface whose bases share a vertex would
   complete to a larger **compact** prism by step2, contradicting the
   stated maximality. Thus their base vertex sets are disjoint. A k-prism
   base uses k of S's d+1 vertices, proving the capacity after summing over
   the n-d cofaces.

The square convention and the possibility of two bases of one selected
prism are handled exactly as in the original proof. The argument does not
assume that the undisplayed vertices of a coface are simple. Its
connectivity input is the original
[Balinski theorem, printed p.432](https://msp.org/pjm/1961/11-2/pjm-v11-n2-p03-s.pdf).

Consequently the same signed incidence-transfer columns are valid when the
model rows represent actual maximal compact prisms and compact simplex
faces of P. This does not supply a charge for noncompact faces, nor a
capacity at an ideal vertex. In particular it does not remove the cusp
term from equation(11) of the
[finite-volume truncation and incidence report](fresh_finite_volume_simplex_truncation_20261017.md).
No new dimension bound or numerical optimization is asserted.

**FINAL/FROZEN.** This is a written proof transfer, with no arithmetic
certificate needed. Frontier_literature independently read the complete
compact proof and checked the three steps above. Construction_audit
separately confirmed the same finite-vertex/Balinski/Boolean-interval
transfer to root. These are two internal conceptual reviews, not external
review or formal proof-assistant validation. No predecessor file was edited.
