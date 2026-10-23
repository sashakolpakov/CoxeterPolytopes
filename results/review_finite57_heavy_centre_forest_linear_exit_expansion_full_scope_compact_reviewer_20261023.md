# Full outside-author review: heavy-centre forest and linear exits

23 October2026. **FULL mathematical SOURCE PASS, no remaining correction.**
This review binds the entire current182-line root-authored source, not an
informal lead. I independently rederived its Schur calculation, actual edge
length, hyperbolic induction and incidence count. No mathematical worker,
import, enumeration, test, compiler, Lean, TeX, Sphinx, plan or Git operation
ran. Bounded file reads and SHA256 metadata authenticated the bindings.

## 1. Frozen object, inputs and independence

| Object | Bytes | SHA256 |
|---|---:|---|
| results/finite57_heavy_centre_forest_and_linear_exit_expansion_candidate_20261023.md | 9102 | 059e3fcf6c80700a0b11778638cbb47e758ca245f6cd9573c747e1ff8caa158f |
| manuscripts/finite_volume57.tex | 156120 | c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0 |
| results/finite57_complete_finite_centre_recipients_and_heavy_component_escape_candidate_20261023.md | 15879 | 994d3b9ed36a32e8698610140ee656b40a88507c7f0a7b6b1978faef21b00186 |
| results/review_finite57_finite_centre_safe_union_and_heavy_escape_root_acceptance_20261023.md | 5601 | 1c2e887918d4277efffd9c568aed5d9e203a61dce2ca7c98eca96665831792f9 |

All three literal input pins match. I read the whole new source, the whole
284-line finite-centre source and97-line acceptance, and the manuscript's
actual all-wall, finite-simplicity and edge interfaces. The complete
finite-centre classification and reverse-port table are accepted inputs
here; this report does not claim transitive regeneration of their charge
tables or of the3322-line manuscript.

I authored related inherited actual-face/allocation foundations and checked
root's proposed long-chain mechanism in discussion. I did not author this
new proof. The pinned acceptance records Replay as the finite-centre source
author and Foundations as its outside reviewer. The initial9024-byte source
45ae91ecfb3d3943a31d4aef2f53f37ec1a7a52464f0ec82a977c5639ca3db04
mistakenly attributed that classification to Composition. Root corrected
only that attribution paragraph during review. Replacing the current
paragraph with its exact predecessor recovers the preceding hash and byte
length; all mathematical bytes are unchanged. This correction has no
effect on the verdict or theorem.

## 2. Literal graph and the actual heavy edge

The graph uses all finite heavy vertices of the ORIGINAL polytope, whose
whole forest is2D17+D23, and only actual fork-centre edges with two heavy
endpoints. There are finitely many vertices. Convex-polytope edges do not
give loops or parallel edges between the same two vertices, so the graph
is finite and simple. Every heavy vertex has exactly three distinct actual
centre edges; this is finite simplicity, not an ideal-vertex assumption.

Deleting a centre leaves the actual56-root support
A_s+2A1+D_j+D_k, with(s,j,k)=(14,17,23) or(20,17,17).
In the complete finite table, row9 is the only heavy completion. The new
wall therefore restores the same path and two actual isolated leaves by
three simple bonds. Either path-end orientation is allowed. It is itself
the fork centre at the other endpoint, so the same edge is a centre move
from both ends. There is no unproved choice of a neighbouring diagram and
no assumption about the product of the released and new walls.

In square2 normalization, the positive contact column has entry1 at the
chosen A_s endpoint and at each isolated leaf. The inverse diagonal at
that endpoint is s/(s+1); each A1 inverse is1/2. The full Schur norm is

    2 - s/(s+1) - 1/2 - 1/2 = 1/(s+1).

The cross term is1+s/(s+1) for matching endpoints and1+1/(s+1)
for opposite endpoints. Inert factors contribute zero. All-wall acuteness
gives a.a'<=0, and consequently the two projected actual endpoint normals
satisfy

    u.u' <= -(s+2)/(s+1).

Their positive norms are both1/(s+1). The actual support perpendicular
space is Lorentzian of dimension2 and contains the real edge. With the two
outward endpoint orientations, normalized product is -cosh L. Thus
cosh L>=s+2>=16. This uses intrinsic projected norms correctly; it does
not apply a native-angle gap to unnormalized projected products. Independent
square2 normalization introduces no relative disconnected-root-scale
assumption. No arithmeticity is used.

At a heavy vertex, the three released centres belong to different
orthogonal whole D components. The corresponding actual tangent edge
directions are dual vectors in those blocks. Their products are zero.
Thus two distinct incident centre edges have actual angle pi/2, including
on a path whose successive segments do not lie in one plane.

## 3. Complete hyperbolic-chain induction

For a chain starting at z0, let the angle between the backward edge and
the direction to z0 at z_k be C_k. At the first step C_1=0. If
C_k<pi/4, orthogonality and the triangle inequality on the tangent unit
sphere put the next triangle's included angle B strictly in
(pi/4,3pi/4); the first such angle is exactly pi/2. This prevents the
next point from equalling z0 and ensures a nondegenerate triangle.

Write R=d(z0,z_k)>0, L for the next edge and D for the new start
distance. The cosine rule gives

    cosh D = cosh R cosh L - sinh R sinh L cos B
            > cosh R (cosh L - sinh L/sqrt2)
            > 4 cosh R.

The first lower bound also holds when cos B<=0: then the subtracted
term is nonnegative. For the second, sinh L<cosh L,1/sqrt2<3/4,
and cosh L>=16 imply the factor is strictly greater than cosh L/4>=4.
There is no hidden restriction on the sign of cos B.

The angle C at the new endpoint obeys, by the sine and cosine rules,

    cot C = (coth R sinh L - cosh L cos B)/sin B.

Its numerator is greater than sinh L-3cosh L/4. Since
sinh L/cosh L>=sqrt255/16>15/16, this is greater than
3cosh L/16>=3. Dividing by0<sin B<=1 preserves the strict positive
lower bound, so cot C>3>1 and C<pi/4. This completes the induction.
Every triangle may have a different hyperbolic2-plane; all statements
are intrinsic ambient metric facts. The resulting quantitative estimate
cosh d(z0,z_m)>=16*4^(m-1), strict for m>=2, follows from the base
edge and these successive inequalities.

A first repeated point yields a subchain starting at that point, with
the same long-edge/right-angle hypotheses, and is excluded by this
induction. In particular a simple graph cycle cannot occur. No global
compactness of P or of a selected face is being inferred or needed.

## 4. Exact expansion and receiver multiplicity

The graph is therefore a forest. For any finite nonempty heavy set S,
the induced graph has |S|-c(S) internal edges. Each internal centre
edge contributes two of the3|S| centre incidences; every remaining
incidence leaves S. Therefore

    b(S)=3|S|-2(|S|-c(S))=|S|+2c(S).

This counts ACTUAL incidences, not distinct target vertices. A full
connected heavy component has H+2 actual exits and each target is
nonheavy finite or ideal. For a general subset S, targets may include
other heavy vertices; already discharged heavy vertices are not silently
treated as free receivers.

If the full component's exits are all finite, the COMPLETE reverse-port
union bounds each nonheavy finite target by at most3 actual centre
preimages. The row1 maximum3, row2 maximum2 and the two-port cross-rank
D38 collision are included. Deleting a chosen actual node determines one
actual edge and at most one other endpoint; automorphisms and equal
component names add no multiplicity. Hence at least ceil((H+2)/3)
distinct finite targets occur. No ideal reverse-port bound is assumed.

## 5. Verdict and precise limit

The forest theorem, exact subset incidence formula, H+2 full-component
exit count and conditional finite-target lower bound all PASS at the
written input scopes. The proof supplies no receiver charge reduction,
allocation or new simultaneous rule. Expensive exits, previously received
debt, source/receiver collisions and ideal capacities remain separate
obligations. Linear exit incidence alone does not discharge a component
or prove a Hall condition.

The original C0 rule and all earlier debts are untouched. No universal
finite cap364, dimension56, new example, volume claim, proof execution or
end-to-end formalization follows. Root acceptance is still required.
All frozen papers and global ceilings remain unchanged.
