# Root review: mixed A1/A2 farthest vertices

20 October 2026. Full independent mathematical read of
`compact_a1_a2_farthest_vertex_energy_exclusion_20261020.md`, SHA256
`3bcf38496cd0aa7ba922dd13101595a479c0c14937d5c08bb368e2df435707e2`.
PASS for its stated pointwise theorem. No execution or dimension ceiling.

I rederived the sign of the farthest-point inequality, the A2 constrained
energy bounds, and the Schur inverse identity proving the distinct-wall
shadow. The distinction between a wall closing several directions and
several distinct walls is retained throughout. Positive definite replacement
subsets give actual endpoints: the inverse Gram has nonnegative entries,
so the projected point satisfies every other wall by acuteness.

For A1-short directions I checked the native-square intervals for supports
of size one through four and the additional A2 costs. The only mixed
possibility is A1_m/A2_3, m>=7. Two such walls can coexist only on opposite
nodes of one A2 block, at distinct A1 nodes. A mixed wall cannot coexist
with an A1-only short wall. In the latter subfamily all shared-wall and
distinct-wall cases give at most two short directions and the stated
half-tanh-square credit. This establishes the debit of less than3/25.

For A2-short directions I checked the double-contact exclusion and all
off-block cost windows. The surviving inventories are exactly isolated
finite/dotted, A2_3/A2_5, and A1_3/A2_m for m=5,7,8,9,10. The affine
null subdiagrams used to discard other native sums are explicit. Releasing
an uncontacted mate cannot introduce a missed profile: the contacted
direction then closes earlier at the same wall.

The finite isolated, F, and N5 mate replacements are positive definite
and actual. For dotted isolated and N7--10 their white Gram is not
positive definite, so their mate has a DISTINCT opposite. The shadow
then supplies its large energy without enumerating that other opposite.
I checked both rational bounds against16/25. Only isolated7 can fall
below the block budget319/2000. At most one such block occurs, and its
wall is incompatible with every A1-short wall. This correctly prevents
adding the two worst debits in the final estimate.

The n>=15 estimate and both n14 parity cases exceed one strictly. At
zero component counts the non-strict zero bound is used, and strictness
comes from a nonempty component. No other vertex is restricted. This
removes the A1/A2 family at the single farthest vertex, without changing
the compact29 or arithmetic compact20 ceilings.
