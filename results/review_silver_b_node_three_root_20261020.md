# Root review: the native-eight star at B node three

20 October 2026. **PASS, complete written derivation.** No scientific
execution, enumeration, build or Lean work.

Read in full:
`results/silver_b_node_three_triple_cap_lattice_extension_20261020.md`,
SHA256 `47ec9726572e9ccc5e87295689849e45a2b23c2499a3457beaa4cd2bffb6ce8f`.
Its unchanged complete flag input is
`silver_b_vertex_invariant_lattice_flags_and_even_neighbor_20261020.md`,
SHA256 `2fed0220943f80acd6747119e98e1edee39d2ea4fe9073acff5ce6b894770264`,
separately fully reviewed in `review_silver_b_vertex_flags_root_20261020.md`.

I independently checked the triple cap's norm and exact reflection
difference, both off-diagonal matrix-unit extraction formulas and the
complete primitive-vector cases. Triple indicators span the full
spatial module over the dyadic DVR: differences give augmentation,
and their total coordinate sum is the unit 3. This is valid for n>=4,
including the smallest case where two other indices are needed.
The resulting p^3 matrix units first give p^3L; the separate spatial-unit
and time-only cases then give p^2L. The argument does not assume a
diagonal invariant lattice.

Modulo p^2, all triple operators have common kernel D and their image
vectors span H. If the nonzero top code is not D, permutation
invariance forces every triple image into the middle code, hence H.
Conversely either top D or middle containing H suffices, for EVERY
equivariant lift. This is exactly the singleton-cap criterion.
Consequently every triple-cap invariant lattice is actually
singleton-cap invariant in the same coordinates; the existing full
flag theorem applies without a new enumeration or equivalence claim.
At odd primes the cap's unit mixing again identifies coordinate
ideals, so the transfer is global.

For the actual native8 star at beta3, scaling the cap by `-1/g`
gives contact -1 and norm `2/(2+sqrt(2))=2-sqrt(2)`. Its spatial
part is e1+e2+e3, leaving time norm `-(1+sqrt(2))`. The resulting
full-rank Lorentz Gram over Q(sqrt(2)) proves equality of the actual
ground field by the nonphysical embedding argument. The finite
left-coset lattice sum supplies a genuinely invariant lattice for
the arithmetic group.

Group containment in O(L), or O(S) for even n, now contradicts the
reviewed literal L obstruction for n>=9 or the even-neighbor facet
obstruction for even n>=10. The full integral chamber is contained
in the supposed compact actual chamber. No reflectivity transfer
from rational commensurability is invoked.

Thus a compact arithmetic first-type polytope with a whole B_n vertex
cannot have an actual isolated native8 cap at node3 when n>=9.
Together with the separately reviewed long-end composition this
excludes contact positions1 and3. Node2, all other contact roles,
and a uniform compact dimension bound remain open. This is one
addition to the silver research packet, not another manuscript or
formalization project.
