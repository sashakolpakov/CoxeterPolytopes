# Fixed-source rank-four coupling classification

The exact classification leaves only two irreducible nontrivial support types: **D4 and B4**. For the fixed seven-facet mixed compact H5 source, facets 0–5 must remain ordinary; facet 6 may remain ordinary or acquire a support weight of squared norm 1. This is a classification of the specified coupling template, not a nonexistence theorem for compact H9 polytopes.

The [producer](../scripts/fresh_h9_alternative_mixed_couplings_20261014.py) and [certificate](fresh_h9_alternative_mixed_couplings_20261014.json) use exact arithmetic in Q(sqrt(2),sqrt(5)). They bind the frozen [mixed-source geometry](fresh_h9_mixed_d4_coupling_20261014.json) and its independently implemented scalar arithmetic. Generation and saved replay pass.

## Hypotheses and reduction

Retain four spherical support facets with roots of norm 2. Use exactly one lifted wall for each of the seven fixed source facets, with its projected normal a positive multiple of that source normal. Normalize the lifted wall to norm 2 and write it as `alpha*r - w`, where the normalized source normal has norm 2, the support contacts are `c_i = 2 cos(pi/m_i)`, and `w = sum_i c_i omega_i`. If J is the support inverse Gram matrix, then

`w² = cᵀ J c < 2`, and `alpha² = 1 - w²/2`.

For two source walls whose normalized pairing is `-s`, their lifted pairing is `-s*alpha*alpha' + cᵀJc'`. The inverse J is entrywise nonnegative. Thus lifting cannot increase the magnitude of a negative finite source pairing.

The finite computation tests all 4,096 support labelings in {2,3,4,5}, reduces the positive definite ones under all 24 node permutations, and obtains 18 spherical types. It then tests 216 positive-projection contact profiles and 5,168 ordered profile pairs at the source's pi/5 edge. Only the ordinary/ordinary pair survives. Each following source pi/3 edge forces the next wall ordinary: a nonzero weight would turn its pairing into a value strictly between -1 and 0. This propagates through facets 0–5. At the remaining pi/4 edge, facet 6 can have only `alpha=1` or `alpha=1/sqrt(2)`, hence weight square 0 or 1. The source's remaining dotted pair stays strictly dotted in every square-1 case, checked exactly from the physical source vectors.

## Unbounded labels are covered

An internal spherical support label at least 6 must be an isolated I2 factor. Otherwise a three-node path already has cosine-square sum at least 1 and cannot be positive definite. The inverse diagonal of that factor is at least 2, so positive projected norm forces every lifted wall to be right to it. Replacing this inactive factor by two isolated A1 nodes embeds its active-profile calculations in the finite census.

A facet contact label at least 6 has contact square at least 3. The census verifies inverse diagonals at least 2/3 on every non-A1 node, so such a contact is possible only at an isolated A1 and must be the wall's sole support contact. Its projected scale is at most 1/2. At a finite source neighbor, either the neighbor is right to this A1, giving a forbidden pairing magnitude strictly between 0 and phi/2, or it also contacts the A1, in which case the positive support contribution makes the lifted pairing obtuse. Here every finite source pairing has magnitude at most phi, and phi<2. Therefore all contact labels at least 6 are excluded without a finite cutoff.

## Remaining constructions

There are 38 square-1 profiles across the 18 canonical types. The irreducible cases are the three outer weights of D4 and one end weight of B4. Every reducible case leaves at least one entire spherical component untouched. This observation alone does not exclude later couplings to additional walls.

For the existing D4 realization, the B4 seed can be obtained by retaining the coupled outer root, the central root, one uncoupled outer root, and the difference of the two uncoupled outer roots. The last root has norm 4 and is a reflection bisector. The bare B4 chamber is one half of the D4 seed; this does not decide which extra polar walls are compatible with the half-chamber. That separate experiment is the next construction step.

No claim is made about changing the H5 source, changing projected directions, adding multiple walls over a source facet, increasing the support rank, or completing one of these seeds by other walls. In particular, this result does not classify arbitrary compact H9 Coxeter polytopes.
