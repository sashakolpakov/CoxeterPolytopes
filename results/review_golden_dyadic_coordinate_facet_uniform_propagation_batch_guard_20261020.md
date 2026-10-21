# Outside review of the literal R_n coordinate-facet propagation

20 October 2026. Paper-only review: **PASS**, conditional on the stated
R5 base obstruction and imported uniform norm pool. No code body,
scientific import, generation/replay, Lean or TeX was executed.

Bound source:
`results/golden_dyadic_coordinate_facet_uniform_propagation_20261020.md`,
SHA256 `e7d593d3e7a3a63cf9a48c7668b73cc8db25da12346fc2ddd351e8ffa0edf732`.
I read it completely and reread the uniform pool and zero-germ argument
in `golden_dyadic_d4_extension_first_edge_20261020.md`, SHA256
`3fc21952a3d119db4370ca46e39aa607666eaea7eaa54d835c3163e7a281035b`.

For an ambient primitive mirror crossing H=e_n^perp, the determinant
of its two-normal Gram gives x^2<q. The other real positive form gives
x'^2<=q. The trace identity in the source is exact. Its smallest
nonzero value is 2; hence the q1 strip and the q4 strip after division
by 2 contain only zero. In the q2 strip, trace<4 leaves exactly
0, +/-1, +/-phi, +/-phi^-1 before individual embedding tests. The last
four fail one embedding, so the surviving coordinate list is precisely
0,+/-1. This uses arbitrary crossing mirrors, without substituting a
fundamental-angle alphabet.

At x=0 the ambient reflection preserves the exact perpendicular
summand R_(n-1). Its normal remains primitive there. At x=+/-1,
q=2 and the projected normal is an integral unit vector of norm 1,
so it is primitive and reflective in R_(n-1). These prove containment
of every induced mirror in the full lower-rank arrangement. Conversely
any primitive reflecting normal of R_(n-1) extends by zero coordinate,
and its reflection preserves the direct orthogonal O*e_n summand.
Thus equality of FULL arrangements, not just prefix containment, is
proved.

The imported full B_n+F4 zero germ makes H an actual facet. Its
relative interior is consequently one chamber of that induced full
arrangement: any crossing mirror would otherwise cut the ambient
chamber interior near H; mirrors disjoint from H cannot add an
intrinsic boundary. The induced zero germ is B_(n-1)+F4, fixing the
claimed chamber. This argument applies at every step n>=6, down to
the literal R5 chamber.

Compactness therefore propagates downward, so a noncompact full R5
chamber excludes compact full R_n chambers for every n>=5. The stated
noncofinite consequence uses a separate arithmetic input, correctly
identified in the source: these admissible forms are K-anisotropic
because the other embedding is positive definite, and their integral
orthogonal groups are uniform arithmetic lattices. A cofinite
reflection subgroup would have finite index and be cocompact; its
full reflection chamber would then be compact. Geometry of a facet
alone would not establish that stronger implication.

No flaw or correction was found. The base e2aeff93 candidate has my
separate paper PASS review but its fixed execution status remains
separate. This record is conditional on that base conclusion and
concerns only the displayed literal lattices, not arbitrary lattices
in their similarity or commensurability classes and not a uniform
dimension bound.
