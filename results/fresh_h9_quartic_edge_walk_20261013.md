# Three finite-edge steps in the quartic H9 host

**Subsequent complete census:** [the first-endpoint report](fresh_h9_quartic_first_endpoint_20261013.md) upgrades the first two steps below to genuine nearest-mirror steps of the full host chamber. At the third edge it finds a previously missed norm-(2+T) mirror before the norm-one proposal, giving a corrected alternative 186-wall cone. The bounded unit-grid artifacts below remain valid subgroup constructions; their third endpoint is not the nearest one.

The construction adds three coherent four-mirror symmetry orbits to the previous 174-wall cone, producing **186 walls in H9 and 76 walls in its forced H5 face**. Each step cuts the preceding H5 escape and creates an actual finite vertex. Both final cones still have exact strict spacelike escapes. This is a constructive extension, not a finite-volume example or a full-host obstruction.

The host is unchanged: over $K=\mathbb Q(T)$, $T=2\cos(\pi/8)$, $T^4-4T^2+2=0$, use
$\langle-(T-1)\rangle+F_4+F_4+\langle1\rangle$.
The starting lattice, complete necessary primitive norm pool, and earlier 174-wall certificate are documented in [the preceding report](fresh_h9_quartic_mixed_20261013.md). The present search changes the construction method to finite-edge equations and explicit relative-unit orbits. It makes no new-lattice or new-group claim.

| Step | Spherical H5 edge | Finite endpoint diagram | New norm | Height $\alpha x^2/q$ | H9 / H5 walls |
|---|---|---|---|---:|---:|
| 1 | $B_2+A_1+A_1$ | $B_4+A_1$ | 2 | 1012.7604 | 178 / 72 |
| 2 | $B_3+A_1$ | $F_4+A_1$ | 2 | 2268.4100 | 182 / 74 |
| 3 | $F_4$ | $F_4+A_1$ | 1 | 3174864.1373 | 186 / 76 |

The authoritative chain uses case 0 in each of [step 1](fresh_h9_quartic_edge_walk_20261013_profiles.json), [step 2](fresh_h9_quartic_edge_walk_20261013_next.json), and [step 3](fresh_h9_quartic_edge_walk_20261013_f4.json). Every previous wall is retained. The four copies come from reversing an F4 block and interchanging the two F4 blocks; only two copies restrict to the chosen H5 face. The saved positive pairings with the successive old H5 escapes are
$-1-2T+T^2/2+T^3/2$, $T-T^2/2$, and $9-11T-2T^2+3T^3$, respectively.

For a representative of the last orbit, the time coordinate is $-200-108T+342T^2+185T^3$, the first F4 coordinates are
$(-232-124T+396T^2+214T^3,\;-408-218T+696T^2+376T^3,\;-412-226T+704T^2+382T^3,\;-248-136T+424T^2+230T^3)$,
the second block is zero, and the scalar coordinate is $-1+2T^2+T^3$.
Its exact norm is 1. This last step introduces a scalar component; the first two new orbits have scalar coordinate zero.

## Exact construction and bounded scope

For each rank-four spherical edge, prescribing its finite pairings determines a component $r_0$ in the edge span. The orthogonal complement is a binary Lorentzian plane. Its exact projection $P$ and skew-adjoint operator $J$ satisfy
$P^2=P$, $JP=J$, and $J^2=\alpha P$, with $\alpha=T-1$.
The explicit identity
\[
y_0=-1+2T^2+T^3=\frac{T}{2-T},\qquad
s_0=2T^2+T^3=\frac2{2-T},\qquad
y_0^2-\alpha s_0^2=1
\]
gives an isometry $H=I+(y_0-1)P+s_0J$. For residual norms that are squares in K, the search tests $r_0\pm\sqrt{q-r_0^2}\,H^k v_0$ for an explicitly reconstructed norm-one vector $v_0$.

The complete finite-pairing census uses all five necessary primitive norm classes and every permitted finite label in this field, then tests positive definiteness at every embedding. There are 18, 8, and 5 finite profiles at the three edges; 15, 6, and 3 have square residual norms. The stated unit grids are $-4\leq k\leq4$ at the first edge and $-2\leq k\leq2$ at the next two. They yield 22, 12, and 4 primitive integral candidate records and 10, 5, and 2 distinct roots compatible with the respective old cones. These are complete counts for these grids only. Nonsquare residual classes and other orbits of the binary norm equation remain open.

An earlier calibration first squared this unit to obtain an integral ambient isometry, then used its rational half-step. The [integral-translation artifact](fresh_h9_quartic_edge_walk_20261013.json) and [half-step artifact](fresh_h9_quartic_edge_walk_20261013_halfstep.json) give alternative 178-wall cones from the original 174 walls. Their larger first caps are not merged into the authoritative chain: one conflicts with the smaller finite-coupling cap. A rational isometry is not asserted to preserve the lattice; each resulting root is separately tested for integrality, primitivity, and integral reflection coefficients.

## Geometry and independent replay

The separate Fraction implementation reconstructs every finite profile, unit-grid root, and saved dictionary. It checks all embeddings and root/reflection integrality, all finite-angle or dotted pair gates, strict timelike interiors, the exact spacelike escapes, and the physical Gram matrices of the displayed finite vertices. The final H9 and H5 checks cover 17,205 and 2,850 pairs. The positive definite incident ranks are nine and five, and the finite points satisfy every wall inequality; these are actual vertices of the displayed cones.

The final H5 escape has exact norm
\[
\frac{-39034+21393T+73052T^2-39610T^3}{2209}>0.
\]
The original H9 escape also survives, with norm
\[
-\frac{56876}{289}+\frac{34263}{289}T+\frac{201861}{578}T^2-\frac{55639}{289}T^3>0.
\]
Both have independently checked strict spacelike perturbations in the future component. Thus the final displayed cones have infinite volume. The next H5 escape exposes another $B_3+A_1$ edge involving the newly introduced scalar cap; it has not been completed here. No conclusion about full reflectivity of the arithmetic host follows.

The [profile and geometry audit](fresh_h9_quartic_edge_walk_20261013_profiles_independent.json) and [initial isometry audit](fresh_h9_quartic_edge_walk_20261013_independent.json) both pass. Replay from the repository root:

```sh
python scripts/fresh_h9_quartic_edge_walk_20261013_independent.py --verify
python scripts/fresh_h9_quartic_edge_walk_20261013_profiles_independent.py --verify
```

Producers are [the first edge construction](../scripts/fresh_h9_quartic_edge_walk_20261013.py), [its finite profiles](../scripts/fresh_h9_quartic_edge_walk_20261013_profiles.py), and [the next two edge constructions](../scripts/fresh_h9_quartic_edge_walk_20261013_next.py). The last accepts `--f4`. These numerical geometry proposals are always followed by exact reconstruction; the independent replay imports none of these producers.

## Relationship to Allcock's edgewalking algorithm

[Allcock's paper, Sections 1 and 9](https://arxiv.org/html/2110.03784v1#S9), proves its algorithm for lattices over Z. It finds the first admissible residual root in each constraint coset, orders candidates by their intersection along the edge, and exhausts a frontier of incident edges. Its binary continued-fraction algorithm does not directly apply to our quartic ring of integers.

This implementation uses the same useful geometric reduction to spherical extensions and a binary orthogonal complement, but is a bounded search for compatible caps. It does not enumerate all residual cosets or all seed orbits, prove that its chosen relative unit generates the integral plane group, or prove that any proposed endpoint is nearest in the full arithmetic chamber. It sorts the saved proposals by height from the original basepoint, which is not an along-edge ordering. The verified finite vertices belong to the displayed reflection-subgroup cones; their identification with corners of the full host chamber is unproved.

A meaningful next step is a complete residual-coset census over a bounded interval between a certified starting corner and a proposed endpoint, including the nonsquare residual classes and all conjugate-embedding bounds. Such a certificate could establish the first endpoint without asserting a generalization of the published integer-lattice algorithm. The starting corner's complete root stabilizer must also be certified. An explicit frontier would then track which outgoing edges remain unresolved.
