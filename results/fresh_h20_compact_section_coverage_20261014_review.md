# Review of the compact-section G2 gate and ten exact A2 samples

**PASS, within the stated scope.** I found no issue in the [coverage report](fresh_h20_compact_section_coverage_20261014.md), its compact-only orthogonal-short-root G2 obstruction, or the actual-vertex proof used for the 162 sampled A2 constructions. The [small gate checker](../scripts/fresh_h20_compact_section_coverage_20261014_gates.py) passed a fresh saved-result replay. This review does not classify all compact sections or repeat all 162 constructions.

## Actual facets and the arbitrary-ratio G2 argument

The proof starts with an **actual** G2 ridge of a compact Coxeter section, represented by two orthogonal old norm-two roots. Its dimension-18 face is compact. A facet of that face meets it in the hyperbolic interior, so the incident three-wall diagram is positive definite. A finite spherical diagram cannot extend G2 by a connected third node. Thus every face boundary is orthogonal to both projected G2 normals, and its intrinsic normal receives no additional projection correction.

For contacts `a,b` and a modified boundary contact `x`, this gives old pairings `xa,xb`. Both are nonzero. Adding back the section-normal rank-one term makes the old three-root Gram positive definite. I independently rebuilt the literal source table over all 756 orthogonal short-root pairs: its only positive-definite row with both pairings nonzero is `(norm,pair,pair)=(2,-1,-1)`. Hence a modified boundary forces `a=b`. This implication does not assume any particular contact ratio, edge position, or other contact support. If no modified boundary exists, the face already has only unchanged rational boundaries.

In the equal-contact case the label-six equation gives `a²=4√3−6`. The modified norm is `N=3/2−√3/3`. Two modified boundaries are divergent, since their negative pairing has magnitude at least `1/a²>N`. A nonright unchanged-short contact has squared cosine `(9+2√3)/23`, strictly between `1/2` and `5/8`, below the next possible value `cos²(π/5)`. A nonright unchanged-long contact is divergent because its squared cosine is `1/N>1`. Therefore a compact modified facet has only unchanged rational boundaries, all perpendicular to its normal. These are statements about boundaries that actually meet the compact face or facet, not raw projected inequalities whose actuality is unknown.

## Why the rational face gives a contradiction

The resulting compact face has dimension 18 if no modified boundary exists, or dimension 17 after choosing a modified facet. Its intrinsic wall normals are literal rational source roots; orthogonality at each face restriction is essential here. Compactness forces those normals to span the face's full Lorentzian vector space. Their rational span therefore defines a rational indefinite form, and its intersection with the old integer lattice is a full lattice.

Each boundary reflection preserves this space and lattice: the source reflection coefficients `2(Qr)/r²` are integral, independently rechecked for all 210 roots. The compact reflection group is consequently a lattice subgroup of the integral orthogonal group. Finite covolume gives finite index in that discrete supergroup. Meyer's isotropy theorem applies in these ranks and makes the integral orthogonal group nonuniform; its finite-index subgroup cannot be cocompact. This proves the contradiction. Merely having rational coordinates, without the invariant lattice and full rational span, would not suffice; the report supplies both premises.

## Sparse support and endpoint scope

At a finite section vertex, choosing one source representative for each of its 20 actual simple facets gives positive-definite `P`, and the old Gram `G=P+aaᵀ` is positive definite too. The old source graph is therefore a crystallographic forest. Three contacts with the same sign give a forbidden negative Coxeter triangle; two of each sign would require an old four-cycle. This proves support at most three. For three contacts, the old graph on them is a star, and at least one old star edge must cancel. The argument does not assume that every old incident inequality survives as a facet.

The report correctly separates vertices in an open old edge from sections passing through an old finite vertex, where restriction can create redundancies or coincident inequalities. It does not claim that a Coxeter-preserving deformation removes those endpoint cases. The G2 actual-face obstruction itself does not require an open-edge parameterization.

## Separate replay of ten A2 witnesses

The [review checker](../scripts/fresh_h20_compact_section_coverage_20261014_review.py) and [certificate](fresh_h20_compact_section_coverage_20261014_review.json) use separate Fraction sparse-radical arithmetic, with rigorous 180-bit rational square-root enclosures for strict signs. They import neither the producer nor its arithmetic engines. SymPy is used for rational matrix inverses and positive-definite decompositions.

The selected cases are `0,27,52,55,80,107,108,135,160,60`. Together they cover all nine source edges, all three contact ratios, all three sampled positions, both residual signs, and all three defect types. The replay checks:

- Unit normal identities, initial timelike section points, and the positive projected-Gram Schur complement.
- The released dual's exact 20 pairings and nonzero pairing with the entering wall, independently reconstructing the saved next point.
- Negative norm, future orientation, the exact 20-wall zero set, and all 2,100 source-wall signs at the ten next points.
- Projected independence and thus actual positive-definite rank-20 vertices, followed by each claimed obtuse or forbidden-angle ridge.

The full producer's source implements these same geometric implications for all 162 cases, but only the ten selected witnesses were independently recomputed in this review. Source and data hashes are bound in the review certificate; no claim of a second exhaustive 162-case replay is made.

```sh
python scripts/fresh_h20_compact_section_coverage_20261014_gates.py --verify
python scripts/fresh_h20_compact_section_coverage_20261014_review.py --verify
```

Both saved-result replays passed. Unsampled continuous A2 profiles, B2 joins, other old source pairs, and general compact H20 sections remain outside this theorem.
