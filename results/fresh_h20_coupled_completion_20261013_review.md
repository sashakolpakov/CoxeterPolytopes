# Source review of the H20 coupled completion

**PASS.** I independently reviewed the [104-wall exact source](../scripts/fresh_h20_coupled_completion_20261013_exact.py), its frame construction, and the [105th-cap source](../scripts/fresh_h20_coupled_completion_20261013_dotted_exact.py). Both frozen replays passed:

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_coupled_completion_20261013_exact.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_coupled_completion_20261013_dotted_exact.py --verify
```

This is an independent mathematical/source review and fresh replay of the existing exact arithmetic, not a second all-root implementation. The previously reviewed 98-wall package is an inherited, hash-anchored dependency. Neither replay runs a numerical search.

The six cusp caps have the form `v+5 delta`, with `delta=r18+r96` future null and `v` perpendicular to it. Their norms are therefore exactly the positive norms of their base vectors. The two highest-root replacements complete the two `D8` components, and the four signed residual directions complete two further affine `A1` factors alongside the original parallel pair. All 603 new pair checks, 642 inherited-cusp strict signs and common-center inequalities pass.

The bounded-link conclusion uses more than rank. Each of the five components has a positive-definite principal minor, a strictly positive kernel and weighted vector sum `mu*delta` with `mu>0`. The lists exhaust every incident wall. All cross-component pairings are nonpositive; pairing a positive kernel sum with any incident wall forces every cross-component pairing to vanish. The Euclidean component ranks therefore add to `8+8+1+1+1=19`, the entire H20 horosphere dimension. A recession vector pairs nonpositively with every normal; the positive relations force all those pairings to be zero, and full quotient rank forces the recession vector itself to vanish. The strict chamber center supplies nonempty horospherical interior. Thus the displayed cusp link is bounded.

For the next cap, the three norm-two leaves and the scaled central root have exact `D4` Gram matrix. With `u=3p+sum(leaves)` and `v=center-p`, every leaf pairing is `−1` for `u+bv`, independently of `b`, because `(v,leaf)=0`. The exact quadratic `A b²+B b+C=0` enforces norm two; the selected negative-square-root branch is isolated in `(119/100,6/5)`. All 104 new pairs pass: 13 right, eight `pi/3`, and 83 dotted. The two intended anchor relations are explicitly dotted. The old strict escape is checked against all 104 retained walls and then cut strictly.

The 105th wall is strict at all 107 inherited cusp rays. It also preserves the newly completed cusp: its pairing with `r18` is `−1` and its pairing with `r96` is nonpositive, hence its pairing with `delta=r18+r96` is strictly negative. Common-center projection and the nonpositive pairings prove that the added inequalities are actual facets.

The 104-wall chamber retains the old strict escape of square `698631401/1887436800`. After the 105th cap cuts it, a new exactly reconstructed strict future escape has square `135725357/41943040 - (49017/40960)sqrt(6)>0`. Both chambers therefore have infinite volume. The audit does not assert a complete cusp or vertex census.

## Later orbit and second cutter

I also source-reviewed and freshly replayed the [175-wall orbit extension](../scripts/fresh_h20_coupled_completion_20261013_orbit_exact.py). Its 72 permutations form an explicitly checked group, and the matrices use the correct row-coordinate action and preserve the full form. Preservation of the old wall Gram determines their action on every old wall because the first twenty normals are a basis. The selected 70 further copies are checked directly against all retained walls and each other, so compatibility does not rely on the entire 104-wall chamber having these symmetries. The saved escape is reconstructed from a rank-20 kernel over the cap field and strictly perturbed. The full 9,765-pair replay **passed**.

The coefficient-encoding uniqueness test is valid for the saved field: its radicand is `D=5291/4-(2299/6)sqrt(6)`. I independently computed `Norm(D)=41700835/48`, which is not a rational square. Hence `D` is not a square in `Q(sqrt(6))`, and the two-coefficient representation is injective. All orbit vectors have the same nonzero center pairing, so distinct vectors cannot describe opposite orientations of one mirror.

The [176th-cap formula](../scripts/fresh_h20_coupled_completion_20261013_second_exact.py) and its nested field embeddings also pass review. The three old normalized roots are mutually orthogonal, the weighted sum `u` has square seven, and the equation for `u+t*p` correctly has constant term six when the target norm is one. The positive-square-root branch is isolated in `(353/100,177/50)`. The prescribed branch, norm, strict center, all 175 new pairs, 107 inherited-cusp signs and strict target cut are checked exactly.

I ran the additional command

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_coupled_completion_20261013_second_exact.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_coupled_completion_20261013_orbit_exact.py --verify
```

Both **passed**. The final cap has 12 right, four `pi/4`, two `pi/6` and 157 dotted relations. The final strict future escape has square `1591730805825/162076295168 - (299252657/79138816)sqrt(6)>0`, so the 176-wall chamber still has infinite volume. These replays do not establish completeness of any wider cap search.

## Obstruction to adding only dotted walls

The [18-neighbor facet obstruction](../scripts/fresh_h20_coupled_completion_20261013_facet_obstruction.py) also passes source review and fresh `--verify`. Wall 175 has exactly 18 intersecting neighbors; every other retained wall is strictly dotted to it. In its orthogonal complement, which has signature `(19,1)` and dimension 20, project any proposed completion's strict center to a timelike point `c`. The conditions of orthogonality to `c` and the 18 finite-neighbor projections leave a nonzero spacelike vector `w`.

At the future null point `c+sqrt(-c²/w²)w`, all 18 neighbor inequalities remain strict. Every dotted wall projects to a timelike normal with the same future orientation, determined by its negative pairing with `c`; hence its inequality is also strict at that null point. This applies to any finite collection of additional walls dotted to wall 175. Moving slightly outward along the line and then inward across wall 175 gives a strict future spacelike escape. The argument uses the proposed completion's own strict center, not an assumption that added walls retain the archived center.

Thus retaining these 176 walls and adding only walls strictly dotted to wall 175 cannot yield finite volume. Added walls intersecting or parallel to that facet, and coherent replacements, remain outside this obstruction. The proof is not limited by a coefficient or angle bound.
