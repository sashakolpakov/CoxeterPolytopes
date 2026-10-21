# The complete affine-ray family cannot repair the fixed lower H22 cusp

The all-real-shift search is closed for one precise family. **No finite addition from the complete transported affine-ray families can make the lowered H22 chamber finite-volume while retaining all 210 original facets and the cap c₄, with every new distinct wall acute-compatible with those retained facets.** The only compatible incident candidates are repeated retained mirrors. This excludes changing the shifts and choosing arbitrary subsets of this family as a repair; it does not exclude replacing retained walls or using other normal directions.

The [exact certificate](fresh_h22_lower_incident_family_20261013.json) and [replay script](../scripts/fresh_h22_lower_incident_family_20261013.py) use the complete, previously verified 120,960-frame orbit. This is an algebraic incidence reduction and exact enumeration, without a parameter grid, numerical tolerance, or upper bound on Coxeter labels. The [earlier cusp diagnosis](fresh_h22_lower_new_cusp_20261013.md) supplies the complete current-wall incidence and the explicit four-dimensional lineality space.

The [independent review](fresh_h22_lower_incident_family_20261013_review.md) includes a fresh successful replay and a separate standard-library `Fraction` check of all three final rejection witnesses; no mathematical or scope issue was found.

## Family and fixed cusp

Let Q be the original 22-coordinate Lorentz form, wᵢ its original facet normals, and let the additional coordinate be represented as its coefficient of √11. The full form in these coordinates is diag(Q,11). The fixed lowered cap and deficient cusp are

```
c₄ = (D₀/6 − w₃₂, 2),
δ₄ = c₄ + w₂₆.
```

Both c₄ and w₂₆ have norm 2 and pairing −2. The future null ray δ₄ has exactly 17 original incident facets and c₄; its Euclidean normal rank is 17 in a 21-dimensional horosphere. Every other facet in the frozen lowered chamber is strict there.

For each of the 120,960 transported cusp frames g, the tested families consist of the seven transported negative highest roots and the four signed residual tangent bases. With each base written as b=B/d, they are

```
r(s) = (B/d,0) + s(Dg/6, σ·2),     s ∈ R, σ ∈ {−1,+1}.
```

Both overall normal orientations are allowed. The seven highest-root bases have d=1 and square 2 or 4. The first residual bases have d=28 and square 279/14; the second have d=5208 and square 55/62496. The script independently checks every base square, every base–cusp orthogonality, and every cusp square. The completeness and transport of the frame/highest-root arrays remain dependencies of the hash-anchored [complete orbit certificate](fresh_h22_new_cusp_completion_20261013_verify.json).

## Reducing all real parameters to exact candidates

A genuinely new incident wall acute-compatible with c₄ and w₂₆ must be perpendicular to both: its pairing with their sum δ₄ is zero and both summands are nonpositive. Positive-proportional copies of c₄ or w₂₆ are trivially redundant and are excluded before this reduction; their self-pairings are not treated as angle conflicts.

Put C=D₀−6w₃₂ and

```
U = B·Q·w₂₆,          V = B·Q·C,
L = Dg·Q·w₂₆,         M = Dg·Q·C + 1584σ.
```

The two incidence equations are

```
6U + sdL = 0,         6V + sdM = 0.
```

For both signs, the exact scan proves that L and M never vanish simultaneously. Thus there is at most one real parameter for each signed base: the cross condition is UM=VL, and the parameter is rational, given by −6U/(dL) or −6V/(dM). No negative or zero shifts are discarded.

The center h is a positive scalar multiple of the sum of all original normals. For any nonzero candidate acute-compatible with all originals, (r,h)<0. Equality would force all original pairings to vanish; their exact span has rank 22, so the old coordinates would vanish, and orthogonality to c₄ would then force the added coordinate to vanish too. The code therefore rejects only zero-center candidates and selects the unique possible orientation for every other candidate, covering both orientations of every mirror.

## Local angle gate and counts

The 17 original incident walls all have square 2. Their pairings with every reduced candidate and the candidate square are rational. Since both normals are perpendicular to δ₄, their two-dimensional Gram is positive semidefinite; a dotted pair is impossible. If the angle is π/m, rationality implies that 4cos²(π/m)−2 is a rational algebraic integer and hence an integer. Consequently all finite labels reduce exactly to m=2,3,4,6, with the parallel endpoint included separately.

After requiring the actual pairings to be rational, the allowed pairings are:

| Candidate square | Allowed pairings with an incident norm-two wall |
|---|---|
|2|0, −1, −2|
|4|0, −2|
|279/14|0|
|55/62496|0|

The last value −2 in the square-two row includes parallel mirrors. Positive-proportional copies of existing old walls are detected separately before applying this angle gate, so a duplicate is not rejected merely for having a positive self-pairing.

| Base family and height sign | Incidence solutions | Zero-center exclusions | Redundant old-wall copies | New literals passing all 17 gates |
|---|---:|---:|---:|---:|
|Highest roots, +|367,336|210|92,160|1|
|Highest roots, −|367,269|162|92,160|0|
|First residual +, +|48,968|32|0|1|
|First residual +, −|48,965|32|0|0|
|First residual −, +|48,968|32|0|1|
|First residual −, −|48,965|32|0|0|
|Second residual +, +|48,964|1|0|1|
|Second residual +, −|48,975|2|0|0|
|Second residual −, +|48,964|1|0|1|
|Second residual −, −|48,975|2|0|0|

In total there are 1,126,349 incidence solutions, 506 zero-center exclusions, and 184,320 literal old-wall copies. The copies represent 16 existing incident walls, each repeated 11,520 times. The five remaining new literals represent three distinct oriented mirrors.

## The final three mirrors fail retained-facet compatibility

All three occur at frame 54,518. Their exact coordinates are saved in the certificate.

| Base | Entry | Shift | Square | Center pairing | Positive retained-wall pairing |
|---|---:|---:|---:|---:|---|
|Highest root|381,632|1/5|4|−32/5|(r,w₈)=1/5|
|First residual +|54,518|3/10|279/14|−108/5|(r,w₂)=333/140|
|Second residual +|54,518|−1/372|55/62496|−11/744|(r,w₀)=11/7812|

The negative residual bases with reversed overall orientation give the same two residual mirrors. The three displayed rational-coordinate vectors have pairwise nonsquare norm ratios in Q, checked exactly, so they are not proportional and define distinct mirrors. In particular, the negative-shift candidate is included. All three fail the necessary nonpositive pairing condition with a retained facet. There are no genuinely new compatible incident normals left.

## Why this excludes finite repairs from the family

Every acute-compatible addition has nonpositive pairing with δ₄=c₄+w₂₆, so it cannot remove this null ray. If it is incident, the exhaustive calculation above shows it is a redundant retained mirror. If it is nonincident, its pairing at δ₄ is strictly negative.

For any finite set of such additions, the certified local escape construction from the earlier diagnosis still works: choose a positive tangent vector v perpendicular to all incident normals and to h, set λ=v²/288, and take

```
X(ε)=δ₄+εv+λε²h.
```

For sufficiently small ε>0 this vector is future spacelike, all incident inequalities are strict, and all finitely many nonincident inequalities remain strict by continuity. Thus the enlarged finite wall system has infinite volume. The existing inherited and lowered nonincident caps also remain strict for sufficiently small ε.

The conclusion requires retained-facet acute compatibility with **all original 210 walls and c₄**. It is not a theorem about arbitrary extra inequalities that might render an obstructive old wall redundant, about deleting/replacing c₄ or w₂₆, or about other bases and cusp directions. It is not a global obstruction to H22 Coxeter polytopes.

Replay:

```sh
python scripts/fresh_h22_lower_incident_family_20261013.py --verify
```

The script recomputes the complete enumeration and compares the resulting certificate exactly. Native integer operations have explicit prebounds; final normals, norms, incidence and positive rejection witnesses use independent `Fraction` arithmetic.
