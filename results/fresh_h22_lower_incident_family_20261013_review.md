# Independent review of the all-shift H22 incident family

**PASS within the stated retained-facet scope.** I read the [source](../scripts/fresh_h22_lower_incident_family_20261013.py), [frozen certificate](fresh_h22_lower_incident_family_20261013.json), and [report](fresh_h22_lower_incident_family_20261013.md), then freshly ran:

```sh
python scripts/fresh_h22_lower_incident_family_20261013.py --verify
```

The replay passed. A separate short calculation using only standard-library `Fraction` arithmetic reconstructed the original form and independently checked the final three normals, their norms, center and incidence identities, every one of the seventeen incident-wall angle gates, and their positive rejection pairings. This is a source/proof review and deterministic census replay, with an independent final-three check; it is not a second implementation of the entire orbit census.

The coverage and geometric implications are sound:

- For a genuinely new acute-compatible mirror incident at `delta4=c4+w26`, both individual pairings must vanish. Trivial positive copies of `c4` and `w26` are redundant and are correctly excluded before imposing this condition.
- The two incidence equations and cross condition `UM=VL` include every real shift. The script proves that `L` and `M` never vanish together, uses the other denominator when either vanishes separately, and includes negative and zero shifts. All remaining parameters are rational. The complete transported arrays cover the seven highest-root and four signed tangent bases, all 120,960 frames, and both added-coordinate signs.
- Selecting the unique orientation with negative center pairing does not impose an extra completion assumption. The center is a positive multiple of the sum of all 210 retained original normals. Acuteness forces its pairing to be nonpositive. Equality forces the old component to vanish by the checked rank-22 span, then orthogonality to `c4` forces the added component to vanish. A nonzero admissible incident normal therefore has strictly negative center pairing.
- The incident-angle gate is complete for all labels. Candidate and old incident normals lie in the null ray's perpendicular space, whose quotient is positive definite; hence a dotted pair is impossible. Their squared cosine is rational, so the cyclotomic algebraic-integer condition leaves precisely labels 2, 3, 4 and 6, with the parallel endpoint included. The subsequent rational-square test correctly narrows the possible literal pairings further.
- Positive-proportional old-wall copies are recognized before their self-pairings could falsely reject them. Zero added coordinate forces zero shift for such a copy. The three final new rational-coordinate normals have pairwise nonsquare rational norm ratios, proving that they represent distinct mirrors.
- Native integer contractions, sums, products, sign changes and final cross multiplications have adequate checked bounds. The inputs exclude the signed-integer minimum before negation. Final normals and rejection witnesses are reconstructed with unbounded rational arithmetic. The full orbit and prior complete cusp-incidence proofs remain explicit source-hash dependencies.

The exact census contains 1,126,349 incidence solutions, 506 zero-center exclusions and 184,320 redundant literal copies of sixteen retained original walls. Five new literal entries reduce to three distinct mirrors. Their rejection witnesses independently reproduce as follows:

| New normal square | Retained wall index | Positive pairing |
|---|---:|---:|
|4|8|1/5|
|279/14|2|333/140|
|55/62496|0|11/7812|

There is therefore no genuinely new acute-compatible incident mirror in the declared family. The earlier exact cusp certificate supplies a positive tangent direction in the four-dimensional lineality space. For any **finite** selection of further nonincident walls, its local strict spacelike perturbation still satisfies every inequality for sufficiently small positive perturbation. This proves infinite volume for finite additions from this family to the frozen lowered chamber.

The scope restriction is essential: the original 210 facets and `c4` are retained, with nonpositive pairings required for genuinely new mirrors. A positive pairing is not used to exclude arbitrary inequalities that could make an original facet redundant. The result does not prohibit coherent wall replacement, different cap directions, or H22 Coxeter polytopes in general. No issue remains in the frozen result.
