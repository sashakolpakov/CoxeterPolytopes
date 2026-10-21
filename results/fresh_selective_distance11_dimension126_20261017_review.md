# Independent review of selective distance 11 and dimension 126

**Mathematical and exact-arithmetic PASS.** No unresolved geometric case or arithmetic discrepancy was identified in the [producer proof](fresh_selective_distance11_dimension126_20261017.md). Its conclusion is dimension at most **126** for every finite-sided, finite-volume real hyperbolic Coxeter polytope, without arithmeticity, compactness or cusp restrictions. The separate compact-facet conclusion is ambient dimension at most **95**.

The [independent checker](../scripts/fresh_selective_distance11_dimension126_20261017_review.py) and [certificate](fresh_selective_distance11_dimension126_20261017_review.json) use only standard-library `Fraction` arithmetic and literal graphs. They import no producer arithmetic, SymPy, numerical linear algebra or optimizer. The finite reconstruction and the mathematical review below are distinct checks. This is internal review, not external refereeing or proof-assistant formalization; no priority claim is made.

## Complete local domain

The inherited actual-face argument gives a positive-definite supporting white diagram W, positive-definite W plus any one side root, and positive-semidefinite W plus an adjacent pair. The two opposite closed sides are disjoint even when the quadrilateral has ideal vertices. Hence each opposite-side diagram is indefinite, and its minimal indefinite core contains both black roots with both black deletions positive definite.

The reviewer reconstructs all 384 marked contexts in the [frozen order-6 through order-10 catalogue](fresh_high_rank_quasilanner_census_20261017.json), including their graph distances and all relevant inverse diagonals. The maximum marked radius is exactly 7. A core of order at most 5 has radius at most 3 by its size; it cannot participate in a disjoint join requiring radius sum at least 11 or an overlap requiring radius sum at least 12. Thus its potentially unbounded finite labels introduce no omitted case here.

At distance at least 4 there are 200 tested ports, of which 55 survive the inverse test. For a label-4 edge, the checker independently constructs a positive diagonal square-root-of-two congruence to a rational Gram, then returns inverse diagonals to the norm-two gauge. An external supporting white has nonpositive contacts with nonzero magnitude at least 1. Since the spherical inverse is entrywise nonnegative, inverse diagonal at least 2 contradicts the required positive definiteness. Equality 2 is correctly excluded: this test concerns an actual edge support, not merely a semidefinite ideal corner.

## Positive witnesses and all disjoint joins

The reviewer independently generates a nonnegative vector whenever a required positive-definite Z-matrix first has a nonpositive elimination pivot. It does so by solving against the positive leading block, rather than importing the producer's accumulated basis transformations. The vector and its quadratic value are checked exactly.

This sign is essential to the proof. Restoring additional nonpositive contacts, or increasing a tested cross-bond magnitude, cannot increase the vector's quadratic value. The rejection therefore applies to the actual Gram matrix, not just its minimal-contact completion. The checker also independently verifies all **204** saved producer witnesses: 119 disjoint-white obstructions, 40 overlapping edge-support obstructions, and 45 additional-cross-edge obstructions.

The 55 allowed ports give 425 ordered marked disjoint joins. Of these, 119 have a nonnegative nonpositive-square witness already in W. The remaining 306 correspond to 306 different ordered marked-context pairs; thus each such pair has exactly one possible cross edge after the short-path, inverse and white-witness exclusions.

For every remaining join the reviewer checks the producer's two wholly white three-node tails directly in the literal native diagrams. With positive tail coefficients (1,2/3,1/3), each tail has square at most 4/3, also when a native edge is labelled 4. A cross bond of magnitude at least sqrt(2) makes their combined square negative. Consequently the only possible nonright contact is simple; parallel and dotted contacts are excluded too. No label cutoff is used.

As a separate check, the reviewer eliminates the full white block on each side and computes the exact upper bound on the squared cross contact; every bound is at most 2 and the positive-definite inequality is strict. Both arguments independently force the same simple bond.

The reviewer reconstructs every resulting Gram and its inertia by rational one-by-one or two-by-two symmetric pivots. All 306 joins have at least two negative directions. It also compares all 13 underlying saved Gram/inertia types. Since every other cross edge is already excluded, these are the actual core Grams, contradicting ambient negative index one.

## Exhaustive overlaps and unknown cross contacts

Every shared white root has its two marked distances summing to at least 12. Both distances are at least 5, and each marked core has at most three such roots. The reviewer enumerates subsets of the left and right node sets and all bijections between them, independently of the producer's filtering of combinations of candidate pairs. It obtains 2,137 candidates and 1,010 injective, Gram-consistent overlaps:

| Outcome | Count |
|---|---:|
| External-neighbor inverse obstruction | 856 |
| Positive witness in W plus one black root | 40 |
| L25/L25 shared far vertex | 60 |
| L25/L25 reversed terminal edge | 27 |
| L25/L25 reversed terminal three-node path | 9 |
| L22/L25 shared far vertex, both orientations | 18 |

The matching is the complete intersection. Each external neighbor used by an inverse rejection is therefore genuinely outside the other core, and its large marked distance ensures it is white. The independently recovered categories agree with every saved producer record.

For each of the 114 surviving contexts, the reviewer tests every unspecified edge between the exclusive parts. Exactly 8,280 tests create a path of length at most 11; the other 45 have explicit nonnegative nonpositive-square witnesses in W with one black root. These witnesses remain valid for stronger or additional acute contacts. Thus all surviving core Grams are fixed before their isolation arguments are applied.

## Supporting-white extensions

The reviewer reconstructs all 114 preferred corners from scratch, including all finite inverse diagonals, affine positive kernels, possible contacts and marked components. It additionally checks positive semidefiniteness of all four adjacent corners in every surviving context.

For the 60 shared-vertex cases, the 19-root Gram has inertia (17,1,1). Its radical is an actual linear relation: its Gram-null root combination has square zero and is orthogonal to a timelike vector in the core span. Positive definiteness perpendicular to that timelike vector makes the combination vanish. Every allowed new contact has a nonzero radical coordinate. The two signs are within distance at most 3 of their respective marked pairs. Acuteness and radical balance force both signs whenever a new contact exists, yielding a path of length at most 8. Thus all further supporting whites are orthogonal, and the marked component has order at most 17.

For the 27 shared-edge cases, the 18-root core is nonsingular Lorentzian. Preferred-corner restrictions allow contacts only in the two end sets of four nodes. The full inverse blocks on those end sets have strictly positive entries and diagonal minimum 4. A contact at both ends gives a path of length at most 6; a nonzero contact vector supported at one end has squared projection at least 4. The orthogonal complement of the core span is positive definite, so a square-two root cannot have the required negative residual. Every further supporting white is orthogonal, and the marked component has order at most 16.

The nine three-node overlaps require a different argument: they do admit an extension. The reviewer independently recomputes the [separate extension lemma](fresh_distance11_extension_review_20261017.md). The one-end inverse block is

```
[[7/2,4,2],[4,6,3],[2,3,2]].
```

Only a sole contact -1 at node 3 or its opposite-copy counterpart is possible. The projection has square 2, forcing zero orthogonal residual and a unique actual root in the core span. The two possible roots have mutual product -2 and cannot both lie in positive-definite W. Every other white orthogonal to the core is automatically orthogonal to the extension. All 18 one-extension corner variants are independently checked; the marked component remains isolated and has order at most 16.

The 18 mixed cases are handled by the [mixed radical lemma](fresh_distance11_extension_review_20261017_mixed.md). The reviewer reconstructs both orientations, including their literal isomorphisms to the frozen lemma. The core has inertia (16,1,1). Its only zero-radical node is excluded by the preferred-corner contact tests, and any nonzero contact forces a short path between opposite signs. Consequently every supporting white is orthogonal and the marked component has order at most 16.

At finite corners, the isolated component is an actual ambient vertex component. At ideal corners, the inherited actual factor-fiber statement puts the entire varying affine factor inside W plus the two side roots. A finite marked component cannot be that factor, and an isolated affine marked component is already the whole factor. Additional incident cusp walls outside this star cannot enlarge the component relevant to the payment. This is why the proof permits nonsimple cusps.

## Payment and the dimension inequalities

The selective rule pays a pair when its distance is at most 11 or it belongs to a connected corner component of order at most 17. The preceding argument pays every quadrilateral. The inherited triangle argument provides two paths of length at most 8. Hence every k-gonal two-face receives at least 5-k.

Small components contribute at most 8 times their order in pairs; the inherited all-rank distance counts give the required bound on larger components. Finite vertices cost at most 11n, and an ideal vertex with factor ranks r_i and link vertex count V costs at most `V sum min(11,r_i/2)`.

The proved convex-order transfer is applied separately below and above cap 22. These transfers preserve the charge, do not increase curvature, and terminate by increasing the sum of squared ranks. No transfer across the cap is presumed. The resulting profiles have copies of 22, at most one smaller part and at most one larger part.

The reviewer independently enumerates all 1,095 reduced profiles for link dimension 126 through 135. Integer convolution and a separate inclusion-exclusion moment formula agree exactly. Every margin is positive; the least is **40603/559682**, at ranks (22,22,22,22,22,16). The finite-vertex slack is positive in ambient dimensions 127 through 136, and the global shape inequality retains the strictly positive term 5(n-1). The [completed dimension-136 theorem](fresh_selective_distance12_dimension136_20261017.md) supplies the higher-dimensional tail. This proves 126 with no new restriction.

At ambient dimension 126, the profile (22,22,22,22,22,15) has margin -238535/12872686; the reviewer confirms this boundary. The same scalar argument does not establish 125. For the compact-facet branch, the independent parity polynomials have coefficients (12,192,4) and (104,196,4), giving ambient dimension at most 95.

## Final validation scope

Final independent generation and fresh saved replay 53606 both passed. The primary source, data and report are bound at their final frozen hashes. Python syntax, all ten direct source bindings, all 25 producer source bindings and seven report links passed. This review source, certificate and report are final/frozen. The published classification, actual incidence statements, all-rank graph budgets, convex-order theorem and shape inequality remain explicit inherited mathematical inputs. They are reviewed geometric arguments, not consequences of the finite verifier alone. The completed package records a complete internally checked computer-assisted proof; it does not assert formal verification or external refereeing.
