# Interior sections, dihedral lattices, and lifts: 13 October 2026

No new compact or finite-volume Coxeter polytope has been established in a new dimension. This batch tested four construction mechanisms and produced exact obstructions, including one that excludes a continuous family of H20 sections. The arithmetic H9 searches remain bounded searches, not proofs that the full lattices are nonreflective.

The targets remain compact H9 and finite-volume H20 or a dimension above 21. The [Felikson–Tumarkin reference page](https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/polytopes.html), checked on 13 October 2026, lists known compact examples through dimension 8 and finite-volume examples in dimensions at most 19 and in dimension 21. None of the calculations below changes that list.

| Construction | Exact result | Scope |
|---|---|---|
| Interior sections of Borcherds's H21 chamber | Every wall-pair pencil has a forbidden actual ridge | All 21,945 pairs, every positive real parameter |
| H3+A6 arithmetic H9 lattices | 45 and 56 walls; both finite prefixes have infinite volume | Complete Vinberg searches through height 30 |
| Cubic-seven dihedral H9 lattices | 18, 25, and 18 walls; all finite prefixes have infinite volume | Complete Vinberg searches through height 30 |
| Single-coordinate lifts of the compact H8 chamber | Two possible nontrivial lifts, both retain an unavoidable ideal ray | No compact completion by any additional acute-compatible walls |

## H20: all sections supported on two wall normals fail

Let P be the frozen 210-wall Borcherds H21 chamber. For two distinct outward normals r_i,r_j, consider the hyperplane with normal

```
e(t) = r_i - t r_j,  t > 0.
```

Every such hyperplane meets the interior of P. Its section is a full-dimensional finite-volume H20 polytope: in the Klein model, a hyperplane section of the finite convex hull of finite and ideal vertices has only finite or ideal vertices. Finite triangulation into such hyperbolic simplices gives finite volume. The question is whether its actual dihedral angles are Coxeter angles.

The 21,945 unordered wall pairs reduce to nine orbits under exactly verified chamber symmetries. Initially the search tested the balanced section through a fixed interior point in each orbit. All nine failed. It then replaced the individual parameter values by exact polynomial calculations:

- For the six elliptic pair orbits, persistent triple faces give actual section ridges for every t>0. Their necessary Coxeter-angle polynomials have greatest common divisor 1 or a power of t. They have no common positive root.
- For two nonelliptic orbits, an actual ridge has squared cosine strictly between 1/4 and 1/2 for every t>0. No Coxeter angle has a squared cosine in this interval.
- For the last orbit, an actual ridge has squared cosine strictly between 0 and 1/4 on 0<t≤1. A verified chamber symmetry interchanges the two defining walls and sends t to 1/t, covering the rest of the parameter range.

All claims that a ridge is actual include a timelike point on precisely its defining hyperplanes with strict inequalities at the remaining walls. Polynomial coefficient signs certify these inequalities over the entire stated parameter intervals. The proof allows arbitrary real t; it is not limited to rational or arithmetic sections.

This excludes every interior section whose normal lies in the span of at most two original wall normals. Coefficients of the same sign cannot define a hyperplane through the interior, and opposite signs give exactly the pencils above. It does **not** exclude normals involving three or more walls, sections of different chambers, or other H20 constructions.

The [full section report](fresh_h20_global_20261013.md) links the certificates. The [independent Fraction-only verifier](../scripts/fresh_h20_global_fraction_verify_20261013.py) checks orbit coverage, all nine central failures, and all nine universal parameter obstructions. Its [saved replay](../data/enumerations/fresh_h20_global_fraction_verify_20261013.json) also rejects forged certificates. An exploratory 28-parameter algebraic census is supplementary; the universal result does not depend on that bounded census.

## Compact H9: change the discriminant with H3+A6

The golden-field search first checked H3³. Its two proposed quadratic spaces are similar over the field to earlier spaces. The displayed transports do not identify their integral reflection groups, and no global group-identification claim is made. The search therefore moved to

```
<-α> + H3 + A6,   α ∈ {φ, 2φ},   φ²=φ+1.
```

The determinant has odd valuation at the inert prime 7. This separates both spaces from the nine explicitly compared recent hosts; their determinant ratio also separates them from one another. It is a limited comparison, not a worldwide novelty claim.

The complete candidate norm alphabets, full zero-height chamber, bounded time controls, and spatial enumeration were established before running Vinberg's algorithm through height 30. The two runs give 45 and 56 walls, respectively. Every accepted wall is a primitive integral reflection and every wall pair passes the exact angle test. Both cones nevertheless have a strict timelike interior and an exact spacelike escape, proving that these finite prefixes have infinite volume. Averaging each escape with its A6-reversed image remains spacelike, so the failure persists in the visible symmetry-fixed subspace.

The [H3+A6 report](fresh_h3_blocks_h9_report_20261013.md) gives formulas, scripts, comparison scope, and [independent audit](fresh_h3_blocks_h9_independent_audit_20261013.json). Seventeen focused regressions passed, including independent block-enumeration comparisons and forged-certificate rejection. The audit does not repeat the entire positive-height spatial enumeration with a second algorithm.

## Compact H9: genuine dihedral-14 mirrors over a cubic field

Write τ=2cos(2π/7), C=τ²+τ−1, and β=2−τ. The new lattice family is

```
<-α> + B^4 + <1>,
B = [[2,-C],[-C,2]],   α ∈ {τ, 2τ, τβ}.
```

Although B initially presents I2(7), an integral coordinate-swap reflection refines its chamber to I2(14). The full zero-height chamber is I2(14)^4+A1, not I2(7)^4+A1. Independent binary enumeration and a trace argument certify its completeness, with 114 oriented primitive zero-height roots.

The primitive norm classes are {1,2,β,2β}. For α=2τ the producer conservatively includes two further imprimitive candidate norms and removes parallel duplicates. The exact bounded enumeration uses all three real embeddings, reconstructs integral vectors from dominant covectors, and processes equal-height packets together.

| α | Walls at height 30 | Exact spacelike escape norm |
|---|---:|---|
| τ | 18 | (-19-37τ+53τ²)/14 > 0 |
| 2τ | 25 | (114-44τ-38τ²)/7 > 0 |
| τβ | 18 | (57-18τ-21τ²)/2 > 0 |

All three finite prefixes have infinite volume. Their escapes use only two spatial blocks, providing a lower-dimensional target for future arithmetic work. H3 and H4 versions were checked as well; none of those finite prefixes closes the volume gap.

The same implementation reproduces the known compact (2,3,14) triangle. Compactness is certified by a positive dependence among spanning affine facet normals together with the exhaustive list of three timelike vertices. This is a calibration, not a new example.

There are valid additional mirrors beyond the chosen cutoff. For α=2τ, the exact Pell root with time and final scalar coordinates

```
X = 6+18τ+8τ²,   Y = 11+28τ+12τ²,
Y²-2τX² = 1
```

is acute-compatible with all 25 saved walls. Its height is 648+1816τ+808τ². Adding it retains the previous escape, so it is an exact additional cap, not a completed search or a finite-volume example. The [Pell certificate](fresh_cubic7_dihedral_pell_cap_20261013.json) demonstrates concretely why the bounded failures must not be upgraded to full-host nonreflectivity.

The π/14 pair also distinguishes these reflection groups from the earlier compared cubic groups whose primitive norm classes are only 1 and 2: the required square identity would have odd valuation at the prime of norm 7. This does not establish a new quadratic space, commensurability class, or literature priority.

The [independent cubic audit](fresh_cubic7_dihedral_independent_audit_20261013.md) records the norm proof and all eleven saved low- and high-dimensional runs. The [producer](../scripts/fresh_cubic7_dihedral_h9_20261013.py) and [Fraction-only verifier](../scripts/fresh_cubic7_dihedral_independent_audit_20261013.py) are separate. The verifier reconstructs the zero-height census and replays every saved root, angle, strict interior, and geometry certificate; it does not independently redo all positive-height candidate enumeration.

## Compact H9: lift the known H8 normals

Keeping the old projected wall directions and adding one Euclidean coordinate gives exactly two nontrivial Coxeter lifts up to symmetry. One changes an endpoint label 5 to 3; the other changes it to 4. Both eleven-wall Grams have rank 10 and signature (9,1,1).

The first contains an affine E7 subdiagram, the second an affine B5 subdiagram. Each supplies an explicit positive combination of wall normals that is a future null vector. Any extra wall acute-compatible with every retained wall has nonpositive pairing with that null vector. Thus every such completion retains an ideal ray and cannot be compact, regardless of the number or height of added caps.

The [lift report](fresh_h8_lift_h9_20261013.md) and its Fraction-only replay give the complete restricted lift classification and obstruction. Changing or removing old projected wall directions falls outside this argument. Noncompact finite-volume completions are not excluded.

The useful next changes of mechanism are sections using at least three original wall normals, H8 extensions that change some old projected directions, or targeted arithmetic cuts in the remaining H9 escape subspaces. The present certificates close the stated families only. There is no discovery claim.
