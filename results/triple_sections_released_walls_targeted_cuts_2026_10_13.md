# New-dimension search: three-wall sections, released H8 walls, and targeted H9 cuts

**No new compact or finite-volume Coxeter polytope was found in this round.** The searches produced new exact candidate mirrors, Coxeter cones with up to 773 walls, and several stronger obstructions for specific construction families. Every cone reported below still has an exact infinite-volume witness; every tested finite-volume H20 section has a forbidden actual ridge angle. These results do not prove nonexistence in H9 or H20.

The intended advances remain a **compact H9** example, a **finite-volume H20** example, or a finite-volume example above H21. The known examples and the distinction between chamber novelty and new algebraic realizations remain documented in the [lineage ledger](higher_dimensional_polytope_lineage_ledger.md). The research frontier is also listed on [Felikson's hyperbolic Coxeter polytope page](https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/polytopes.html).

This round changed the construction mechanisms rather than merely raising earlier cutoffs: it allowed continuous three-wall section parameters, deleted an inherited wall before extending an H8 example, broke block symmetries in arithmetic H9 hosts, and replaced premature trial mirrors when they excluded better cuts.

| Direction | New result | Exact limitation |
|---|---|---|
| H20 sections of the frozen Borcherds H21 chamber | All central sections supported on at most three original normals fail; 64 noncentral angle repairs also fail | Central theorem is continuous and complete; noncentral search covers only its declared pencils and labels |
| H9 extension after releasing one H8 wall | Rank-ten seed and 452 accepted individual caps, including quartic free-distance caps | Every subset of this finite cap bank fails; other caps and seeds remain open |
| Two golden H3+A6 arithmetic H9 hosts | 1,134 new compatible primitive mirrors and ten Coxeter cones | All cones escape; two exact slice obstructions apply only to their stated subspaces and retained walls |
| Cubic-seven arithmetic H9 host | Symmetry-breaking and scalar cuts, with replacement of conflicting earlier caps | All tested cones escape; bounded support dictionaries are not full ambient enumerations |

## H20: all central three-wall sections, then exact noncentral angle repair

The [H20 report](fresh_h20_triple_sections_20261013.md) proves that every interior H20 section through the fixed center `h`, with normal in the real span of at most three original walls, is non-Coxeter. The 210-wall input has 1,521,520 unordered triples, decomposing into 49 exact symmetry orbits. Rescaling the walls to have equal pairing with `h` reduces all central coefficients to 294 ordered continuous half-pencils. Exact actual-ridge certificates exclude every real parameter in them. The earlier all-parameter two-wall result handles boundary supports.

The proof checks actual local facets, strict timelike ridge points, polynomial signs, and necessary angle equations. It uses polynomial gcds and exact Sturm counts; it is not a coefficient grid. Interpolating timelike points from two ambient faces supplies ridge witnesses where orthogonal projection of the center leaves the chamber.

A separate experiment deliberately repairs a selected ridge angle in 17 noncentral symmetric pencils. Labels 2–6 give exactly 64 genuine parameter choices: nine rational, four quadratic, and 51 quartic. Each finite-volume section has another exact forbidden actual ridge. This bounded noncentral experiment does not exclude unequal positive coefficient ratios, higher labels, or arbitrary three-wall normals away from `h`.

Both independent standard-library verifiers pass, including orbit and root-denominator reconstruction, forged-certificate rejections, and seven Sturm calibration cases. The report also corrects a potentially misleading shortcut: rational projected normals imply quasi-arithmeticity, but do not by themselves imply arithmeticity. Arithmetic H20 exclusion results cannot replace the angle tests for rational sections.

## H9: release an inherited wall and solve a new dotted distance

The [released-wall report](fresh_h8_release_h9_20261013.md) starts from the known compact H8 example, applies a label-5-to-3 lift, and deletes original wall 4. The ten remaining normals span a genuine H9 Lorentz space and have no affine subdiagram. This escapes the earlier obstruction that depended on retaining all eleven old projected wall directions.

A direct pairing construction gives 434 norm-two caps over Q(φ), two additional caps involving new square roots, and sixteen caps obtained by solving a free dotted pairing. One particularly simple cap is orthogonal to eight seed walls, has label 5 with a ninth, and has a dotted distance given by a quadratic equation over Q(φ) with nonsquare discriminant. It cuts the common escape of the first cap bank.

All 120 pairs among the sixteen free-distance caps conflict. Exact signs handle 103 pairs; the other seventeen have pairings in an interval that permits only labels 2, 3, 4, all then excluded exactly. This remains valid even when two quartic cap fields have a degree-eight compositum.

For each possible free-distance cap, a strict spacelike witness survives an overinclusive set of every potentially compatible earlier cap. If none is selected, the old escape survives. Thus no subset of these 452 archived normals completes the ten retained seed walls to finite volume. The standalone Fraction verifier reconstructs the seed and verifies every cap, conflict, and escape. This excludes the finite bank, not all possible cap constructions.

## H9: targeted cuts in two H3+A6 hosts

The [golden-field report](fresh_h3_a6_target_report_20261013.md) considers the hosts `<−α>+H3+A6`, with `α=φ` or `2φ`. A saturated integral slice involving the two middle A6 weights contains no primitive reflection normal: reduction modulo the inert prime 3 gives an impossible square equation over F9. This is an arithmetic obstruction at every bound.

Adding a neighboring H3 weight produces 772 and 362 new compatible primitive normals in a declared coefficient box. Five maximal compatible selections in each host produce ten cones, with wall counts ranging from 153 to 773. All have strict timelike interiors, exact Coxeter angles, and exact spacelike escapes.

Positive linear identities between restricted wall pairings prove that, with the original 45 or 56 walls retained, no real cuts confined to this neighboring slice can eliminate the protected escape at any height. These are precise slice obstructions; the full hosts and cuts outside those slices remain unresolved. An independent exact audit checks the 1,134 accepted roots, all conflict edges, the ten selections, geometry, saturated bases, and both obstructions.

## H9: break a dihedral block symmetry and reconsider premature mirrors

The [cubic-seven report](fresh_cubic7_target_cuts_20261013.md) concerns `<−2τ>+B⁴+<1>`, where `τ=2cos(2π/7)` and each binary block has a complete `I2(14)` zero-height chamber. An exact projection identity proves that any mirror cutting the old escape must break one specific block-swap symmetry. The result applies to all real compatible normals fixed by that swap, independently of field or height.

Allowing an asymmetric block produces 48 new compatible mirrors and a 73-wall cone. A complete bounded two-block dictionary then exposes useful earlier mirrors excluded by some trial caps. Rebuilding from the original 25-wall prefix gives a different 122-wall cone; retaining all previous caps gives 158 walls. Both include a separately verified Pell cap and still escape.

A scalar-plus-equal-block search produces 25 further dictionary mirrors. Three completions have 182 or 146 walls and a common new escape. Rebuilding in height order correctly rejects the old Pell cap when it conflicts with a lower-height scalar cut. Exact finite-angle gates distinguish admissible `π/7` from the arithmetic but forbidden `2π/7` angle.

The larger slice with a full block, an equal block, and the scalar then produces 663 compatible dictionary mirrors, 76 cutting that previous escape. Four completions or replacements have 742 or 706 walls. A final exact witness satisfies all 855 mirrors in the combined dictionary, the original 25 walls, and the separate Pell cap. Thus every compatible subset of this bounded bank still has infinite volume. This union is not a complete ambient enumeration. Even the successful removal of several successive gaps therefore does not supply finite volume.

The independent Fraction replay verifies the saved mirrors, selections, all cone angles, strict interiors, and strict spacelike witnesses. A separate [enumeration review](fresh_cubic7_two_block_completeness_audit_20261013.md) proves the support-specific embedding bounds and checks an independent height-30 calibration against the full host enumerator: all 20 roots agree. It does not claim a second exhaustive height-1000 run.

## Reproducibility and remaining scope

Each linked report supplies exact JSON artifacts, producers, and independent replay commands. Negative volume certificates are algebraic identities and strict inequalities; numerical solvers only propose candidates. The cubic regression suite has six passing tests, the golden targeted suite has four, and the H20 and released-wall reports record their independent mutation checks.

The strongest completed results close specific mechanisms. They do not close general noncentral H20 sections, arbitrary released-wall H9 caps, or either full arithmetic H9 host family. A new-dimensional claim still requires a complete valid Coxeter angle check, a finite-volume or compactness certificate, and a separate prior-art comparison. None is claimed here.
