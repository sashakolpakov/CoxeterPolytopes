# Independent composition review of the dimension ten argument

20 October 2026. **PASS, within the explicit written-geometric and exact-scalar scope below.** No checker, Lean compiler, build or numerical job was run in this review lane. This is a full mathematical/source review of the assembly, rather than a new claim of formalized hyperbolic geometry. The reviewer authored the balanced-count and ten-direction reductions; the assembly and the other scalar proofs received separate authorship and review.

## Source binding and scope

The assembly reviewed in full is `right_angled_mixed_dimension10_assembly_candidate_20261020.md`, SHA256 `70ef7de40e9fbd62dbe763d64d258c40ab777a42d7550e1608074e4c29aa0326`. Its historical pending-check language is superseded by the completed checks and this review, without changing its pinned bytes.

The critical written inputs were read in full: `right_angled_mixed_eleven_three_face_corner_balance_20261020.md` (`73b8073843171d4bf11ad0fc11829b0855df36b2a78307fc024cd6d26dd450df`); the one-cusp proof (`19590da4197ee04940df141a492ff89d69aea44b49bf87fe3c63699c252b56d4`); the general subset/moment note (`8c2051b3be0f692e5097eedd81cf494f15313037d36b0f6625ca78cd01aabd34`); the seven/eight/eleven-direction note (`04d3b8fa433365210600234709894d777803e7db8591cbb32876a0796e87a1ce`); the nine-direction global proof (`bf288985ade5acf93e13faafb25544b0642b38f13e422ad38a3f1223e9706d15`); and the ten-direction rational envelopes (`7fdcad68b9b4ae113c482e9e17266a6f5e9aea0d21409eb7db332ec915627e45`). Their standard-library exact checkers were also read. The older nine-direction interval proof was independently read and supplies additional support, but the global analytic proof suffices.

The actual-wall projection, canonical P3 recognition, bad-triangle attachment and opposition formulas were checked against the previously fully reviewed focused dimension eleven manuscript, including the complete all-wall feasibility argument. Published inputs remain the compact right-angled ceiling four, the entirely ideal ceiling six, and the facet minimum 82 in dimension eight used in the weighted face row. Standard acute-angle face/link geometry, Euler and Dehn--Sommerville identities are explicit mathematical inputs. These are not proved by coefficient certificates.

## Actual geometry and ordinary counting

The two-normal projection argument uses the correct sign: a product in the interval (-1,0) would produce an actual intersecting two-wall face with a non-right angle. The M-matrix/inverse-sign argument ensures feasibility against all other walls. The six canonical P3 walls are actual walls, and convexity of the closure of their feasible vertex hull forces the whole three-face. The bad-triangle attachment lemma uses the native product gap and closes the six-wall chamber; it is not merely a local drawing.

The canonical opposition move changes exactly the three selected normals. Direct substitution preserves `h=3q−Σn_i`, whose square is two in dimension eleven. A single future null vector perpendicular to this spacelike vector exists and is fixed on each finite opposition component. Its pairing with every future unit vertex is strictly negative. Thus all denominators are positive, each actual reverse ratio is the reciprocal of the forward ratio, and the same vector is used at both ends of every undirected edge. No maximum-height vertex or weighted stationary average is substituted for ordinary vertex counting.

The finite vertices and ideal vertices have the expected orthant and cube links. The finite-corner multiplicity is 165, the ideal three-face incidence multiplicity is 11520, and each polygon has nine three-cofaces. The displayed Euler and truncated Dehn--Sommerville identities therefore refer to the same actual counts. The weighted inequality is strict; its use is justified by the published eight-dimensional facet minimum and elementary ridge incidences, not by an unproved mixed bound.

The one-cusp ten-facet exclusion has a complete actual-facet argument. In the equality pattern, the four cusp quadrilaterals have eight distinct boundary vertices; repetition of opposite cusp facets is excluded by their asymptotic geometry, and convex facet intersections exclude repeated pentagonal neighbors. Repetition two steps apart would close a three-edge boundary of a pentagon. The resulting pentagon adjacency count is 13 versus 11. This is a weaker own proof of a known stronger one-cusp minimum, not a novelty claim about that three-dimensional minimum.

The small three-face inventory and the R7 degree proof are sufficient for `D≥−30E`. The R7 contribution is zero. For ν≥11, the one-cusp exclusion and Euler bound imply `p≤10(ν−10)`, so the negative contribution is bounded by three times this quantity. The corner polynomial row is checked for every k from zero to eleven, not sampled. Combining it with `X≤18E+4R`, the exact Euler identities and the strict weighted row gives `3722B>159195f+2794E` with the same canonical face count B used by transport.

## A direct proof of the required H budget

The assembly's `H≤2E` can be proved without metric recognition of D3 or a classification of the ν=11 three-faces. Here is the full weaker budget needed by the new argument.

For a noncanonical right-angled three-face Y, let p be its number of finite vertices, s its cusps, F its facets, ν=F+s, and A its finite vertices whose three edges all end ideally. Then `p+4s=2ν−4`, `p≥A`, and `3A≤4s`.

If A>0, take such a corner. Each incident polygon that is not a bad triangle contributes at least one polygon-excess unit. Each bad triangle has an opposite actual polygon across its ideal--ideal edge. That polygon is not a bad triangle, by the attachment lemma, and contributes at least one excess unit. The opposite polygons are distinct and avoid the chosen corner (their normals are the distinct `q−n_i−n_j`). Hence there are three excess units in all, so `F+2s≥15`.

For ν≤10, the inequalities `s≥15−ν` and `4s≤2ν−5` contradict each other. For ν=11 they force s=4, p=2, so A≤2. For ν=12, A≥5 would give p≥5 and `4s≥3A≥15`; integrality gives 4s≥16 and hence `p+4s≥21>20`. Thus A≤4. Finally, for ν≥13, `4A≤p+4s=2ν−4`, whence `A≤(ν−2)/2≤2(ν−10)`. Summing these bounds over actual noncanonical three-faces proves `H≤2E`. It also confirms E≥0. No global hypothesis about finite--finite edges enters.

## Scalar contract and exact certificates

The frame equations are exactly `Σp_i=3`, `Σp_i²=1`, with real coordinates of either sign. Their Cauchy range makes every formal triple ratio positive. Formal triples are used only to enlarge a sum of nonnegative terms; they are not asserted to be actual vertices or faces.

All subset sizes are covered: sizes zero through six by the elementary paired-triple argument; seven and eight by the full rational moment envelopes; nine by the global analytic proof; ten by the joint-moment rational envelopes; and eleven by its separate rational moment bound. The ten-direction proof correctly replaces the entire nonnegative quartic bracket when bounding its remainder coefficient. Its negative-mean and positive-mean envelopes have the required signs and domains. The 840 Bernstein intervals cover all three polynomial domains, and all 12550 recorded coefficients are strictly positive. This is a finite polynomial-positivity supplement, not a geometric census. The seven/eight interval certificates and the nine-direction global polynomial identity have separate actual root-owned generation/replay ledgers. The independent all-subset review `review_mixed_eleven_all_subset_contracts_resource_20261020.md` is an additional source-level PASS.

Actual P3 triples lie within the ideal-ending subset. Positive-part deletion therefore gives the displayed pointwise row with constants 9/25 and 27/100. Summing over ordinary finite vertices pairs every actual opposition edge in its two directions and gives `r²+r⁻²−2≥0`. Consequently `4B≤168f+3E`. Twice the strict lower row and 1861 times this upper row give the impossible `5742f+5E<0`, since a mixed polytope has f>0 and E≥0. Every count and strictness needed for this final comparison is present.

## Universal dimensional tail and verification boundary

The no-mixed-eleven theorem alone yields the universal ceiling ten. Any finite-volume right-angled polytope of dimension at least eleven has an eleven-dimensional face, obtained by a face flag. Right angles are inherited in the induced Lorentz subspace. The finite and ideal face links are faces of the corresponding orthant and compact cube links, so this face has finite volume. An eleven-face cannot be compact or entirely ideal by the published ceilings four and six. It would therefore be mixed, contrary to the theorem. This avoids dependence on the earlier internal ceiling eleven or the old unrestricted ceiling twelve.

The complete written composition passes. Scalar checks certify the stated rational identities and positivity obligations; they do not formalize actual-wall geometry or prove that every formal subset is realizable. The new conditional Lean module is source-prepared at this review date and remains pending the parent's capped compilation and axiom audit. No new Lean verification is claimed here.
