# Complete small-core domain for distance-nine quadrilateral payments

**FINAL / FROZEN. Generation and fresh saved replay PASS; narrow internal mathematical/source review PASS.** This is a local classification and isolation lemma. It does not establish a new ambient dimension bound or complete the remaining high-rank distance-nine analysis.

Consider an actual quadrilateral face star in a finite-volume hyperbolic Coxeter polytope, with supporting white roots W and opposite black pairs U,V. Normalize roots to square 2. Assume that no U-to-V path of length at most nine has white internal nodes. If a minimal indefinite core of W+U or W+V has order at most five, the configuration is either impossible or has an adjacent marked corner whose whole connected component is isolated and has **order at most 12**.

The [source](../scripts/fresh_distance9_low_core_20261017.py) and [certificate](fresh_distance9_low_core_20261017.json) cover every finite label and every additional supporting white under these hypotheses. They bind the [high-rank catalogue](fresh_high_rank_quasilanner_census_20261017.md) and the actual-face, positive-kernel, and cusp-factor arguments of the [completed distance-ten proof](fresh_selective_distance10_dimension116_20261017.md). No frozen predecessor is changed.

## Why smaller cores enter here

Each minimal indefinite core contains both of its black roots: omitting either black leaves a positive-definite actual-edge support. Its proper principal subdiagrams are positive semidefinite. The inherited complete classification bounds its order by ten, and the high-order marked radius by seven.

A connected core of order k with two marked roots has radius at most k−2. Thus cores of order at most three cannot participate in an unpaid disjoint join: their radius plus seven plus the joining edge is at most nine. Two cores of order at most five are excluded similarly. An overlap with a small core requires a shared node with distance sum at least ten. Consequently only order-five radius-three nodes can overlap a high-order radius-seven node.

Disjoint joins are different. Order-four radius-two ports can meet radius-seven ports. Order-five radius-two or radius-three ports can meet the appropriate radius-six or radius-seven ports. The order-five radius-three tree list alone would miss this domain.

## Complete all-label enumeration of orders four and five

Every edge of a connected four-node graph belongs to a connected proper three-node subdiagram. If its norm-two magnitude exceeds sqrt(3), adjoining any nonright edge already gives negative determinant, even when the third possible bond is removed. Hence every bond in an order-four minimal indefinite core has label at most six. Equality [6,3] can be affine, so label six must be retained. Parallel and dotted contacts are excluded by the same proper-principal test.

For order five, extend an edge to a connected proper four-node subdiagram. A label at least six gives an indefinite connected triple, or the affine [6,3] equality. In the latter case the connected fourth node has a negative product with its positive kernel, so the four-node diagram is indefinite. Thus order-five labels are at most five. These bounds are proved before enumeration; neither is a numerical cutoff.

The checker enumerates every weighted graph on four nodes with edge alphabet {right,3,4,5,6}, and every graph on five nodes with alphabet {right,3,4,5}. It retains exactly the connected indefinite graphs whose one-node deletions are all positive semidefinite. Every smaller proper principal subdiagram lies in such a deletion. Black pairs are selected only from strictly positive-definite deletions.

The positive-semidefinite recognizer is elementary and exact. For a tree component, leaf Schur elimination uses squared bond magnitudes 1,2,(3+sqrt(5))/2,3 in Q(sqrt(5)); a zero pivot with a remaining nonzero incident bond is indefinite. A component containing a cycle can be positive semidefinite only when it is precisely a simple affine A cycle. Indeed, the all-ones vector on the cycle has nonpositive square. Equality requires all its bonds simple, and its positive kernel excludes every chord or attachment. This handles disconnected finite/affine deletions without assuming a selected diagram list.

Canonicalization tests every vertex permutation, preserving all bond labels. The resulting census is:

| Order | Labelled minimal-indefinite graphs | Canonical types | Relevant marked ports | Inverse-admissible ports |
|---|---:|---:|---:|---:|
| 4 | 392 | 32 | 17 | 11 |
| 5 | 1,050 | 14 | 68 | 18 |

All marked pairs and ports are retained within each canonical graph. A relabelling therefore cannot omit a marked case. The exact far-port inverse diagonal is computed from spherical-deletion determinants. A native simple contact requires both inverse diagonals strictly below two; a doubled contact requires both strictly below one. The strict inequalities come from actual-edge positive definiteness, not an ideal-corner semidefinite test.

The high-order inverse filters leave seven relevant ports, all in L25: three radius-seven cases at node 9, three radius-six cases at node 9, and one radius-six case at node 8. Their literal contexts and inverse values are saved.

## Disjoint joins and every possible joining label

Each low port has a two-node white path; each relevant high port has a three-node white path. Coefficients (1,1/2) give square at most 3/2 on the first path, and (1,2/3,1/3) give at most 4/3 on the second. With native cross magnitude beta, their combined square is at most

    17/6 − 2 beta.

Supporting-white positive definiteness forces beta<17/12. This permits precisely labels 3 and 4; label 5, parallel contacts, and dotted contacts are excluded. Any other acute bonds only strengthen the nonpositive-vector test.

For a proposed cross label, the inverse tests at both endpoints are applied. Eliminating the two spherical white blocks gives the additional strict condition

    beta² (G_Wlow inverse)xx (G_Whigh inverse)yy < 1.

The complete low-first count after the endpoint tests is:

| Low core order | Simple joins with two negative directions | Simple white-PD rejections | Label-four white-PD rejections |
|---|---:|---:|---:|
| 4 | 33 | 0 | 2 |
| 5 | 60 | 2 | 7 |

Thus all nine possible doubled joins fail the white test. The 93 surviving simple joins have 23 underlying types, with inertia (12,2,0) for order-four low cores and (13,2,0) for order-five low cores. Tree parts are eliminated exactly using squared bonds. At most five cyclic nodes remain, and their exact symbolic congruences are recorded. No floating-point signature is used.

One must also consider multiple cross edges. Six low-first marked-context pairs admit two candidate simple edges. In three, using both produces a wholly white simple four-cycle; in the other three it produces a wholly white simple triangle. The certificate stores the vertex set and its all-ones affine kernel. These contradict positive definiteness even after extra acute contacts are restored. Every nonempty subset of allowed cross edges is therefore excluded. Reversing low/high roles gives the corresponding ordered cases without changing any obstruction.

## The two low/high overlap cores

For an overlap, only an order-five radius-three node can be shared. Every such low marked context has at most one eligible node, so larger or different identifications are impossible. After the high inverse gate, there are 30 low-first candidates. Twenty-four fail the low inverse test. The other six pass actual-edge positive definiteness and must be handled geometrically; they are not discarded as numerical near-misses.

There are two literal 14-root cores, each with three high marked pairs. Use low nodes 0 through 4. The right L25 copy maps its node 9 to the shared low node, and its nodes 0 through 8 to 5 through 13.

* The low fork has edges 0–4, 1–4, 2–3, 3–4, with labels 3,4,3,3. The low black pair is {0,1}; node 2 is shared. Its Gram has signature (12,1,1).
* The low path has edges 0–4, 1–3, 1–4, 2–3, with labels 3,3,3,5. The low black pair is {2,3}; node 0 is shared. Its Gram has signature (13,1,0).

The high black pairs are {1,2}, {1,3}, and {2,3}, transported into the right copy. Reversing the two cores gives 12 ordered surviving contexts.

Every unspecified contact between exclusive core nodes would be a marked path of length at most nine. There are 36 such pairs in each low-first context, all checked exactly. Thus no additional negative cross bond can change either core Gram while the no-short-path hypothesis holds.

For every adjacent corner, the checker decomposes the principal Gram into finite and affine components. A further supporting white is orthogonal to any affine component by its positive kernel. On a finite component, a nonzero norm-two contact can occur only at inverse diagonal at most two. Omitted black nodes are included separately in the possible contact domain. These necessary domains include contacts to every core node; they do not assume that new whites meet only marked endpoints.

## Isolation of the fork core

In the literal normalization above, the fork core has radical

    (-1,-sqrt(2),0,-1,-2,6,3,4,2,5,4,3,2,1).

This is an actual root relation in any Lorentzian realization. Its root combination is orthogonal to a timelike vector in the indefinite core span and has square zero; the timelike perpendicular space is positive definite, so the combination vanishes.

Each chosen adjacent-corner contact domain omits the zero coordinate. Its negative coordinates are on the low end and its positive coordinates on the high end. Acuteness and the relation force any nonzero external contact vector to touch both signs. The saved distances then give a marked path of length at most four through that new white. Therefore every additional supporting white is orthogonal to the core. The chosen marked component has order 12.

## Isolation of the label-five path core

This core is nonsingular Lorentzian. The selected corner permits contacts only in two end sets. Touching both ends gives a marked path of length at most four. Within either end, the full inverse-Gram block is entrywise positive, and every diagonal is strictly larger than two. The smallest relevant diagonal is 8+4sqrt(5).

For a norm-two external root with contact vector b, projection onto the Lorentzian core gives

    b^T G^{-1} b <= 2,

because the perpendicular complement is positive definite. A nonzero one-end contact has magnitude at least one and therefore already exceeds two; extra contacts at that end only increase its cost. A two-end contact gives the forbidden short path. Thus every extra supporting white is orthogonal to this core as well. The selected connected marked component has order at most 12.

This is complete isolation, not merely a bounded search for single-node extensions. Arbitrarily many indirectly connected supporting whites cannot enlarge a component to which none of them is incident. At an ideal corner the inherited factor-fiber lemma puts the entire varying affine factor in W plus the marked pair, so additional cusp walls outside the face star cannot enlarge the paid component either. No simple-cusp assumption is made.

## Validation and scope

Generation 59182 and fresh saved replay 5697 both terminated successfully. The parent independently read the all-label bounds, forest/cycle PSD recognizer, complete permutation/port domain, joining-label exhaustion, white-cycle obstructions, and all six overlap contact-domain and isolation arguments, and found no gap. That review did not duplicate the million-graph census or its exact arithmetic. The certificate contains all canonical low graphs, all relevant marked ports, 104 low-first joining-label candidates, the six multiple-edge obstructions, all 30 low-first overlap candidates, all adjacent-corner decompositions for the six survivors, and 216 extra-cross shortcut checks.

The conclusion is restricted to the complete small-core portion of a distance-nine face-payment argument. It supplies neither the remaining high-order extension classification nor an ambient ceiling. The source uses the frozen independent distance-ten review's exact quadratic-field class, and the small cyclic remainders and extension matrices use exact SymPy arithmetic. A fresh replay is a recomputation by this source, not an independent implementation. No external refereeing, proof-assistant formalization, construction, or publication-priority claim is made.
