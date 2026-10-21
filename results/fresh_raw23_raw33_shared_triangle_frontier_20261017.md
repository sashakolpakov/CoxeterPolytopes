# A bounded raw23/raw33 shared-triangle frontier

This is a completed local diagnostic in the compact arithmetic H13 model over K=Q(sqrt5). It finds an actual event equality but no improvement to either current fixed-mass witness. It does not assert physical realization of a graph cover, target12 feasibility, or a dimension theorem.

The source is [the exact producer](../scripts/fresh_raw23_raw33_shared_triangle_frontier_20261017.py), with [saved data](fresh_raw23_raw33_shared_triangle_frontier_20261017.json). Its baseline is the frozen [shared-tetrahedron packet](fresh_shared_prism_tetrahedron_incidence_20261017.md). The frozen [1325 packet](fresh_forest1325_cross_vertex_incidence_20261017.md) supplies the strict-center finite-support criterion and actual full-star uniqueness used below. All pair owners, including both owned unknowns at1763, remain unchanged.

## Why this pair was selected

Forest1806 has type2A1+2A4+H3. The literal pair stars32/61 comprises raw23 face26410 (provenance483) and raw33 face27350 (its inherited unique actual provenance). Its prior restricted attaining-cover mass was0.0063635749 at an almost tight vertex budget, large enough to warrant this bounded test. The earlier candidate982 raw23/raw33 event had only0.0009802684 in the current cover; its larger historical mass was stale. No repeated scan of the completed13-source family was performed.

There were two genuinely changed small local LP proposals, one at each fixed face-mass choice, followed by exact rational recovery. There was no global model optimization.

## Exact native configuration and all finite supports

The passive forest is A4+H3, of rank7. Active base roots0,1 are isolated; roots2–3–4–5 form a simple A4 chain. External roots6=A,7=B,8=C have the following nonzero norm2 products, all with the acute negative sign:

- 0–A has squared product (7+sqrt5)/2;
- 1–B has squared product20+8sqrt5;
- 5–A and2–B are simple;
- 2–C has label5;
- A–B has squared product9+4sqrt5, hence product−(2+sqrt5).

All other products outside the base chain are zero. Positive rescaling of0 and1 puts the complete9-root Gram matrix in K. Its inertia is(6,1,2). Each of28 complete native embeddings independently rechecks all three positive external heights, all three future-side cross products, both actual star contact pairs, their native omitted products, provenance, and occupied-port complements. No coefficient vector is moved between different pair fields.

The source exhausts all C(9,6)=84 active supports. Exactly13 are positive definite. Adding the passive rank7 gives actual finite vertices: for an actual wall set with positive-definite acute Gram G and a strict future timelike chamber center h, the orthogonal projection x=h−N G^{-1}(N^T h) is future timelike. Since G^{-1} is entrywise nonnegative and all other actual wall products are nonpositive, every additional wall remains strictly inactive at x. Rank13 leaves the unique vertex line. Thus arbitrary additional chamber facets cannot invalidate these vertices.

The raw23 prism fixes active root1 and the passive7; it omits B from its full star and has10 finite vertices. The raw33 prism fixes active roots0,4,5 and the passive7; it omits A and has6 finite vertices. Their union contains all13 positive supports; their intersection contains exactly three. There is no new finite vertex among these nine normals outside both prisms. This is a bounded support conclusion, not a statement about other unknown facets or other configurations. The dotted A–B product also directly precludes any positive-definite support containing both A and B.

## The actual triangle and its event bijection

Both actual prisms intersect in the common11 independent fixed walls, passive7 plus{0,1,4,5}. The remaining sides are{2,3,8}. Their Schur Gram has inertia(2,1,0), and its inverse is entrywise strictly negative. These sides bound a compact hyperbolic triangle. Its three vertices are the already verified actual finite supports:

| Forest | Active support | Complete event orientations |
| --- | --- | ---: |
| 1806 | 0,1,2,3,4,5 | 8 |
| 1505 | 0,1,2,4,5,8 | 8 |
| 944 | 0,1,3,4,5,8 | 12 |

The chamber face lies inside this three-side simplex; all three simplex vertices belong to the chamber. Convexity gives the reverse inclusion, so this is the actual intersection face even in the presence of arbitrary additional facets.

Each counted event is the same unordered pair of the two actual prism types, pointed at one of the three distinct forest types. Complete component embeddings reconstruct every8/8/12 native pair; distinct embeddings give distinct literal pair keys. The occupied-port complement assertions identify the actual supporting walls of each prism, not merely its abstract diagram. The native future-height equations determine the same external actual normals in both directions. Each compact prism is the convex hull of its finite vertices; at a finite simple ambient vertex its support and side walls exhaust the incident walls. Consequently the reverse transport recovers the same actual prism and side ownership. Each intersection triangle has exactly one vertex of each of the three distinct types, with no automorphism multiplier. Hence

    X1806 = X1505 = X944.

This is a new equality of local event totals. No target-vertex quota or1325 pointed injection is assumed for this different configuration.

## Complete local cover test and exact zero gain

All inherited compatibility graphs are reconstructed before the event test. Unknown pairs remain compatible. The three graphs have336,316,7114 maximal cliques. Every nonempty clique subset is included; grouping by actual corner-orbit counts and the new event indicator yields35,87,195 vectors, respectively. There is no matching shortcut. Orbit coefficients count actual prism corners, not embedding multiplicities.

For each source f, let lambda_fj be the mass of saved vector(v_fj,e_fj), let d_f be the current actual face-derived demand, and let mu_f be the existing vertex multiplier. The local proposal minimizes sum(delta_f), subject to

    sum_j lambda_fj v_fj >= d_f,
    sum_j lambda_fj <= mu_f + delta_f,
    sum_j lambda_fj e_fj = X,
    lambda_fj, delta_f, X >= 0.

An exact dual has nonnegative orbit charges w_f, source costs0<=a_f<=1 and signed event coefficients z_f with sum(z_f)<=0; every vector obeys w_f·v<=a_f+z_f e. Summing and canceling the common event gives a valid lower bound on the total increase. Here the all-zero dual gives lower bound zero. The substantive certificates are exact nonnegative rational covers meeting all demands, staying within the old source budgets, and sharing a common event mass. Both achieve zero increase:

| Fixed masses | Common event X | Total multiplier increase |
| --- | ---: | ---: |
| Full | 7943653/1250000000 | 0 |
| Named248 restriction | 64197493/10000000000 | 0 |

No capacity column is exported. Every alpha and mu remains identical to the baseline, so all439410 inherited model columns, signed-face columns, previous incidence columns and owners remain preserved. The rest of the571-forest data is untouched. The ratios are unchanged:

- full:111428537994904130/9288556982227791 =11.99632388627268;
- named248 restriction:18670750401222348/1549226494433359 =12.05165963034431.

The surviving above12 restriction means this triangle equality does not remove the current obstruction. The exact covers certify the stated conservative local relaxation, not a globally realizable chamber. No additional geometry bank was opened after this bounded result.

## Verification scope

Producer generation23253 and separate fresh saved-certificate replay26709 both completed terminal exit0/PASS (21.09s and26.59s). Fresh replay reconstructs the84 supports,28 native embeddings and complete three graph/event domains, then checks saved rational covers without invoking an optimizer. The separate [stdlib saved checker](../scripts/check_raw23_raw33_shared_triangle_saved_20261017.py) independently checks saved support positivity and triangle matrices, enumerates saved graph cliques/subsets, reconstructs demands, and checks the exact common-event covers and unchanged witnesses; it does not independently reconstruct native projections or inherited pair ownership. The final stdlib checker run completed terminal exit0/PASS:84 exact support positivity decisions, the triangle Schur matrix and inverse,7,766 saved graph maximal cliques,317 count/event vectors,634 charge inequalities,78 positive mixture terms, six common-event equalities,13 direct source hashes and five report links. Its lower dual charges are all zero; the positive rational covers are the nontrivial evidence.

Root and construction peer each completed a full producer source and mathematical/report read, including actual face identification, reverse prism ownership, complete embedding/subset domain, exact zero-gain cover logic and unchanged-witness scope. Neither review duplicated the84-support,28-native-map, graph or LP census. The root has not yet separately reviewed the saved checker; no such review is claimed here. No live job or additional bank remains.
