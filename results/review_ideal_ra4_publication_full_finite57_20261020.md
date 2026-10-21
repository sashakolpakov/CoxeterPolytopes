# Independent full mathematical source audit: ideal RA4 volume and chains

20 October 2026. **PASS.** I found no mathematical gap in the held publication text or in the scope of its new fixed-model algorithm statements. This is a fresh source audit, not a new computation, external peer review, or end-to-end formal verification.

Source: `manuscripts/ideal_ra4_volume_and_structure_publication.tex`, SHA256
`e0c2ca24cec24ffcd736e201d3c7e06d75b28ef7be8093decd321ae3ac54bda1`.
I read the complete 1102-line source and derived the geometric/counting arguments before consulting supporting sources for the declared finite inputs. I did not author this integration or edit the primary. No script, imported module, solver, TeX build, or Lean process was executed.

## Volume gap and equality spectrum

The box-link counts give `e=4v`, `r=3v+f`, and total ridge-vertex incidence `12v`. Consequently `W=3(v-f)` and the stated volume formula agree. For a facet, four-valence gives `d=v_F-6=sum(k-3)` and total facet defect `2W`. Defect zero is exactly the octahedron; positive defect is at least two by the antiprism bound.

The opposite-cusp injection has the necessary global distinctness. Opposite box faces give congruent rectangles in asymptotic facet walls, and those walls have exactly one common ideal point. They cannot repeat at another cusp or be ridge-neighbors. Thus the inequality with `q+s` counts distinct positive facets, rather than merely local incidences.

I checked the antiprism model and width calculation: `r^2+h^2=1`, `bc=r`, and `mu=2h/r`; the two opposite-unit-normal sums have coefficients `2/mu` and `h/r+1/mu`, whose ratio is `2 cos(pi/k)`. The minimum-defect split is exhaustive and yields `2W>=20`; divisibility by three then gives `W>=12`.

At equality, the exclusions of minimum defect two, at least five, and four are sound. In particular, a square-antiprism cusp requires positive triangular neighbors in addition to its two quadrilateral neighbors and eight cusp opposites; these lists are disjoint. For minimum four, the four-quad matching bounds every neighbor's available positive neighbors and forces the defect sum at least28. Four pairwise adjacent right-angled facets would meet at a finite point by the stated timelike projection, so the clique contradiction is geometric.

For seven positive facets, the three excess-defect partitions are complete. A defect-four facet cannot have a D9 quadrilateral neighbor, using the paired quadrilateral and the absence of positive triangular neighbors of D9. The remaining defect-six facet has six nonadjacent quadrilaterals and exactly two alternating quadrilaterals at every cusp. Their incidence graph is simple four-regular on six vertices. The embedded medial recovery identifies the cuboctahedral face lattice; it does not assume that an abstract graph alone fixes an arbitrary embedding. In the displayed metric model the unit triangle-normal sums are sqrt(2) times the square-normal sums, so all twelve cusps are nonsquare. Opposite-cusp injection then contradicts seven positive facets. Hence equality really implies exactly eight D9 facets **before** the chain theorem is invoked.

The published minimum, antiprism bound/equality, admissible nine-vertex case, ten-vertex Figure7 case and its matching, and three-dimensional right-angled rigidity remain explicit external inputs. I did not reprove their historical classifications in this audit.

## Restricted Oct/D structure and actual sections

The cusp-box classification is correct under the stated metric facet spectrum. The `2:2:1` case would give a square end with four quadrilateral sides, impossible for either permitted facet; the `2:1:1` case has four doubled belt facets and octahedral ends. The distinguished triangles therefore form a closed two-manifold after adjoining ideal points. Finite edge points lie inside actual ridges, and each ideal link is the single four-segment middle square. In particular no unseen triangle sheet or joining of unrelated belt components is omitted.

The Euler argument does not presume orientability or a sphere. A primal spanning tree has a disjoint dual spanning tree because a separating face boundary would be a nonempty even-degree subgraph of the primal tree. Thus `chi<=2`; together with `chi=n/4` and the seven distinct same-component facets supplied by cusp injection, this forces exactly eight triangles. The six-cusp graph has no parallel edges: otherwise at least three belt facets would contain the same two cusps, including an opposite asymptotic pair. Its complement is a matching and its eight distinct triangle cycles give precisely the cube indexing used next.

Equal short widths give the opposite-normal sum equations with **unit** normals. Their sign-character expansion has only constant and linear terms. The unit and orthogonality equations give `a^2=-1/2`, `b_i^2=1/2` with all cross products zero. Hence their span is Lorentzian and defines a genuine hyperplane. Its six null rays are the actual six cusp rays; convexity contains their hull, while the eight wall inequalities give the reverse inclusion. This proves `P intersect H=S`, not merely an abstract octahedral incidence pattern.

The no-outside-facet argument covers every finite point of S: interior points, triangular-face interiors, and diagonal-edge interiors. The latter lie inside actual quadrilateral ridges, whose only two incident facets are the chosen doubled facets. Thus cutting introduces only the six old ideal vertices, has right new angles and nonempty interior on both sides, and leaves every outside facet uncut. Other belt components and their convex hulls stay on one side; cusp-link uniqueness prevents shared ideal endpoints. Iterating decreases the doubled-facet count by eight and gives disjoint cuts.

All-octahedral recognition is valid. A shared triangle and orthogonality fix the neighboring wall. A regular ideal octahedron on its specified side is fixed by that triangle: every face symmetry extends, and the pointwise stabilizer of three distinct ideal points in H3 is the identity or their plane reflection. Propagation matches every model 24-cell facet, and the halfspace/hull inclusions give equality with the whole model cell.

Finally, the D4-root opposition argument excludes two nonopposite disjoint cuts in a cell. For spatial root product minus one their root sum gives a facet orthogonal to both; its two reflected neighbors would force three interior-disjoint octahedra inside one original facet, exceeding both allowed intrinsic volumes. The placement used here is justified by the subsequent flag-transitive no-twist argument. Therefore intermediate cuts are opposite, the tree is a path, and its existing isometry type is fixed by cell count. The argument retains the Oct/D hypothesis and makes no unrestricted decomposition or existence-for-every-length claim.

## Fixed-model algorithms and their proof bridges

I read the relevant independent saved-data checker sources without importing them. The new descriptions match their finite domains:

* The original 24-cell test considers all10626 four-wall subsets. At an extreme point the active spatial coefficients span R4, and every feasible independent four-wall intersection is extreme. Boundedness and a strict origin are proved directly. Thus the retained list is complete for these24 halfspaces.
* For the double, the old/reflected halfspace inequalities identify the entire intersection with the convex union on the two sides of the gluing wall. This supplies the missing completeness bridge that a saved feasible vertex list alone would not give. The42 distinct null points are all extreme by strict convexity of the ideal sphere.
* Repeated intersection with the facet vertex sets halts within the finite power set. Every generated set is a face, and every face is an intersection of containing facets. Affine rank supplies its dimension. This is complete for the specified model, not a census of ideal four-polytopes.
* The cuboctahedron's364 active triples use the same argument after the rational coordinate rescaling. The455 triples of missing graph edges retain exactly the15 labeled perfect matchings, covering precisely the simple four-regular six-vertex graphs. Cube row reduction concerns the six displayed square equations; its four-dimensional kernel is the constant plus three sign characters.
* The276 D4-root pairs and1152 reflected-octahedral vertex incidences have the stated fixed scope. They support the model facts; the disjoint-interior, actual-cut and gluing conclusions remain written geometric proofs.

The reported face vectors, cusp ratios, graph counts and exact ranks remain inherited recorded outputs, not fresh execution results from this audit. Their declared domains match the claims for which the paper uses them. The theorem's f-vectors and volumes also agree independently with the incidence identities and the loss of six cusps and ten facets per gluing. Equality at `8 pi^2/3` therefore gives precisely the reflection double and `(42,168,164,38)` under the explicitly stated published inputs.
