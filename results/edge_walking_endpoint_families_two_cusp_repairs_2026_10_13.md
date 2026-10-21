# Endpoint families, finite-edge constructions, and two-cusp repairs

This round constructs further exact Coxeter chambers with finite vertices in H9 and tests a new two-parameter section family in H20. **No new compact or finite-volume Coxeter dimension has been established.** The H9 chambers retain certified strict spacelike escapes; every tested finite-volume H20 section has an actual forbidden ridge.

The user's proposed comparison of Vinberg wall enumeration and Allcock edgewalking is checked against the primary paper in [the algorithm audit](vinberg_allcock_search_strategy_20261013.md). Allcock's published algorithm is itself arithmetic, over integer lattices. Its geometric reduction suggests a broader diagram and metric search, but our current local constructions do not exhaust the edge frontier or certify the full chamber of an arithmetic lattice.

## Split-seed endpoint families

The previous exclusion of every completion from the [3,034-mirror bank](fresh_h9_split_leaf_endpoint_cover_20261013.md) motivated changing the normal parameterization. Allowing two independent opposite-block distance totals, and imposing a finite endpoint diagram, produces **728 exact distinct mirrors** compatible individually with the sixteen-wall seed. Independent exact comparison finds 32 in the old bank and **696 outside that literal bank**. Ten cut the previous bank's universal escape. These are local mirror constructions; they are not a claim of 696 previously unknown polytopes. See [the independent certificate](fresh_h9_split_endpoint_families_independent_20261013.md).

A numerical selection from the combined banks gives a 72-wall chamber with four finite vertices. The first edge-extension scout considered only seed-wall anchors and produced 852 proposals; coherent replacements made more finite vertices but reopened gaps when old walls were removed. Allowing every existing wall to supply the extra angle condition gives 7,138 numerical proposals. Of these, 3,857 retain every old wall, and 306 both retain every old wall and cut the selected leading escape. These counts describe the restricted numerical scout, not an exhaustive search.

Ten of the no-removal cutters use the same quartic anchor field and admit exact reconstruction. Their compatibility graph is K6,4: all 24 cross-pairs are strictly dotted and all 21 within-part pairs have positive inner product. Each maximal pair produces **74 genuine facets and at least six certified actual finite vertices**, retaining all 72 old walls. The complete affine-subdiagram check is empty. However, a single strict spacelike vector satisfies all 82 old-and-new inequalities, with square

\[
\frac{11764468488123570049-6234856568308221881\varphi}
     {16000000000000000000}>0.
\]

Thus no selection from these ten caps completes the chamber to finite volume. The construction and independent replay are recorded in [the exact endpoint-walk report](fresh_h9_split_endpoint_walk_exact_20261013.md). This is a scoped exclusion of these ten additions, not of the full numerical scout or arbitrary completion of the seed.

## Arithmetic finite-edge construction

Three successive finite-edge constructions add coherent four-mirror orbits to the previous quartic arithmetic chamber, taking it from 174 to **186 H9 walls**, with **76 walls in its forced H5 face**. The successive H5 edge diagrams are B2+A1+A1, B3+A1, and F4. Each addition cuts the preceding H5 escape and creates an actual finite endpoint. Every old wall is retained. Exact independent replay checks the roots, reflection integrality, all embeddings, all chamber pairs, actual vertices, and surviving escapes. See [the arithmetic edge-walk report](fresh_h9_quartic_edge_walk_20261013.md).

The unit searches are complete only for their stated finite profiles and bounded grids. They do not exhaust all residual norm-equation orbits or prove that each new endpoint is the nearest one in the full lattice chamber. Both final cones still have infinite volume.

## A separate 73-wall walk and H20 sections

The [other finite-vertex walk](fresh_h9_vertex_walk_20261013.md) makes four coherent replacements while retaining the original ten seed facets and the required polar wall. The final 73-wall chamber has four simultaneously certified new finite vertices and passes its exact affine checks. Its smaller remaining numerical gap is still an exact infinite-volume escape. A [second geometry audit](fresh_h9_vertex_walk_independent_geometry_20261013.md) independently checks the final configuration.

In H20, a new family uses three null cusp vectors, with section normal E=7u0−x u3−12y u4 and square 336(x+y−xy). This gives two parameters while retaining seventeen prescribed orthogonal walls. Nine selected exact sections have finite volume and eighteen actual prescribed angle repairs; each also has an actual forbidden ridge. Seven irrational section normals require at least five original wall normals and have nondegenerate rational coefficient planes. A separate six-chamber star is obstructed throughout its specified positive A2 family. See [the H20 construction and independent replay](fresh_h20_two_cusp_repairs_20261013.md).

The [explicit edge-frontier diagnostic](fresh_h9_split_edge_frontier_20261013.md) numerically finds 32 unresolved outgoing directions at the original four finite vertices. With cap 2747 added, it discovers the fifth finite vertex and records 39 unresolved directions. It enumerates edges of the specified chambers, not missing roots of a full lattice group. The productive continuation is to branch on these unresolved endpoints. Additional walls or finite vertices alone do not certify progress toward finite volume; the exact surviving escapes remain the decisive checks.
