# Full source audit of the compact right-angled four-volume constraints

20 October 2026. **PASS for the stated constraints.** This is an internal mathematical reading, with the external facts below assumed. It is not a new computation, formal verification, novelty assessment or outside refereeing. No source change is requested.

Source: `manuscripts/compact_ra4_volume_constraints.tex`, SHA256 `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`. I read the complete source and reconstructed its deductions before consulting the two three-dimensional metric references. I did not use a prior project review as a proof, run a scientific checker, or change the manuscript. As the coordinating project agent, I am not independent of its development.

## Normal geometry and the capacity count

The positive-principal projection has the right signs: a positive-definite two-wall Z-matrix has nonnegative inverse, so projection of a strict interior point adds nonnegative outward-normal multiples and preserves all other strict inequalities. Its squared norm is more negative and its product with the original future point is negative. Thus a product in `(-1,0)` gives an actual nonright ridge, rather than merely intersecting supporting planes. The `-1` case produces a nonzero feasible future null ray and contradicts compactness. The same calculation for orthogonal collections supplies the actual intersections repeatedly used later. Finite vertex simplicity follows.

For a disjoint pair of facets, their normal plane has one negative direction, so its orthogonal complement is positive definite. Two common neighbors cannot have product below `-1`, and equality would make their normals opposite, incompatible with an interior point. Their product is therefore zero. This proves the clique footprint assertion with its actual-face interpretation.

The perpendicular-edge endpoint formula gives exactly the displayed opposite normal. Its product with a neighboring normal vanishes precisely at vertices on that neighbor. Projection to the original facet proves distinctness of all opposites. A second actual normal on the same perpendicular would have positive product `cosh(L-L_x)` with the first; hence the uniqueness assertion used in both dimension-three and dimension-four counts is justified.

I checked both complementary-disc arguments in the polygon lemma. With no extra face, the known cyclic links force the boundary of a disc with no interior vertices; an ear has total degree four. With one extra face meeting a neighbor, the inserted extra vertex also lies on that boundary and needs degree at least four there, while every opposite vertex needs degree at least three. An ear is again impossible. In the remaining equality case a chord cuts off an empty disc; an ear away from the chord endpoints excludes the chord. This gives the claimed Loebell lattice, not an unproved classification input.

In dimension four, all neighboring and opposite facets counted inside each neighbor are exact lists. An extra facet serves at most two neighbors: three would force it to be a previously counted perpendicular opposite. A neighbor needing just one extra consumes a singleton footprint, by the polygon equality case. The final singleton/double-incidence count therefore requires at least `M` extras and proves `f>=1+V+2M=2V+5`. No presumed disjointness of excess contributions enters this count.

## Excess, volume and rigidity

Simplicity gives `e=2v`, Euler gives `r=v+f`, and the tetrahedral link gives six ridge corners per vertex. Hence `W=v-5f`. Inside a facet, trivalence and Euler give `d(F)=(V-20)/2=M-12`; summing gives `sum d(F)=2W`. The compact reflection-orbifold face weights give `chi=(v-4f+16)/16`. Chern--Gauss--Bonnet consequently yields the stated volume formula.

For a `k`-gonal ridge the neighboring facet has excess at least `2(k-5)`. Summing over the distinct neighbors, together with the original facet, gives `2W>=3d(F)`. The thresholds at `V=48` and `V=98` are respectively `f+W>=122` and `>=260`; parity and monotonicity give the asserted 46 and 96 cutoffs. The 120-cell double has `v=1160`, `f=226`, `W=30`; the two comparison volumes and strictness at the second threshold are consistent.

Zero excess makes every facet a regular right-angled dodecahedron, using the declared three-dimensional input. After matching one whole facet with the regular 120-cell, its ridge plane and the right angle determine the next supporting hyperplane. The fixed whole regular pentagon determines the adjacent regular dodecahedron on the required side. Connectedness of the two facet graphs proves equality of the complete halfspace lists. This argument requires every facet to be regular and is not applied to a single regular facet with arbitrary neighbors.

## Metric and small-graph arguments

The edge-weight bound `d(F)>=2w` excludes degree one. At degree two both weights agree and the equality polygon lemma gives the Loebell facet. Marked three-dimensional uniqueness realizes its rotations, so its base is a regular right-angled polygon. Its edge is longer than a regular pentagon edge for `k>=6`; all `2k` lateral neighbors are therefore positive. Summing the distinct selected facets gives `W>=5k-15`.

For `L_6`, the pentagon identity gives the displayed length cosines `2,sqrt(2),sqrt(3)`, with a unique longest base edge. An excess-two lateral neighbor must itself be `L_6`; its longest edge identifies the same base parent using exactly three facets at an edge interior. The two rows are disjoint selected facet sets. Their excess increments are each at least five, proving the degree-two bound `W>=20` without counting any selected facet twice.

The three-consecutive-regular-edges argument fixes the two closing supporting lines. Their orthogonality and the actual-face projection prevent another wall truncating their intersection. The two-adjacent-edge pentagon argument instead fixes two ultraparallel lines and their unique common perpendicular. In the hexagon lemma, the absence of a three-edge block forces the two remaining positions to be opposite; the displayed two equations give `cosh(a)=cosh(d)=phi^3=cosh(2L)`. These distinctions between supporting-line geometry and actual polygon closure are respected.

For each positive facet, a nonpentagonal neighbor must be eligible in its induced positive link: zero or one positive carrier, or two adjacent carriers, would leave a forbidden consecutive regular-edge block. These link graphs have no `K4` and no induced square by the normal geometry. The four-vertex graph lemma is exhaustive, giving minimum degree five in the positive adjacency graph once weighted support degree is at least three.

I reconstructed the nine-vertex complement lemma. In the triangle case the saturated three outside slots force the stated isolated vertices and independent pair. In the triangle-free nonbipartite case the shortest odd cycle is `C5`; an outside vertex has either zero or two cycle contacts, and the five contact slots leave at least two zero-contact vertices. Such vertices are isolated, since two cycle contacts cannot cover all five cycle edges. The displayed independent five follows. The smaller complement cases at seven and eight vertices also hold, giving at least ten positive facets.

For an excess-three facet with five positive neighbors I checked the full marked five-vertex analysis. A universal vertex leaves a four-vertex forest, allowing at most two marks. With maximum degree three, the triangle cases allow at most two marks, including the disjoint/nested outside-contact argument. In the triangle-free case the only configurations with three admissible marks are the internal vertices of `P5` or three vertices of `C5`; in each case a nonmarked pentagon must be regular and shares a forced `2L` edge. Both are impossible. No claim is made that the weaker abstract graph conditions alone are contradictory.

For the ten-vertex complement, deleting a degree-four vertex works unless there are two nonadjacent such vertices. In that case their common-neighbor count is 2, 3 or 4; the spare-degree and induced-disjoint-edge tests give the independent five in every case. Finally `0<W<=16` forces exactly ten positive facets and at most two excesses above three. The marked-link exclusion makes all other adjacency degrees at least six. The ten-vertex complement lemma then contradicts the absence of five mutually orthogonal wall normals in signature `(4,1)`. This proves exactly `W=0` or `W>=17`.

## External interfaces and scope

The acute-chamber all-wall inequality, standard Lorentz and convex-polytope topology, and compact orbifold Chern--Gauss--Bonnet are declared inputs. The two specifically cited three-dimensional metric inputs were checked in primary sources after reconstructing the proof:

- [Alexandrov--Bogachev--Vesnin--Egorov, Section 1.2](https://www.mathnet.ru/php/archive.phtml?jrnid=sm&option_lang=eng&paperid=9740&wshow=paper) states that 20 vertices characterize the regular right-angled dodecahedron. Trivalence converts 12 faces to 20 vertices here.
- [Roeder--Hubbard--Dunbar, Theorem 1.4, printed pp.832--833](https://aif.centre-mersenne.org/article/AIF_2007__57_3_825_0.pdf) gives compact nonobtuse uniqueness. Every polyhedron to which it is applied has more than four faces; the marked combinatorial symmetries preserve all the specified right angles.

The paper uses no finite census or search output. It explicitly leaves `f+W>=120`, a universal compact minimum-volume theorem, and the second-volume/equality problem open. Its recognition and constraints are not a classification of all compact right-angled four-polytopes. Within those stated limits, all arguments pass this source audit. No fresh TeX, Lean or scientific execution is claimed.
