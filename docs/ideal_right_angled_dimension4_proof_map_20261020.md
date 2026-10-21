# Sharp entirely ideal right-angled dimension bound

This is the preserved predecessor map. The current primary is the
[unified sharp8/4 manuscript and reproduction map](right_angled_sharp_dimensions_proof_map_20261020.md),
which shortens the ideal proof using one finite-point projection.
The old graph-propagation Lean source was reviewed but not compiled;
its deferred core run is not a Lean verification claim.

20 October 2026. The complete written proof and independent mathematical
reviews establish dimension at most four; the regular ideal 24-cell
shows sharpness. This page separates the actual proof from the finite
model checker and the conditional Lean interface.

## Result and short manuscript

An entirely ideal finite-volume right-angled hyperbolic polytope has
dimension at most four. There is no arithmetic or restricted facet-type
assumption. The new argument excludes dimension five directly. Every
higher-dimensional example would have an entirely ideal right-angled
five-face, giving the same contradiction.

The focused [manuscript](../manuscripts/ideal_right_angled_dimension4.tex)
has a complete [authoritative proof note](../results/ideal_right_angled_dimension4_proof_20261020.md).
The previous ideal5/6 size and cusp-compatibility papers remain preserved
supporting history. None is needed for this proof. The separate ideal4
volume/structure paper keeps its existing scope and frozen evidence.

The complete [root review](../results/review_ideal_right_angled_maximum_height_kissing_root_20261020.md)
and two independent adversarial reviews
([first](../results/review_ideal_right_angled_maximum_height_kissing_batch_guard_20261020.md),
[second](../results/review_ideal_maximum_cusp_height_finite57_20261020.md))
pass. The initial candidate's Kolpakov theorem number is corrected in
the authoritative note and manuscript: the journal statement is Theorem4.

## Published inputs

|Input|Exact statement used|Primary reference|
|---|---|---|
|Four-dimensional kissing bound|A finite subset of S3 with distinct-pair scalar products at most1/2 has at most24 points.|[Musin, Annals168(2008), Main theorem, printed p4](https://annals.math.princeton.edu/wp-content/uploads/annals-v168-n1-p01.pdf)|
|Ideal4 minimum and equality|An entirely ideal finite-volume right-angled four-polytope has at least24 facets; equality implies isometry to the regular ideal right-angled24-cell.|[Kolpakov, AGT12(2012), Theorem4, printed p1945](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf)|
|All-wall acuteness|Distinct actual outward unit normals have nonpositive Lorentz product, including nonadjacent facets.|[Vinberg, Hyperbolic reflection groups(1985), Section2.1](https://www.mathnet.ru/links/8684a695fdf254b4efd6aa83de004ad9/rm2140_eng.pdf)|

The other face and graph facts are standard convex geometry: actual
facets inherit ideal vertices and finite volume, orthogonal walls
restrict to unit intrinsic walls, right angles are inherited in
dimensions at least three, and the actual ridge-adjacency graph of
facets is connected.

## Internal proof chain

1. Choose a cusp vector v and heights H(m)=-v.m. They are nonnegative,
   some are positive, and cusp-incident walls have height zero.
2. At a positive globally maximal-height wall n, let q=v/H(n)+n.
   Lorentz algebra and all-wall acuteness prove q is an actual future
   unit point in the facet Q supported by n.
3. Every intrinsic wall m_i of Q is an actual orthogonal neighbor of n.
   Relative to q its time coordinate is t_i=H(m_i)/H(n) in [0,1].
   Its spatial unit direction u_i satisfies u_i.u_j<=1/2.
4. Musin gives at most24 walls. Kolpakov gives at least24 and identifies
   Q with the regular ideal24-cell.
5. Opposite standard normals (1,r),(1,-r) force q to be its center:
   both t<=1 imply a<=1 for q=(a,w), whereas a²-|w|²=1 gives a>=1.
   Hence all neighboring heights equal H(n).
6. The same conclusion holds at every positive maximum-height facet.
   Propagate along the connected actual facet graph to a cusp-incident
   height-zero facet, a contradiction.

No combinatorial enumeration or kissing-code uniqueness is assumed.
Kolpakov supplies metric rigidity of Q; this is essential in step5.
Maximality is a choice among finitely many actual walls, not a
restriction on a selected subset of vertices.

## Verification boundaries

The exact [model checker](../scripts/check_ideal_right_angled_dimension4_model_20261020.py)
uses integer and rational arithmetic. It checks the24 standard unit
normals, all276 Gram pairs,12 opposite pairs,24 displayed null vertices,
their cubical links, actual facet witnesses and a scalar polynomial
identity for the code reduction. It does not enumerate all possible
ideal polytopes or prove the external theorems or universal geometry.

The conditional [Std Lean kernel](../formal/CoxeterProofs/IdealRightAngled4.lean)
checks the cardinality squeeze, propagation along finite paths and the
final dimension implication. Geometric and published input interfaces
remain explicit. A compilation pass will not be described as an
end-to-end Lean formalization of hyperbolic geometry. Actual execution
ledgers and the one-command reproduction entry will be recorded after
source reviews and resource-gated batches.
