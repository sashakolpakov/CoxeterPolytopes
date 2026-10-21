# Finite-vertex and coupled-cap completion in H9

This round reaches **two exact 73-facet Coxeter chambers in H9**, with new finite vertices and finite angles between newly constructed facets. Both still have infinite volume, certified by strict future spacelike escapes. No new finite-volume example or dimension record is claimed.

The substantial change is to release the old requirement that every added cap have a finite angle with polar wall `old348`. Caps orthogonal to eight retained walls can have a dotted pairing with that polar wall. This creates actual finite vertices, which the preceding cap additions had not done. A second construction imposes finite angles between these new vertex caps.

## The coherent starting progression

Start with case 0 of the independently checked [rebuilt 45-wall chamber](fresh_h9_replaced_core_20261013.md), already outside the earlier fixed-core cusp obstruction.

The newest retained cap, `round4_caps[2]` from that round, supplies a different quartic anchor:

`Q(phi, sqrt(24(phi−1)))`.

Three free seed distances, a finite angle with this anchor, and a finite angle with `old348` leave a quadratic norm equation. Six exact caps pass all 45 retained-wall pair gates. Four cut all twelve saved escapes. Two coexist, giving a 47-wall chamber. Eight tested configurations have exact strict escapes.

The seed has the exact symmetry exchanging coordinates 6 and 7. Transporting old `new6` by this symmetry and replacing that one wall cuts the leading 47-wall gap while preserving all other 46 walls. This replacement passes the full pair and affine gates. A different proposed replacement using transported `old496` creates an affine subdiagram and is rejected.

A fresh use of the retained left quartic anchor, over `Q(phi,sqrt(132−60phi))`, adds a cap to the changed 47-wall configuration. Four of six tested 48-wall configurations are affine-free; two are rejected for positive affine diagrams. All four accepted chambers still have exact strict escapes. This establishes the authoritative 48-wall base used below.

Sources and exact replays:

- [Initial gap geometry](../data/enumerations/fresh_h9_gap_completion_20261013_geometry.json): eight configurations, 8,510 pair checks, 5,583 affine inertia checks, 96 strict escapes.
- [Symmetry replacements](../data/enumerations/fresh_h9_gap_completion_20261013_swap.json): three configurations, 3,197 pair checks, 2,057 affine inertia checks, 24 strict escapes.
- [Changed left-anchor geometry](../data/enumerations/fresh_h9_gap_completion_20261013_left_geometry.json): six configurations, 6,768 pair checks, 4,452 affine inertia checks, 48 strict escapes.
- [Producer-free replay](../scripts/fresh_h9_gap_completion_20261013_verify.py), with corresponding `*_verify.json` outputs.

The numerical vertex census retains only fifteen finite vertices through this progression. Thus cutting the saved escapes was not, by itself, creating the finite vertices needed for a compact chamber. That observation motivated the following change of construction.

## Creating finite vertices directly

Choose eight retained normals whose Gram matrix is positive definite. Their orthogonal complement has Lorentzian signature `(1,1)`. A spacelike normal in this two-dimensional complement gives eight right angles. Its norm-two hyperbola still has one free real parameter.

The numerical proposal scan finds 147 rank-eight positive-definite blocks. Of these, 129 admit sampled caps dotted to every other old wall. These are proposal counts, not an exhaustive numerical nonexistence theorem. The [exact cap archive](../data/enumerations/fresh_h9_gap_completion_20261013_vertex_caps.json) certifies thirty representatives whose eight defining normals are all over `Q(phi)`. Every one has:

- exact norm two and a strict common timelike interior;
- exactly eight old right-angle relations, with all other old pairings dotted;
- a positive-definite rank-eight defining Gram matrix;
- an explicit future timelike point on the eight old walls and the new cap, strict on every other old inequality.

The last point is an actual finite vertex. It is stronger than merely recognizing a positive Gram subdiagram numerically.

For example, the first cap is the normalization of the raw seed-pairing vector with

`p7=−173/702, p8=−1`,

and all other entries zero. Its eight orthogonal walls are seeds `0,1,2,3,4,5,6,9`, an `A7+A1` block. If `I=G⁻¹` and `q=pᵀIp`, multiply the raw normal by `sqrt(2/q)`. The saved normalization radicand is

`(1655743084164+2686162737492 phi)/21732829031`.

Its finite-vertex witness has exact negative norm

`−(124155089757429+200960248200337 phi)/173862632248`.

The full point, future-orientation check, and every strict/zero inequality are in the archive. In this example `old348` is dotted; the earlier finite-polar-angle families could not express it.

All 435 pairs among the thirty caps are checked exactly. There are 421 compatible pairs, all dotted; five others have positive pairing and nine have non-Coxeter finite-angle values. For this particular bank, every normalized squared pairing lies in `Q(phi)`, so a possible finite label must have real cyclotomic degree at most two. The gate includes all such labels (`2,3,4,5,6,8,10,12`), making this thirty-cap compatibility graph exact. Twenty-two caps can coexist. Four tested 70-wall chambers retain the twenty-two new finite vertices but still have exact strict escapes. Their numerical finite-vertex count is 37.

The [vertex geometry replay](../scripts/fresh_h9_gap_completion_20261013_vertex_verify.py) checks all thirty individual finite vertices, all cap pairs, and all four selected chambers: 9,660 chamber pair checks, 3,012 affine inertia checks, 48 strict escapes, and 88 checks that the selected finite vertices survive the other selected caps.

## Coupling new facets at finite angles

Generic caps chosen independently on those hyperbolas meet only by dotted relations in the accepted graph. To introduce a stronger geometric constraint, hold one new vertex cap fixed and choose a second rank-eight block sharing seven old normals with the first.

The second cap is orthogonal to its eight old normals. Its norm equation and its prescribed angle `pi/2`, `pi/3`, or `pi/5` with the first cap give a quadratic equation over the first cap's quartic field. The shared seven old normals and the two new caps then form a positive-definite rank-nine diagram whenever the new pair has a finite Coxeter label. This creates a finite vertex involving **both** newly constructed facets.

The [coupled-cap archive](../data/enumerations/fresh_h9_gap_completion_20261013_coupled_vertices.json) contains fifteen exact caps on two anchor stars. All have eight old right angles and are dotted to every other old wall. The finite-angle equation with the new anchor is verified separately.

Joint selection then gives [two 73-wall chambers](../data/enumerations/fresh_h9_gap_completion_20261013_coupled_geometry.json). Each retains the changed 48-wall base, twenty-one generic vertex caps, and four coupled caps. They include actual new/new `pi/3` relations as well as right-angle attachments.

Case 0 saves three joint finite vertices, and case 1 saves two. Each witness has exactly nine zero inequalities, all others strict, negative norm, correct future orientation, and an exactly positive-definite incident Gram matrix. These are explicit actual-vertex certificates. The numerical total finite-vertex counts are 43 and 42 respectively; a complete exact enumeration of every finite vertex is not claimed.

Every chamber pair is exactly Coxeter-compatible and nonparallel. The complete connected affine census finds no affine principal subdiagram, testing 757 connected subsets in each final chamber. As in the preceding round, the census is exhaustive because it extends every connected positive-definite subset and every proper principal subset of an irreducible affine diagram is positive definite.

The [coupled producer-free replay](../scripts/fresh_h9_gap_completion_20261013_coupled_verify.py) passes:

- fifteen coupled-root equations and all old-wall/anchor relations;
- 5,256 pair checks across the two final chambers;
- 1,514 affine inertia/frontier checks;
- all five actual joint finite vertices;
- twenty-four strict future spacelike escapes.

All claimed facets are genuine: a common strict timelike interior and nonpositive pairings allow its Lorentz-orthogonal projection onto any one wall while keeping every other inequality strict.

A [second arithmetic audit](fresh_h9_gap_completion_independent_geometry_20261013.json), using a separately implemented recursive radical representation, checks all 73 normal/Gram identities, the strict center, and the first strict escape in case 0. It also checks the actual joint vertex at `vertex12` and `coupled13`: all 73 inequalities, exactly nine active walls, its physical positive-definite rank-nine Gram matrix, and the new `pi/3` pairing. Its [replay](../scripts/fresh_h9_gap_completion_independent_geometry_20261013.py) uses a frozen exact input bundle and imports no construction producer. The complete pair-angle and affine claims remain covered by the separate producer-free replay above.

## The remaining gap

Neither final chamber has finite volume. Their first common strict escape has exact positive norm

`−6343656658040413429/5000000000000000000
+ (8190710171640288811/10000000000000000000) phi`.

Its unperturbed numerical proposal has seed support `[2,3,6,8,9]`. The actual certificate is a rational strict perturbation, checked against all 73 walls and against the future center. This remains the next constructive target; it is not a proof that further completion is impossible.

The numerical maximum normalized vertex norm falls from about `0.0633932` at the starting 45 walls to `0.0567845` at the changed 48-wall base. The finite-vertex and coupled additions do not reduce that leading value further. Their contribution is certified finite-vertex creation and new inter-cap angle constraints, not volume closure.

## Reproduction and exact scope

```sh
PYTHONPATH=. python scripts/fresh_h9_gap_completion_20261013_verify.py --verify
PYTHONPATH=. python scripts/fresh_h9_gap_completion_20261013_verify.py --swap --verify
PYTHONPATH=. python scripts/fresh_h9_gap_completion_20261013_verify.py --source data/enumerations/fresh_h9_gap_completion_20261013_left_geometry.json --verify
PYTHONPATH=. python scripts/fresh_h9_gap_completion_20261013_vertex_verify.py --verify
PYTHONPATH=. python scripts/fresh_h9_gap_completion_20261013_coupled_verify.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_gap_completion_independent_geometry_20261013.py --verify
```

These exact replays do not import this round's numerical search or geometry producers. The bounded searches, numerical discovery filters, rational parameter choices, and selected compatible subsets are not exhaustive over their continuous families. The finite-angle pool is an acceptance alphabet for the general composita, not a universal exclusion of every possible label. No comparison with all historical cap banks and no claim of global mirror novelty is made.
