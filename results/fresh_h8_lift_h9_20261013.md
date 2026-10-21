# Exact single-coordinate lifts of the compact Bugaenko H8 chamber

There are two nontrivial Coxeter lifts of the eleven old wall directions into one extra spatial coordinate, up to the diagram symmetry and reversal of the new coordinate. Both have the desired H9 Gram rank, but both contain explicit affine null-root combinations. Consequently **no completion by additional walls acute-compatible with these retained walls can be compact**, regardless of the number or heights of the added walls.

This differs from the earlier orthogonal prism-cap search: the old wall Gram changes and its rank rises from 9 to 10. The obstruction below also differs from literal-facet tilt rigidity: the original H8 section is not required to be a facet of the new polytope.

## Classification of the lifts

Normalize the archived H8 normals `r_i` to square 2. Consider `R_i=(r_i,h_i)` with one extra positive Euclidean coordinate. Positive rescaling of each whole normal does not change this formulation, provided its old projected direction is retained.

The Bugaenko diagram has finite labels 3 and 5, with no label-4 edge. Its graph of negative pairings has no triangle. If two nonzero heights have the same sign, acuteness forces their old pairing to be negative. Hence each sign-support is a clique and has at most two vertices. At most four heights can be nonzero.

If an old label-3 edge has one endpoint of height zero, the other endpoint must also have height zero: otherwise its new nonzero cosine is strictly between 0 and 1/2, which is not `cos(π/m)` for any integer `m≥3`. The label-3 graph has components

```
{0}, {1,2,3,4,5,6,7,9,10}, {8}.
```

The nine-node component cannot have all heights nonzero, so this zero-height condition propagates throughout it. Only endpoints 0 and 8 can lift.

Each endpoint has exactly one nonright old edge, of label 5. Increasing its squared norm decreases that cosine, so the only nontrivial possibilities are labels 3 and 4. The required added squared lengths are

\[
5\to3:\quad h^2=2\phi,
\qquad
5\to4:\quad h^2=\phi-1.
\]

Both endpoints cannot lift together. Their heights would need opposite signs, because their old normals are orthogonal. The resulting values of `4cos²` for the new mutual pair, for endpoint labels `(3,3),(3,4),(4,4)`, are respectively

\[
8-4\phi\in(1,2),\qquad
20-12\phi\in(0,1),\qquad
52-32\phi\in(0,1).
\]

None is an allowed finite Coxeter value. Thus precisely one endpoint lifts. Reversing the diagram exchanges endpoints 0 and 8, leaving two seeds.

## The two exact H9 seeds

For each seed, append a coordinate of squared length `k` to the archived H8 form, append coordinate 1 to old root 0, append coordinate 0 to every other root, then divide the whole lifted root 0 by `φ`.

| New label on edge 0–1 | `k` | Lifted root-0 square | Gram signature |
| --- | --- | --- | --- |
| 3 | `2φ` | 2 | `(9,1,1)` |
| 4 | `φ−1` | 1 | `(9,1,1)` |

Both eleven-wall Grams have exact rank 10. The old strict timelike interior point, extended by a zero final coordinate, remains strict. The replay verifies the host signature using the positive principal minors of the E8 spatial block and the positive new coefficient.

For the label-3 seed, the walls `0,1,2,3,4,5,6,9` have an affine E7 Gram. In that order, their positive null coefficients are

\[
(1,2,3,4,3,2,1,2).
\]

Their combination is the future null vector

\[
\delta_3=(1,2,3,4,6,5,4,3,2,\phi-1).
\]

Its only nonzero pairings with the eleven seed walls are `−1` at wall 7 and `−2−4φ` at wall 10.

For the label-4 seed, the walls `0,1,2,3,4,9` have an affine B5 Gram. Their positive null coefficients are

\[
(1,1,1,1,1/2,1/2).
\]

Their combination is

\[
\delta_4=(1,3,9/2,6,9,15/2,6,4,2,\phi-1).
\]

Its only nonzero seed-wall pairings are `−1/2` at wall 5 and `−φ` at wall 10. Both displayed null vectors have negative pairing with the strict timelike center.

## Why extra compatible caps cannot make these compact

Each `δ` is a positive combination of the displayed seed normals. Every additional normal `s` acute-compatible with all retained seed walls therefore satisfies `(s,δ)≤0`. Hence the same future null ray lies in every such enlarged system of inequalities. If the final intersection is full-dimensional, a path from an interior point toward that ray reaches arbitrarily far into hyperbolic space. The polytope cannot be compact.

This conclusion allows any number of new cap or truncation walls and places no bound on their fields, coefficients, or angles. It assumes the eleven old projected directions remain the retained Coxeter walls. Removing those walls, allowing them to become redundant while changing their incompatible pairings, or changing their old projected directions is outside the argument. Noncompact finite-volume completions are not excluded here and were not pursued further in this batch.

## Exact artifacts

- `scripts/fresh_h8_lift_h9_20261013.py`
- `data/enumerations/fresh_h8_lift_h9_20261013.json`

The replay uses standard-library `Fraction` arithmetic over `Q(φ)`. It checks the graph restrictions, both full Gram realizations, signatures, strict centers, affine ranks, all positive coefficients, and every null-vector pairing. The following command passed:

```sh
PYTHONPATH=. python3 scripts/fresh_h8_lift_h9_20261013.py --verify
```
