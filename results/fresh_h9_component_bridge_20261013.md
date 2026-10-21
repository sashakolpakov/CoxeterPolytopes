# Trivalent endpoints outside the closed symmetry family

Two new exact caps retain all 102 preceding walls and cut both saved escapes. They close two old open edges at actual finite vertices. The resulting **104-facet H9 chamber still has infinite volume**, proved by a new strict future spacelike escape. This is not a new compact or finite-volume Coxeter polytope.

The construction changes the endpoint diagram search. Earlier local scouts used at most two nonright relations to the eight walls of a spherical edge. Here the new node has three relations of label 3. A numerical scout considered the 272 spherical rank-eight subsets it found in the old chamber, then fixed one additional old-wall angle. It produced proposals, not an exhaustive all-angle classification. Two compatible proposals were reconstructed exactly; no old wall was removed.

| New cap | Eight-wall edge | Its three label-3 relations | Extra finite anchor |
|---|---|---|---|
| `trivalentRemote` | seeds 1,3,4,5,7; wall22; walk2720; second33 | seed3, seed5, second33 | wall54, label 3 |
| `trivalentLeading` | seeds 0,1,2,3,6,7; coupled0; coupled1 | seed1, seed3, coupled1 | wall17, label 5 |

All other relations to those eight edge walls are right angles. Each cap has exactly nine finite old neighbors and 93 dotted old relations. The new caps are mutually dotted. Thus the 205 new pairs comprise 10 right angles, seven label-3 angles, one label-5 angle, and 187 dotted pairs.

The first cap cuts the exact common escape of the entire previous 13,824-image family. Consequently its mirror lies outside that family. The second cuts the leading escape of the previous 102-wall chamber. The first does not cut that second target; the second does not cut the symmetry-family target. Their joint use is deliberate.

## Exact construction and geometry

Let (A) be the eight retained edge normals, (M) their positive-definite Gram matrix, and (p) the prescribed vector of zero and minus-one pairings. Set (u=A M^{-1}p). Project the extra anchor (a) to (E^\perp), where (E=\operatorname{span}(A)), obtaining (b). If the desired anchor pairing is (c), define

\[
\gamma=c-(u,a),\qquad w=u+\frac{\gamma}{(b,b)}b.
\]

Project the common strict center to (E^\perp\cap b^\perp), obtaining (z\ne0). The two possible normals are

\[
r=w\pm\sqrt{\frac{2-(w,w)}{(z,z)}}\,z.
\]

All denominators are certified nonzero and both selected radicands positive. These projection identities prove norm two and all nine prescribed pairings exactly. Outward integer intervals certify branch orientation, all remaining dotted relations, target cuts and geometry. Merely enclosing zero is used only as a consistency check, never as proof of a finite-angle equality.

The common strict timelike center and nonpositive pairings imply every retained inequality is an actual facet: project the center orthogonally onto that wall; the point remains timelike and satisfies every other inequality strictly.

The complete new affine gate checks all 1,024 subsets drawn from either new wall and its nine finite neighbors. It excludes 352 subsets by a dotted pair, finds 564 positive-definite subsets and 108 indefinite subsets, and finds **no positive-semidefinite singular subset**. A subset containing both new walls is already indefinite because their pair is dotted. The previous 102-wall affine census is an explicit dependency.

Exactly two new rank-nine positive-definite systems occur, of types E7+A1² and D4+5A1. Projection of the common center gives actual finite vertices, with every nonincident wall checked strictly. The new endpoints close these old edges:

- the edge of `trivalentRemote` previously had finite end wall `seed6`;
- the edge of `trivalentLeading` previously had finite end wall `seed5`.

Every ninth-wall possibility on the 18 edges incident to these two new vertices is checked. Two edges are now finite; 16 still run to the light cone. This advances two inherited vertex components; it does **not** join those components. Existing label-10 relations in the endpoint gate are handled by rigorous nonzero LDL pivots over their actual radical values, not by treating them as elements of Q(phi).

## Remaining volume gap

The numerical affine-section census finds 80,102 projective vertices, 80,065 spacelike, and maximum squared norm approximately 0.04645. These counts are diagnostics, not a compactness proof. A rational perturbation of one such vertex is independently verified against all 104 inequalities. Its exact norm is

\[
\frac{4552626030126814721772817-2354371445870741771052801\phi}
 {16000000000000000000000000}>0,
\]

and its center pairing is (-535537826939/250000000000<0). Every wall inequality is strict. This proves infinite volume. The unperturbed proposal was incident to seeds 2,5,6,7,11; wall17, wall22, wall57; and bridge0. Those data are a new target, not an assertion that one further cap suffices.

The exact [recipe](../data/enumerations/fresh_h9_component_bridge_20261013_caps.json) and [certificate](../data/enumerations/fresh_h9_component_bridge_20261013_verify.json) are replayed by the [producer-free verifier](../scripts/fresh_h9_component_bridge_20261013_verify.py). It replays the preceding exact dependencies, with deep-copy memoization to avoid repeated identical nested checks. The [producer](../scripts/fresh_h9_component_bridge_20261013.py) uses numerical geometry only to propose the final exact escape.

```sh
python scripts/fresh_h9_component_bridge_20261013.py
python scripts/fresh_h9_component_bridge_20261013_verify.py --verify
```

The geometric search remains unrestricted outside this particular construction. Neither the numerical trivalent scout nor this escaped chamber establishes a general obstruction to compact H9.

A subsequent [all-angle frontier certificate](fresh_h9_component_bridge_20261013_continuous.md) excludes a cap that both cuts the current escape and creates a finite endpoint on one of the 16 new open edges while retaining all 104 walls. It covers all 940 profiles with edge labels 2,3,4,5 and separately covers every higher-label isolated dihedral case by 51 exact two-wall dual bounds. It does not exclude other frontier components, replacement, or noncutting intermediate caps.
