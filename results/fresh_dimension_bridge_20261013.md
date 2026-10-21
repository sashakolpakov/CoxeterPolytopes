# A coupled lift from the Vinberg H19 chamber into H20

This round constructs an exact **85-facet Coxeter chamber in H20**, with all 107 inherited cusp links complete, but a strict spacelike escape proves that its volume is infinite. It therefore does **not** supply the requested new finite-volume dimension. The construction remains open to additional compatible facets: the saved escape rejects these particular 70-, 82-, and 85-wall stages, not every completion.

The mechanism changes the ambient dimension of the known H19 chamber, rather than taking another section of the Borcherds H21 chamber. It first replaces each of twelve old long walls by two coupled lifts, then completes eleven deficient cusp links simultaneously, then adds twelve exact polar truncations and three caps satisfying two prescribed nonright angles. All calculations below use the primary real square root.

## Input and explicit normals

The input is [the archived Vinberg X3 H19 chamber](../data/enumerations/vinberg_x3_h19.json), with 36 walls and 107 ideal vertices. Its lattice and chamber construction are discussed in §1 of [Looijenga, *The period map for cubic fourfolds*](https://arxiv.org/pdf/0705.0951). Here “short” and “long” refer to the unscaled squared norms 2 and 6; Looijenga uses a different scaling/naming convention for the second root class.

Negate the input's wall Gram matrix to obtain signature (19,1). The first twenty short normals form a basis. Write the original normals as r_0,...,r_35 in that basis, and adjoin an orthogonal coordinate z with z²=3. The new ambient form has signature (20,1). Inequalities have convention (x,r)≤0.

The graph on the twelve norm-six walls with edges (r_i,r_j)=−6 is connected and bipartite. Let ε_i=±1 be its bipartition, with ε_24=+1. Keep the 24 short normals at height zero and replace every long wall by

```
A_i = (r_i,  2 ε_i),   A_i² = 18,
B_i = (r_i, −  ε_i),   B_i² =  9.
```

These paired walls are orthogonal, and `(A_i+2 B_i)/3=(r_i,0)`. All 48 normals satisfy the complete Coxeter pair gate. The center h=Σ_i(r_i,0) has h²=−1518 and is strictly inside every wall.

Ninety-six inherited cusp links already gain the required extra rank through these split long walls. The remaining eleven original cusps involve short walls only: nine have type A11~+D7~ and two have type E6~³. For each archived primitive future null vector δ, add

```
C_δ,+ = ((2+√6/2) δ,  1),
C_δ,− = ((2+√6/2) δ, −1).
```

Both have square 3, their pairing is −3, and they are orthogonal to the old walls incident to δ. These 22 additions give the **70-wall stage** over Q(√6).

## Exact angle, facet, and cusp certificates

The [70-wall archive](../data/enumerations/fresh_dimension_bridge_20261013.json) freezes the form, all normals and their names, every pair classification, the strict center, all 107 cusp rays and positive null relations, and a strict escape. Its [producer](../scripts/fresh_dimension_bridge_20261013.py) uses numerical optimization only to propose an escape; the acceptance gate is exact.

The [separate Fraction verifier](../scripts/fresh_dimension_bridge_20261013_verify.py) imports no producer, SymPy, NumPy, or LP implementation. It reconstructs the original coordinates directly from the H19 Gram matrix and checks all of the following:

- The ambient inertia is (20,1), all normals are spacelike, and the center is strict and timelike.
- Every one of the 2,415 unordered pairs is valid: 1,287 right, 147 of angle π/3, 60 of angle π/4, 832 dotted, and 89 parallel.
- Each inherited null ray is feasible; its incident normals have actual span rank 20 and have a strictly positive combination equal to a positive multiple of that ray. In the horospherical quotient they span dimension 19 with a positive zero relation, so the link is bounded.
- The saved future-oriented point is strict on all 70 inequalities and has square `121580092279441/10000000000000000 > 0`.

Every wall is an actual facet, not a redundant inequality. Indeed, project the strict timelike center orthogonally to one spacelike wall. Its square stays negative, and global nonpositive pairings make all other inequalities strictly negative. The analogous projection to any positive-definite two-wall span certifies the finite-angle ridges. This argument relies on the complete acute pair gate; it would not be valid for an arbitrary nonacute list of inequalities.

The strict spacelike witness and the strict timelike center give an open set of ideal directions in the chamber, proving infinite volume. Closing all **inherited** cusp links is consequently not a complete cusp or frontier census. The replay result is [archived separately](../data/enumerations/fresh_dimension_bridge_20261013_verify.json), including rejection of a forged normal and a timelike “escape.”

## Twelve polar extensions

The main search added twelve exact polars in Q(√6). Each new normal is orthogonal to twenty independent active walls, has positive square and negative center pairing, and passes every pair check against all preceding walls. It cuts its own feasible spacelike polar ray; the first also cuts the saved 70-wall escape.

The [82-wall archive](../data/enumerations/fresh_dimension_bridge_polar_20261013.json) was replayed with the separate arithmetic in [the polar verifier](../scripts/fresh_dimension_bridge_20261013_polar_verify.py). It checks all 906 new-versus-previous pairs, the rank and orthogonality of every active set, and all 107 inherited cusps. Every new polar is strictly negative at every inherited cusp, so these bounded links are preserved.

| Pair type | 70 walls | 82 walls | 85 walls |
|---|---:|---:|---:|
| Right | 1,287 | 1,529 | 1,583 |
| π/3 | 147 | 147 | 147 |
| π/4 | 60 | 63 | 63 |
| π/6 | 0 | 17 | 26 |
| Dotted | 832 | 1,473 | 1,659 |
| Parallel | 89 | 92 | 92 |
| Total | 2,415 | 3,321 | 3,570 |

All 82 walls remain actual facets by the projection argument. However, the final saved strict escape has square **3849/32768** and center pairing **−14839/128**. The 82-wall chamber therefore still has infinite volume. The [independent replay result](../data/enumerations/fresh_dimension_bridge_20261013_polar_verify.json) also records rejection of a forged polar norm and a timelike escape. No assertion is made that all possible polar or nonpolar extensions have been exhausted.

A subsequent bank of **107 exact cap entries** solves 40 rank-20 systems: nineteen prescribed right-angle equations and one finite-angle anchor, with labels 2,3,4,5,6,8,10,12. The [separate recursive-field audit](fresh_dimension_bridge_angle_repairs_20261013_independent.md) verifies all 8,774 cap–base relations, every norm and strict center, and cutting of all fourteen distinct targeted feasible polar rays. Nevertheless the same base-82 escape is strictly feasible for **all 107 entries simultaneously**, hence no subset of this bank can yield finite volume while retaining the 82 walls. This is a finite-bank obstruction, not a universal no-completion theorem. Entries may repeat mirrors and were not required to be mutually compatible. See the [bank](../data/enumerations/fresh_dimension_bridge_angle_repairs_20261013_exact.json), [separate checker](../scripts/fresh_dimension_bridge_angle_repairs_20261013_independent.py), and [audit result](fresh_dimension_bridge_angle_repairs_20261013_independent.json).

## Three coherent caps with two prescribed angles

Allowing eighteen orthogonality equations and two prescribed finite angles produced 933 numerical proposals on three further target rays. Exactly **three selected proposals** were reconstructed and checked: each has square 2/27, eighteen right angles and two prescribed π/6 anchor angles. Each also has a third, unprescribed π/6 relation. They cut the three intended feasible hyperideal rays and are mutually dotted, so they can be added simultaneously to obtain 85 actual facets. The full 933-entry numerical scout is not an exact construction census.

The [three-cap archive](../data/enumerations/fresh_dimension_bridge_two_angle_20261013_exact.json) has an [independent recursive-field audit](fresh_dimension_bridge_two_angle_20261013_independent.md). It checks all 249 new pairs, all three rank-20 prescribed systems, strict center, actual target feasibility and cutting, and the surviving strict escape. Cross pairs are recalculated in the reverse order of radical adjunction. The [checker](../scripts/fresh_dimension_bridge_two_angle_20261013_independent.py) and [saved result](fresh_dimension_bridge_two_angle_20261013_independent.json) also reject forged normals, norms, and target pairings.

The 107 inherited bounded cusp links remain unchanged under these acute-compatible additions. To see this without an additional cusp enumeration, write a positive incident-normal relation as `Σ a_i r_i = μδ`, with every a_i>0. A new acute-compatible normal r has `(r,δ)≤0`. Equality would force `(r,r_i)=0` for every incident wall. Their span is the full 20-dimensional δ-perpendicular subspace, whose orthogonal complement is the null line spanned by δ. This contradicts r being spacelike. Hence `(r,δ)<0`.

The same exact escape of square **3849/32768** remains strict on all 85 walls. These successful local cuts therefore still do not provide finite volume.

Any finite-volume completion retaining these facets would be nonarithmetic. Already at the 70-wall stage, the normalized length-two cyclic product for walls 6 and 48 is

```
4 (r_6,r_48)² / (r_6² r_48²) = (11+4√6)/3.
```

Its field trace is 22/3, so it is not an algebraic integer. This violates the integrality condition in Vinberg's arithmeticity criterion, as stated in Theorem 2.2 of [Bogachev–Yorov, *Quasi-arithmetic hyperbolic Coxeter prisms*](https://nvbogachev.github.io/Quasi_arithmetic_Coxeter_prisms-2024.pdf). This is a conditional statement about a future finite-volume completion, not a claim that the present infinite-covolume reflection group is a nonarithmetic lattice.

## Complete classification of the underlying one-coordinate lifts

The separate [classification script](../scripts/fresh_dimension_bridge_20261013_lift_classification.py) and [exact archive](../data/enumerations/fresh_dimension_bridge_20261013_lift_classification.json) classify **all real height assignments** `(r_i,h_i)` that retain all 36 original projected wall directions as Coxeter facets, with z²=3. Positive rescaling of an individual normal does not change this scope.

The 24 short heights must vanish. Their old negative-pair graph is connected and triangle-free. Nonzero heights of one sign must be supported on a clique, so at most two can be positive and two negative. A zero height therefore exists. An old π/3 edge with one zero height and one nonzero height would have cosine strictly between 0 and 1/2, which is impossible for a Coxeter angle. Connectivity propagates zero to all 24 short walls.

Each long wall has an old π/6 edge to a short wall. Its height must consequently belong to

```
{0, ±2, ±1, ±√2(2−φ)},   φ=(1+√5)/2.
```

The seven possibilities account for all angles from π/6 through π/3, with zero recovering π/6. The script checks all 147 long-pair cases in exact Q(√2,√5) arithmetic and enumerates the complete constraint system. Its finite-angle alphabet is justified by algebraic degree, with exact cosine-polynomial and Sturm checks; it is not an arbitrary label cutoff.

There are exactly **8,195 labelled assignments**:

- 8,192 rational assignments: choose a global phase ±1 and independently choose each long height from `{2 phase ε_i, −phase ε_i}`. The paired 48-wall construction contains all choices for one phase simultaneously; the other phase is its z-reflection.
- Two golden assignments: `h_i=phase ε_i √2(2−φ)`.
- The all-zero assignment.

These are labelled assignments, not a count modulo diagram isometries. This classification does not constrain later facets with new projected directions.

Both golden assignments have a stronger obstruction. An archived positive affine D16 null relation persists, and its forced cusp is incident to an actual π/5 pair. The product of those two reflections on the 19-dimensional horosphere has trace `16+φ`, which is irrational. A bounded Euclidean reflection chamber gives a cocompact Euclidean group preserving a full translation lattice, so every such linear trace is an integer; see Bieberbach's theorem in [Gelander, *Lectures on lattices*, Theorem 2.9](https://www.math.utah.edu/pcmi12/lecture_notes/gelander.pdf). Hence **neither golden lift admits any finite-volume completion by adding acute-compatible walls while retaining its walls**.

For this universal statement the strict point used in the horospherical facet argument is a strict interior point of the proposed completion, not the saved original center. New acute-compatible walls cannot remove the forced null ray because it is a positive combination of retained normals. The rational branches, the 70/82/85-wall construction, and arbitrary H20 chambers are not excluded by this golden-branch obstruction.

## Reproduction

Run from the repository root; these five checks require only the standard library:

```sh
python scripts/fresh_dimension_bridge_20261013_verify.py --verify
python scripts/fresh_dimension_bridge_20261013_lift_classification.py --verify
python scripts/fresh_dimension_bridge_20261013_polar_verify.py --verify
python scripts/fresh_dimension_bridge_angle_repairs_20261013_independent.py --verify
python scripts/fresh_dimension_bridge_two_angle_20261013_independent.py --verify
```

All five passed, with the last two independently executed by the corresponding audit agents. The 70-, 82-, and 85-wall checks include meaningful forged-certificate rejection. The files provide exact finite input for the next search: it must add or replace walls to cut the saved new escape while preserving the full angle gate, rather than merely recount the inherited cusp links. No global novelty claim is made for these infinite-volume chambers.
