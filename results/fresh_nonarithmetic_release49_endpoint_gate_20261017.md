# Two mandatory endpoint gates for a released wall49

The proposed three-pentagon replacement of wall49 is impossible while the other63 walls of the exact64-wall control are retained. This is an all-label obstruction to that particular hexagonal patch. It does not exclude different patches, further wall replacements, or walls missing the face S6.

The [checker](../scripts/fresh_nonarithmetic_release49_endpoint_gate_20261017.py) and [saved certificate](fresh_nonarithmetic_release49_endpoint_gate_20261017.json) use the frozen [complete intrinsic boundary](fresh_nonarithmetic_replace47_close12_boundary_20261017.md) and its exact coordinates. The underlying chamber has198 finite ambient vertices and15 open intrinsic ends. No chamber is changed here.

## Mandatory endpoints and complete contact domain

Write S6={0,1,2,3,5,6}=H4+A2. The side49 has cyclic neighbor list [7,33,30,59,60,58]. A three-pentagon disk with each new wall meeting three consecutive old sides has two alternating placements. The first necessarily has one wall making both edges [60,58] and [58,7] finite; the second necessarily has one wall making both [58,7] and [7,33] finite.

Any actual S6 facet is orthogonal to its H4 component and has one of seven contacts to its A2: both right, or a single bond3,4,5 to either A2 node. This is the finite-extension classification already used in the frozen endpoint domain. Both specified endpoints must be finite in the proposed compact patch; their positive-definite diagrams consist of S6, the two old side roots and the new root.

There is no unexamined large-label tail. Each of the three anchors has an existing nonright neighbor in a required finite corner:60 has a5-bond to6,7 has a3-bond to6, and58 and33 have a3-bond to7. If the new wall had label m>=6 to an anchor, the three relevant norm2 roots would have a positive principal determinant requirement but determinant

`8−2(a²+b²+c²)−2abc <= 8−2(3+1)=0`,

where a>=sqrt(3), b>=1 and c>=0 are the nonnegative contact magnitudes. This is impossible. Thus every new-anchor label is2,3,4 or5, in every coefficient field.

For each frame the checker examines all7×4³=448 assignments. After removing the orthogonal H4, each required corner has rank5. It recognizes every finite connected component by exact graph isomorphism to A1–A5, B2–B5, D4/D5, F4, I2(5), H3 or H4, including arbitrary orthogonal products. This uses the standard complete finite Coxeter classification, not a numerical positive-definiteness test.

## Exact exclusion of both branches

In the residual Lorentz four-space, three prescribed native contacts determine an affine line. If G is the Gram matrix of its three projected anchors, the particular solution is their linear combination with coefficients G^-1 times the contact vector. Projection of an independent old root gives a positive-square kernel k. The norm2 equation, after subtracting the A2 component, is

`r = p ± sqrt((q−p²)/k²) k`.

The saved interval operations enclose these exact expressions. Their algebraic projection definitions establish the contact and norm identities; overlapping interval evaluations are consistency checks, not substitutes for those identities. Each real branch has a strictly positive product with a retained old wall. Such a normal cannot be part of any acute Coxeter extension, regardless of a common-center choice.

| Mandatory frame | Finite-corner recipes | Real branches | Strict obtuse blockers used |
|---|---:|---:|---|
| [60,58,7] | 5 | 10 |41,61,62|
| [58,7,33] |13|26|41,62|

All878 other assignments fail the exact finite-diagram test. The18 surviving recipes have strictly positive norm residuals, so each has exactly two real branches; there are no unresolved tangent or degenerate cases. Every blocker is retained when49 is removed.

The associated [bounded scout and deletion ledger](fresh_nonarithmetic_release49_hexagon_20261017.md) supplied the proposed topology. It is not used as acceptance evidence in this exact gate. Deleting49 loses one old finite vertex and no certified compact face of dimension at least2. The deletion-only reduction in the number of recorded open ends enlarges the chamber and is not construction progress.

## Verification and scope

**Final/frozen.** Generation and a fresh saved replay both terminate with PASS. The checker exhausts896 finite-label assignments,18 exact root recipes and36 exact obtuse-branch certificates. The high-label isolation argument received a separate narrow conceptual check; no independent interval replay is claimed. The physical interval arithmetic and exact parent geometry are inherited explicitly through source hashes. Source syntax, direct source hashes and local report links pass.

This proves failure of the two alternating three-pentagon hexagon disks retaining the specified63 walls. It is not a no-completion theorem for the underlying H9 chamber. No new finite-volume example, ambient cone census or dimension assertion is made.
