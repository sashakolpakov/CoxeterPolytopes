# Exact arithmetic gate for the 196 selected triangle families

Exactly **195 of the 196 selected literal triangle families pass** the local physical and conjugate Gram tests. The sole excluded family is **core51_attach_D5_root3**, the short-leaf D5 attachment with corner types H4+D5, A1+E8, and H3+E6. Its full ten-root Gram matrix has inertia (9,1,0) at both real embeddings, so it is not arithmetic over Q(√5).

Deleting white node 2 from its conjugate matrix leaves an affine nine-root principal matrix of inertia (8,0,1). In the node order [0,1,3,4,5,6,7,8,9], an exact null vector is

\[
(4/3,\ 2/3,\ 1/3,\ 1,\ 2,\ 4/3,\ 2/3,\ 5/3,\ 1).
\]

Its product with the deleted root is (φ−1)/3, which is nonzero. A positive semidefinite extension of a singular positive semidefinite principal matrix must annihilate its nullspace; this provides a short independent reason for exclusion. The first finite-corner Schur remainders are −5/4−2φ physically and −13/4+2φ at the conjugate embedding, both negative. Orthogonal spectator components cannot remove that negative direction.

The domain is precisely the existing 196-family selection: 76 minimal cores, 68 chosen A-tail extensions, 19 I2(10)-tail extensions, and 33 additional one-component attachments of core51. All three spherical corners of every literal family are checked at both embeddings by exact positive Schur pivots. Squareclass rescaling puts the Gram entries in Q(φ); cycle-field consistency and integral scaled entries are checked by the inherited exact engine. The final one-dimensional Schur remainder determines the full-star inertia. The other 195 families pass these necessary local Gram tests; this does not itself prove their global chamber realizability.

The [source](../scripts/fresh_sqrt5_marked_triangle_arithmetic_gate_20261016.py) and [certificate](fresh_sqrt5_marked_triangle_arithmetic_gate_20261016.json) passed generation and a fresh saved replay. The certificate stores every literal row, all three finite-corner pivot sequences, both full-star inertia results, and the explicit affine-principal witness for the removed row.

Earlier payment-model obstructions remain correct for their explicitly broader, overcomplete domains. They must not be reinterpreted as geometric obstructions when they use the excluded family. Conversely, a feasible payment rule on all 196 families also covers the retained 195. Arbitrary other white extensions and formerly paid quadrilaterals still require complete coverage; this gate makes no dimension-bound claim.
