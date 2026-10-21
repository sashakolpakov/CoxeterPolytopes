# An actual local B6 quadrilateral underpaid by the global proposal

The [standalone exact witness](../scripts/fresh_sqrt5_global_b6_discount_witness_20261016.py) reconstructs the eight-root diagram from finite first-extension class 39 using Fraction arithmetic in Q(sqrt(5)). Its [data](fresh_sqrt5_global_b6_discount_witness_20261016.json) record the full scaled Gram matrix, coordinates, strict center, root relation and every endpoint test.

The physical and conjugate inertias are respectively (6,1,1) and (7,0,1). Quotienting the one relation gives an admissible seven-dimensional quadratic space, hence local hyperbolic dimension six. The four white roots span a positive-definite space. In their orthogonal three-dimensional Lorentz space, all six choices of two black active walls and both ray orientations are checked: precisely the four intended adjacent pairs are feasible, future-pointing and timelike. The remaining two wall products are strictly negative at every corner. Together with the strict common center, this gives an actual compact quadrilateral in the local chamber.

An orthogonal spectator forest of eight A1, two I2(5) and one A3 has rank 15. Appending it gives 23 root vectors of rank 22, physical inertia (21,1,1), and conjugate inertia (22,0,1). All 529 Cartan coefficients are algebraic integers. The Z[phi] span of the roots in the quotient space is a full lattice stable under their reflections. A strict center for the enlarged chamber is checked explicitly. This provides arithmetic-admissible rank-21 local face data; a compact ambient 21-dimensional chamber is not asserted.

A literal recount of all four full rank-21 corner forests gives

`payment = 1 − delta_adjacent3 + (w14894 + w21728)/10`.

For the [frozen global40 proposal](fresh_sqrt5_global_discount_20261016.md), the two bonus coefficients are zero and the exact adjacent discount is approximately 0.2296366581. The payment is therefore approximately 0.7703633419, below the required one. The exact fractions are saved, and the geometric calculation does not use the optimizer's Gram arithmetic.

This establishes a real local coverage failure of that restricted proposal. Repair must impose `(w14894+w21728)/10 >= delta_adjacent3`, or change the eligible discount pattern. The full repair also needs every other admissible stopping case, not just this witness. The proposal's exact selected-domain certificate remains valid with its stated scope, and the separate compact arithmetic dimension-21 proof candidate is unchanged. No new finite-volume polytope or dimension bound is claimed.
