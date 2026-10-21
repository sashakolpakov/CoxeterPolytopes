# Independent review of the D4+2E8 local edge

**PASS, generation and fresh saved replay.** The
[independent checker](../scripts/fresh_sqrt5_d4_two_e8_self_edge_20261016_review.py)
and [certificate](fresh_sqrt5_d4_two_e8_self_edge_20261016_review.json) rebuild
the literal22-by22 Gram matrix in project Golden arithmetic, without importing
the [producer](fresh_sqrt5_d4_two_e8_self_edge_20261016.md).

It checks determinant `−4(1+3φ)`, physical signature(21,1), definite conjugate,
integral reflection coefficients and the field-generating cycle−φ. It computes
both endpoint rays from the inverse Gram, checks their complete incidence and
positive-definite active diagrams, and obtains `cosh(edge length)=1+φ`.
The common edge diagram is D4+2E8 and both endpoints are D5+2E8.

The infinite-volume witness is reconstructed independently: a strict timelike
core vector u, an orthogonal positive E8 vector h, and
`t²=(169φ−176)/62000>0` give `u+t h` null for the positive real square root t.
Its core inequalities are strictly negative and all sixteen E8 inequalities
equal−t. A strict timelike interior point is also checked. Thus the chamber
has an open ideal-boundary patch. No K-rational null vector is asserted;
that would contradict the definite conjugate.

This shows that the recurring D4+2E8 edge with D5+2E8 endpoints has actual
V-arithmetic compatibility. Eliminating that vertex type requires more than
declaring this edge impossible. It does not settle all incidences around a
vertex or their compatibility in a finite-volume chamber. It supplies no new
finite-volume example, arithmetic reflection lattice, or dimension bound.

```sh
python3 scripts/fresh_sqrt5_d4_two_e8_self_edge_20261016_review.py
python3 scripts/fresh_sqrt5_d4_two_e8_self_edge_20261016_review.py --verify
```
