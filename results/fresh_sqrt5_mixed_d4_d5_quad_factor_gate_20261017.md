# Arithmetic factor gate for mixed D4/D5 quadrilaterals

The free-product determinant equations in the [complete finite-pattern census](fresh_sqrt5_mixed_d4_d5_quad_patterns_20261017.md) give the following exact reduction. All rational univariate branches are impossible; the remaining branches are coupled equations, not certified quadrilaterals.

| Vertex forest | Original colored patterns | Patterns with a coupled branch | Remaining literal patterns | Remaining possible corners |
|---|---:|---:|---:|---:|
| 3D4 | 5 | 2 | 18 | 18 |
| D5+2D4 | 38 | 19 | 90 | 30 |
| 2D5+D4 | 68 | 19 | 114 | 34 |
| 3D5+A2 | 60 | 18 | 127 | 37 |

At an actual compact edge, the d finite normals at one endpoint span a positive-definite space. The new wall at the other endpoint does not contain the first endpoint, so its normal lies outside that space. These d+1 normals therefore have a full-rank Lorentzian Gram matrix. In the present patterns every entry of this edge Gram is rational except the opposite framing product. If that product were rational, the entire normalized Gram would be rational.

This is incompatible with arithmeticity over Q(sqrt5). On each connected block the rational norm-two entries give a rational Gram gauge. Conjugation in a field-valued root gauge gives a matrix diagonally congruent to that same rational matrix: the common norm factors have positive conjugates, and possible normalization sign changes are a diagonal sign congruence. Thus its physical negative direction persists, contradicting positivity at the nonidentity real embedding. This argument applies to either framing edge and does not assume the two free products individually belong to the ground field in every pattern.

The [checker](../scripts/fresh_sqrt5_mixed_d4_d5_quad_factor_gate_20261017.py) independently rebuilds every four-root Schur matrix from the literal labels using SymPy rational arithmetic. It compares the complete polynomial matrix and determinant with the frozen input, factors the determinant exactly over Q, and checks the factor product. Every factor involving only one variable is linear with a rational root. Each such branch is excluded by the preceding edge argument. A pattern is retained precisely when it has a factor involving both variables. Literal corner counts are then recomputed from all saved attachment maps.

Generation and fresh saved replay pass, with [exact factor records](fresh_sqrt5_mixed_d4_d5_quad_factor_gate_20261017.json):

```sh
python scripts/fresh_sqrt5_mixed_d4_d5_quad_factor_gate_20261017.py
python scripts/fresh_sqrt5_mixed_d4_d5_quad_factor_gate_20261017.py --verify
```

No numerical parameter cutoff or minimality assumption is used. The retained coupled factors must still satisfy physical signature, conjugate positivity, Coxeter angles, integrality and face coherence. This package gives an upper list of possible all-special corners, not a realization, a mixed incidence inequality or a dimension bound.
