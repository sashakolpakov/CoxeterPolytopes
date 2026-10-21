# Independent exclusion of four orthogonal edge-cutting families

The [primary four-axis gate](fresh_h20_h4_edge_escape_gate_20261016.md)
passes a separate [standard-library interval verifier](../scripts/fresh_h20_h4_edge_escape_gate_20261016_review.py)
with its own [certificate](fresh_h20_h4_edge_escape_gate_20261016_review.json).
All geometric bounds come from the
[independent 145-wall reconstruction](fresh_h20_residual_h4_corner_20261016_review.md),
not the primary gate or producer physical coordinates. Generation and fresh
saved replay pass.

A spacelike additional normal right to the sixteen support roots and three
of the four orthogonal residual caps has the form

    alpha*(-n_i + lambda*C), alpha > 0,

where C is the finite corner. The rank-19 orthogonal support leaves precisely
this Lorentz two-plane. Acuteness to n_i makes its coefficient nonpositive;
a zero coefficient would be timelike and is excluded.

For three choices of n_i, acuteness to old wall 133 forces lambda strictly
above the upper limit allowed by nonnegative pairing with E4. On the fourth
axis, acuteness to chain wall 132 leaves a closed interval. Throughout its
interior and upper endpoint the normalized cyclic invariant
4*(r.n132)^2/(r^2*n132^2) is strictly between zero and one. It cannot be a
Coxeter contact: the smallest positive finite value is one, for pi/3.
At the remaining endpoint, where wall 132 is right, wall 31 instead has a
cyclic invariant strictly between zero and one. Positive norm and all strict
margins are independently verified.

Thus every compatible Coxeter wall in these four entire families has
strictly negative E4 pairing. There is no height, norm, or angle-label cutoff.
Nonorthogonal contacts to the other caps and walls outside this residual
space remain open.

```sh
python3 scripts/fresh_h20_h4_edge_escape_gate_20261016_review.py --verify
```
