# Different root norm classes cannot meet obliquely

In the discriminant-eleven lattice, a norm-two reflecting hyperplane and a norm-`2β` reflecting hyperplane are **orthogonal or strictly disjoint**, where `β=3+φ`. They cannot meet at any nonright angle or be parallel. This holds at all heights and does not require a classification of Coxeter angles.

Let the normals be `r,s`, with `r²=2`, `s²=2β`. Reflection integrality for the longer root gives `(r,s)=βz`, with `z∈O=Z[φ]`. If the product is zero, the walls are orthogonal. Otherwise put

`a=(r,s)²/β=βz²`.

This is a totally positive algebraic integer with `Norm(a)=11 Norm(z)²`. The two root lines are distinct because `β` is not a square in the field: its norm eleven is not a rational square. At the conjugate embedding the full form is positive definite, so its two-root Gram matrix gives `0<a′<4`.

If the hyperplanes intersect or are parallel at the hyperbolic embedding, `0<a≤4`. Therefore `0<Norm(a)<16`, forcing `Norm(a)=11`. The integer trace is positive and less than eight. Its quadratic discriminant is `Tr(a)²−44≥0`, so the trace must be seven. Consequently `a=(7±√5)/2`; one of its embeddings is greater than four, a contradiction.

The [certificate](fresh_h9_beta11_mixed_norm_gap_20261014.json) and [checker](../scripts/fresh_h9_beta11_mixed_norm_gap_20261014.py) verify the finite algebraic inputs to this argument, binding the independently reviewed norm alphabet and reflection filter. The universal conclusion follows from the displayed norm-and-trace deduction. It constrains mixed-class contacts in the full arithmetic search but does not establish reflectivity.

```sh
python scripts/fresh_h9_beta11_mixed_norm_gap_20261014.py --verify
```
