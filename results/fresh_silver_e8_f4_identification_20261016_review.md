# Independent check: the E8 host duplicates the F4² host

The two full Lorentz lattices are integrally isometric over Z[√2]. The
[separate verifier](../scripts/fresh_silver_e8_f4_identification_20261016_review.py)
and [certificate](fresh_silver_e8_f4_identification_20261016_review.json)
pass generation and saved replay using the older standalone Fraction engine.

The explicit matrix in the [primary identification](fresh_silver_e8_f4_identification_20261016.md)
has integral entries and integral inverse, determinant one, and satisfies
`Pᵀ G_E8 P = G_F4²`. It preserves the future cone. Both Gram matrices are
independently reconstructed from their literal block definitions. These
identities suffice to prove an integral lattice isomorphism; the finite
search that discovered the matrix is not needed for verification.

The change of basis mixes the time coordinate with the spatial coordinates.
Different positive blocks therefore did not make these distinct Lorentz
hosts. Their full integral reflection groups are conjugate. The new E8
symmetry proof remains valid, but excludes no new integral isometry class
beyond the previously excluded F4² host. Future candidate hosts should be
screened for such identifications before extending long root searches.

```sh
python3 scripts/fresh_silver_e8_f4_identification_20261016_review.py
python3 scripts/fresh_silver_e8_f4_identification_20261016_review.py --verify
```
