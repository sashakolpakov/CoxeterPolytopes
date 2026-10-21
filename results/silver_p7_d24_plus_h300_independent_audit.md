# Independent audit: silver prime-seven `D24+` H=300 ray

## Verdict

The direct replacement-ray certificate in the schema-0.3 H=300 archive is
valid.  The complete 348-wall prefix has an exact spacelike extreme ray and
therefore infinite volume.  This conclusion is only about the H=300 prefix:
later Vinberg walls may cut the ray, and the eventual chamber remains
unresolved.

No precision or scope defect was found in the backend, generator, frozen
artifact, tests, or primary report.

## Clean-room reconstruction

`scripts/independent_audit_silver_p7_d24_plus_h300.py` imports only the Python
standard library.  It implements `Q(sqrt(2))` arithmetic with `Fraction`,
constructs the 24 negative `D24` simple walls in orthonormal coordinates, and
then reads the 324 accepted positive-height walls from the frozen archive.  It
checks, independently of the search backend:

- all 348 wall norms and the 324 archived height/target identities;
- `D24+` membership, the reduced glue-coset labels, and reflection
  integrality against an explicit integral `D24+` basis;
- chronological wall orientation, including same-height walls;
- the Lorentzian coordinate convention
  `-(2sqrt(2)-1)x_0y_0 + sum(x_i y_i)` and the future orientation `x_0>0`.

For prime-seven norms, the archive's `coset` field labels the reduced spatial
vector before multiplication by `delta=3+sqrt(2)`.  Multiplication by `delta`
permutes the three nonzero spinor residues modulo 2.  The audit accounts for
this and separately verifies that every displayed orthonormal vector itself
lies in `D24+`.

## Replacement ray

For

```text
c = (172-120sqrt(2))/7,
v = (1; c^19, 0^5),
```

clean-room exact pairing gives 324 negative, 24 zero, and no positive
pairings against the reconstructed 348 walls.  The zero set is

```text
0,...,17, 19,...,23, 59.
```

Thus 23 active facets are height-zero `D24` walls.  Facet 59 is the
zero-based fourth accepted wall in the `62+44sqrt(2)` packet.  Exact Gaussian
elimination over `Q(sqrt(2))`, with no modular specialization, gives rank 24
for these active wall functionals.  Hence their common nullspace is the line
spanned by `v`, and feasibility makes it an extreme ray.

Its exact Lorentzian norm is

```text
(1109345-784418sqrt(2))/49 > 0.
```

The field norm of this element is `67439/7 > 0`, while its rational
coefficient is positive, which independently certifies the displayed sign.

## Reproducibility

The clean-room command is

```text
python3 scripts/independent_audit_silver_p7_d24_plus_h300.py
```

Regenerating schema 0.3 with the project generator produced a byte-identical
archive.  Both files had SHA-256

```text
4ea87341148489598dc85b3cc8271965a036e5ffba2f24052784705e4af614cf
```

The primary report correctly limits its conclusion to infinite volume of the
H=300 prefix and explicitly leaves the eventual chamber and reflectivity
question unresolved.
