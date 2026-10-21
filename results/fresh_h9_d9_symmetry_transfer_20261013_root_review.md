# Root source review and fresh replay of D9 nonreflectivity

PASS. The [global checker](../scripts/fresh_h9_d9_symmetry_transfer_20261013.py)
was read and freshly run with `--verify` after the construction agent's
certificate production. The [saved certificate](fresh_h9_d9_symmetry_transfer_20261013.json)
reproduces exactly, with status
`PASS_EXACT_INFINITE_FULL_D9_CHAMBER_STABILIZER`.

The proof uses four complete first-edge calculations: three steps along the
original chain, and the new initial A2+B6 edge. Their separate arithmetic
audits and recursive source hashes are explicit dependencies. They connect
two full B3+B6 corners of the same reflection chamber. The checker verifies
the joins: eight identical oriented edge normals, one released and one
arriving normal, the appropriate strict signs, and timelike endpoints.

The transport is constructed by mapping nine incident primitive normals and
one scaled timelike corner vector. I checked that the chosen permutation
preserves the full Gram, that the endpoint scale is positive, and that the
checker verifies exact form preservation and the inverse. Conjugation into
the actual time-plus-D9 basis makes both matrices integral. This checks the
sum-in-2O lattice condition, not merely coordinate integrality. The generated
finite root closures contain 45 mirrors at each corner; completeness of these
corner sectors rests on the complete first-contact censuses, not on the finite
closure count alone.

The oriented simple-system mapping takes the local sector of one corner to
the other. The certified paths identify these sectors as parts of the same
full chamber. Consequently the lattice isometry stabilizes that chamber.
Its exact trace is

```text
-7378736 - 3993344*T + 12596300*T^2 + 6817068*T^3 > 10.
```

A finite-order real ten-dimensional matrix has trace at most ten, so this
isometry has infinite order. The arithmetic lattice isometry group is
discrete; an infinite chamber stabilizer forces infinite chamber volume.
This excludes finite-covolume reflection subgroups using this lattice's
mirrors, since their chambers contain full reflection chambers. It does not
exclude other lattices, other quadratic spaces, or H9 examples in general.

The checker reuses the existing independent Fraction field engine. The
number-field unit and norm-orbit censuses retain their certified PARI
dependency. This source review and replay do not claim a new implementation
of those number-field algorithms. The supplementary overclosure calculation
and failed eight-element transport group are unnecessary to the successful
global proof.

```sh
python scripts/fresh_h9_d9_symmetry_transfer_20261013.py --verify
```
