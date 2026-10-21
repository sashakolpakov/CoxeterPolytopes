# Independent review of the new tripled-D9 chamber symmetry

**PASS.** Separate Fraction polynomial arithmetic independently reconstructs the transport in the new host
\(\langle-3(T-1)\rangle+D_9\), with \(T^4-4T^2+2=0\) and \(T=2\cos(\pi/8)\). No matrix or nonreflectivity conclusion is transferred from the preceding host.

The [independent checker](../scripts/fresh_h9_tripled_d9_20261013_transport_independent.py) reconstructs the matrix from the two new complete oriented B3+B6 corners. It checks the form, inverse, positive corner scale, future orientation, all nine primitive reflective simple-normal maps, and integrality of both maps on the actual D9 basis. Each corner generates exactly 45 finite mirrors. The exact trace exceeds 10, so the matrix has infinite order.

It also reconnects all four required paths from the original endpoint artifacts. Each step retains eight oriented normals, has the correct arrival/release signs and future timelike endpoints, and lies in the same retained prefix. Completeness and absence of an intervening mirror depend explicitly on the separately replayed initial and chain first-edge censuses: 136 profiles, 1,161 interval points, and all 26 safe norm classes. Their nested source manifests and the certified PARI norm/unit interface are hash-bound dependencies; this small matrix audit does not independently replace those arithmetic enumerations.

These premises are enough for the global conclusion. A lattice isometry permutes the full reflection arrangement. Transporting the full oriented corner sector between points connected in one chamber therefore stabilizes that chamber, rather than merely matching a finite collection of walls. A finite-volume hyperbolic polytope has finite isometry stabilizer; the infinite-order transport excludes a finite-volume full reflection chamber in this specific host. The trace test is elementary: a finite-order real 10-dimensional matrix has trace at most 10.

The implementation uses the preceding **independent Fraction field engine**, with the form replaced by the new coefficient. It imports no producer field-class arithmetic. The frozen [replay result](fresh_h9_tripled_d9_20261013_transport_independent.json) records both D9 module matrices, the four path checks and a rejected forged matrix.

```sh
python scripts/fresh_h9_tripled_d9_20261013_transport_independent.py --verify
```

The conclusion is nonreflectivity of this explicit lattice. It is not nonexistence of compact H9 polytopes, and does not rely on commensurability or identification with an older F4 host.
