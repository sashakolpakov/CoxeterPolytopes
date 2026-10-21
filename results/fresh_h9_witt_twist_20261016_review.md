# Independent review of the new local-Witt H9 setup

**PASS.** The [independent checker](../scripts/fresh_h9_witt_twist_20261016_review.py)
and [certificate](fresh_h9_witt_twist_20261016_review.json) confirm the exact
setup for `L = <−φ> ⊥ I7 ⊥ β I2`, with β=4−φ. Generation and fresh saved
replay pass. They use project Golden arithmetic and no primary-code imports.
The [primary proof](fresh_h9_witt_twist_20261016.md) gives the complete local
and integral arguments; they were reviewed separately from the calculations.

The physical signature is (9,1), the conjugate is definite, and the determinant
is −φβ² of norm−121. Since βG⁻¹ is integral, primitive reflection integrality
implies q divides2β. The PID and unit-square facts for Z[φ] then leave exactly
the four totally positive norm representatives1,2,β,2β. The reviewer checks
the dual multiplier explicitly, the inertness of2, and N(β)=11. There is no
height truncation in this norm-alphabet argument.

At (β), the residue of φ is4. Hensel lifting makes φ a local square. The
residue binary form I2 is anisotropic; all121 residue pairs are checked.
The nonsingular isotropic vector (1,1,3,0) for I4 lifts, and its square
discriminant makes the remaining binary complement hyperbolic. Thus I4 is
split. The form `I2 ⊥ β I2` is anisotropic by two successive reductions
moduloβ: a hypothetical isotropic vector normalized with a unit coordinate
would have every coordinate divisible byβ. Hence L has local decomposition
`H³ ⊥ (I2 ⊥ β I2)` and Witt index3.

An explicit rational basis for E8 is verified directly against its literal
Cartan matrix. It identifies its rational form with I8, locally H⁴. Every
E8-plus-binary form, including every scalar similarity, therefore has local
Witt index at least4. The distinction covers the entire family, not just the
three tested discriminants. It also rules out a rational spherical E8
component anywhere in this host. This local argument does not assert
reflectivity or nonreflectivity.

The maximal-integral claim uses the same anisotropic residue plane. Outside
(β) the lattice is self-dual; at (β), an integral vector in a larger lattice
would give a nonzero isotropic vector in the residue I2. There is none.
Here integrality is for the bilinear Gram form; this is not a claim about
every possible commensurable reflection lattice.

The zero-sector proof excludes mixed blocks using
`11 N(z)² ≤ 4`, impossible for a nonzero algebraic integer. The reviewer then
enumerates all424 sparse coefficient vectors allowed by the identity
`sum(aᵢ²+bᵢ²) ≤ 2`, rather than importing the primary root list or its dynamic
counts. It recovers exactly106 normalized oriented roots:14,84,4,4 in the
four norm classes. All reflection-integrality tests, nine simple walls and36
pair relations agree. The complete zero sector is B7+B2, with53 mirrors.

The subsequent [independent positive-height pilot](fresh_h9_witt_twist_pilot_20261016_review.md)
is a separate calculation. These setup checks do not substitute for a
complete root search or a finite-volume test.

```sh
python3 scripts/fresh_h9_witt_twist_20261016_review.py
python3 scripts/fresh_h9_witt_twist_20261016_review.py --verify
```
