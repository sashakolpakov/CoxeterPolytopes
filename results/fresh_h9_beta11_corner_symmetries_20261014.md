# An infinite-order chamber symmetry for the beta11 lattice

The lattice `L = E8 ⊥ [[−2φ,−1],[−1,2]]`, where `φ=(1+√5)/2`, is **nonreflective**. An exact integral future-preserving isometry stabilizes its full integral-reflection chamber and has trace `4+8φ>10`. Consequently the chamber has infinite volume. No finite-volume Coxeter polytope can have all its wall reflections preserving this literal lattice.

The [certificate](fresh_h9_beta11_corner_symmetries_20261014.json) is generated and replayed by the [exact checker](../scripts/fresh_h9_beta11_corner_symmetries_20261014.py). It uses the frozen standalone Fraction arithmetic of the [independent complete height300 replay](fresh_h9_beta11_vinberg_h300_20261014_review.md). The two complete primitive norm classes and the zero sector are inherited from the [norm-alphabet proof](fresh_h9_h8_cap_height_20261014_norm_pool_review.md) and [zero-root certificate](fresh_h9_beta11_zero_sector_20261014.md). No numerical equality or matrix-integrality test is used in acceptance.

## The literal witness

Use the basis `(x,y1,…,y8,t)` of the height300 prefix. Its initial finite corner is `H=(2,0,…,0,1)`, incident to outward normals `0,…,8`. The target is finite vertex735 of the [35-wall geometry](fresh_h9_beta11_geometry_h300_20261014.md), with incident set

`[2,4,5,6,7,11,20,29,30]`.

The exact map sends the ordered initial normals to

`[7,4,6,5,20,29,30,11,2]`

and sends `H` to `(8+14φ)P`, where `P` is the saved future timelike target point. The positive scale is reconstructed from the determinant ratio of the two root-plus-apex frames. The JSON records the complete 10×10 matrix, its integral inverse, both frames and all checks. In particular:

- `TᵀGT=G`, `det(T)=1`, and every entry of `T` and `T⁻¹` lies in `Z[φ]`;
- the nine outward-root identities hold literally, without sign reversal or rescaling;
- both corner points are future timelike, vanish on exactly their nine incident roots, and are strict at every other prefix wall;
- the common spherical Gram is positive definite, with entrywise nonnegative inverse;
- `tr(T)=4+8φ>10`.

This is one verified witness. The exploratory corner-map search is not claimed to classify all symmetries.

## Why a finite prefix corner is already a full-chamber corner

Let `v` be a finite corner of a genuine acute Vinberg prefix, with incident positive-definite Gram `G_S`. The matrix has nonpositive off-diagonal entries, hence `G_S⁻¹` is entrywise nonnegative. Decompose the future Vinberg basepoint as

`H = a v + Σ b_i n_i`, with `a=(H,v)/(v,v)>0` and `b=G_S⁻¹((H,n_i))≤0`.

Every later simple wall `r` is oriented by `(H,r)<0` and satisfies `(n_i,r)≤0`. Therefore

`(H,r) = a(v,r) + Σ b_i(n_i,r) ≥ a(v,r)`,

so `(v,r)<0`. No later wall cuts or touches `v`. The checker verifies this decomposition physically at the target corner. The initial corner is also covered by the exact zero-sector census.

There is no missing distance threshold here. The complete two-norm alphabet and positive-definite conjugate form make the integral mirror arrangement locally finite: bounded physical mirror distance bounds the physical root coordinates, while the fixed conjugate norm bounds the conjugate coordinates. The Minkowski embedding then leaves only finitely many lattice vectors. Thus the nine incident walls describe an actual neighborhood of each corner in the full chamber.

The height300 prefix used here is the complete Vinberg prefix for **both** primitive norm classes, independently replayed. The older norm-two-only host search is not used to infer this premise.

## From the local map to infinite volume

The integral inverse makes `T` an automorphism of `L`; it therefore permutes all integral reflection mirrors. Its future orientation and exact outward-root map carry a small interior sector at `H` onto the corresponding sector at `P`. Two chambers of the same locally finite mirror arrangement cannot share an interior open set. Hence `T` stabilizes the entire reflection chamber, including walls beyond height300.

A finite-order real 10×10 matrix has all eigenvalues on the unit circle, so its trace has absolute value at most10. The verified trace proves that `T` has infinite order.

The future integral isometry group is discrete. Indeed its other real embedding preserves a positive-definite form, so all conjugate matrices lie in a compact set; compact physical matrix bounds together with these conjugate bounds leave finitely many integral matrices. An infinite discrete chamber stabilizer supplies infinitely many disjoint translates of a sufficiently small interior ball, all inside that chamber. Its volume is therefore infinite, and the integral reflection subgroup has infinite index in the full lattice isometry group.

The conclusion applies to every Coxeter polytope whose wall reflections preserve this `L`, not just to the displayed prefix: its generated reflection group would be a subgroup of the full integral reflection group. A finite-covolume subgroup would force that larger discrete group to have finite covolume, contradicting the result above.

## Validation and scope

Generation and saved `--verify` replay pass. The certificate checks 100 isometry entries, 200 matrix-and-inverse integrality entries, nine literal root images, 70 corner-wall signs, the exact positive scale and the preservation decomposition. It also verifies positive definiteness of the conjugate host form.

This closes the reflection-lattice route for this host, including its earlier equivalent signed/permuted presentation. It is not a nonexistence result for arbitrary compact hyperbolic9 Coxeter polytopes or nonarithmetic hyperbolic20 constructions. A larger geometric construction using reflections that do not preserve this lattice is outside the theorem.
