# An infinite-order symmetry closes the discriminant55 integral host

The lattice `L55=E8 ⊥ [[−2φ,−1],[−1,4]]`, where `φ=(1+√5)/2`, is **nonreflective**. Its full integral reflection chamber has infinite volume. An exact integral future-preserving isometry maps one full finite corner to another and has trace `5+8φ>10`, hence infinite order.

Consequently, no finite-volume Coxeter polytope can have all its wall reflections preserving this literal lattice. Other lattices, including different lattices in a related quadratic space, and arbitrary nonintegral geometric walls are outside this result.

The [producer](../scripts/fresh_h9_beta55_corner_symmetries_20261014.py) and [exact certificate](fresh_h9_beta55_corner_symmetries_20261014.json) reconstruct one map from literal root images and corner vectors. They do not claim an exhaustive symmetry search. Generation and saved replay pass. The arithmetic is the frozen standalone Fraction implementation in the basis `1,φ`; the map computation does not use floating-point acceptance.

The [independent review](fresh_h9_beta55_corner_symmetries_20261014_review.md) also passes generation and saved replay. It uses separate Golden arithmetic and vertex-first frames to reconstruct the matrix and inverse, both corner-preservation decompositions, signatures and trace. It independently binds the complete height100 root replay.

## The two literal corners and the map

The source is vertex3 of the [height100 geometry](fresh_h9_beta55_geometry_h100_20261014.json), with outward incident wall indices

`[0,1,2,3,4,5,6,9,10]`.

Their ordered images are

`[25,0,35,2,3,4,5,24,7]`,

which form the complete incident set of vertex42. Both points satisfy every one of the44 prefix inequalities, with equality exactly at their nine incident walls. Their norms are respectively `7−9φ/2` and `97/8−15φ/2`, both negative, and both are future-directed relative to the Vinberg basepoint `B=(4,0,…,0,1)`.

The two root-plus-apex frames determine a unique Lorentz map after scaling the target apex by the positive factor `2+2φ`. The certificate records the full matrix and its inverse, checks all200 entries lie in `Z[φ]`, and verifies

`TᵀGT=G`, `det(T)=−1`, `T(v3)=(2+2φ)v42`.

All nine outward simple roots map exactly to the listed outward roots, without a sign reversal. The determinant−1 is allowed: the positive apex relation proves preservation of the future sheet. The trace is exactly `5+8φ`. A finite-order real10-dimensional matrix has all eigenvalues on the unit circle and therefore trace at most10; this map has infinite order.

## Why these prefix corners are full chamber corners

The [complete height100 enumeration](fresh_h9_beta55_vinberg_h100_20261014.md) uses all four possible primitive reflection norm classes. Its walls are actual walls of the full integral reflection chamber. The [height20 proof](fresh_h9_beta55_pilot_20261014.md) supplies the complete norm alphabet and zero sector; the zero sector is E8 only. In particular B is not a finite corner, and the symmetry is not asserted to map B to a corner.

For either checked finite corner v, let S be its nine outward normals and J the inverse incident Gram. The certificate checks positive definiteness, nonpositive off-diagonal entries, and entrywise nonnegativity of J. It independently verifies the physical decomposition

`B=a v+Σ b_i n_i`, with `a=(B,v)/v²>0` and `b=J((B,n_i))≤0`.

Every later simple wall r has `(B,r)<0` and `(n_i,r)≤0`. Therefore

`(B,r)=a(v,r)+Σ b_i(n_i,r) ≥ a(v,r)`,

so `(v,r)<0`. No later wall can cut or touch either corner, regardless of its height. The exact all44 zero sets thus remain the complete outward incident sets of the full chamber.

The full integral mirror arrangement is locally finite. After normalization there are only four root norms; the conjugate form is positive definite, and roots whose hyperplanes meet a compact physical set are bounded at both embeddings. Integrality then gives finitely many such roots. A sufficiently small neighborhood of either finite corner consequently has precisely the local chamber sector determined by its incident roots. Since T preserves the lattice, form and future sheet and maps those outward roots positively, it maps one local sector into the other. Both belong to the same full reflection chamber, so T stabilizes that chamber.

## Infinite covolume and the host-wide consequence

The conjugate form is positive definite, with explicit positive elimination pivots in the certificate. Integral isometries are discrete in the physical factor: a bounded physical sequence has bounded conjugate matrices in a compact orthogonal group, and the two embeddings of their integral entries form a discrete lattice.

The infinite cyclic chamber stabilizer therefore supplies infinitely many disjoint translates of a sufficiently small ball around a generic strict interior point. Every translate remains inside the chamber and has the same positive volume. The chamber and the full integral reflection group have infinite covolume.

If a finite-volume Coxeter polytope had all reflections preserving L55, its reflection group would be a subgroup of the full integral reflection group. A discrete supergroup of a finite-covolume group also has finite covolume, contradicting the preceding result. This is a statement about all integral reflection walls of L55, not merely the44 enumerated walls or the finite height100 prefix.

The certificate depends on the complete bounded Vinberg enumeration, its norm/zero-sector proof, and the literal two corner witnesses. A complete classification of all129 prefix corners or all possible corner maps is not needed for the obstruction. No higher-height continuation was run after this witness was found.
