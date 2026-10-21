# Independent proof that the discriminant-55 host is nonreflective

**The full integral reflection group of `E8 ⊥ [[−2φ,−1],[−1,4]]` has infinite covolume.** No finite-volume Coxeter polytope can have all wall reflections preserving this literal lattice. This conclusion is specific to the lattice; arbitrary geometric walls and other lattices remain outside it.

The independent [checker](../scripts/fresh_h9_beta55_corner_symmetries_20261014_review.py) and [certificate](fresh_h9_beta55_corner_symmetries_20261014_review.json) reconstruct the witness using separate `Golden` arithmetic and freshly ordered frames. They import neither the symmetry producer nor its arithmetic helper. The [independent complete height-100 enumeration](fresh_h9_beta55_vinberg_h100_20261014_review.md) supplies the complete root-prefix dependency, including all four norm classes.

## The exact map

Finite corner 3 has incident walls `[0,1,2,3,4,5,6,9,10]`. The map sends these ordered outward normals to `[25,0,35,2,3,4,5,24,7]`, the complete incident set at finite corner 42. Its action on the vertex ray has positive scale `2+2φ`. Reconstructing the ten-dimensional matrix from these vectors gives:

- integral matrix and integral inverse;
- `TᵀGT=G`, determinant `−1`, and preservation of the future timelike sheet;
- trace `5+8φ > 10`;
- physical inertia `(9,1)` and positive-definite conjugate form.

All entries agree with the producer's matrix. The two incident Grams are positive definite, and their inverses have nonnegative entries. Every coordinate of both basepoint decompositions is checked independently.

## Why these are full chamber corners

Let B be the Vinberg basepoint, v either finite corner, and n_i its nine outward incident roots. The checker verifies

`B = a v + Σ b_i n_i`, with `a>0` and every `b_i≤0`.

For any later simple root r, completeness of the zero sector gives `(B,r)<0`, while acuteness gives `(n_i,r)≤0`. Therefore

`a(v,r) = (B,r) − Σ b_i(n_i,r) < 0`.

No later wall cuts or touches either corner. Local finiteness of the integral reflection arrangement then makes the checked outward sectors the complete local sectors of the full chamber. This argument needs no additional height cutoff beyond the verified prefix containing the two corners.

The integral future isometry T preserves the complete reflection arrangement. Since it maps one complete outward corner sector to the other, it stabilizes their full chamber. A finite-order ten-dimensional real matrix has trace at most ten in absolute value, so the checked trace proves that T has infinite order. The admissible integral isometry group is discrete; a sufficiently small ball about a generic chamber-interior point has infinitely many disjoint equal-volume translates under powers of T. Thus the full chamber has infinite volume.

Any finite-volume Coxeter polytope whose wall reflections preserve this lattice would contain a chamber of its full integral reflection arrangement. That chamber already has infinite volume, which is impossible. This also excludes finite-index geometric variants whose wall reflections preserve the same literal lattice; it does not transfer the theorem to merely similar or commensurable lattices.

Generation and saved independent replay pass. No classification of all chamber symmetries or all arithmetic H9 hosts is claimed.
