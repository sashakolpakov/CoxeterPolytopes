# Independent exclusion of the local-Witt H9 host

**PASS, generation and fresh saved replay.** The literal lattice
`<−φ> ⊥ I7 ⊥ (4−φ)I2` is nonreflective. The
[independent checker](../scripts/fresh_h9_witt_twist_symmetry_20261016_review.py)
and [certificate](fresh_h9_witt_twist_symmetry_20261016_review.json) reconstruct
two integral symmetries of its full reflection chamber. Their product has
trace `50+76φ>10`, and therefore infinite order.

This is a different quadratic space from every E8-plus-binary host, as proved
by the [independent local invariant review](fresh_h9_witt_twist_20261016_review.md).
Its exclusion does not contradict that distinction: it supplies a separate
obstruction in the new space. No new finite-volume example is produced.

The proof uses the independently reconstructed complete
[height20 prefix](fresh_h9_witt_twist_pilot_20261016_review.md), with all13
controls and16 walls. The all-height primitive reflection norms1,2,β,2β and
complete B7+B2 zero sector were checked independently. No higher-height
enumeration is needed for the obstruction.

The two literal source/target frame pairs, in zero-based root indices, are:

| Generator | Source frame | Ordered target frame |
|---|---|---|
| 14 | 0,1,2,4,6,7,8,9,11 | 4,3,2,0,6,7,8,11,9 |
| 3 | 0,1,2,3,4,5,6,8,10 | 10,14,2,3,4,5,6,8,0 |

For each of the four frames, the reviewer computes the apex directly by
projecting H=e₀ off the nine normals. It checks the positive-definite corner
Gram, inverse-Gram nonnegativity, a future timelike apex, all16 inequalities
and exactly the nine claimed incident walls. This yields
`H = V + Σ bᵢNᵢ` with every bᵢ≤0.

The two maps are reconstructed from their root-plus-apex frames in project
Golden arithmetic. The only imported implementation is the independently
written linear algebra from the norm131 reviewer; no producer or producer
arithmetic helper is imported. Both matrices and inverses are integral,
preserve the literal form, and carry each oriented source normal to its
specified target. Positive apex scaling proves preservation of the future
sheet. Their determinants are−1 and1, and their traces are4 and6. The product
is recomputed and compared with the saved discovery witness; its trace is
exactly50+76φ and its determinant is−1.

The full-chamber argument is the same directly checked geometric argument
used in the [norm131 review](fresh_h9_beta131_h300_20261016_review.md).
For any later simple root r, `(H,r)<0` and `(Nᵢ,r)≤0`, hence
`(V,r)=(H,r)−Σ bᵢ(Nᵢ,r)<0`. Neither finite corner can be cut or reached by a
later wall. Local finiteness of lattice mirrors preserves the whole corner
sector. An integral isometry maps the lattice mirror arrangement to itself;
agreement of its image sector with a sector of the same chamber therefore
makes it a full-chamber stabilizer.

Each of the two verified generators stabilizes the full chamber, so their
infinite-order product does too. The full reflection subgroup has infinite
index in the arithmetic lattice and infinite covolume. Every group generated
by reflections preserving this literal lattice lies in that full reflection
subgroup, so none can have finite covolume.

The producer's37-corner cone census helped discover the maps; its completeness
is not used by this independent proof. Likewise the discovery screen's44
integral maps and149 tested maps were only a bounded search, not a full
stabilizer classification. The certificate needs precisely the two displayed
maps and their four actual corners. Other lattices and arbitrary geometric
walls remain outside its scope.

```sh
python3 scripts/fresh_h9_witt_twist_symmetry_20261016_review.py
python3 scripts/fresh_h9_witt_twist_symmetry_20261016_review.py --verify
```
