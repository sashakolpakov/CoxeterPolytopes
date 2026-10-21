# Complete two-norm Vinberg prefix through height 100

The exact enumeration through Vinberg height **100** gives **21 actual walls**: the nine height-zero walls and twelve positive-height walls. It exhausts both primitive norm classes, `2` and `2β`, where `β=3+φ`. There are 46 candidate levels and 68,242 recursion nodes. All 210 wall pairs pass the Coxeter gate, and explicit strict facet points certify all 21 walls.

The [producer](../scripts/fresh_h9_beta11_vinberg_20261014.py) and [transcript](fresh_h9_beta11_vinberg_20261014.json) do not presume that the earlier 15-wall chamber, or the old norm-two-only Vinberg prefix, is a full-group prefix. This is a bounded enumeration; termination, finite volume and reflectivity are not asserted here.

## Host identification and zero sector

In the [certified free basis](fresh_h9_h8_cap_height_20261014_host_basis.md), the form is

`L = E8 ⊥ [[−2φ,−1],[−1,2]]`, with `δ=1+4φ=φβ`.

This is the **same discriminant-eleven host** previously named `H11` in the repository, despite its hyperbolic dimension being nine. The transcript certifies the literal signed-permutation identity between the two Gram matrices. The new E8 nodes map to old nodes `[0,7,1,2,3,4,5,6]`, and the old second binary coordinate is minus the new last coordinate. All enumeration matrices are reconstructed from the literal new host; no E8 ordering is silently assumed.

The earlier [H11 search](h11_vinberg_prefix.md) treated norm-two roots only. The new [complete norm alphabet](fresh_h9_h8_cap_height_20261014_norm_pool.md) requires the additional norm-`2β` class. Its first wall occurs before the former first short wall and rejects that wall. Thus the old accepted list cannot be reused as a full-group prefix without a new selection.

Take `H=(2,0⁸,1)`, so `H²=−2δ` and `−H·r=δx` for `r=(x,y,t)`. The [independent zero-sector certificate](fresh_h9_beta11_zero_sector_20261014.md) proves that the complete height-zero arrangement consists of the 240 ordinary E8 roots and the two A1 roots. The chosen simple walls are the eight negative E8 basis roots and the last binary basis root. Their sector is exactly `E8+A1`.

## Exhaustive finite enumeration

Let `C` be the literal E8 Cartan matrix and put

`λ=C y`, `l=x−2t`, `q=r²`.

The zero-sector inequalities are `λ_i≥0` and `l≥0`. Completing the square gives

`2λᵀ C⁻¹ λ+l² = 2q+δx²`.

The geometric height is `h=δx²/(2q)`, the squared hyperbolic sine of the wall's distance from `H`. The bound `h≤100` controls `x` at the hyperbolic embedding. Positivity at the conjugate embedding gives `2q′+δ′x′²≥0`, and hence controls `x′`. Exact coefficient boxes in the basis `1,φ` exhaust all such positive `x`: ten values for norm two and 36 for norm `2β`. The two classes are merged by **exact height**, not by control coordinate. All 46 levels are distinct.

For each level, the positive form is `diag(1,2C⁻¹)` in the coordinates `(l,λ)`. Its inverse has every diagonal entry equal to one. Consequently every coordinate's square is bounded by `2q+δx²` at both embeddings, giving a complete finite ring alphabet. Exact rational prefix Schur complements give lower bounds for every partial assignment.

The integrality condition `t=(x−l)/2∈O` is enforced coefficientwise. For norm `2β`, the reflection residue filter becomes the particularly simple condition that **every** coordinate `l,λ_1,…,λ_8` is divisible by `β`. The reconstructed vector's norm and its reflection integrality are checked again directly.

For an earlier wall `R=(X,Y,T)` with `L=X−2T`, acuteness is

`2λ·Y+lL ≤ δxX`.

All coefficients on the left are nonnegative in the chosen sector. A partial sum exceeding the right side can therefore be rejected safely. Schur bounds are tested at both embeddings; these wall bounds are tested at the hyperbolic embedding. Every comparison in the search uses exact integer-pair arithmetic. No numerical tolerance or floating-point pruning is used.

The alphabet and recursion helpers are reused from the old exact implementation, after checking their hypotheses on the literal new form. The new norm-class merge, residue filter and accepted wall list are reconstructed here. No old positive-height root list is imported.

## Nonempty levels and actual walls

| Norm | Control coordinate `x` | Height | New walls |
| --- | --- | --- | ---: |
| `2β` | `φ−1` | `(φ−1)/4` | 1 |
| `2` | `1` | `1/4+φ` | 1 |
| `2` | `φ` | `5/4+9φ/4` | 1 |
| `2` | `2φ` | `5+9φ` | 1 |
| `2` | `1+2φ` | `37/4+15φ` | 2 |
| `2` | `2+2φ` | `14+23φ` | 3 |
| `2` | `2+3φ` | `97/4+157φ/4` | 3 |

The remaining 39 candidate levels are exhausted and empty after earlier-wall tests. There is one long wall and twenty short walls in the resulting prefix. Its pair counts are 137 right, 32 at `π/3`, fifteen at `π/5`, and 26 dotted; no parallel pair occurs.

The transcript constructs a strict future center by moving `H` a dyadic distance `1/8` into its zero sector. For each normal `N_i`, the point

`P_i = center − (center·N_i/N_i²)N_i`

is checked to be future timelike, zero at exactly that wall, and strictly negative at every other wall. This verifies all 21 actual facets. The full Gram, Cartan matrix, pair classifications, center and facet points are stored in the transcript.

## Verification and limits

Generation and saved replay pass. Source hashes bind the exact host, complete two-norm proof and independent review, complete zero sector, arithmetic helpers and producer. The transcript records every candidate control value, coordinate alphabet, per-depth pruning counter and accepted root.

```sh
python scripts/fresh_h9_beta11_vinberg_20261014.py --height 100
python scripts/fresh_h9_beta11_vinberg_20261014.py --verify
```

The saved replay reads the bound from the saved transcript when `--height` is omitted. No larger height is included implicitly. This proves completeness through the stated height for the integral reflection group of this lattice; it does not restrict arbitrary nonintegral geometric caps. The global corner, edge and volume diagnosis is a separate task. No frozen predecessor was changed.
