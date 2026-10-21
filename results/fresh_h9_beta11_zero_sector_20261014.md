# Complete zero-height sector for the discriminant-eleven H9 host

The lattice `L=E8 ⊥ [[−2φ,−1],[-1,2]]` has a completely determined starting sector for Vinberg's algorithm. In coordinates `(x,y,t)`, choose `H=(2,0,1)` and write `δ=1+4φ=φ(3+φ)`. Then

`H²=−2δ`, `−(r,H)=δx`, and `height(r)=δx²/(2q)`.

Every height-zero root has `x=0`. Its remaining form is the ordinary rational even positive-definite lattice `M=E8⊥A1`, extended to `Z[φ]`. The independently established [norm alphabet](fresh_h9_h8_cap_height_20261014_norm_pool_review.md) leaves norms 2 and `2β`, where `β=3+φ`.

For norm two, write a vector as `u+φv`, with `u,v` in the rational integral lattice `M`. Comparing coefficients gives `u²+v²=2` and `2(u,v)+v²=0`. Since `M` is even and positive definite, the first equation permits only `u²=2,v=0` or `u=0,v²=2`. The latter violates the second equation. Thus every such root is an ordinary integral norm-two root of `M`.

The [checker](../scripts/fresh_h9_beta11_zero_sector_20261014.py) enumerates that entire integer ellipsoid using exact rational square completion. It finds 242 roots: 240 from E8 and the two A1 roots. The recursion visits 770 nodes. Every root is a combination of the proposed nine simple normals with coefficients all of one sign, so none of these mirrors cuts the chosen sector.

For norm `2β`, the reflection-integrality filter at `x=0` forces every remaining coordinate to be divisible by `β`. Dividing the vector by `β` would give an integral vector of norm `2/β`, impossible in an even lattice because `1/β=(4−φ)/11` is not integral. There are no zero-height roots in this class.

Consequently the complete zero-height sector is **E8+A1**, with simple normals the eight negative E8 coordinate roots and the final binary basis vector. The [certificate](fresh_h9_beta11_zero_sector_20261014.json) stores the basepoint, all 242 roots and the simple normals. Generation and saved replay pass. This starts the positive-height enumeration without assuming that the earlier fifteen-wall chamber is a complete prefix.

```sh
python scripts/fresh_h9_beta11_zero_sector_20261014.py --verify
```
