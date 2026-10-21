# Golden balanced cusp: finite-angle radial chains

A new cap can cross the golden branch’s radial gap while remaining connected to the central cap by prescribed finite angles. The exact construction retains the 22 seed facets and all 392 endpoint caps of radius `sqrt(28)`. A two-cap chain has successive angles `π/3`; its outer cap cuts the probe `7z/20`. Taking the full outer layer and eight compatible inner bridges gives **815 actual facets and cuts all 392 transported probes**.

This is not a finite-volume construction. A different strict spacelike escape survives, the spherical-cap contact graph remains disconnected, and the new caps create explicitly deficient cusp links. The result identifies a way past the single central-contact radius barrier, while making the remaining global problems concrete.

The [exact checker](../scripts/fresh_h20_balanced_golden_20261014.py) and [certificate](fresh_h20_balanced_golden_20261014.json) use `Fraction` arithmetic in `Q(φ)(sqrt(12N−12))`, with exact recursive signs. The source [balanced geometry](../scripts/fresh_h20_balanced_20261014.py) is hash-bound and its factor metric formulas are reconstructed. No floating-point decision appears in the certificate.

## The base and the coupled recipe

The cusp link is the product of unit-inradius `A6`, unit-inradius `A6`, and inradius-`φ` `A7` alcoves, where `φ²=φ+1`. Write `z_abc` for their 392 product vertices. Every vertex has

`N=|z_abc|²=112+84φ²=196+84φ`.

Their distance squares are reconstructed exactly as

`14d_A(7−d_A)+14d_B(7−d_B)+16φ²d_C(8−d_C)`.

The 22 vertical facets, the central sphere of radius two, and all 392 spheres centered at `z_abc` of radius `sqrt(28)` are coherent. The central seed angles are `π/3` in both `A6` factors and `π/5` in `A7`. Nearest endpoint pairs in either `A6` direction have angle `π/3`; all other endpoint pairs are dotted. Each endpoint has a positive-definite `A6+A6+A7+A1` incident Gram, and is an actual finite vertex.

Define

`k=(N+2−sqrt(12N−12))/(N−4)`,

so `0<k<1` and

`(N−4)k²−2(N+2)k+(N−4)=0`.

On one vertex ray, let `C_j` have center `(1−k^j)z` and radius `2k^j`; `C_0` is the central sphere. Each `C_j`, `j≥1`, has twelve retained seed contacts at `π/3`, seven at `π/5`, and three dotted omitted-seed relations. For two such spheres on the same ray, their normalized negative pairing is

`[(N−4)(k^(i−j)+k^(j−i))−2N]/8`.

The defining quadratic makes adjacent levels exactly `π/3`. Every nonadjacent pair used here is checked strictly dotted. Thus the finite angles are algebraic identities, not matches to numerical tolerances.

The first two levels have approximately

| cap | center multiple | radius |
|---|---:|---:|
| `C_1` | 0.173873485 | 1.652253031 |
| `C_2` | 0.317514981 | 1.364970039 |

The outer cap has power about `−1.512881266` at `7z/20`, so it cuts that old gap. Both caps retain every radius-`sqrt(28)` endpoint sphere. Their finite contacts are actual codimension-two faces: the relevant pair Gram is positive definite, and projecting the strict common center to that pair gives a timelike point strict on every other wall. This does not assert a new full rank-20 vertex.

## Three exact checkpoints

The eight selected rays use omissions

`(a mod 7, floor(a/7), a)`, for `a=0,…,7`.

Their least vertex distance square is `84+112φ²`. All caps on different selected chains are dotted, as verified using the full distance formula.

| checkpoint | actual facets | transported `7z/20` probes cut | explicit finite-contact paths |
|---|---:|---:|---|
| Two new caps, endpoint radius squared 28 | 417 | 1 | central → `C_1` → `C_2` |
| Eight chains through level eight, changed endpoint radius | 479 | 8 | eight paths from central to old endpoint facets |
| All 392 level-two caps plus eight level-one bridges, endpoint radius squared 28 | 815 | 392 | eight central → inner → outer paths |

The middle checkpoint deliberately replaces the endpoint radius by

`R_e²=(N−4)k^16 ≈15.434980142`.

This is below 21, so its complete endpoint layer is mutually dotted. On each selected ray, `C_8` is exactly right to that ray’s endpoint sphere. Earlier levels and every other endpoint pair are dotted. All 392 original finite endpoints remain, with their new heights. This is a coherent wall replacement, not an addition to the radius-`sqrt(28)` layer.

In the largest checkpoint, all 392 copies of `C_2` are mutually dotted; the minimum normalized negative pairing is approximately 1.272644437. Each selected `C_1` has angle `π/3` with its own `C_2` and is dotted to every other copy; the different-ray minimum is approximately 1.528128015. These decimal values are illustrative only: exact squared pair gates supply the certificate.

The 815-wall Gram has 5,636 `π/3` pairs, 2,808 `π/5` pairs, 7,657 right pairs, and 315,604 dotted pairs. All 331,705 pairs are covered by literal enumeration with cached exact distance classes. Its 245,800 pairs involving a new wall are all checked. The common point `(0,200,1)` has norm `−400` in the Lorentz form `|x|²−2DU` and is strict. Orthogonal projection to any individual wall therefore proves that all reported facets are actual. The complete pair gates also preserve all 392 finite endpoint points, with strict nonincident inequalities.

## Global geometry still fails

For the 815-wall chamber, put

`x=z_(3,3,0)/8`,

`X=(x,(|x|²−1)/2,1)`.

Then `X²=1`, every seed gap is at least `7/8`, and every spherical ball power is at least

`(21φ−15)/16 ≈1.186169610 >1`.

Thus every wall inequality at `X` is strict. It is future-oriented relative to the saved timelike center and proves infinite volume. Separate strict escapes are saved for the 417- and 479-wall checkpoints as well.

The spherical-cap finite-contact graph also shows the remaining separation. It excludes the vertical seed facets. The 417-, 479-, and 815-wall configurations have respectively **9, 385, and 393 connected components** in that graph. The largest construction covers all the tested radial probes but does not establish coverage of the alcove. The eight long chains connect selected endpoint facets, while the full second layer leaves many isolated cap components.

There is a second concrete deficiency. Every new cap has `π/3` contact to both members of each adjacent retained `A6` seed pair. The resulting three unit normals form affine `A2`: their Gram has positive kernel `(1,1,1)`. There are ten such pairs per new cap, explicitly listed. Pairwise acuteness places the resulting future null ray in the chamber.

The complete Gram gates also identify every incident wall at these particular null rays. Except for a level-eight cap’s own right endpoint mate, no other wall is orthogonal to all three positive summands. Consequently these cusp links have quotient rank two, or rank three when that endpoint mate is present; the required rank is nineteen. The checkpoints contain respectively 20, 640, and 4,000 such explicit deficient cusps. For the 479-wall case, 80 have quotient rank three and 560 have rank two. This is not a census of all cusps. Any finite-volume continuation retaining these walls must also repair these links.

## Why the outer step is a different kind of bridge

For any cap whose center lies strictly inside this golden alcove, the seven gaps in either `A6` factor sum to seven. One positive gap is at most one. Coxeter compatibility makes the corresponding positive cosine at least `1/2`, so the cap radius is at most two. If the cap has finite or parallel contact with the central radius-two sphere, its center length is at most `R+2`; its entire Euclidean ball therefore lies within radius six.

But `|7z/20|>6`. Hence no cap in that direct central-contact family can cut this probe. The outer `C_2` is dotted to central and extends beyond radius six, while retaining the exact `π/3` contact to `C_1`. The coupled construction crosses this restricted barrier. Boundary-center caps, other mixed profiles and arbitrary interior caps dotted to central remain separate possibilities.

The next useful work is a profile change that joins different chains while controlling the newly created links. More copies of the same radial recipe cannot by themselves be treated as a global closure result.

Verification:

```sh
python scripts/fresh_h20_balanced_golden_20261014.py --verify
```

This checker is a new exact implementation of the formulas above using the shared rational seed geometry. It does not claim an independent second replay of itself, a complete vertex census, or a finite-volume example.
