# Mixed face caps with coupled endpoint radii

A new mixed-factor construction gives **1,202 exact acute Coxeter facets**. Its 392 new caps meet neighboring new caps at `pi/3` and a resized endpoint layer at `pi/8`. They cut the square-centered escape of the preceding face hierarchy. However, they create **no new finite vertices**, introduce deficient parallel cusps, and leave another exact strict spacelike escape. This is a coherent local metric construction, not a finite-volume example.

The [checker](../scripts/fresh_h20_endpoint_replacement_20261014_mixed.py) and [certificate](fresh_h20_endpoint_replacement_20261014_mixed.json) pass generation and saved replay. A separate [FLINT matrix and proof review](fresh_h20_endpoint_replacement_20261014_mixed_review.md), with no producer imports, also passes generation and saved replay; it checks the transported pairs, actual geometry, both cusp representatives and escape. The preceding [continuous hierarchy obstruction](fresh_h20_endpoint_replacement_20261014.md) remains unchanged.

## The mixed orbit and its finite contacts

Use the unit-inradius `A6 × A6 × A7` cusp product. Write its 392 vertices as `z_abc`, with cyclic indices of lengths `7,7,8`; every vertex has squared norm 196. The factor vertex-distance identity is

`|v_i-v_j|² = 2n*d*(n-d)`, where `d=(i-j) mod n`.

Put

`R = 14-4sqrt(7)`, `Q=R²`, `t=R/7`,

and transport the center

`w = z_000 + t*(z_100-z_000) + t*(z_010-z_000)`

over all 392 cyclic translations. Its squared norm is `916-288sqrt(7)`, approximately `154.0236224`. Each cap has 17 right seed contacts, two parallel seed contacts and three dotted seed contacts. The nontrivial seed ratio is `b=(1+2sqrt(7))/6>1`, with `R*(1+b)=7`.

The identity `R²-28R+84=0` makes the squared distance between neighboring mixed centers equal to `3R²`. Consequently, all 76,636 distinct mixed-cap pairs are Coxeter: 784 are `pi/3`, and 75,852 are dotted. The complete pair check transfers from the 391 reference displacements by the exact cyclic distance identities.

Replace the entire endpoint layer by radius

`rho = R*(sqrt(119+49c²)/7-c)`, where `c=cos(pi/8)`.

Here `rho≈3.0335533327` and `rho²≈9.2024458224`. The defining radical identity gives

`rho²+2R*rho*c+R² = (24/7)R²`.

This proves the `pi/8` relation with the unique nearest endpoint cap. All other mixed–endpoint relations are dotted: 392 finite contacts and 153,272 dotted pairs. The equality is proved algebraically; rational interval enclosures certify the strict inequalities. The endpoint caps are mutually dotted because `rho²<21` while their minimum squared center distance is 84.

The mixed centers lie just off the old square circumcenters. For the previous target `x=99*(z_000+z_100+z_010+z_110)/400`, the mixed ball power is

`2270077/5000-176sqrt(7) < 0`.

Thus this orbit leaves the previously excluded single-factor face formulas and cuts their saved square-centered escape.

## Retained walls and actual facets

The final set comprises the 22 vertical seed walls, the central radius-two cap, 392 resized endpoint caps, all 395 old middle/motif caps from the [810-wall checkpoint](fresh_h20_balanced_integer_exact_20261014.json), and the 392 mixed caps.

Every retained middle/motif ball has center norm less than six and radius at most two, so lies strictly inside the radius-eight ball about the origin. The mixed balls and resized endpoint balls lie strictly outside that ball. These bounds certify all cross pairs as dotted without a floating-point pair test. The old retained subset inherits its exact pair certificate through checked source hashes.

In the Lorentz convention `|(x,D,U)|²=|x|²-2DU`, the common point `(0,200,1)` is strict and timelike. Projecting it onto any wall, using the nonpositive pairings, gives a timelike point strict on every other wall. Thus all 1,202 walls are actual facets.

All 392 old vertical endpoints survive, now at height `rho`: the corresponding point is

`P_z = (z,(196+rho²)/2,1)`, with squared norm `-rho²`.

Its incident diagram is `A6 + A6 + A7 + A1`; every other inequality is strict. No complete finite-vertex census is claimed.

## Why the new contacts do not close the frontier

A mixed cap has only 17 right seed neighbors, four mixed neighbors at `pi/3`, and its unique endpoint neighbor at `pi/8`. The last five neighbors are pairwise dotted. Hence a finite-pair clique containing a mixed cap has at most 19 members. A finite vertex in `H20` requires 20 independent incident normals with positive-definite Gram matrix. There are therefore **no finite vertices incident to a mixed cap**.

Each mixed cap also has two parallel seed partners. For one reference partner, the complete incident diagram is

`affine A1 + spherical A4 + A5 + A7`.

Its Euclidean quotient rank is 17 rather than 19, so its cusp link is unbounded. Cyclic transport gives 784 distinct listed cusps with this diagnosis; this is not a census of every cusp. As explained in the [810-wall review](fresh_h20_balanced_integer_exact_20261014_review.md), deficient links of this form can admit further local completion. Their deficiency is not a proof against arbitrary later additions.

## The surviving global escape

Let `x=99*(z_000+z_001)/200`, an inward perturbation of an A7 edge midpoint. Then

`|x|²=205821/1250`, `X=(x,204571/2500,1)`, `X²=1`.

Every vertical seed gap is at least `1/100`. The exact or rigorously enclosed minimum ball powers at `x` are:

| Family | Minimum or certified lower bound |
|---|---:|
| Mixed orbit | `965821/1250-(6776/25)sqrt(7)` |
| Resized endpoints | greater than `18.81435417` |
| Old 395 middle/motif caps | `164217/2500` |
| Central cap | `200821/1250` |
| Every permitted real vertex-radial cap | `702863/30000` |

Every bound exceeds one. The point `X` is future-oriented and strict against all 1,202 walls. It also survives the continuous vertex-radial relaxation from the preceding report. Consequently, this chamber has infinite volume, and adding caps only from that radial family cannot repair the remaining gap.

A matched three-factor variant would additionally shift the A7 center by `R/8`. Its squared distance from the current mixed center is `(7/4)R²`, giving positive mutual pairing `R²/8`. It therefore requires replacing the matched current cap; it cannot be appended to this retained orbit. No further variant orbit is constructed here.

```sh
python scripts/fresh_h20_endpoint_replacement_20261014_mixed.py --verify
```

The exact scope is the stated cyclic orbit, resized endpoint layer, retained old subset, finite-neighbor argument, listed cusp links and strict escape. No arithmetic-lattice, global finite-volume or complete cusp-census claim is made.
