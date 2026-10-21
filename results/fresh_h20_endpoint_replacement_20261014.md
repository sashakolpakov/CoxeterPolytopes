# Changed endpoints and a face hierarchy: a continuous global gap

The endpoint replacement opens a real local construction mechanism, but the tested face hierarchy does **not** cover the cusp chamber. A single exact strict escape survives the entire continuous hierarchy, every vertex-radial cap permitted by its seed contacts, and all 395 previously saved middle/motif caps. This rules out finite volume for every finite compatible selection from this stated union, including selections which delete optional walls.

The [standard-library exact checker](../scripts/fresh_h20_endpoint_replacement_20261014.py) and [certificate](fresh_h20_endpoint_replacement_20261014.json) pass generation and saved replay. The separate [metric scout](../scripts/fresh_h20_endpoint_replacement_20261014_metric_scout.py) and [numerical archive](fresh_h20_endpoint_replacement_20261014_metric_scout.json) record a compatible numerical proposal; they are not an exact chamber certificate. All predecessor files remain unchanged.

## Why longer vertex-radial chains do not suffice

The unit-inradius `A6 × A6 × A7` product has 392 vertices `z`, each with squared norm 196. A cap centered at `gamma*z`, with `gamma<1`, has 19 equal positive seed gaps `1-gamma`. Coxeter acuteness therefore forces `R <= 2(1-gamma)`, without any finite-label cutoff. A center with `gamma>1` violates a seed inequality. The case `gamma=1` is the replacement endpoint layer.

For a test point `x`, relaxing even to every real `gamma` gives the lower bound

`power(x) >= |x|² - 4 - (x·z - 4)²/192`.

Let `m` be the midpoint of a nearest A6 edge: its two endpoints have squared distance 84 and `|m|²=175`. At `x=99m/100`, every radial ball has power at least `1407119/76800`, and every endpoint ball with `R²<=21` has power at least `7/400`. All seed gaps are positive. Taking the Lorentz norm to be `1/100` therefore gives a strict escape for this complete radial-plus-endpoint relaxation. Extending a chain indefinitely cannot remove this gap while remaining in these families.

## A different local mechanism, with one metric coupling

Write `r` for the common squared endpoint radius. For nearest vertices `z0,z1` in a factor with `n` affine nodes (`n=7` or `8`), put a cap at their midpoint with

`R_edge² = n(n-1)/2 - r`.

It is right to both endpoint caps. In particular, `r=16` gives radius `sqrt(5)` on an A6 edge. Its 18 retained seed contacts are right and its other seed contacts are dotted.

For three consecutive vertices of the same factor, use the shifted center and radius

`w = (r*z0 + (n²-2r)*z1 + r*z2)/n²`,

`R_triangle² = 4r²/n² - r`.

This cap is right to the middle endpoint and the two adjacent midpoint caps. The exact checker verifies these norm/contact identities. The shift is essential: the true A6 triangle circumcenter would require squared radius 20 to be right to radius-four endpoints, but squared radius 10 to be right to the two radius-`sqrt(5)` midpoint caps.

At `r=16`, adjacent shifted A6 triangle caps have cosine `109/120`, so the whole orbit fails the Coxeter gate. The radius can instead be coupled by the equation

`(4-4sqrt(5))*r² + (49sqrt(5)-539)*r + 8232 = 0`.

Its unique positive root is approximately `16.1620249977362`, strictly between 16 and `49/3`. The adjacent A6 triangle relation is then exactly `pi/5`. The numerical complete-orbit test uses the central cap, 392 endpoints, 1,176 nearest-edge caps, and 1,176 shifted triangle caps. All 3,766,140 cap pairs pass: 5,880 right, 784 at `pi/5`, and 3,759,476 dotted. All seed gates and 1,084,275 pair gates against the old 395 middle/motif caps also pass numerically. The exact contact equations make this a reproducible metric proposal; no full exact pair replay or complete vertex census is claimed.

## The global test which stops this hierarchy

Let

`m = (z_000 + z_100 + z_010 + z_110)/4`, `x = 99m/100`.

This is an inward perturbation of a square center mixing the two A6 factors. Its squared norm is `754677/5000`. The certificate considers every family member for **every** `16<=r<=49/3`; different selected caps may even use different parameter values. Zero-radius A7 triangle limits at `r=16` can be omitted and cannot affect the result.

| Family | Uniform lower bound on ball power at `x` |
|---|---:|
| Central radius-two cap | `734677/5000` |
| All 392 endpoints | `385231/15000` |
| All real vertex-radial caps | `15428063/480000` |
| All 1,176 nearest-edge caps | `80077/5000` |
| All 1,176 shifted triangle caps | `381731/15000` |
| All 395 old middle/motif caps | `295677/5000` |

The edge powers are affine in `r`. For the triangle powers, the quadratic center-norm term cancels the quadratic radius term, leaving another affine function of `r`; exact endpoint evaluation therefore covers the whole interval. All 392 transports are checked for the radial bound. Compatibility of the relaxed union is unnecessary for this exclusion.

Every listed power exceeds one, while the minimum vertical seed gap is `1/100`. Hence

`X = (x, (|x|²-1)/2, 1)`

has norm one, is future-oriented, and is strict against the whole union. This is a global infinite-volume certificate, not merely an open-edge diagnostic.

The remaining constructive direction must leave these formulas. Genuinely mixed-factor face caps, or differently positioned caps even within a single factor, are not excluded. The metric coupling above repairs a local orbit conflict, but it does not repair the cross-factor covering gap.

```sh
python scripts/fresh_h20_endpoint_replacement_20261014.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_endpoint_replacement_20261014_metric_scout.py --verify
```
