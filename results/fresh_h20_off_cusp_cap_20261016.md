# An exact off-cusp cutter: 141 actual H20 facets

The algebraic wall

```
r = 8*x - (3/2)*e - (77/10)*f - (14/5)*omega17
```

is compatible with all 140 walls of the
[completed G2+A17 cusp chamber](fresh_h20_affine_a17_completion_20261016.md)
and cuts its strict escape E1. It has positive square **73/50**; there is
no need to force a square-2 normalization or enlarge the field for this
recipe. The resulting Coxeter chamber has **141 actual facets** and
retains all four previously certified proper cusps. A new strict escape
remains, so its volume is still infinite.

The [source](../scripts/fresh_h20_off_cusp_cap_20261016.py) and
[certificate](fresh_h20_off_cusp_cap_20261016.json) reconstruct all old
physical normals and use rational outward intervals for strict signs.
The norm and named escape pairings below are exact rational identities.
Generation and fresh saved replay both pass. Independent review of this
new cap is not included in this primary package.

The frame x,e,f is the exact orthonormal Lorentz frame of the original
133-wall chain construction. The vector omega17 is the old positive
support dual root, with square1/2 and perpendicular to that frame.
Consequently

```
r^2 = -8^2 + (3/2)^2 + (77/10)^2 + (14/5)^2/2 = 73/50.
```

The new wall is right to the 17 support roots other than root17. Every
one of its other 123 contacts is strictly dotted. The checker verifies
`r.n<0` and `(r.n)^2>2*r^2` for every such old square-2 root n, including
all 26 rotational-chain walls, both recent caps133/134 and all five A17
connectors. No finite-label approximation is used. The weakest
normalized dotted product is about -2.29077, at old cap133.

The retained strict future center H pairs negatively with r. For each
wall n, `H-(H.n/n^2)*n` has negative square and is strict for every other
wall by acuteness. Thus all 141 walls are actual facets. The usual
positive-definite incident-Gram decomposition of H proves that every
old finite vertex persists. There are no new finite vertices: only 17
old walls are finite neighbors of r, which is insufficient for a
20-dimensional finite corner. The already certified lower bound of170
finite vertices remains a lower bound, not a complete ambient census.

The pairing with the completed cusp ray nu is strictly negative. Its
bounded G2+A17 horospherical link therefore stays unchanged. The three
other previously certified cusp points also pair strictly negatively
with r; all four proper cusps persist.

The original escape is cut by the exact amount

```
r.E1 = 961/1000 > 0.
```

This wall necessarily lies outside the
[through-nu obstruction](fresh_h20_null_link_escape_gate_20261016.md):
it is nonright to the retained affine G2 triple and has `r.nu<0`.

A new exact strict escape is

```
E2 = x-e/100-(49/50)*f-(3/10)*omega17-(2/5)*omega31
     -(1/100000)*sum(other16 support dual roots),
E2^2 = 427500061/5000000000 > 0,
r.E2 = -19/1000 < 0.
```

All 141 pairings are strictly negative, and E2 lies in the future affine
chart of the retained center. This certifies infinite volume. The
construction removes the particular open region witnessed by E1, but
does not establish a finite-volume H20 example or a new dimension.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_off_cusp_cap_20261016.py --verify
```
