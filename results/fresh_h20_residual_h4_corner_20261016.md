# A residual H4 cluster creates a new finite H20 corner

Three exact rational caps extend the
[142-wall mixed-support chamber](fresh_h20_mixed_isolated_cap_20261016.md)
to **145 actual facets** and create **exactly one new finite vertex**.
Together with the existing mixed cap s (wall141), the new caps u,v,w
are pairwise right. Their common finite corner has diagram
`A10+A4+A2+4A1`. This closes a genuine finite corner in the residual H4
geometry, although a continuous family of strict escapes survives and
the full H20 chamber still has infinite volume.

The [source](../scripts/fresh_h20_residual_h4_corner_20261016.py) and
[certificate](fresh_h20_residual_h4_corner_20261016.json) use exact rational
Gram identities and outward rational intervals for the old physical
normals. Generation and fresh saved replay both pass. Independent
review is not included in this primary package. No frozen predecessor
is changed.

In the chart `(x,e,f,omega17,omega31)`, whose Gram matrix is
`diag(-1,1,1,1/2,1/2)`, the four normals are

| Wall | Chart coordinates | Square |
|---|---|---|
|141, s|`(17,-3/2,-67/4,-2,-3)`|`5/16`|
|142, u|`(13,0,-64/5,-21/10,-3)`|`309/200`|
|143, v|`(2727/100,-773/300,-267/10,-5,-5)`|`19717/22500`|
|144, w|`(276851/39846,-8739/13282,-135452/19923,-1,-3/2)`|`6567095/1058469144`|

All six mutual products are exactly zero. The arbitrary positive norms
are retained. Every new normal is right to the 16 old support roots
`S16=S\{17,31}`. The wall u is also right to old95 because its e
coefficient is zero. These right relations follow algebraically;
interval enclosures containing zero are only consistency checks.

The full new-old and new-new census consists of **429 pairs**, with
**55 right** and **374 strictly dotted** relations. Every dotted pair
is certified by `g<0` and `g^2>q_i*q_j`. This includes all26 chain walls,
all five A17 cusp connectors and the previous off-cusp wall140.
In particular u.wall140 is exactly `-5/2`. The same strict future center
lies inside every halfspace. The acute projection formula then supplies
a strict interior point for each of the145 facets.

The unique new finite corner has incident set

```
S16 union {141,142,143,144}.
```

Its exact ray, with x coefficient1, is

```
C=(1,-5971425/64223509,-314433802/321117545,
     -46580508/321117545,-67243924/321117545),
C^2=-10669433158276/103116477706827025<0.
```

The checker verifies all20 incident equations, strict negativity for
every other wall, the future chart, and all20 positive rational LDL
pivots of the incident Gram matrix. The old S16 block is positive and
the other four blocks are orthogonal A1 roots.

The new-corner census is complete even though the old ambient finite
census is not. In the final145-wall chamber the finite neighbors of the
new caps are

```
u: S16 + {95,s,v,w},
v: S16 + {s,u,w},
w: S16 + {s,u,v}.
```

Any finite corner containing a new cap must contain20 independent
pairwise finite normals. The lists for v and w force the displayed
corner. For u, any alternative19-neighbor choice includes95 together
with at least one of s,v,w, which is a dotted pair. The displayed corner
is therefore the only new finite vertex. All old finite vertices remain
strict for the new inequalities, so the total certified lower bound
increases from170 to171. No old finite vertex shares19 walls with this
new corner, hence no finite edge joins it to an old finite vertex. Ideal
endpoints of its edge directions have not been classified.

Each new cap is strictly negative at all four previously certified
proper cusp points, including nu with its completed affineG2+affineA17
link. These cusps persist unchanged. Both u and v cut the previous
escape E3:

```
u.E3=21/250,
v.E3=1403/30000.
```

A new strict spacelike escape is

```
E4=x-e/100-(49/50)*f-(3/10)*omega17-(1/20)*omega31
   -(1/100000)*sum(other16 support dual roots),
E4^2=33750061/5000000000>0.
```

More generally, replace its isolated coefficients by `-a,-b`. Every
point in the entire rectangle

```
29/100 <= a <= 31/100,
3/100 <= b <= 3/50
```

is a strict spacelike escape. Every old or new wall pairing and the
center pairing is affine in a,b and is strictly negative at all four
rectangle corners. The squared norm is a constant plus `(a^2+b^2)/2`,
so its minimum on this positive rectangle is exactly
`15000061/5000000000>0`. Thus the remaining obstruction is a continuous
escape region, not only a single chosen ray.

The numerical discovery searches were bounded and are not proof inputs.
The exact conclusion is the stated145-wall Coxeter extension, its one
new finite corner, preserved cusps, and infinite volume. It is not a
finite-volume construction or a new-dimension example.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_residual_h4_corner_20261016.py --verify
```
