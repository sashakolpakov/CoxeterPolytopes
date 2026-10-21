# A finite/parallel residual partner cuts E4

An exact norm-2 wall extends the
[145-wall chamber](fresh_h20_residual_h4_corner_20261016.md) to **146 actual
facets** and cuts its strict escape E4. It meets wall140 at pi/3, is
parallel to141, is right to142 and143, and is strictly dotted to144.
Every other nonright old contact is also dotted. The four proper cusps
and all old finite vertices persist. A new ideal vertex is certified,
but its link is unbounded; a strict spacelike escape E5 remains.

The [source](../scripts/fresh_h20_h4_parallel_partner_20261016.py) and
[certificate](fresh_h20_h4_parallel_partner_20261016.json) give the exact
radical recipe and all strict sign bounds. Generation and fresh
saved replay both pass. This primary package includes no independent
review of the new wall and changes no frozen predecessor.

Let s=wall141, w=wall144, and C the exact finite corner of the four
orthogonal residual caps. Set k=-C^2>0 and q_i=n_i^2. Define

```
A=-sqrt(2/q_141)*s,
B=sqrt(2/k)*C-sqrt(2/q_144)*w,
z=(-sqrt(q_140/2)-A.n_140)/(B.n_140),
r=A+z*B.
```

The denominator is strictly negative and z>1 (approximately5.023193).
The orthogonal residual basis gives exact identities

```
A^2=2, B^2=0, A.B=0, r^2=2,
r.n_140=-sqrt(q_140/2),
r.s=-sqrt(2*q_141),
r.n_142=r.n_143=0,
r.w=-z*sqrt(2*q_144).
```

Thus the specified pi/3, parallel, right and dotted relations are exact,
not rounded angle tests. The remaining16 support-right relations follow
from the dual-root construction. The full145 new-old checks give
**18 right, one pi/3, one parallel and125 dotted** pairs. All26 chain
contacts and all five earlier A17 connectors are dotted. The retained
strict future center has negative pairing with r; the acute projection
formula verifies every one of the146 facets is actual.

The E4 cut has certified positive pairing approximately0.00188526.
This wall is in the class singled out by the
[all-real residual cutter gate](fresh_h20_residual_cutter_domain_20261016.md):
a residual E4 cutter must be dotted to144. The gate is context, not an
arithmetic dependency of this construction.

There are no new finite vertices. The only19 finite old neighbors are
`S16 + {140,142,143}`. Their sole candidate20-root support with r contains
the dotted pair140/142 and cannot be positive definite. Every old finite
vertex persists, including the residual A1^4 corner. The certified total
remains at least171; the old ambient finite census remains incomplete.

There is an explicit new ideal vertex

```
nu2=C/sqrt(k)-w/sqrt(q_144),
r/sqrt(2)+s/sqrt(q_141)=z*nu2.
```

Its norm is zero and it is future-oriented. Its complete incidence is
`S16 + {141,142,143,145}`. All other old wall pairings are strictly
negative. The positive rank18 block S16+{142,143} is orthogonal to the
two independent parallel-pair normals; hence the20 incident normals
have actual rank20 and Gram rank19. This proves that the null ray is an
actual ideal vertex.

Its link diagram is

```
affine A1 + A10 + A4 + A2 + 2A1.
```

The spherical factors make this Euclidean link unbounded. It is not a
new proper cusp. This also identifies a concrete local completion
problem: the five spherical components would each need affine closure.
The present certificate does not add those further walls.

All four previously certified proper cusp rays have strictly negative
pairing with r, so their bounded links remain unchanged. A new exact
strict escape is

```
E5=x-e/100-(49/50)*f-(3/10)*omega17-(1/25)*omega31
   -(1/100000)*sum(other16 support dual roots),
E5^2=31500061/5000000000>0.
```

Every one of the146 pairings is strictly negative, with E5 in the future
chart of the common center. Thus volume is still infinite. No complete
ideal-vertex census, finite-volume construction or new-dimension claim
is made.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_h4_parallel_partner_20261016.py --verify
```
