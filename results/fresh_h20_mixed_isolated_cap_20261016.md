# A mixed-support cutter: 142 actual H20 facets

The exact normal

```
s = 17*x-(3/2)*e-(67/4)*f-2*omega17-3*omega31
```

has square **5/16** and is compatible with every wall of the
[reviewed 141-wall chamber](fresh_h20_off_cusp_cap_20261016_review.md).
It cuts that chamber's strict escape E2 by **s.E2=33/100**. The resulting
chamber has **142 actual facets**. Its four certified proper cusps and
all old finite vertices persist. A new strict escape E3 survives, so the
volume remains infinite.

The [source](../scripts/fresh_h20_mixed_isolated_cap_20261016.py) and
[certificate](fresh_h20_mixed_isolated_cap_20261016.json) reconstruct the
physical old normals and verify strict inequalities by outward rational
intervals. The recipe, norms and named escape pairings are exact rational
identities. Generation and fresh saved replay both pass. This primary
package does not include an independent review of the new wall.

The two support dual roots are orthogonal, each with square1/2, and are
perpendicular to the orthonormal Lorentz frame x,e,f. Therefore

```
s^2 = -17^2+(3/2)^2+(67/4)^2+2^2/2+3^2/2 = 5/16.
```

The positive norm is retained without a square-2 rescaling. The wall is
right to the 16 roots of the original support other than17 and31. This
follows exactly from the dual-root definition; interval zero enclosures
are consistency checks. Its other **125** old contacts are strictly
dotted, including both isolated roots, all26 chain walls, both earlier
escape-cutting caps and all five A17 connectors. For every nonright old
normal n the certificate checks

```
s.n<0, (s.n)^2>s^2*n^2.
```

In particular, with the previous wall
`r=8*x-(3/2)*e-(77/10)*f-(14/5)*omega17`,

```
s.r=-79/40.
```

The same strict future center H lies on the negative side of all142
walls. Acuteness makes each projected point
`H-(H.n/n^2)*n` strict for all other inequalities, with negative square;
all walls remain actual. The positive-definite incident-Gram argument
preserves every old finite vertex. There are no new finite vertices,
since s has only16 finite old neighbors and a finite H20 vertex needs20
independent incident walls. The inherited count at least170 remains a
lower bound; no complete ambient vertex census is asserted.

All four proper cusp rays pair strictly negatively with s, including
nu, whose completed affineG2+affineA17 link is unchanged. The cut escape
E2 and a new surviving escape satisfy

```
s.E2=33/100,
E3=x-e/100-(49/50)*f-(3/10)*omega17-(3/20)*omega31
   -(1/100000)*sum(other16 support dual roots),
E3^2=83750061/5000000000>0,
s.E3=-9/200,
r.E3=-19/1000.
```

Every one of the142 E3 pairings is strictly negative, and E3 has negative
pairing with H. This gives a strict spacelike escape and certifies
infinite volume.

A bounded numerical search over121 fixed pairs of isolated contact
coefficients suggested this recipe. That search is not a proof dependency
and does not exhaust mixed-support caps. The exact certificate concerns
this one wall and its stated consequences. It neither closes the entire
continuous escape region nor constructs a finite-volume example.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_mixed_isolated_cap_20261016.py --verify
```
