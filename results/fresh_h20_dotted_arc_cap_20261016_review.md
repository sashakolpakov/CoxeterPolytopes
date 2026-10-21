# Independent review of the H20 dotted-arc cap

**PASS.** The cap coherently extends the retained 133-wall chamber to 134
actual facets. It has exactly two new finite vertices, cuts the stronger
847-bank escape, and leaves a strict escape proving **infinite volume**.
No discrepancy was found in the producer's geometric claims.

Artifacts: [independent source](../scripts/fresh_h20_dotted_arc_cap_20261016_review.py),
[literal review certificate](fresh_h20_dotted_arc_cap_20261016_review.json),
and [producer report](fresh_h20_dotted_arc_cap_20261016.md).

The review imports no producer arithmetic. It uses the frozen independent
sparse-radical `Field` engine, rational support inversion, the independent
complex-power implementation of the chain, and outward 512-bit rational
intervals. All 107 old physical contact covectors are reconstructed from
literal normal coordinates. The unique chain parameter and old 133-wall
acuteness are explicitly inherited from the frozen independent chain review.

## Contacts and one coherent normal

The retained support has rank 18. Wall 31 is an isolated norm 2 root there, so
`omega31=n31/2`, with square `1/2`, is its dual weight. The independently
checked Lorentz3 frame is orthogonal to the support.

Write the central chain normal as `n129=a*x+rho*v`, where
`rho^2=a^2+2`, `x^2=-1`, `v^2=1`, and `x.v=0`. With
`h=sqrt2-1`, define

`t=(h*a^2+2sqrt2)/(2a)`, `y=t*x+(h*rho/2)*v`,
`c=sqrt(5+4/a^2)`, and `r=y-c*omega31`.

Exact identities in `Q(sqrt2)` establish
`y^2=-1/2-2/a^2` and `r^2=2`. The central contact is
`y.n129=h-sqrt2=-1`. The two adjacent chain contacts are exactly `-sqrt2`,
using their already certified orthogonality to `n129`; the 17 remaining
support contacts are exactly zero. These equalities are symbolic, not
inferred from intervals containing zero.

Independent strict interval checks complete all 133 contacts:

| Contact | Count |
| --- | ---: |
| Right | 17 |
| Angle `pi/4` | 2 |
| Angle `pi/3` | 1 |
| Dotted | 113 |

The new normal is a single vector in the same ambient Lorentz space, rather
than a proposed Gram row with unresolved realization.

## Actual facets and the complete new-corner census

The inherited center recipe is checked against all 134 walls, with strictly
negative pairings and negative square. For each norm 2 wall `n`, its projection
`C-(C.n)n/2` has negative square, is future, and lies strictly inside every
other wall halfspace. Thus every listed wall is an actual facet.

There are 20 old walls having finite contact with the new wall: the 17 support
roots other than 31, and 128, 129, 130. A new finite vertex requires 19 of them.
Among all 20 candidate 19-subsets, 18 contain the dotted pair 128/130 and are
therefore not spherical. The two remaining sets use chain pairs 128/129 and
129/130. Their final Schur pivot is exactly `1/2`.

Both vertices are solved independently in the Lorentz3 frame plus the
`omega31` direction. Their negative squares, future orientation, all 133 old
wall inequalities and defining zero contacts are checked. Exactly two new
finite corners therefore occur. Their incidence with the two inherited
finite chain vertices gives three actual finite edges; each shared support
has rank 19 and the connecting segment satisfies every other inequality
strictly in its interior.

The acute-extension argument strictly retains every old finite vertex. The
inherited lower bound 168 hence becomes at least 170. This does not establish
a complete ambient finite-vertex census.

## The cap cuts one escape and leaves another

The stronger 847-bank escape has strictly positive pairing with the new wall
and is cut. The older 133-wall escape has strictly negative pairing and
survives. A separate escape is also independently reconstructed:

`E=x-e/100-omega5-(1/100000) sum_{i in S, i!=5} omega_i`.

Its squared norm is exactly `190038001053/110000000000`, computed from the
orthogonal Lorentz3/support decomposition. All 134 wall pairings are strictly
negative, and its pairing with the strict future center is negative. The ray
`C+tE` reaches the future null cone while retaining every strict wall
inequality. This supplies an open ideal-boundary patch and proves infinite
volume. The three inherited proper cusps also have strictly negative new
wall pairings.

## Validation and limits

Generation and a fresh saved-output replay pass. Source hashes are checked
recursively. Replay command:

```sh
python3 scripts/fresh_h20_dotted_arc_cap_20261016_review.py --verify
```

The review establishes this exact partial construction, its actual facets,
complete new-corner census and surviving escape. The old ambient census and
ray census remain incomplete. It does not certify a finite-volume polytope,
an arithmetic lattice, or an example in a new dimension.
