# Multiple basepoints and targeted H9 escape cuts

This batch finds additional **exact Coxeter facets**, including cuts that
remove previously certified escape directions. It does not find a
finite-volume or compact polytope. The final odd-A1 cone has 372 facets
and a new independently verified spacelike escape.

No uniform Vinberg height cutoff was increased. The work combines
previous bounded searches through explicit lattice isometries, then
solves small arithmetic problems aimed directly at surviving escapes.

## Compatible wall unions

The H4-to-E8 isometries use the Bugaenko E8 basis, while the earlier
searches use the standard E8 basis. The exact ambient-coordinate map
from Bugaenko indices to standard indices is

```
(0,1,2,3,4,5,6,7,8,9) -> (0,1,8,2,3,4,5,6,7,9).
```

After this permutation, the translated H4 height-30 walls are compatible
with both earlier cones. Identical normals are removed. Every remaining
cross-pairing is nonpositive, and both archived strict interior points
remain strictly feasible for the union. No reflection word was needed.
The independent verifier checks the entire Gram congruence and all these
inequalities.

A further exact symmetry swaps the two H4 blocks. If `T` is the transport
to standard E8 coordinates and `S` is that swap, then `U=T*S*T^-1` is an
integral isometry with `U^2=1`. Applying it to the earlier wall lists
adds substantially more facets. The resulting unions are invariant
under `U`.

| Host `<-alpha>+E8+<1>` | Old complete height | Old facets | Direct H4 additions | Further swap additions | Merged facets |
|---|---:|---:|---:|---:|---:|
| alpha=2phi | 300 | 184 | 10 | 170 | 364 |
| alpha=sqrt5 | 100 | 43 | 10 | 32 | 85 |

Every added wall has old-basepoint height above the old cutoff. Their
height ranges are, respectively,

```
119+193phi <= h <= 20736+33552phi,
38+(123/2)phi <= h <= (19251+31149phi)/2.
```

These are isolated additional walls from different basepoints, not an
exhaustive search through either maximum height.

Acuteness and a shared interior do not by themselves prove that two
prefixes belong to the same full reflection chamber. The claim here is
the directly verified compatible Coxeter wall union. For every distinct
pair, the independent verifier checks nonpositivity and checks

```
4*(r,s)^2 / ((r,r)*(s,s))
```

against `{0,1,2,3,1+phi,2+phi}` whenever it is below four. These correspond
to angles `pi/m` for `m=2,3,4,6,5,10`; values at least four give parallel
or disjoint supporting hyperplanes. In particular, the non-Coxeter acute
angle `2pi/5` is rejected. This gate covers **all** pairs in both merged
cones, including pairs from different basepoints.

Each distinct normal actually defines a facet. Orthogonally projecting
the strict timelike interior onto its hyperplane remains timelike and
strictly satisfies every other inequality, because all cross-pairings
are nonpositive. Thus the additions are not redundant inequalities.

Both old escape vectors survive these unions. Their norms remain
`(61-37phi)/49` and `(38486-23762phi)/841`. Neither escape is fixed by
`U`; its symmetry image is another feasible spacelike direction.

## Directly cutting an escape

For an escape `v=(1,v_spatial)`, saturate its spatial line in the integral
lattice to obtain a primitive vector `z=c*v_spatial`, with `c>0`. Put
`beta=(z,z)`. Vectors in the rank-two slice through the original timelike
basepoint `p` and `v` have the form

```
r=(x,t*z),       beta*t^2-alpha*x^2=q,       x,t in Z[phi].
```

If `t>0` and `x>=c*t`, then `r=(x-c*t)*p+c*t*v` pairs nonpositively
with every existing wall. An exact solution is accepted only after
checking integrality, reflection coefficients, its positive pairing
with the targeted escape, its negative pairing with the common
interior, and every finite Coxeter angle. This is a nonchronological
construction of compatible reflection cuts, not a claimed continuation
of a complete Vinberg transcript.

For the first odd-A1 escape, saturation gives

```
z = (14,28,40,32,24,16,8,20,sqrt5),       beta=61.
```

The norm-one equation has the exact solution

```
t=11+18phi,       x=48+78phi,
61*t^2-2phi*x^2=1.
```

Its pairing with the targeted escape is `(6+5phi)/7>0`, and its
old-basepoint height is `27144+43920phi`. Only five conjugate-admissible
values of `t` were visited before this cut was found. Its symmetry mate
is compatible too: the two roots have pairing `-6071-9824phi`.

Three further bounded adaptive rounds find fresh exact escape vertices,
try their simplest saturated rank-two slices, and add three more cut
orbits. Each round considers at most four planes, with `0<t<=100000`;
the source saves every attempted plane. A numerical linear program
proposes an extreme-ray support, and exact field arithmetic verifies
the resulting ray before any arithmetic search uses it.

| Cut orbit | Root norm | beta | x | t | Old-basepoint height |
|---|---|---|---|---|---|
| Initial escape | 1 | 61 | `48+78phi` | `11+18phi` | `27144+43920phi` |
| Adaptive 0 | 1 | `40-11phi` | `5+8phi` | `2+3phi` | `288+466phi` |
| Adaptive 1 | 2 | `9+phi` | `4+7phi` | `2+4phi` | `105+170phi` |
| Adaptive 2 | 1 | 65 | `2182+3532phi` | `487+788phi` | `55777344+90249640phi` |

Each orbit supplies two walls. The final cone therefore has **372
facets**, compared with 184 in the original odd-A1 prefix. All finite
angles introduced by the eight targeted walls are right angles or
`pi/4`. The independent checker nevertheless checks the full angle gate
for every pair of the 372 walls.

The former escapes have been cut off, but a new exact escape remains:

```
(v_final,v_final) = (424-262phi)/175 > 0.
```

Its active facets include the first newly added cut. This explains the
geometric outcome: removing one escape can create new hyperideal
intersections on the added facet. It is not a global obstruction to
further completion.

The sqrt5 escape was also targeted. Its saturated spatial norm is
`beta=73-8phi`. Norms one and two are excluded in that plane by exact
congruences at primes 31 and 151. The two remaining norm classes
`2+phi` and `4+2phi` yield no cut for `0<t<=1000000`, after testing
57,504 reflection-admissible values. This is only a failure in the
specified rank-two slice. It does not exclude a compatible cut elsewhere
in the lattice.

## Independent checks and artifacts

The standalone verifiers use only `Fraction` arithmetic and their own
field and matrix routines. They replay both merges, all four cut orbits,
the final 372-wall cone, all finite angles, and all six saved local
obstructions. The square-root and bounded plane searches themselves are
not independently rerun by these certificate checkers.

All six regression tests pass, including an actual integral H4 root pair
at angle `2pi/5` that must be rejected, a forged cut norm, and a forged
local obstruction.

- [Merge generator](../scripts/fresh_multibasepoint_h9_merge.py)
- [Direct escape-plane search](../scripts/fresh_multibasepoint_h9_escape_plane.py)
- [Three adaptive rounds](../scripts/fresh_multibasepoint_h9_targeted_rounds.py)
- [Merge verifier](../scripts/fresh_multibasepoint_h9_verify.py)
- [Targeted-cut verifier](../scripts/fresh_multibasepoint_h9_targeted_verify.py)
- [364-wall odd-A1 union](fresh_multibasepoint_h9_2phi_20261013.json)
- [85-wall sqrt5 union](fresh_multibasepoint_h9_sqrt5_20261013.json)
- [Final 372-wall construction](fresh_multibasepoint_h9_2phi_targeted_rounds_20261013.json)
- [Merge audit](fresh_multibasepoint_h9_independent_audit_20261013.json)
- [Targeted audit](fresh_multibasepoint_h9_targeted_independent_audit_20261013.json)
- [Regression tests](../tests/test_fresh_multibasepoint_h9.py)

```sh
PYTHONPATH=. python3 scripts/fresh_multibasepoint_h9_verify.py \
  results/fresh_multibasepoint_h9_2phi_20261013.json \
  results/fresh_multibasepoint_h9_sqrt5_20261013.json
PYTHONPATH=. python3 scripts/fresh_multibasepoint_h9_targeted_verify.py \
  results/fresh_multibasepoint_h9_2phi_targeted_rounds_20261013.json
PYTHONPATH=src:. python3 -m pytest -q tests/test_fresh_multibasepoint_h9.py
```

No new finite-volume example, new compact dimension, global
nonreflectivity theorem, or completeness beyond the original uniform
height bounds is claimed.
