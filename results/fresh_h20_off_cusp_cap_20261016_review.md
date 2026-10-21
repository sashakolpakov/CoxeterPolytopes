# Independent review of the rational off-cusp H20 cap

**PASS:** the exact wall

```
r = 8*x - (3/2)*e - (77/10)*f - (14/5)*omega17
```

cuts E1 and coherently extends the 140-wall chamber to 141 actual facets.
It is strictly away from the completed cusp and preserves all four
proper cusps. A new strict spacelike escape E2 survives, so the chamber
still has infinite volume.

The [independent checker](../scripts/fresh_h20_off_cusp_cap_20261016_review.py)
and [certificate](fresh_h20_off_cusp_cap_20261016_review.json) pass
generation and fresh saved replay. The
[primary construction](fresh_h20_off_cusp_cap_20261016.md) is separate.
The reviewer independently reconstructs all 140 retained normals,
including the five cusp connectors, using the earlier independent field
and interval engines. It does not import the producer's arithmetic.

In the exact Lorentz frame, the old support is perpendicular to x,e,f,
and omega17 is its Gram-dual vector. Direct rational calculation gives

```
r^2 = 73/50,
r.E1 = 961/1000 > 0.
```

There are 17 exact right contacts, namely the old 18-root support with
root 17 removed. These equalities follow from the frame and dual-vector
identities. Every other old wall, including all five cusp connectors,
is strictly dotted. The reviewer checks all 123 negative pairings and
their squared inequalities `(r.n_i)^2 > 2*r^2` with rational bounds.
It also compares each contact with the primary certificate.

The inherited timelike center pairs strictly negatively with r. As in
the old chamber, projection of that center onto any retained or new
wall is timelike and strict for every other wall. Thus all 141 facets
are actual. The new wall is strict at each of the four previously
certified proper cusps, so their bounded links persist.

There are no new finite vertices: a finite vertex in H20 would need a
positive-definite 20-root support, but r has only 17 finite neighbors.
Old finite vertices persist by acuteness and the retained common center,
preserving the certified lower bound of 170. This does not make the old
partial vertex census complete.

The independent remaining escape is

```
E2 = x - e/100 - (49/50)*f
     - (3/10)*omega17 - (2/5)*omega31
     - (1/100000)*sum(other16 support duals).
```

Its norm is exactly 427500061/5000000000>0, and its new-wall pairing is
exactly -19/1000. Every one of its 141 wall pairings is strictly negative,
as is its pairing with the timelike center. This proves infinite volume.

The [through-cusp obstruction](fresh_h20_null_link_escape_gate_20261016_review.md)
does not apply to r: its pairing with the persistent null direction nu
is strictly negative. The construction supplies a concrete cutter in
the class left open by that theorem. No complete ambient cone, new
finite-volume example, or new construction dimension is asserted.

```sh
python3 scripts/fresh_h20_off_cusp_cap_20261016_review.py --verify
```
