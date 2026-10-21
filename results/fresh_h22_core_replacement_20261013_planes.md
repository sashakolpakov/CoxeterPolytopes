# Two earlier changed-core planes, and why the cusp-based replacement differs

These exact experiments preceded the [813-facet one-wall replacement](fresh_h22_core_replacement_20261013.md). They show that deleting old facets really can permit spacelike old projections, but two particularly symmetric projection planes still cannot cut an inherited old cusp. The obstructions are confined to those planes, not to general facet replacement.

## Five-short-wall deletion

Choose old short wall `w=W0` and delete its five short neighbors, indices `21,22,23,48,49`. Retain the other 205 old walls and all 602 signed parent caps. In the form `diag(Q,2)`, insert

```
t=(sqrt(31)-1)/2,
r±=(h-t*w, ±1/sqrt(2)).
```

The old projection has square 1, the caps have square 2, and their mutual pairing is zero. Their complete new relations comprise **72 π/3 pairs, 1 right pair and 1,542 dotted pairs**. A strict common center proves that this is an exact 809-facet Coxeter chamber. An irrational cyclic product `56-9sqrt(31)` and a retained rational short-root block with negative all-ones value `−86` supply the same conditional non-quasi-arithmetic obstruction described in the main replacement report.

Nevertheless neither cap cuts any old cusp. The complete old cusp census gives

```
0 <= -(w,δ)/(-(h,δ)) <= 1/3
```

for every short wall and every old cusp. The checker verifies this on all 42 short walls and the 18 representative cusps. The previously independently verified symmetry/corner census supplies coverage of their full orbits. Since the chosen `t<3`, the new old projection is strictly negative on every old cusp.

There is an all-real obstruction in this entire projection plane. Consider a new normal `(s*h-t*w,a)` acute and Coxeter-compatible with the retained walls. A short wall orthogonal to `w` forces `s>=0`. If `s=0`, retaining `w` and a long wall paired `−2` with it forces `t=0`; the opposite parent heights then exclude a nonzero pure-vertical normal. Otherwise scale to `s=1`. The pairing `−1` with a retained short wall implies norm `q<=2`: every nonzero finite Coxeter cosine has absolute value at least `1/2`, and parallel/dotted relations give an even stronger bound. Hence

```
2t²+2t-14 <= q <= 2,
```

so `t<3`. Together with the complete cusp ratio bound this rules out a positive pairing with any old cusp. This includes every finite label and every dotted or parallel relation, not just a chosen label alphabet.

## Deleting one whole incidence part

Let the incidence-center vectors satisfy

```
v²=vprime²=-2,  v·vprime=-5,  h=v+vprime.
```

Delete all 21 point walls but retain the 21 line walls, all 168 long walls and all 602 parent caps. A potential normal with old projection in `span(v,vprime)` either cannot cut an old cusp by its signs, or can be positively rescaled to

```
u=-s*v+vprime.
```

Retained long walls force `s<=1`. The complete cusp census gives ratios between `1/2` and `2` for the negative pairings with `v` and `vprime`, including the possible interchange of the two incidence parts. Therefore cutting any old cusp requires `s>1/2`.

Every retained line has pairing `−1`, so again `q<=2`. But on `1/2<s<=1`,

```
u²=-2s²+10s-2 >= 5/2 > 2,
```

which is impossible. Thus this second three-dimensional plane, including the added vertical direction and arbitrary real coefficients, cannot supply an acute Coxeter wall cutting an inherited old cusp after the specified deletion.

The source checker verifies that the 42 short walls span the old 22-dimensional space and their incidence graph is connected and bipartite. Their pairing patterns uniquely determine `v,vprime`; every source-diagram symmetry fixes or interchanges the two. This justifies reducing the second ratio check to the 18 cusp representatives without assuming a larger continuous symmetry.

## Reproduction and scope

- [Exact checker](../scripts/fresh_h22_core_replacement_20261013_planes.py)
- [Frozen output](fresh_h22_core_replacement_20261013_planes.json)

```
python scripts/fresh_h22_core_replacement_20261013_planes.py --verify
```

The producer and fresh replay pass. Rational input reconstruction is shared with the standalone 813 checker; the full old cusp coverage is a hash-bound dependency on its independently replayed closed corner census. No new global cusp enumeration is claimed here. The later successful cusp-cutting pair uses the projection perpendicular to 20 walls at a one-component affine cusp and lies outside these two tested planes.
