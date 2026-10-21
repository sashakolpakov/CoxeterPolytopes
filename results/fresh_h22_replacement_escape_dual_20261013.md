# Exact restrictions on walls that can cut the H22 escape

The strict escape of the 813-facet chamber cannot be removed by varying the
height or angle in the previous three-dimensional replacement plane. More
generally, every acute new wall that cuts this escape must pair strictly
negatively with both old wall 7 and old wall 168. These are necessary
conditions for a successful addition, not a global nonexistence theorem.

The [rational checker](../scripts/fresh_h22_replacement_escape_dual_20261013.py)
and [certificate](fresh_h22_replacement_escape_dual_20261013.json) use the
frozen 813-facet source. Producer and saved replay pass. A numerical linear
program proposed the last two identities below; the certificate checks literal
rational vector equalities and signs, without an optimization dependency.

## The entire old replacement plane misses the escape

Use the source vectors p and delta, and write a completely general normal
in their extended plane as `r=(a*p+b*delta,z)`. There is no restriction on
its norm, height, field or angle labels. The pairings with three retained
walls and the saved escape E are the following linear functions:

| Vector | Coefficient of a | Coefficient of b | Coefficient of z |
| --- | ---: | ---: | ---: |
| Old wall 5 | -4 | -1 | 0 |
| Old wall 188 | -192 | -2 | 0 |
| Parent wall 229 | -423 | -8 | 2 |
| E | -128984 | -2346 | 80 |

The last row is exactly `20608/23` times the first, `12995/23` times the
second, plus 40 times the third. All coefficients are positive, and the
first three rows have determinant -368. Therefore every nonzero r acute
to these retained walls has `(r,E)<0`. In particular, no choice of the
mutual angle between the two symmetric caps can make this family cut E.

The same identity lifted to the full ambient space has a remainder in
the span of the 20 retained incident old walls. Fourteen coefficients
are negative and six are positive. Consequently every nonzero acute new
normal perpendicular to those fourteen negative-coefficient walls also
has strictly negative pairing with E. The certificate records all
coefficients and the full-rank argument for strictness.

## Two sharper necessary pairings

Two further exact decompositions are

```text
E = positive combination of retained normals - (161/2) W7,
E = positive combination of retained normals - 50 W168.
```

The first has 19 positive terms and the second 21; every coefficient and
normal index is included in the certificate. All terms belong to the
retained 811-parent set. If r is acute to those normals and `(r,E)>0`,
each identity separately forces `(r,W7)<0` and `(r,W168)<0`.

For a general wall these pairings may be finite-angle, parallel or dotted relations.
If either old wall belongs to an endpoint's spherical or affine support,
however, its negative pairing with r must be finite-angle or parallel:
positive semidefiniteness excludes a strictly dotted relation. For
a new norm-two normal the smallest permitted magnitude is
`sqrt(2*W_i²)/2`, corresponding to label three. This gives an exact lower
bound for the Schur complement in endpoint searches.

The [separate endpoint analysis](fresh_h22_replacement_edges_20261013.md)
applies this necessary condition to all 22 edges of the new finite vertex.

Any finite proposed completion whose additional walls all leave E feasible
still has infinite volume: combine E with a sufficiently small positive
multiple of a common strict timelike interior point to get a strict
spacelike point. This applies even when some added inequalities vanish at E.
The restrictions here do not exclude successive preparatory additions,
walls outside these orthogonality classes, or further old-wall deletions.

```sh
python scripts/fresh_h22_replacement_escape_dual_20261013.py --verify
```
