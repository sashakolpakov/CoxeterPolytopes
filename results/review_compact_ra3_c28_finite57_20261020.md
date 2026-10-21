# Independent review of the C28 regular-pentagon exclusion

Status: **PASS**, 20 October 2026. This is an outside-author written
mathematical review; no scientific computation, worker import, Lean or TeX
execution was performed. The C26 packet and all reviewed sources were left
unchanged.

Reviewed source:
`results/compact_ra3_twenty_eight_vertices_no_regular_pentagon_20261020.md`,
SHA256 `076864a5d836456e56ac377dffbcaa581a4738ec977dbf58a7905580d9838afc`.
I read the full source and independently rederived its metric calculations,
graph-incidence arguments and stated RA4 consequences. I did not coauthor
the C28 candidate. Earlier reviews of some inherited geometric inputs are
separate from this review.

## Exhaustion and the three metric cases

For 28 vertices the cubic Euler identities give 16 faces. Polygon capacity
excludes faces of size above seven. In the heptagonal equality case, the
full equality-pattern argument in the pinned capacity manuscript identifies
L7, rather than merely producing a disjoint extra face. Its two bases are
regular by cyclic symmetry and metric uniqueness. I checked
`cosh L_k = 1 + 2 cos(2 pi/k)`: every lateral pentagon has an L7 edge,
strictly longer than the regular-pentagon edge. This excludes regular
pentagons in that case.

Otherwise Euler gives twelve pentagons and four hexagons. I checked the
primary published count: Table 2, row `(n_v,n_f)=(28,16)`, of Brinkmann,
Goedgebeur and McKay, [The Generation of Fullerenes, arXiv:1207.7010v2](https://arxiv.org/pdf/1207.7010v2)
(PDF page 15) has exactly two fullerenes. Both displayed dual complexes
have 16 vertices, 42 edges and 28 triangular faces, cyclic vertex links
and each edge incident with two triangles. Their degree-six induced
graphs differ, so they are the two distinct types. This uses the published
finite classification; it is not a new exhaustive graph enumeration.

For the tetrahedral type, the S4 action is a genuine combinatorial action
and hence acts isometrically by compact right-angled metric uniqueness.
The common fixed interior point and the stabilizer of an ordered pair
justify the two-parameter spatial ansatz. I independently recovered
`B^2=3/4`, `A=2B`, `a^2=7/4`, `h^2=3/4`, and the actual distinguished
edge value `cosh length=4/3 < phi`. Every pentagon has such an edge.

For the second type, I checked both displayed permutations on all 28
triangles, their commutation and the three pentagon orbits. The propagation
from either outer orbit uses two equal regular-pentagon edge lengths;
the nonadjacent-edge case follows from the stated right-angled pentagon
identity. It then forces five consecutive such edges on H0, contradicting
the actual-wall three-edge closure lemma.

For the remaining orbit I rederived the full normal calculation. In
particular, with `z=k^2`, elimination gives

`(33-20 phi) z^2 - (8-2 phi) z + phi^2 = 0`.

The roots are `phi^4` and `(6+7 phi)/29`; the latter is below one. Thus
`k=phi^2`, `e=-phi`. The four H normals span a Lorentz three-space,
and A, C and the distinct facet G have identical products with all four.
The remaining positive one-dimensional fiber permits at most two unit
normals, giving the stated contradiction. The inference that the displayed
linear relation is zero is valid because it is null and perpendicular to
a timelike vector. As an optional shorter observation, the same equations
already give `C.H1=0`, contrary to their nonadjacency via the actual
orthogonal-wall projection lemma; the longer argument needs no change.

## RA4 consequences and boundaries

A positive facet adjacent to a dodecahedral facet contains an entire
regular right-angled pentagon. The inherited exclusions at 22, 24 and 26
vertices, followed by this C28 result, therefore give at least 30 vertices
for such a positive neighbor. The inherited facet capacity yields
`f >= 2*30+5 = 65`. Existence of a dodecahedral facet and connected facet
adjacency supply such a neighbor whenever the polytope is not all
dodecahedral. With the inherited `W>=17` and exact volume identity,

`Vol = pi^2 (f+W+16)/12 >= 49 pi^2/6`.

The all-dodecahedral 120-cell separately satisfies this lower bound.
These statements do not assert that every positive facet has excess at
least five, nor prove the unrestricted 120-cell minimum or a second-volume
classification. I found no unresolved mathematical step within the stated
scope. Metric uniqueness, the fullerene count and the earlier capacity,
radius and excess theorems remain explicit imported inputs, not conclusions
of an executable or a formalized geometry check.

Critical inherited sources read and rehashed:

| Source | SHA256 |
| --- | --- |
| `manuscripts/compact_ra4_volume_constraints.tex` | `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f` |
| `results/compact_ra3_twenty_six_vertices_no_regular_pentagon_20261020.md` | `1bbe6adc630ff2f44424f6ef3c573b3c71dcef4f25496def9b13fd35b2e4a03f` |
| `results/compact_ra4_farthest_dodecahedron_wall_radius_20261020.md` | `d082e4a918d4504a59ac200f32bbb141f1f2cb52d701006eaea7ba5b0e5bcc22` |
