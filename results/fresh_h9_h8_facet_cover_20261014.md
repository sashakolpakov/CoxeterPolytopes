# The fixed H8 source cannot attach to the literal dotted cap

The compact Bugaenko H8 chamber with **11 facets and 41 finite vertices** cannot be placed perpendicular to the literal dotted cap `t=65323/65536` while retaining acute Coxeter relations to the twelve old mixed-D4 walls. This excludes every isometric placement of this source, with no cutoff on finite labels. Ordinary diameter and sampled geometric placement tests alone did not establish this; the decisive obstruction is the first nonzero Coxeter-contact gap.

The [exact checker](../scripts/fresh_h9_h8_facet_cover_20261014_gap_gate.py) and [certificate](fresh_h9_h8_facet_cover_20261014_gap_gate.json) use the frozen [H8 source](../data/known_examples/bugaenko_compact_h8.json), [mixed-D4 chamber](fresh_h9_mixed_d4_coupling_20261014.md), and [literal dotted cap](fresh_h9_mixed_d4_global_caps_20261014.md). The source really has dimension eight: its form is `<−2φ>⊕E8`, its normal list has eleven nine-dimensional vectors, and its complete compact vertex list has 41 entries. The JSON scalars `{"a":a,"b":b}` mean `a+bφ`, not `a+b√5`.

This is a restriction on a particular compact face source and perpendicular attachment. It does not exclude other H8 sources, nonperpendicular sidewalls or changed old walls. The separate [exact parameter-line theorem](fresh_h9_h8_facet_cover_20261014_line_exclusion.md) uses the uniform source lemma below to exclude every all-dotted cap on this particular H–E segment, not just the literal parameter.

## The all-label source lemma

Normalize the eleven H8 source facet normals `u_j` to square two. Suppose a future timelike vector `p` lies in its closed chamber cone and has

`h_j=−(p,u_j) ∈ {0} ∪ [1,∞)`.

Put `a=(√5−1)/4`. Then either

* `−p² ≥ 4a = √5−1`, or
* `p` lies on one of exactly three source vertex rays, numbered **0, 32, 33** in the frozen source.

On each exceptional ray, write `p=s v`, where the smallest nonzero contact of v is one. Then `−v²=a` and `s≥1`. Their nonzero contact triples are, in source-facet order:

| Source vertex | Nonzero contacts |
|---|---|
| 0 | `h_4=φ`, `h_9=h_10=1` |
| 32 | `h_0=φ`, `h_8=1`, `h_10=(7+3√5)/2` |
| 33 | `h_0=1`, `h_8=φ`, `h_9=(7+3√5)/2` |

Here is why a finite linear enumeration proves the lemma over all real p. Fix the set Z of zero contacts. The remaining contacts obey `h_j≥1`; this defines a polyhedron in the nine-dimensional Lorentz space. Its recession cone is contained in the source chamber cone. Since that source is compact, every nonzero recession vector is future timelike. The polyhedron is pointed and, when nonempty, is the convex hull of its vertices plus its recession cone.

On the future cone, `f(x)=√(−x²)` is homogeneous and satisfies the reverse triangle inequality. Consequently, for a convex combination of vertices plus a nonnegative combination of recession rays, f is at least the smallest vertex value. Thus its minimum is attained at a polyhedron vertex. Every such vertex has nine independent active equations, each a contact equal to zero or one.

The checker exhausts all **55** choices of nine source roots. Exactly **43** are independent; all **512** binary right-hand sides are tested for each. The remaining two contacts must also be zero or at least one. There are **3,666 feasible occurrences**, giving **3,657 distinct vectors** in **1,411 zero patterns**. Every vertex except the three displayed ray generators has `−p²≥√5−1`. The three exceptional zero patterns have eight independent zero equations, so each entire pattern is confined to its single vertex ray. This proves the stated all-real lemma, including zero contacts and boundary patterns. No grid of p or angle labels is used.

## Application to the literal cap

Let z be the saved raw dotted cap, `q=z²>0`, and n an old wall normal. Its projection to the proposed H8 facet is

`p=n−(n,z)z/q`, with `κ²=−p²=(n,z)²/q−n²`.

It is future timelike by the frozen dotted-cap certificate. A perpendicular source sidewall u has `(p,u)=(n,u)`. For an old norm-two wall, every acute Coxeter relation therefore gives exactly the source lemma's necessary condition `−(p,u)=0` or `≥1`: finite labels start with `2cos(π/3)=1`, while parallel and dotted relations are at least two. There is no upper-label assumption.

For old walls **0, 1, 2, 3**, the exact projected norms are strictly below `√5−1`; their numerical values are approximately `1.1911500, 1.2009636, 1.1882550, 1.1799527`. All four normalized projected centers are distinct. The certificate checks each of the six strict squared-distance margins

`(p_i,p_j)²−p_i² p_j² > 0`.

The source lemma would put these four distinct centers on just three source vertices. That is impossible.

There is also an independent one-wall contradiction. Old wall **7** has

`a < κ_7² < 2a`, approximately `κ_7²=0.5154439887`.

It must therefore lie on an exceptional ray with `1<s<√2`. One of its source contacts is exactly s, strictly between the consecutive allowed magnitudes `2cos(π/3)=1` and `2cos(π/4)=√2`. It cannot be right, any other finite label, parallel, or dotted. Thus wall 7 alone excludes the literal source attachment.

## Preliminary numerical diagnostics

These remain numerical necessary tests, superseded for the literal cap by the exact result above.

The [high-precision geometry helper](../scripts/fresh_h9_h8_facet_cover_20261014_geometry.py) evaluates exact input data at 80 decimal digits. The source vertex diameter is approximately `4.24510024762`, attained by vertices 32 and 33. The literal twelve-ball union has diameter approximately `4.47152651941`; the center diameter is `2.06987037163`. The basic size test therefore passes. The source diameter endpoints can only be assigned to ball 7 and one of balls 0–3 under the necessary pair-distance bound.

The [assignment diagnostic](../scripts/fresh_h9_h8_facet_cover_20261014_csp.py) and [saved results](fresh_h9_h8_facet_cover_20261014_csp.json) test the literal cap and 21 parameters in the normalized null-to-first-parallel interval. They impose both lower and upper possible distances between points in each pair of balls. The literal case has a full assignment after 41 search nodes. Twenty sampled parameters also pass this necessary relaxation; the lowest sampled point fails already by the diameter condition. An assignment is not an isometric placement, and no numerical UNSAT result is promoted to an exact certificate.

The [continuous placement scout](../scripts/fresh_h9_h8_facet_cover_20261014_placement_scout.py) and [output](fresh_h9_h8_facet_cover_20261014_placement_scout.json) use 40 deterministic local starts: eight each for the literal parameter and `u=0.85,0.90,0.95,0.98`. Each is limited to 250 function evaluations. The 36 isometry parameters are constrained by the 132 old-center/source-facet acuteness inequalities and 41 source-vertex/ball-union inequalities. No feasible relaxation was found. The best squared positive-violation residual is approximately `0.157923` at `u=0.85`; no discrete Coxeter pair checks or full base-cover acceptance followed. Local optimization failure is not a nonexistence proof. Even successful vertex coverage would not establish that the entire source polytope is covered by the ball union.

The old twelve finite vertices give additional necessary projected containment points. They do not individually give an uncovered-point obstruction: for an incident old wall n, the raw outside-ball expression at the projection of an old finite vertex v reduces to `(z,n)² v²<0`. Their projected mutual diameter is only approximately `0.05544` at the literal parameter. No broader projected-vertex hull exclusion is claimed.

Exact generation and saved replay pass. The numerical checks have their own saved replays and remain explicitly scoped as diagnostics. No walls were added and no old frozen artifact was changed.

```sh
python scripts/fresh_h9_h8_facet_cover_20261014_gap_gate.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h9_h8_facet_cover_20261014_csp.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h9_h8_facet_cover_20261014_placement_scout.py --verify
```
