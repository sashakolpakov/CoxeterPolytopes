# Two exceptional-corner anchors are impossible

For **every all-dotted cap direction** over the frozen twelve-wall mixed-D4 chamber, a coherent perpendicular placement of the fixed Bugaenko H8 source can have **at most one** projected old center below the source norm threshold `√5−1`. This is a necessary condition for a fresh cap direction, not an exclusion of all attachments. Single-corner anchors and configurations above the threshold remain open.

The [exact checker](../scripts/fresh_h9_h8_corner_anchors_20261014_pair_gate.py) and [certificate](fresh_h9_h8_corner_anchors_20261014_pair_gate.json) test the old Gram equations before constructing a new cap. Generation and saved replay pass. All finite angle labels are included; there is no label cutoff or numerical candidate acceptance.

Normalize the cap r and every old normal n to square two. If the unnormalized old normal is N and the raw cap is z, the relevant projected norm is

`κ²=2(z,N)²/(z² N²)−2`.

This normalization matters for old polar wall 11, whose stored square is not two. The [source gap theorem](fresh_h9_h8_facet_cover_20261014.md) forces every `κ²<√5−1` center onto source vertex 0, 32 or 33. Put `a=(√5−1)/4` and normalize their future ray generators v so that their smallest nonzero source contact is one. Their Gram is

```
             0                 32                33
0           −a                −b_A              −b_A
32          −b_A              −a                −b_B
33          −b_A              −b_B              −a

b_A=(11+5√5)/4,  b_B=(23+9√5)/4.
```

The checker reconstructs this matrix by an exact SymPy solve from the literal source Gram and the frozen contact vectors. Each low-norm projection is `p=s v`, with

`s ∈ {1,√2} ∪ [φ,2)`.

This superset includes every finite Coxeter contact magnitude `2cos(π/m)`, `m≥3`, below two. The upper limit follows from `κ²=a s²<4a`.

Two old normals with projections `s v_A`, `t v_B` have product

`(n_i,n_j)=−b st+√((2+a s²)(2+a t²))`.

The positive square-root sign is forced: both old normals have negative pairing with the all-dotted cap. Two projections cannot lie on the same ray, since then their old product is at least two, contradicting old acuteness.

For distinct rays, define

`F_b(s,t)=b st−√((2+a s²)(2+a t²))`.

This is strictly increasing in either variable for `s,t≥1`. For example,

`∂F/∂s > bt−√a√(2+a t²) ≥ t(b−√(a²+2a)) > 0`.

Both exact source constants satisfy `b²>a²+2a`. In particular, every distinct-ray anchor pair has normalized negative old product at least `F_A(1,1)=2φ`.

Of all 66 old pairs, **64** have a strictly smaller normalized negative product, so they cannot furnish two low-norm anchors. The two remaining pairs are `(4,11)` and `(8,11)`, with magnitudes approximately `8.81060214` and `4.77259356`. Their complete all-label exclusion uses only the first two contact gaps:

* For pair `(8,11)` and type A, its magnitude lies strictly between `F_A(1,1)` and `F_A(1,√2)`. Both scales equal to one are too small; increasing either scale to the next allowed value is already too large. Type B fails even at its minimum.
* For pair `(4,11)` and type A, its magnitude exceeds `F_A(1,2)`, so neither scale can be one. It lies strictly between `F_A(√2,√2)` and `F_A(√2,φ)`. Thus both scales would have to be √2, but their product is too small. Increasing either to any later finite label is too large.
* For pair `(4,11)` and type B, its magnitude lies strictly between `F_B(1,1)` and `F_B(1,√2)`, giving the first-gap contradiction again.

All eight strict comparisons use separate 192-bit rational outward bounds on the physical quartic coefficients and normalization radicals. The polar normal is never treated as a norm-two stored vector. The checker also reproduces all 64 ordinary pair bounds.

This stops the proposed two-exceptional-corner anchoring branch before any expensive placement. It leaves a concrete smaller search: choose at most one exceptional projected center, enforce its discrete ray scale, and require every other normalized projected norm to be at least `√5−1`. A source placement must still satisfy all remaining old/new Coxeter relations and the full ball-cover condition. No extra walls or finite-volume example are claimed here.

```sh
python scripts/fresh_h9_h8_corner_anchors_20261014_pair_gate.py --verify
```
