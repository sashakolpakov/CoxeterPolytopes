# Short-face incidence at vertices consisting of D4 copies

Let P be a compact arithmetic Coxeter polytope of simplest type over Q(sqrt5), of dimension d=4k. Call a vertex special when its entire finite diagram is kD4. If V_D and V_other count its special and other vertices, respectively, the following scoped incidence argument gives

`(4k−9) V_D < (12k+2) V_other` for k≥3.

For example, in dimension12 it gives `3 V_D <38 V_other`. In particular, all vertices cannot be 3D4. This is a restriction on vertex frequencies, not a general dimension bound. It does not cover mixed D5/D4 vertex forests and does not justify a universal discount on D4 corners.

## Forced same-type edges

Deleting a leaf of D4 leaves A3; deleting its center leaves three A1 components. An endpoint whose full diagram is again kD4 must restore precisely that D4 with the same contacts. An added root cannot join an intact D4 component and still belong to a D4 component. Thus every edge between special vertices copies the finite attachments of its deleted root. The checker enumerates all two-valued reattachments on the retained three roots and finds exactly one for each port.

Use norm-two roots and write g for the mutual product of the old and new roots. Their common neighbor gives a triangular cyclic product equal to g, so g is an algebraic integer in Q(sqrt5). Conjugate positivity applies to the resulting local Gram matrix.

## Triangles

If such an edge lies on a triangle, its two opposite root normals occur together at the third finite vertex. Since g belongs to Q(sqrt5), its finite Coxeter label can only be2,3 or5: the other real quadratic cosine values are in Q(sqrt2) or Q(sqrt3). Exact physical and conjugate signatures exclude2 and3, and exclude5 when the replaced root is the D4 center.

For a leaf replacement, label5 is possible only when the other released root is the center. Otherwise that center remains in the third corner, producing a nonspherical3–3–5 triangle. The third finite diagram is then `(k−1)D4+2A1+I2(5)`.

Consequently there is no triangle with three special vertices. Every triangle with two special vertices uses a center–leaf corner at both of them. There are exactly3k such corners at each special vertex.

## Quadrilaterals

Suppose all four vertices are special. Each edge is then a forced copy as above. Eliminate the common spherical support, and write the two-root Schur matrix at one vertex as

`H=[[a,c],[c,b]]`.

For framing products g_i,g_j put λ_i=2−g_i and λ_j=2−g_j. The four-root residual determinant is exactly

`λ_i λ_j ((2a−λ_i)(2b−λ_j)−4c²)`.

The source verifies this polynomial identity over rational coefficients for all six pairs in D4 and all16 pairs between two D4 copies.

For ports in different copies, c=0. Each actual hyperbolic edge requires λ_i>2a and λ_j>2b, so the determinant cannot vanish. Such a quadrilateral would have too high a rank.

For two leaves in one copy, `a=b=4/3` and `c=−2/3`. Rank and the physical/conjugate signs would force

`A=−3g_i−2, B=−3g_j−2, AB=16`,

with A,B positive physically and negative at the conjugate embedding. Both are integral and congruent to1 modulo3. Thus `Norm(A)≡1 mod3` and `Norm(A)<0`.

But2 is inert in Z[phi]: its norm form `x²+xy−y²` vanishes modulo2 only when both coordinates vanish. Since Norm(A) divides Norm(16)=256, its negative possibilities are `−1,−4,−16,−64,−256`. They are all2 modulo3, a contradiction.

Hence an all-special quadrilateral also uses only center–leaf corners. This restriction is sharp locally: the checker supplies a compact arithmetic local quadrilateral with all four corners D4, both framing products `−2phi`, and physical/conjugate inertias `(4,1,1)` and `(5,0,1)`. It checks the complete four-ray face cone and actual facet witnesses. It similarly supplies a compact local triangle with corners D4,D4,2A1+I2(5). Neither local construction is claimed to have a compact ambient completion.

## Counting

Let T_2 count triangles with exactly two special vertices and Q_D quadrilaterals with four. The common corner restriction gives

`2T_2+4Q_D ≤3k V_D`.

Pay two units to each triangle and one to each quadrilateral. A triangle counted by T_2 receives one from its other vertex and one half from each special corner. Every other triangle has at least two other vertices and needs at most one unit per such corner. Any quadrilateral with an other vertex can be paid entirely at an other corner; an all-special quadrilateral needs only one quarter per special corner. Thus

`2T+Q ≤(3k/2)V_D+binom(d,2)V_other`.

The strict h-vector inequality gives

`(2k²−3k)(V_D+V_other)<2T+Q`.

Combining the inequalities and cancelling k/2 proves the stated frequency bound. The h-vector formula is the one checked in the [structural-floor package](fresh_sqrt5_discount_structural_floor_20261016.md).

The [source](../scripts/fresh_sqrt5_d4_short_face_incidence_20261017.py) and [saved exact output](fresh_sqrt5_d4_short_face_incidence_20261017.json) use rational Q(sqrt5) arithmetic, symbolic rational determinant coefficients, complete small reattachment tests, and explicit finite local face cones. Generation and a fresh saved replay pass:

```sh
python scripts/fresh_sqrt5_d4_short_face_incidence_20261017.py
python scripts/fresh_sqrt5_d4_short_face_incidence_20261017.py --verify
```

The norm/sign and incidence deductions above are mathematical arguments recorded by the package, rather than claims that finite numerical samples prove them. No minimality14 assumption or sampled framing-distance cutoff is used. The source checks its exact identities and finite reductions; independent mathematical review remains separate.
