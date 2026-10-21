# A common-space constraint on the fixed uniform self-edges

The frozen [universal self-edge theorem](fresh_sqrt5_universal_self_edges_20261016.md)
already prevents exclusion of a rank21 finite vertex diagram merely by
asking whether each root port has an arithmetic neighboring vertex.
Every port has such local data. The additional requirement that the edges
live in **one quadratic space over Q(phi)** is stronger.

For the particular uniform parameter `lambda=phi^8`, this note proves a
common-space incompatibility in302 of the306 maximizing marked forests
of the selected weight proposal. It does **not** exclude these vertex
forests, arbitrary endpoint diagrams, varying edge parameters or compact
dimension21. The [standalone Fraction source](../scripts/fresh_sqrt5_vertex_common_space_20261016.py)
and [certificate](fresh_sqrt5_vertex_common_space_20261016.json) use no
producer arithmetic or numerical square tests. Generation and a fresh
saved replay pass.

The subsequent [common-lattice translation argument](fresh_sqrt5_common_lattice_self_edges_20261016.md)
shows that varying the endpoints supplies same-type local edges at every
port in one fixed admissible lattice. Accordingly this diagnostic does
not support a common-form vertex-exclusion strategy.

## The discriminant condition

Fix the old positive-definite finite Gram A and let `d_i=A_ii`,
`c_i=(A^-1)_ii`. The old uniform rule duplicates root i, retaining all
its other contacts and replacing its self contact by `d_i(1-lambda)`.
The full edge Gram M_i has Schur remainder

```
p_i = 2*d_i*lambda-d_i^2*lambda^2*c_i,
det(M_i)=det(A)*p_i.
```

Every M_i is nonsingular of the full ambient rank. If two such edge
configurations share one ambient quadratic form over F=Q(phi), their
Gram determinants differ by the square of a basis determinant over F.
Consequently **p_i/p_j must be a square in F**, irrespective of any
mutual product chosen between the two new roots. This requirement is
absent from an endpoint-type flow model.

The test here uses two ports within a single connected finite component.
Rescaling that entire component multiplies both p values by the same
factor. Individual F-valued root rescalings multiply them by squares.
The obstruction therefore persists under these changes of finite-root
normalization; it is not an artifact of the chosen integral symmetrization.
Arbitrary rescalings of different orthogonal components are not restricted
by this within-component test.

## An explicit D5 witness

In the saved norm2 D5 Gram, one port has `c_i=1`, another `c_j=5/4`.
Since `lambda=13+21phi`, its norm is1 and its trace is47. Thus

```
p_i=4lambda(1-lambda),       Norm(p_i)=-720,
p_j=4lambda(1-5lambda/4),    Norm(p_j)=-899.
```

The rational number `Norm(p_i/p_j)=720/899` is not a square. Hence
`p_i/p_j` is not a square in F. These two exact uniform edges cannot
coexist in a single F-quadratic space, even if the new-to-new product is
left unrestricted. The earlier real partial edge configurations remain
valid; their use of additional real square roots did not assert a common
F-realization.

## Finite diagnostic on the maximizing forests

The source independently inverts all68 saved finite component Grams and
reconstructs every p_i. It tests every within-component ratio exactly.
For `x=a+b*phi`, write `x=A+B*sqrt5`. A square root requires rational
`sqrt(Norm(x))`; the possible rational first-coordinate squares are
`(A±sqrt(Norm(x)))/2`. Both alternatives, including the zero-coordinate
case, are checked exactly. Every successful result includes its literal
field square root.

All component types except A1,A2 and I2(5) have an incompatible port pair
for this fixed lambda. Applying the result to the
[306 exact maximizing witnesses](fresh_sqrt5_h3_center_selected_diagnostic_20261016.json)
detects302 conflicts. The four undetected marked forests are

```
5A1+8I2(5),   7A1+7I2(5),   9A1+6I2(5),   11A1+5I2(5).
```

No common-space existence is claimed for those four: this particular
within-component invariant simply does not obstruct them. Conversely,
none of the other302 vertex diagrams is excluded. Choosing a different
parameter at each port changes p_i, and different endpoint diagrams need
not use the duplicate-root construction at all.

A single shared squareclass remains necessary for specific full-rank
edge Grams, but it can be satisfied by the variable endpoint construction
in the linked note. The present certificate covers only the fixed uniform
parameter; it is not an all-parameter obstruction or a new dimension
theorem. Simultaneous compact two-face closure remains a different issue.

```sh
python scripts/fresh_sqrt5_vertex_common_space_20261016.py --verify
```
