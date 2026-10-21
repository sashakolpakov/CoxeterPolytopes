# Global incidence restrictions for the small D4/D5 forests

Let P be a compact arithmetic hyperbolic Coxeter polytope of simplest type over Q(sqrt5). In each row below, call a vertex special when its entire finite diagram is the indicated forest. Write V_special and V_other for the numbers of special and remaining vertices. The following strict inequalities hold separately:

| Dimension | Special forest | Frequency restriction |
|---:|---|---|
| 12 | 3D4 | `3 V_special < 38 V_other` |
| 13 | D5+2D4 | `13 V_special < 264 V_other` |
| 14 | 2D5+D4 | `3 V_special < 77 V_other` |
| 17 | 3D5+A2 | `43 V_special < 448 V_other` |

In particular, none of these four forests can be the vertex type at every vertex. This addresses the specific small-component forests that obstructed the previous weight family. It is a restriction on their frequencies, not a dimension bound or an exclusion of these forests at individual vertices. The first row reproduces the previously independently reviewed D4 result.

The proof uses the [complete endpoint and triangle census](fresh_sqrt5_mixed_d4_d5_migration_20261017.md), the [complete quadrilateral finite patterns](fresh_sqrt5_mixed_d4_d5_quad_patterns_20261017.md), and the [arithmetic rational-factor gate](fresh_sqrt5_mixed_d4_d5_quad_factor_gate_20261017.md). These retain every same-type endpoint reattachment, including movement of a D5 component to another D4 and exchange with an A2 component. Opposite quadrilateral products are free; no chosen edge length, minimality bound or copied-attachment assumption is imposed on the whole census.

For a special vertex, let I be the set of root pairs that may support a triangle with exactly two special vertices, and J the set that may support an all-special quadrilateral. The exact lists satisfy I contained in J. Inside a D4 component, a pair of leaves is excluded from J by the [independently reviewed D4 arithmetic obstruction](fresh_sqrt5_d4_short_face_incidence_review_20261017.md). The new checker verifies that every affected literal attachment is exactly the copied leaf attachment with new–new product zero, so this application does not omit migrating alternatives.

| Special forest | Size of I | Size of J | Maximum special-vertex payment |
|---|---:|---:|---:|
| 3D4 | 9 | 9 | 9/2 |
| D5+2D4 | 11 | 24 | 35/4 |
| 2D5+D4 | 13 | 31 | 11 |
| 3D5+A2 | 16 | 37 | 53/4 |

There is no triangle with three special vertices. Every root pair in its full star occurs in a finite corner, so the entire norm-two Gram would be rational. Its physical negative direction contradicts arithmetic conjugate positivity, as in the rational-edge argument.

Now pay one unit at every nonspecial corner of a triangle or quadrilateral. A triangle with exactly two special vertices additionally receives one half at each special corner. An all-special quadrilateral receives one quarter at each of its four corners. Every triangle receives at least two units and every quadrilateral at least one.

Each root pair at a simple special vertex determines one incident two-face. It therefore cannot simultaneously incur the half-unit triangle charge and the quarter-unit quadrilateral charge. Since I is contained in J, the special-vertex payment is at most

`alpha = |I|/2 + (|J|−|I|)/4 = (|I|+|J|)/4`.

Every other vertex contributes at most `binom(d,2)`. If T and Q count triangles and quadrilaterals, this proves

`2T+Q <= alpha V_special + binom(d,2) V_other`.

The strict h-vector inequality, in the normalization checked in the [structural-floor package](fresh_sqrt5_discount_structural_floor_20261016.md), is

`[5 floor((d−1)^2/4)/2 − binom(d,2)] (V_special+V_other) < 2T+Q`.

Its coefficient is9,12,14,24 in the four dimensions. Subtracting the payment bound gives the displayed frequency inequalities.

The [source](../scripts/fresh_sqrt5_mixed_d4_d5_incidence_20261017.py) and [saved certificate](fresh_sqrt5_mixed_d4_d5_incidence_20261017.json) record all corner sets, the literal D4 maps removed, every corner charge, and the exact rational counting. Generation and a fresh saved replay pass:

```sh
python scripts/fresh_sqrt5_mixed_d4_d5_incidence_20261017.py
python scripts/fresh_sqrt5_mixed_d4_d5_incidence_20261017.py --verify
```

This is a new composition of the stated frozen classifications and arithmetic arguments, not an independent reimplementation of their graph engines. The three mixed inequalities require separate external review. The nonspecial vertices remain unrestricted, and these separate estimates do not yet control a union of several competing high-cost vertex types. No general compact dimension bound follows at this stage.
