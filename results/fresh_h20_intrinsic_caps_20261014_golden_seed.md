# An intrinsic golden H20 cusp chamber

The [exact certificate](fresh_h20_intrinsic_caps_20261014_golden_seed.json) constructs **275 actual Coxeter facets in H20**, with 252 certified finite vertices that terminate every original cusp-edge germ. It still has infinite volume: a saved strict spacelike point satisfies all 275 inequalities. This is a new intrinsic metric construction, not a section of the old 210-wall chamber or a rational-ambient cap extension of it.

The [stdlib Fraction checker](../scripts/fresh_h20_intrinsic_caps_20261014_golden_seed.py) works in `Q(phi)`, where `phi²=phi+1` and `phi>1`. Its only geometric input is the frozen [affine endpoint diagram](fresh_h20_intrinsic_caps_20261014_endpoint.md), whose accepted records and transports have a [separate arithmetic review](fresh_h20_intrinsic_caps_20261014_review.md). The inherited diagram is affine `A11 × A2 × E6`; the last factor is E6, not D6. The construction chooses inradii `1, phi, 1`, so the 22 seed walls bound a compact Euclidean product alcove at the common ideal vertex.

## Exact construction

Use 19 spatial coordinates with a rational positive-definite unit-root Gram `B`, and two null coordinates `D,U` with `D.U=-1`. A vector is written `(z,d,u)` and has square `z^T B z-2du`. The resulting signature is `(20,1)`.

The seed normals are `(v_i,b_i,0)`, of norm1, with `b_i=phi` on the triangle and `b_i=1` otherwise. They have only right and label-3 mutual angles. Their product alcove has exactly `12*3*7=252` vertices `z_v`, computed by rational linear equations with the triangle coordinates scaled by phi.

Add the central hemisphere and a radius-one hemisphere at every alcove vertex:

```text
central cap:      (0,-2,1),                  norm4;
vertex cap v:    (z_v,(z_v²-1)/2,1),         norm1.
```

The central cap meets the A11 and E6 walls at angle pi/3 and the triangle walls at pi/5. Every vertex cap is right to its nineteen incident seed walls and dotted to the three omitted walls, the central cap, and every other vertex cap. The exact full pair census is:

| Relation | Pairs |
|---|---:|
| Right | 4,998 |
| pi/3 | 40 |
| pi/5 | 3 |
| Dotted | 32,634 |
| Total | 37,675 |

The minimum squared norm of an alcove vertex is `302+4phi`; the minimum squared separation of distinct vertices is `12+12phi`. These are exact minima, not numerical estimates.

The common strict timelike point `(0,9/2,1)` has square−9. Each original cusp edge ends at

```text
X_v=(z_v,(z_v²+1)/2,1),   X_v²=-1.
```

The checker verifies all 275 inequalities at all 252 points, exactly nineteen old zeros and one cap zero, and positive definiteness of each rank20 incident Gram. The incident index sets are distinct. Every old facet occurs at one of these vertices; the central facet has the separate strict point `(0,2,1)`. Thus all 275 inequalities define actual facets. The 252 vertices are a certified subset, not a complete global vertex census.

## Global volume and arithmetic scope

The central cap cuts the old incenter escape `(0,-1/2,1)`, with pairing `5/2`. It does not close the chamber globally. For the first alcove vertex, put `z=z_v/2`; then

```text
E=(z,(z²-1)/2,1),   E²=1.
```

All 275 pairings are strictly negative, and its pairing with the strict center is `-167/4-phi/2`. This is an exact infinite-volume certificate. Closing all old cusp edges therefore does not establish global closure; the [power-cell covering test](fresh_h20_intrinsic_caps_20261014_review.md) explains the missing condition.

The archive also records an irrational Cartan two-cycle `phi²` and a rational indefinite principal block (the twelve A11 walls together with the central cap). Their restriction to the sum of the A11 roots and the central cap is `[[0,-12],[-12,4]]`. Under Vinberg's necessary conjugate-positivity criterion, any cofinite Coxeter completion retaining these facets would consequently be non-quasi-arithmetic: an embedding moving phi extends to a totally real Gram field if one exists, but leaves the rational indefinite block unchanged. This is conditional; it does not assert that such a completion exists. The criterion is stated in [Bogachev–Yorov, Theorem 2.2](https://nvbogachev.github.io/Quasi_arithmetic_Coxeter_prisms-2024.pdf).

The [separate half-vertex extension](fresh_h20_intrinsic_caps_20261014_golden_extension.md) tests a coherent second layer without modifying this certificate.

```sh
python3 scripts/fresh_h20_intrinsic_caps_20261014_golden_seed.py --verify
```

Generation and saved replay pass. The [independent geometry review](fresh_h20_intrinsic_caps_20261014_golden_review.md) checks the formulas and actuality arguments; the [separate matrix arithmetic audit](fresh_h20_intrinsic_caps_20261014_matrix_review.md) also verifies the larger 587-wall checkpoint. The construction does not provide a finite-volume H20 example.
