# Complete endpoint profiles for the balanced H20 cusp

The intrinsic affine A6 × A6 × A7 seed has **37 finite and four positive-semidefinite endpoint profiles**, with no label-five contact, H3, H4, or unbounded dihedral tail among the nineteen retained edge walls. The [exact census](fresh_h20_balanced_endpoint_20261014.json) transports these 41 profiles to all 392 original cusp edges, giving 14,504 finite and 1,568 PSD occurrences. This classifies endpoint diagrams for arbitrary real cap normals; it imposes no arithmetic root restriction.

Unlike the preceding A11 × A2 × E6 seed, these endpoint families can cut the half-vertex point. A cap through another face is therefore not required merely to cut that particular gap. Endpoint caps alone still miss the incenter, at both unit and golden-third-factor scales; the central cap being constructed separately removes that obstruction.

The [standalone producer](../scripts/fresh_h20_balanced_endpoint_20261014.py) uses the norm-two Gram `B=A6⊕A6⊕A7`, exact rational inverses, and the radical basis `1,√2,√3,√6,√5,√10,√15,√30`. All accepted records in fact need only √2. Its independent index construction agrees with the [shared intrinsic geometry](../scripts/fresh_h20_balanced_20261014.py): the affine factors have indices `0..6`, `7..13`, and `14..21`; the representative edge omits `0,7,14`.

## Completeness without an angle cutoff

Normalize a proposed cap to square two. Write its pairings with the retained nineteen norm-two walls as `−c_i`, where `c_i=2cos(π/m_i)` or zero. The endpoint Gram is positive semidefinite precisely when `cᵀB⁻¹c≤2`, and positive definite when the inequality is strict. Its inverse is entrywise nonnegative, with smallest diagonal `6/7`.

Every nonzero finite contact has `c_i≥1`, so three contacts cost at least `18/7>2`. Every label at least five costs at least `(6/7)phi²>2` even by itself. A parallel contact costs at least `24/7>2`, and a dotted pair cannot occur in a PSD Gram. Therefore **support at most two and labels three or four are exhaustive**, including ideal endpoints. All 723 assignments in that finite domain are checked; no numerical proposal filter is used.

The accepted profiles are:

| Profile | Finite | PSD |
|---|---:|---:|
| All right | 1 | 0 |
| One pi/3 contact | 18 | 1 |
| One pi/4 contact at a chain end | 6 | 0 |
| Pi/3 contacts at ends of two distinct factors | 12 | 0 |
| Pi/3 contacts at both ends of one factor | 0 | 3 |
| Total | 37 | 4 |

The exceptional single contact is the middle node of A7 (wall18). The other three PSD profiles join the two endpoints of an A6, A6, or A7 chain. These are PSD endpoint diagrams, not assertions that the resulting ideal links are bounded. Every finite profile has actual finite endpoint realizations when its cap radius is sufficiently small: its three omitted contacts then become dotted and the endpoint lies above the product vertex.

All affine marks are one. Independent cyclic rotations of the three factors give 7·7·8 explicit Gram-preserving transports, and the checker verifies every transported support. Thus the one profile list covers all original cusp edges.

## Radius equations and the incenter bound

Let the three inradii be `a,b,c`, let `z` be a product vertex, and solve `n_i·y=cos(π/m_i)` on its nineteen unit-normal incident facets. Then `p=|y|²≤1`, and the cap sphere with radius `R` and center `w=z−Ry` has its prescribed endpoint profile. The three omitted-wall ratios are exactly

```text
7a/R − B_A,    7b/R − B_B,    8c/R − B_C,
```

where each `B` is the sum of that factor's retained contact cosines. These functions are archived profile by profile. Their allowed Coxeter values are zero, a finite cosine, or at least one for parallel/dotted walls; labels on these three omitted walls are not bounded by the endpoint census.

This sphere orientation loses no new acute facet. In null coordinates write its last coordinate as `U`. Summing its nonpositive seed pairings with the positive affine marks gives `U≥0`. If `U=0`, each seed pairing vanishes, so its spatial part is orthogonal to a spanning rank-nineteen system and is zero; the normal would then be null. A new spacelike wall therefore has `U>0` and can be positively scaled to the exterior-of-sphere form used above.

The [separate radius checker](../scripts/fresh_h20_balanced_endpoint_20261014_bounds.py) and [certificate](fresh_h20_balanced_endpoint_20261014_bounds.json) use the universal necessary set `{0,1/2} ∪ [1/√2,∞)`. Define the continuous bound `R0=min A_f/(B_f+1/√2)`. Every admissible radius above `R0` must be one of the finitely many exact values `A_f/B_f` or `A_f/(B_f+1/2)`. The checker retains every such value passing all three necessary ratio gates, including equality and isolated points; it never treats this relaxation as sufficient Coxeter compatibility.

For each profile it checks the largest relaxed radius exactly. The incenter ball power is

```text
|z−Ry|²−R² = |z|² − 2R(z·y) + (p−1)R².
```

Here `z·y≥0` and `p≤1`, so the expression is nonincreasing for `R≥0`. It is at least68 at unit scales and strictly above75 at scales `(1,1,phi)`. The maximum relaxed radius is at most14. Factor rotations fix the incenter, so these 82 representative checks cover all 16,072 profile occurrences at each scale and every real radius.

In coordinates with `D·U=−1`, the point `(0,128,1)` is a common strict timelike center for these entire endpoint families. The point `(0,−1/2,1)` has norm one, future pairing `−255/2` with that center, and is strict on every vertical wall and every endpoint cap. Consequently no finite-volume completion using only these endpoint caps exists. A central or other interior cap is outside this statement; there is no obstruction to the full intrinsic seed.

## Concrete half-vertex cutters

At unit scales every product vertex satisfies `|z|²=196`. The all-right cap centered at `z` with radius `R=7/cos(π/7)` has two omitted pi/7 contacts and a dotted omitted A7 contact. The exact bound `7/8<cos(π/7)<1` gives `7<R<8`, so it strictly cuts `z/2`, whose distance from `z` is seven. The checker isolates the largest root of `8x³−4x²−4x+1`, with its Chebyshev identity, to certify these bounds.

At scales `(1,1,phi)`, the all-right radius `R=7√2` instead gives omitted labels `(4,4,dotted)`. The squared cut margin at `z/2` is `49−21phi>0`. These are genuine individual endpoint caps retaining all22 seed walls; their compatibility with a jointly selected cap bank is a separate task. No finite-volume or complete new-vertex census is claimed.

```sh
python scripts/fresh_h20_balanced_endpoint_20261014.py --verify
python scripts/fresh_h20_balanced_endpoint_20261014_bounds.py --verify
```

The separate [central-plus-endpoint corollary](fresh_h20_balanced_endpoint_20261014_central_gap.md) strengthens the incenter result: every endpoint cap is strictly dotted to the radius-two central cap, and an explicit nearby strict escape survives that central cap as well. Thus interior-face bridges remain necessary even after including the central cap.
