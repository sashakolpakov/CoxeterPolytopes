# Three finite cube types and three infinite product families in `H9`

## Outcome

Four of the five parked rank-ten Alexeev sections have enough doubling-wall
geometry to produce substantially more than four isolated chambers.  The
exact result is:

1. `A14_10_1_C01` has three mutually orthogonal doubling walls.  Its seven
   nonempty marked subsets collapse to exactly three labelled diagram types,
   with 14, 15, and 16 facets and volume multipliers 2, 4, and 8.
2. `A14_8_1_C02` has a parallel pair and one orthogonal doubling wall.  For
   every `N>=1` it gives a finite-volume noncompact polytope with
   `4N+13` facets and volume multiplier `2N`.
3. `A14_10_1_C02` has a parallel pair and two orthogonal doubling walls.  For
   every `N>=1` it gives one with `3N+14` facets and volume multiplier `4N`.
4. `A16_14_1_C01` has two parallel pairs, with every wall in one pair
   orthogonal to both walls in the other.  For every `m,n>=1`, the finite
   rectangular union has

   \[
   \operatorname{facets} P_{m,n}=mn+m+n+12,
   \qquad
   [W(P):W(P_{m,n})]=mn,
   \qquad
   \operatorname{vol}P_{m,n}=mn\operatorname{vol}P.
   \]

These are three infinite families, not three long finite searches.  Taking
`n=1` in the rectangular family already gives pairwise distinct volumes for
all `m`.  All walls remain primitive integral reflection divisors on the same
rank-ten K3 carrier lattices.

The frozen construction is
[`certify_alexeev_parked_h9_doubling_products.py`](../scripts/certify_alexeev_parked_h9_doubling_products.py),
with certificate
[`alexeev_parked_h9_doubling_products.json`](../data/enumerations/alexeev_parked_h9_doubling_products.json).
Its independent all-parameter companion is the
[`hostile product audit`](alexeev_parked_h9_doubling_products_hostile_audit.md),
with a separately reconstructed
[`certificate`](../data/enumerations/alexeev_parked_h9_doubling_products_hostile_audit.json).

## Complete finite cube

The three one-based doubling walls are `3,5,11`.  They are pairwise
orthogonal, their integral reflections commute, and all seven nonempty
subsets have been constructed with strict facet witnesses.  Exact labelled
isomorphism gives:

| copies | facets | marked subsets | relation histogram | finite + ideal + hyperideal |
|---:|---:|---|---|---:|
| 2 | 14 | each singleton | `73m2+12m3+4m4+2 parallel` | `94+8+0` |
| 4 | 15 | each two-wall subset | `83m2+16m3+2m4+4 parallel` | `121+10+0` |
| 8 | 16 | `3,5,11` | `93m2+21m3+6 parallel` | `148+14+0` |

No non-power-of-two copy count occurs in this cube: a selected set of `k`
walls glues exactly `2^k` seed chambers.

## Why the rectangular formula is exact

In the oriented seed, walls `(1,14)` and `(8,15)` are the two parallel
pairs.  Each normal has square `-2`; the inner product within each pair is
`2`; all four cross-pair products are zero.  More strongly,

```text
alpha_1 + alpha_14 = alpha_8 + alpha_15
```

is the same primitive isotropic cusp direction.  The four reflections
therefore generate the rectangular affine factor
`D_infinity x D_infinity`.  All finite rectangles are obtained by applying
the ordinary finite-path construction successively in the two commuting
directions.

The remaining eleven seed walls split exactly as follows:

- seven walls, `3,4,5,6,7,10,11`, are orthogonal to all four doubling walls
  and stay fixed;
- four walls, `2,9,12,13`, are the corners of one fundamental rectangle.

An `m` by `n` union therefore has four outer endpoint mirrors, seven fixed
mirrors, and one corner mirror at each vertex of an `(m+1)` by `(n+1)` grid:

\[
4+7+(m+1)(n+1)=mn+m+n+12.
\]

The certificate replays this formula exactly for every `1<=m,n<=5` and the
two path-times-cube formulas for `1<=N<=8`.  The all-parameter theorem is the
displayed reflection-orbit argument; it is not inferred by polynomial
interpolation from those finite checks.  The hostile companion closes this
point independently: it proves the infinite-dihedral and direct-product
normal forms are faithful, solves every wall-stabilizer and collision equation
for arbitrary integer translation exponents, and derives all three component
counts symbolically.

The doubling and finite-path mechanism is classical: see Daniel Allcock,
[*Infinitely many hyperbolic Coxeter groups through dimension
19*](https://doi.org/10.2140/gt.2006.10.737), especially the doubling lemma
and finite-path/index consequences.  The project-level addition is the exact
Alexeev/K3 seed extraction, its product decomposition, the facet formulas,
and the archived cone/Vinberg replays.

## Direct finite-volume checks

Six representatives were rebuilt independently from their boundary roots.
GMP rational double description, the time sheet, and Vinberg's literal
elliptic-edge and complete parabolic-neighbourhood tests agree:

| ID | copies | facets | finite + ideal + hyperideal | elliptic edges | affine components |
|---|---:|---:|---:|---:|---:|
| `A14_10_1_C01_D2_W3` | 2 | 14 | `94+8+0` | 515 | 17 |
| `A14_10_1_C01_D4_W3_5` | 4 | 15 | `121+10+0` | 661 | 22 |
| `A14_10_1_C01_D8_W3_5_11` | 8 | 16 | `148+14+0` | 828 | 31 |
| `A14_8_1_C02_P2_X_C2` | 4 | 21 | `334+60+0` | 2,407 | 137 |
| `A14_10_1_C02_P2_X_C2SQ` | 8 | 20 | `312+49+0` | 2,096 | 109 |
| `A16_14_1_C01_RECT_2X2` | 4 | 20 | `210+24+0` | 1,294 | 58 |

Every row has zero hyperideal rays, all extreme rays in one time cone, zero
bad edges, and zero bad parabolic neighbourhoods.  General finite volume is
also immediate from convex gluing: the interiors of the seed chambers are
disjoint and the displayed index is the number of copies.  Noncompactness
persists under the finite-index passage from the nonuniform seed group.

## Algebraic and novelty scope

Each family inherits the seed's weak full-Picard K3 divisor carrier.  This is
a genuine `surface carrier => infinitely many polytopes` theorem.  It does
not say that the enlarged chambers are nef cones, that their reflections are
surface automorphisms, or that the carrier surface changes with `N,m,n`.

The six directly audited diagrams are pairwise nonisomorphic and do not
match the three prior factory diagrams, the nine earlier doubling-cube
types, or the five parked seed diagrams.  This proves that the three finite
cube rows are new to the project ledger and supplies new project
representatives in all three infinite families.  It is not a worldwide
priority theorem.  A separate hostile literature/corpus audit records the
bounded lineage denominator and finds zero exact encoded labelled-graph
matches in the pinned CoxIter snapshots.

## Replay

```bash
PYTHONPATH=src:. python3 scripts/certify_alexeev_parked_h9_doubling_products.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_alexeev_parked_h9_doubling_products.py
PYTHONPATH=.:src python3 -u scripts/audit_alexeev_parked_h9_doubling_products_hostile_audit.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_alexeev_parked_h9_doubling_products_hostile_audit.py
```

The JSON archive SHA-256 is
`10e6323895bacd2d18f897f33ea5f9d4bb1ee4eedad061d545aa5302b447c694`;
its internal payload digest is
`3aebc89f3e91e1dd2f8f95569bc0d5ca0d01b81a5b76b95820e837b63ae9c0ea`.
The hostile archive SHA-256 is
`29366fe58e52aa2492693704134b78967e4340ca972e5f56d7d7f69972d646bd`.
