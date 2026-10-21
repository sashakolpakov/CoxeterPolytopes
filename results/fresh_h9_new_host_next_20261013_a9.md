# A changed A9 host and its first finite edge endpoint

**Superseded status.** This historical checkpoint is superseded by the [complete edge-walk chamber-stabilizer certificate](fresh_h9_a9_outer_20261013.md): the specified quartic A9 host is now proved nonreflective, with an explicit certified PARI dependency for first-edge completeness. Its earlier open-host statements describe the status at this checkpoint.

Subsequent status: the [next complete walk](fresh_h9_a9_next_20261013.md)
resolves all eight central directions described below. Its 19-wall prefix
has 36 further unresolved reversal orbits and still has a strict escape.
This report records the earlier first-step certificate.

A genuinely changed arithmetic host now has an **11-wall exact H9 Coxeter cone and a certified first finite edge step** to an `A4+A4+A1` corner. The cone still has an exact strict spacelike escape; eight outgoing directions remain open against its existing walls. This is progress toward a new compact example, not a finite-volume polytope or a reflectivity theorem.

The first-step completeness uses a **certified PARI absolute-norm and fundamental-unit interface**, followed by a separate exact Fraction replay. It is not represented as a CAS-free coefficient census. The [preceding F4 host](fresh_h9_new_host_next_20261013.md) was proved globally nonreflective and is not being searched further.

## A different quadratic space

Use `K=Q(T)`, with `T=2cos(pi/8)` and `T⁴−4T²+2=0`, and

```
alpha = T−1,
L = <-alpha> + A9.
```

The A9 Gram has diagonal 2 and adjacent entries −1. The form has signature `(9,1)` at the distinguished embedding and is positive definite at the other three. Its determinant is `−10 alpha`. The prime 5 is inert in K, as checked directly modulo 5. Its odd determinant valuation distinguishes this space from both compared `F4²+<1>` spaces, whose determinants have no factor at 5. Scalar similarity in even rank ten preserves determinant squareclass. This is a comparison to these searched hosts, not a worldwide novelty assertion.

The dual exponent is 10. The previously certified PID and unit structure therefore gives a safe complete primitive norm pool dividing 20:

```
{1, pi, T², 2/pi, 2, 2pi, 2T², 4/pi, 4} × {1,5},  pi=2+T.
```

Some of these eighteen classes need not occur. All are included in the searches and first-endpoint interface.

The complete zero-height reflection chamber is A9. The integral isometry group of the totally definite positive block is finite, because matrices are bounded at every embedding. Its reflection subgroup contains irreducible `W(A9)`, so it is itself irreducible. The rank-nine finite real reflection classification leaves A9, B9, or D9; the order `10!` divides neither `2⁹·9!` nor `2⁸·9!`. Thus no extra zero-height reflections enlarge A9. The frontier agent independently reviewed this argument and the norm pool.

## Initial census and a constructive cap

The complete producer census through height `alpha*x²/q ≤ 650` has 80 time controls and one positive-height wall, giving a 10-wall prefix. Reflection divisibility is imposed before the dominant spatial recursion. A separate Fraction implementation repeats all spatial enumeration through height 150: 46 controls and the unique first positive root at height

```
−16−8T+26T²+14T³ ≈ 146.3084969.
```

It also replays every saved H650 root, all 45 pairs, the strict interior, and the exact escape. This is **not** an independent repetition of all spatial controls through 650. The initial H30 census is repeated completely as well. See the [census audit](fresh_h9_new_host_next_20261013_a9_independent.json).

The requested H1000 run was stopped after complete height groups beyond 650; its archive is explicitly marked incomplete. The [H650 checkpoint](fresh_h9_new_host_next_20261013_a9_h650.json) regenerates and matches all time controls through 650 before extracting the completed portion. No live or completed-H1000 claim is left behind.

For the middle outgoing edge, put

```
D = (1,2,3,4,5,4,3,2,1),   (D,D)=10.
```

An integral solution of `z²−10 alpha x²=20` gives the primitive norm-two root `r=(x,(z/10)D)`. The chosen solution is

```
x = −5160636−2792914T+8809758T²+4767806T³,
z = −15025870−8131940T+25650760T²+13882090T³.
```

Both are positive at the distinguished embedding. This cap is compatible with the entire 10-wall prefix. The point

```
v = (1,(alpha*x/z)D),     (v,v) = −20 alpha/z² < 0
```

lies on eight initial walls and the cap, satisfies every remaining inequality, and has positive-definite rank-nine incident Gram of type `A4+A4+A1`. The [endpoint artifact](fresh_h9_new_host_next_20261013_a9_endpoint.json) records all eleven normals, the actual finite vertex, a strict escape, and every outgoing edge's existing-wall lookahead. Its [separate replay](fresh_h9_new_host_next_20261013_a9_endpoint_independent.json) checks all 55 pairs, all eleven root/integrality conditions, the finite vertex, and all nine edge directions.

## A complete first-mirror interface

There are fifty possible finite Coxeter coupling profiles from a normalized primitive root to the eight walls of the initial middle `A4+A4` edge. The census uses the entire norm pool, every finite angle whose normalized square lies in K, reflection integrality, and positivity at all embeddings. All profiles reduce to

```
z²−10 alpha x² = 10R,
```

where R is the residual norm after projection to the edge walls. There are twenty-four target values.

The common degree-eight extension has polynomial

```
X⁸+40X⁶+200X⁴−4000X²−10000,
```

signature `(2,3)`, and class number two. The archived PARI run uses `bnfcertify=1`, a full fundamental-unit basis, and `bnfisintnorm` for all twenty-four absolute target norms. The latter supplies complete integral solutions modulo units of positive absolute norm; its completeness is an explicit software dependency. No class-number-one assumption or merely particular `rnfisnorm` solution is used.

The [official PARI documentation](https://pari.math.u-bordeaux.fr/dochtml/ref-stable/General_number_fields.html#bnfisintnorm) specifies this completeness and its dependence on certified number-field data. The archived script calls `bnfcertify` without the optional flag that would give only partial certification. This software contract was checked against the primary documentation on 2026-10-13.

The first three fundamental units lie in K, the first is alpha, and the fourth has relative norm `−alpha`. Consequently the full relative norm-one unit group is generated by a sign and `epsilon=u4²/alpha`. Every surviving absolute representative is normalized by a checked square unit in K. Exact predecessor/successor inequalities then exhaust its unit orbit inside the finite segment bounded by the displayed cap.

Only five orbit points occur in that interval: two are nonintegral, two are imprimitive, and **the sole primitive reflective mirror is the displayed norm-two cap**. No prefix-compatibility rejection is needed for this uniqueness. The [independent first-endpoint certificate](fresh_h9_new_host_next_20261013_a9_first_independent.json) rebuilds all fifty profiles in reverse coordinate order, uses separate Fraction arithmetic, replays the complete relative unit intervals, and checks the initial/edge/endpoint reflection counts `45 / 20 / 21`.

Why this proves firstness: the initial full A9 corner and its middle edge are complete. The explicit cap bounds a compact timelike segment. Any earlier first mirror of the full arrangement would give a finite vertex. Its local spherical chamber is simplicial, so its new simple facet has finite Coxeter angles to all eight existing edge facets and belongs to the profile census. Its integral norm solution must occur in one of the complete orbit intervals. The uniqueness above rules out an earlier mirror and identifies the full endpoint stabilizer as `A4+A4+A1`.

The frontier agent independently reviewed this first-contact reduction. The root agent reran the separate first-endpoint checker successfully. These reviews do not replace the stated certified PARI interface with a second number-theory implementation.

## Live frontier and the symmetry check

The cap closes the middle initial edge. At its new corner, releasing the cap returns to the initial A9 vertex. Releasing any of the other eight incident walls reaches the light cone before an existing wall. These are **open search directions**, not proofs that the full chamber escapes there. Their exact starts, directions, and current wall intersections are in the endpoint artifact.

The short symmetry mechanism that excluded both F4 hosts does not repeat at these two corners. The [exact diagram-symmetry census](fresh_h9_new_host_next_20261013_a9_symmetry.json) tests all eight automorphisms of `A4+A4+A1`. Each fixes its unique A1 cap and the endpoint, hence the initial timelike axis. Only the identity and the initial A9 reversal are integral; the other six have denominator 5. The tested integral group therefore has order two. This is a scoped symmetry test, not evidence excluding a later infinite stabilizer.

## Reproduction

```sh
PYTHONPATH=src:. python scripts/fresh_h9_new_host_next_20261013_a9_first_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_new_host_next_20261013_a9_endpoint_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_new_host_next_20261013_a9_symmetry.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_new_host_next_20261013_a9_independent.py \
  results/fresh_h9_new_host_next_20261013_a9_h30.json \
  results/fresh_h9_new_host_next_20261013_a9_h650.json \
  --output results/fresh_h9_new_host_next_20261013_a9_independent.json --verify
```

The raw norm interface can be regenerated with

```sh
gp -q scripts/fresh_h9_new_host_next_20261013_a9_orbits.gp
```

Its archived [transcript](fresh_h9_new_host_next_20261013_a9_orbits_raw.txt), [exact orbit replay](fresh_h9_new_host_next_20261013_a9_orbits.json), [profile census](fresh_h9_new_host_next_20261013_a9_profiles.json), and source hashes make the dependency reviewable. All jobs in this lane are complete.
