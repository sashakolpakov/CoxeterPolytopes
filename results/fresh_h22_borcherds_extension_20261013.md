# A complete cap orbit above the Borcherds H21 chamber

This tests a new construction in H22. It gives an exactly verified family
with **121,172 actual Coxeter facets**, but one exact strict escape proves
infinite volume for every member. No new finite-volume dimension or
literature novelty is claimed. The useful outputs are a complete compatible
cap orbit, an explicit remaining target, and an all-real obstruction to
changing only the original wall heights.

## Original directions cannot tilt in one extra coordinate

The input is the repository's exact 210-wall Borcherds H21 chamber. Write
its normals as r_i and adjoin an orthogonal coordinate of square one.
Suppose each old direction is retained as a Coxeter wall `(r_i,t_i)`,
allowing every t_i to be real and individual positive rescaling.

The 42 norm-two normals have a connected, triangle-free label-three graph;
all other short-short pairings are zero. Each same-sign nonzero height
support must occupy a clique, since equal-sign heights on old orthogonal normals
give positive pairing. Thus at most four short heights are nonzero. A zero
height exists, and propagates along every label-three edge: a zero/nonzero
height pair would have cosine strictly between zero and 1/2. Consequently
all 42 short heights vanish.

Every norm-four long wall has pairing −2 with a short wall. Its height must
therefore be 0 or ±2, corresponding to labels four or three. Every long
wall also has a label-three long neighbor. The complete nine-case height
table shows that such a pair requires both heights zero: a single nonzero
height gives squared cosine 1/8, equal nonzero signs give positive pairing,
and opposite signs give squared cosine 9/16, between the label-four and
label-five values. All 210 heights are zero.

This excludes nontrivial one-coordinate lifts retaining the old projected
directions. It leaves additional walls and changed projections open. The
[standard-library certificate](fresh_h22_borcherds_lift_rigidity_20261013.json)
checks all graph and pair premises and the complete finite height table.

## A symmetric pair and a complete orbit of additional caps

Let h be the primitive sum of the old normals. The exact input gives

```
h² = −14,
(h,r_i) = −1 for the42 short walls and −8 for the168 long walls.
```

For every integer m≥3, put

```
q_m = 1/(2cos²(pi/m)),
C_± = (h, ±sqrt(14+q_m)).
```

Each center cap has norm q_m, angle π/m with every short wall, and dotted
relations with every long wall and with its opposite center cap.

An old-edge search supplies a further cap `(v,sqrt(5/2))`, where v is a
half-integral old vector with v²=−3/2 and (v,h)=−33/2. Its exact coordinates
are archived in the [orbit certificate](fresh_h22_symmetry_orbit_verify_20261013.json).
The six saved chamber symmetries generate **60,480 distinct projections**.
The audit verifies integral generator matrices, preservation of the form,
center and old walls, closure of the orbit, and connectivity from the
chosen representative. It does not assume a claimed group order to infer
completeness.

All 120,960 signed caps `(g(v), ±sqrt(5/2))` coexist. Their norm is one.
One representative has 21 right angles, one π/3 angle, two π/4 angles and
186 dotted relations with the old walls. Within its same-sign orbit it has
four π/3 relations, four parallel relations and 60,471 dotted relations.
Every opposite-sign pair is dotted. Every orbit cap is dotted with both
center caps, uniformly for m≥3, because

```
−33/2 + 2sqrt(10) < −2,   q_m ≤ 2.
```

Symmetry transfers these checks to all **7,341,266,206 unordered pairs**.
The center `(h,0)` is strict timelike. Projecting it to any one wall gives
a timelike point strict on every other wall, proving that every inequality
defines an actual facet. Exact rational elimination also certifies ambient
signature (22,1).

## Exact remaining escape

The certificate stores a rational point w with

```
w² = 6604718592837147 / 10000000000000000 > 0,
(w,(h,0)) = −187152873 / 50000000 < 0.
```

Every old wall pairing is strictly negative. Both signs of every orbit cap
are checked exactly by a squared integer inequality. Since the center-cap
height is at most four, the uniform upper bound for either center pairing
is `−70003/50000000 < 0`. Thus the same strict escape works for every m≥3.
Together with the strict timelike center, it proves infinite volume.

The [verifier](../scripts/fresh_h22_symmetry_orbit_verify_20261013.py)
imports no scout or project arithmetic implementation. It uses rational
arithmetic and native integer matrix products with explicit overflow
bounds. Two separate source reviews and fresh replays passed, including
the actual-facet argument and the all-m escape. These reviews replay this
arithmetic; they are not additional arithmetic implementations. The orbit's
[compressed integer coordinates](../data/enumerations/fresh_h22_symmetry_orbit_20261013.npz)
are hash-anchored by the certificate.

## Numerical exploration and remaining scope

A preliminary 1,000-objective LP scout reconstructed 1,000 old vertices
exactly: 773 finite and 227 ideal. This is a sample, not the old vertex
census. A bounded angle search produced 26,894 cap proposals over those
directions. Seven center-label choices all retained numerical escapes even
after imposing their entire compatible-with-base banks simultaneously.

Allowing projections to move in twenty-wall residual spaces produced
1,933 two-anchor entries, from 39,103,560 quadratic profiles. Deduplication
leaves 1,813 numerical proposals. None cuts the saved m=5 whole-orbit gap.
That last statement is numerical and confined to those saved proposals;
it is not an obstruction to arbitrary additional walls. The angle alphabet,
support limits and tolerances are recorded in the
[two-anchor scout](../data/enumerations/fresh_h22_two_anchor_20261013.json)
and [extension screen](../data/enumerations/fresh_h22_orbit_extension_scout_20261013.json).

The next constructive change would need directions outside that bank,
different orbit choices, or replacement of retained walls. Merely adding
more members of the already complete orbit cannot alter its escape.

```sh
python scripts/fresh_h22_borcherds_lift_rigidity_20261013.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_symmetry_orbit_verify_20261013.py --verify
```
