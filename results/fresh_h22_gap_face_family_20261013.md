# H22: a complete face-based cap orbit with an unrestricted angle label

A new wall family cuts the preceding exact H22 escape and is compatible with all **2,252,588 retained facets**. Completing its 241,920-element symmetry orbit and taking both height signs gives **2,736,428 actual facets**. Its angle label is any integer `k≥6`; the old center label remains any integer `m≥3`.

The family nevertheless has **infinite volume for every such pair of labels**. At `k=6`, a new parallel pair exposes a deficient cusp link. An explicit strict spacelike perturbation near that cusp remains feasible for every `k≥6`, proving the uniform conclusion. This is a constructive new mechanism and an exact family-level exclusion, not a new finite-volume dimension or an obstruction to adding other walls.

## The face and the exact new normal

The preceding [cusp-orbit construction](fresh_h22_cusp_orbits_20261013.md) left an escape above a two-dimensional face of the old Borcherds H21 chamber. Its 19 active old facets are

```text
1,2,3,9,10,14,21,22,23,28,32,34,39,45,46,48,50,63,198.
```

Their old-space orthogonal complement has dimension three. Prescribe a rational vector `p` perpendicular to those 19 normals and with pairing `−1` to old normals 5, 12, and 26. These 22 independent equations determine it exactly. The [certificate](fresh_h22_gap_face_family_20261013.json) stores its coordinates and proves

`p²=−31/3`, `(h,p)=−36`.

For an integer `k≥6`, set

```text
q_k = 1 / (2 cos²(π/k)),
H_k = sqrt(31/3 + q_k),
r_g,± = (g p, ±H_k).
```

Here `g p` runs through the full orbit of the six integral old chamber isometries. Every new wall has norm `q_k`, 19 right angles and three angles `π/k` with the old walls. The other 188 old-wall relations are dotted. This uses an exact face-based linear system, rather than restricting old projections to sampled old vertices or cusps.

## Why the complete orbit is compatible

The integer array stores `6 g p`. The [exact checker](../scripts/fresh_h22_gap_face_family_20261013.py) verifies closure under every generator, connectivity by earlier predecessors, distinctness, and preservation of the form, center, and old wall set. The orbit has 241,920 vectors. Pairing one representative against the full orbit gives

```text
p² = −31/3,
max (p,gp), gp ≠ p = −35/3,
all remaining distinct pairings ≤ −17.
```

Exactly one other orbit vector achieves the maximum. Since `1/2<q_k≤2/3`, every distinct same-sign new-wall pairing is at most

`−35/3 + H_k² = q_k−4/3 ≤ −q_k`.

Thus these pairs are parallel at equality or dotted. Opposite-sign pairs are dotted, including the two signs over the same old vector. Orbit transfer proves all pairs without quadratic enumeration.

The previous timelike caps have old projections `u/2` with squared norm `−3/2`. The exact representative maximum is `(p,u/2)=−8`. For all prior cusp caps `(δ,±t)`, the exact bound is `(p,δ)/t≤−5`. Together with `H_k≤sqrt11`, these make every new/prior-cap relation dotted:

```text
−8 + sqrt(55/2) < −1 < −sqrt(q_k),
−5 + sqrt11 < −1 < −sqrt(q_k).
```

The old center caps have height at most four and norm at most two, so `(h,p)=−36` gives an ample strict dotted bound uniformly for `m≥3`. The preceding complete-orbit certificates transfer these representative inequalities to every new cap.

The same strict timelike center works throughout. Its projection onto each wall is timelike and strict for all other inequalities, proving that all retained and new walls are actual facets. Each previously completed inherited cusp link remains unchanged because the new caps have timelike old projections and are strict at every old cusp.

The small labels also have a precise explanation. In this same full-orbit ansatz, `k=3` makes the nearest distinct pair positive; `k=4` gives forbidden normalized magnitude `1/3`; and `k=5` gives `sqrt5/3`, strictly between `cos(π/4)` and `cos(π/5)`. These are actual facets even for the small labels: projecting the center onto a new wall gives `x=h+(36/q_k)r`. Every distinct new pair is strictly less than `q_k`, so every other new inequality at `x` is strict; all new/old pairings remain nonpositive because `H_k≤sqrt(37/3)<4`. Thus projection also satisfies all old inequalities strictly. The small-label failures are not artifacts of redundant inequalities. No omitted small integer label repairs this whole-orbit family. This observation does not exclude selecting a proper suborbit or changing the old projection.

## A new deficient cusp and a uniform escape

At `k=6`, each new signed cap has exactly one parallel partner. For the saved representative and its partner `p'`,

`δ_new=(p+p',2sqrt11)`

is future null. The complete incident old set has 18 independent normals:

```text
1,2,3,9,10,21,22,23,28,32,34,39,45,46,48,50,63,198.
```

The parallel pair is orthogonal to these. Every other new or previously added cap is strict at this ray. The horospherical normal rank is therefore `18+1=19`, leaving at least two free directions in the H22 horosphere of dimension 21. Its link is unbounded.

The checker constructs an exact tangent `v` in these missing directions, with old-space norm `279/14`, perpendicular to both `p,p'` and all 18 common old normals. With `ε=1/1000`, the explicit point

`x = δ_new + ε(v,0) + (ε²/100)(h,0)`

is strictly feasible for every wall and has square

`647099999951 / 35000000000000000 > 0`.

All old coordinates are rational; the added coordinate is `2sqrt11`. The largest new-cap pairing at `k=6` is `−9/25000000`. Exact integer inequalities check every old, timelike-orbit, cusp-orbit, and new-orbit wall, with explicit overflow bounds. Squared comparisons handle the remaining square roots.

Increasing `k` decreases `H_k`, so the paired new-wall constraints become weaker. The same strict escape therefore works for every `k≥6`; the center-cap bound works for every `m≥3`. Conversely, the new representative cuts the preceding exact escape even at the limiting smallest height `sqrt(65/6)`, so the target cut is also uniform in `k`.

## Search scope and replay

A [numerical polar diagnosis](../data/enumerations/fresh_h22_gap_polar_20261013.json) found 91 pairs outside its selected finite-angle alphabet, explaining why simply adding the gap polar failed. A preliminary broader anchor search located the useful triple and was interrupted; it carries no census claim. The archived [27-profile face scout](../data/enumerations/fresh_h22_gap_face_20261013.json) tests only anchors 5,12,26 with labels 3,4,6. It finds 15 individually compatible numerical cutters. The complete exact family above is derived and checked separately; compatibility among those 15 proposals is not asserted.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_gap_face_family_20261013.py --verify
# Optional reconstruction of the orbit array, followed by certificate production:
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_gap_face_family_20261013.py --rebuild-orbit
```

The [independent source/proof review and exact replay](fresh_h22_gap_face_family_20261013_review.md) pass, including the small-label facet premise and the uniform escape bounds. This review reuses the existing arithmetic implementation.

The [orbit array](../data/enumerations/fresh_h22_gap_face_orbit_20261013.npz) specifies every old projection. The checker imports no numerical proposal code. Existing compatibility and Lorentz signature are explicit hash-anchored dependencies of the preceding certificate. A further completion would have to address these new ends or change the orbit; extending the integer label alone cannot yield finite volume.
