# Targeted H20 caps and the old-vertex direction census

The previous 85-wall H20 construction extends to **88 mutually compatible, actual Coxeter facets**. Three exact caps successively cut identified spacelike gaps, and the first and third create certified finite vertices. The final 88-wall chamber still has an exact strict future spacelike escape. **This is an infinite-volume construction, not a new finite-volume example in dimension20.** A separate exact audit excludes finite-volume completion of the87-wall stage by the entire8346-direction old-vertex height family, with arbitrary real heights and unrestricted Coxeter labels. The88th cap changes that family. The earlier [85-wall construction](fresh_dimension_bridge_20261013.md) and its source files remain unchanged.

All indices below start at zero. The old form has signature `(20,1)`, with the added coordinate having square3. An entry `[[a,b],[c,d]]` in a quadratic cap represents

`a + b sqrt(6) + (c + d sqrt(6)) sqrt(D)`.

The indicated square roots are positive at the geometric embedding. Cross-field signs are evaluated in an iterated real quadratic extension; raw coefficient equality is not used to infer equality when radicals may be dependent.

## Exact constructive checkpoints

| Walls | New cap | Relations with preceding walls | Exact global result |
|---|---|---|---|
| 86 | Old projection of rejected polar27; height `-4+sqrt(10)` | 19 right, 3 of angle `pi/4`, 63 dotted | Cuts the strict 85-wall escape; creates a finite vertex; another strict escape remains |
| 87 | Projection from the next open direction; height `-sqrt(268/3+36sqrt(6))` after fixing its short-wall pairing | 18 right, 1 of angle `pi/3`, 67 dotted | Norm2; cuts the 86-wall escape; another strict escape remains |
| 88 | Two-angle solution in the complement of 18 short walls | 18 right, 2 of angle `pi/3`, 67 dotted | Norm2; cuts the exact universal escape of the old finite-vertex direction family; creates a finite vertex; another strict escape remains |

For each extension, the preceding walls are retained, the saved timelike center remains strictly inside, every new normal is spacelike, and every old/new pairing passes the exact Coxeter gate. These acute inequalities define true facets: orthogonally projecting the strict center to any one wall gives a timelike point strict on every other wall. A complete compatible Gram gate is therefore meaningful here; it is not a test of unverified inequalities from a nonacute section.

The exact strict escape in the final chamber has

```
(w,w) = 478875483478515721 / 32000000000000000000 > 0,
(w,h) = -48815755199 / 4000000000 < 0,
(w,r_i) < 0 for all 88 walls.
```

Together with the strict timelike center, such a point gives an open set on the ideal boundary of the chamber, and hence infinite volume. Numerical LPs only proposed these vectors; the archived signs and norms are checked exactly.

### The first cap

Write the rejected polar as `(v,-1)`, where `(v,v)=-2`. Replacing its last coordinate by `-4+sqrt(10)` gives

```
(r,r) = 76-24sqrt(10) > 0,
(r,h) = -110,
(r,w_85) = 1225/128-3sqrt(10) > 0.
```

At `t=4-sqrt(10)`, the vector `(3t v,-2)` is the saved actual finite vertex. It lies on 19 old short walls and the new wall, with all other inequalities strict. All107 inherited ideal vertices have strictly negative pairing with this cap.

Artifacts: [cap](../data/enumerations/fresh_h20_targeted_escape_20261013_exact.json), [escape](../data/enumerations/fresh_h20_targeted_escape_20261013_volume.json), [standalone independent audit](fresh_h20_targeted_escape_20261013_independent.json). The independent audit uses separate Fraction arithmetic and rejects forged normal, norm and cut data.

### The second cap

The next open direction is normalized so its pairing with short wall21 is `-1`. Its last coordinate is `-sqrt(268/3+36sqrt(6))`; its norm is exactly2. It cuts the preceding escape and is strictly negative on all107 inherited cusp rays. No finite vertex is asserted for this particular cap.

Artifacts: [cap](../data/enumerations/fresh_h20_targeted_escape_20261013_second.json), [escape](../data/enumerations/fresh_h20_targeted_escape_20261013_second_volume.json), [standalone independent audit](fresh_h20_targeted_escape_20261013_second_independent.json). Its strict escape has square `73993/8388608`.

### The third cap changes the direction family

The previous finite-vertex height family left a common numerical gap. At that point17 original short walls were active. Extending them to18 short walls leaves a three-dimensional Lorentzian complement, in which two independent angle equations and the norm equation can be solved.

The selected exact cap is scout entry1. Its18 orthogonal short walls are

```
0,2,3,4,5,7,9,10,11,12,13,14,15,16,17,19,20,21.
```

It has norm2 and pairings `-1` with old wall8 (norm2) and `-3` with old wall42 (norm18), so both angles are `pi/3`. All67 other nonorthogonal relations are strictly dotted. Its coordinates lie in `Q(sqrt(7))`; the archive retains the common tower format with `D=63` over `Q(sqrt(6))`.

The actual finite vertex is incident to the18 listed walls, **old wall42**, and new wall87. Its20-by-20 Gram matrix is positive definite by exact elimination; projection of the strict center gives a future timelike point with precisely those20 zero pairings and every other pairing negative. The alternative set using wall8 is not asserted to be a finite vertex.

The saved cut target is a rational perturbation of the common numerical bank point, checked to be strictly inside all87 previous inequalities and spacelike. This certificate does not by itself certify its relation to every numerical bank entry.

The separate all-height direction audit also checks that this cap cuts its
exact universal strict escape. The positive pairing is

```
(-78644145/2097152 + 45sqrt(10)/4)
  + (-4718647/1048576 + 3sqrt(10)/2)sqrt(63) > 0.
```

Its old spatial projection is orthogonal to exactly18 short walls and no
long walls, explicitly placing it outside all8346 excluded finite-vertex
directions. This comparison is exact and does not rely on the numerical
bank's completeness.

Artifacts: [cap and finite vertex](../data/enumerations/fresh_h20_targeted_escape_20261013_third.json), [88-wall escape](../data/enumerations/fresh_h20_targeted_escape_20261013_third_volume.json), [standalone independent audit](fresh_h20_targeted_escape_20261013_third_independent.json). The independent audit checks all87 new pairings, all107 inherited cusp signs, the actual finite vertex, target cut and remaining88-wall escape, and rejects forged normal, norm and cut data. The producer intentionally retains its original `VOLUME_PENDING` status: the separate immutable volume artifact supplies the subsequent negative volume result.

## Complete original finite-vertex direction census

The original H19 chamber has24 norm-two short walls and12 norm-six long walls. A finite vertex can contain at most one long wall, because every pair of long walls already has nonpositive two-by-two determinant. The exact census finds:

| Direction class | Candidate supports checked | Positive-definite supports |
|---|---:|---:|
| 19 short walls | all `C(24,19)=42504` | 2046 |
| 18 short walls and one long wall | 22116 after exact local exclusions | 6300 |
| Total | | 8346 |

For every accepted support, exact rational projection of the strict H19 center gives a timelike point, zero exactly on its19 incident walls and strictly negative on every other original wall. Positive definiteness is decided by sparse rational symmetric elimination, not floating eigenvalues.

The one-long support reduction is exact: two short neighbors of a long wall have a negative three-by-three principal determinant; one such neighbor connected to another selected short wall gives a zero or negative principal determinant. Thus a long wall has no selected neighbor, or just one short neighbor isolated in the selected short diagram. Every remaining support is enumerated.

Artifacts: [2046 short directions](../data/enumerations/fresh_h20_targeted_escape_20261013_short_vertices.json), [6300 one-long directions](../data/enumerations/fresh_h20_targeted_escape_20261013_long_vertices.json).

There is an immediate **all-real-height obstruction for the6300 one-long directions**. At such a projected vertex `v`, an incident old long wall `L` satisfies `(v,L)=0`. Its retained lifts have last coordinates `2 epsilon` and `-epsilon`. For a proposed new normal `(v,t)`, acuteness with both requires

```
6 epsilon t <= 0,   -3 epsilon t <= 0,
```

so `t=0`. But then the normal is timelike, since `(v,v)<0`. This excludes every real height in these directions, independently of any finite-angle alphabet or later cap.

## Exact obstruction for all real heights in the8346 directions

The [independent direction-family certificate](fresh_h20_direction_family_audit_20261013.json) re-enumerates the short supports by cycle rejection and rational leaf Schur elimination, verifies the long support reduction, and checks every saved projection. It then reconstructs a common feasible spacelike point in `Q(sqrt(10))`, without trusting the numerical bank coordinates.

For all2046 short directions, acuteness against the first48 retained walls supplies exact rational upper and lower height bounds. The target pairing is affine in height. Maximizing it over these bounds leaves only directions1400 and1511 as possible cutters. The two remaining continuous problems are excluded exactly:

* Direction1400 is scaled to old norm `-2`. Cutting requires `-1 <= t < sqrt(10)-4`. Against wall26 its squared cosine is `2(1+t)^2/(-2+3t^2)`, which lies in `[0,1/2)`. Thus only a right angle or `pi/3` is possible, giving `t=-1` or `t=(-8+sqrt(14))/5`. Wall49 rejects both: the first has squared cosine strictly between0 and1/4, and the second lies strictly between `cos^2(pi/7)` and `cos^2(pi/8)`.
* Direction1511 is scaled to old norm `-6`. Cutting implies `-3/2 <= t < -29/20`. Its wall26 squared cosine lies in `[0,1/4)`, forcing `t=-3/2`. The resulting squared cosine with short wall3 is `2/3`, which cannot equal a finite Coxeter angle: `4cos^2(pi/m)` is an algebraic integer, whereas `8/3` is not.

The6300 one-long directions are excluded by the opposite-height lift argument above. Therefore every allowable cap in this entire direction family has nonpositive pairing with the same exact spacelike point. A fixed positive rational multiple of the strict center makes the point strict on the87 retained walls and every such cap simultaneously, while preserving its positive norm and future orientation.

**No collection of caps from this family can give a finite-volume completion retaining the87-wall stage.** This includes arbitrary real heights, all finite Coxeter labels, parallel relations and strictly dotted relations. It does not exclude the changed18-short direction used for the88th cap, other spatial directions, or replacement of old walls. The1593-entry bank and239 interval scouts motivated this theorem but are not used as completeness assumptions in its proof.

The independent checker uses its own rational census and `Q(sqrt(10))` operations; it reuses the separate Fraction quadratic-extension checker for the87-wall signs. Source review found no issue with the continuous argument.

## Bounded and continuous numerical searches

The exact direction census is separate from the following numerical cap searches.

* For each of the2046 short directions, vary the last coordinate and solve finite-angle equations against each of the first82 walls for labels `2,3,4,5,6,8,10,12`. Test proposals against all87 walls using labels2 through120 or dotted relations. The search skips near-zero quadratic leading coefficients and deduplicates numerically. It produced1593 entries, six of which cut the then-current87-wall escape. It is not an all-label classification.
* Twenty-four greedy compatible selections each contained238 new caps, giving325 walls with the retained87. Every selection had a numerical spacelike escape. Even the relaxation imposing all1593 entries had a numerical spacelike escape, so no exactification of the selected325-wall cones was pursued.
* A separate continuous-height scout partitions each height line at norm, sign and dotted-boundary events. It records239 open intervals with19 right relations and all other relations strictly dotted. None of the sampled interior heights cuts the common bank point. Endpoints and interval ordering remain numerical in this archive; this is not an exact universal interval theorem.
* The changed18-short family uses the17 short walls active at the common bank escape, extended by one additional short wall. Four of the seven supports are positive definite numerically. For every pair of remaining anchors it solves the norm-two equations with labels `2,3,4,5,6,8,10,12` and the parallel case. Both quadratic branches and the linear case are included. Among753867 profiles it found76 distinct numerical compatible cutters. Only scout entry1 is promoted here to an exact cap and exact volume checkpoint.

Artifacts: [1593-entry bank](../data/enumerations/fresh_h20_targeted_escape_20261013_vertex_bank.json), [whole-bank and24 selection scouts](../data/enumerations/fresh_h20_targeted_escape_20261013_bank_geometry.json), [continuous intervals](../data/enumerations/fresh_h20_targeted_escape_20261013_dotted_intervals.json), [18-short two-anchor scout](../data/enumerations/fresh_h20_targeted_escape_20261013_edge_directions.json).

## Replay and remaining scope

Run from the repository root. Exact cap and saved-escape replays do not rerun LP searches:

```sh
python scripts/fresh_h20_targeted_escape_20261013_exact.py --verify
python scripts/fresh_h20_targeted_escape_20261013_volume.py --verify
python scripts/fresh_h20_targeted_escape_20261013_second.py --verify
python scripts/fresh_h20_targeted_escape_20261013_second_volume.py --verify
python scripts/fresh_h20_targeted_escape_20261013_third.py --verify
python scripts/fresh_h20_targeted_escape_20261013_third_volume.py --verify
python scripts/fresh_h20_targeted_escape_20261013_independent.py --verify
python scripts/fresh_h20_targeted_escape_20261013_second_independent.py --verify
python scripts/fresh_h20_targeted_escape_20261013_third_independent.py --verify
PYTHONPATH=. python scripts/fresh_h20_direction_family_audit_20261013.py --verify
```

Rebuild the exact direction censuses and numerical scouts:

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_targeted_escape_20261013_vertex_bank.py --vertices --bank
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_targeted_escape_20261013_long_vertex_bank.py --vertices --bank
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_targeted_escape_20261013_bank_geometry.py
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_targeted_escape_20261013_dotted_intervals.py
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_targeted_escape_20261013_edge_directions.py
```

The remaining exact88-wall escape is an explicit target for a later construction. None of these results excludes additional mirrors outside the tested families, replacing old walls, or a different H20 construction. The76-entry18-short scout is a proposal list, not a proven complete reflection arrangement. There is no claim of global novelty for the88-wall infinite-volume chamber.
