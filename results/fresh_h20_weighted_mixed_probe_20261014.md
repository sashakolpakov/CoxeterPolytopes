# Mixed weighted-box caps: an all-label covering obstruction

The retained weighted-box core has **no finite-volume completion using only radius-two spherical caps and seed-compatible caps centered at old box vertices, with any radius**. All finite Coxeter labels, right angles, parallel pairs and dotted pairs are allowed. Any pairwise Coxeter-compatible collection from this mixed family leaves at least **84,480 of the 524,288 half-vertex probes** outside every closed cap ball.

This closes the mixed-family question left open in the [frozen weighted-box report](fresh_h20_seed_deformation_20261014.md). It does not exclude arbitrary nonvertex radii, a changed central cap, deletion of old box walls, or other H20 constructions. No arithmetic hypothesis is used, and no finite-volume example is claimed.

Reproduce with:

```sh
python3 scripts/fresh_h20_weighted_mixed_probe_20261014.py --verify
```

The [new exact certificate](fresh_h20_weighted_mixed_probe_20261014.json) binds the [checker](../scripts/fresh_h20_weighted_mixed_probe_20261014.py) and both frozen precursor files. It freshly replays the precursor's 88 cutter recipes, all 7,744 adjacent-probe pair checks, and the radius-four compatibility bounds before checking the new claims. The [independent review](fresh_h20_weighted_mixed_probe_20261014_review.md) separately checks the new rational bounds and covering counts; it states its inherited-premise scope.

## Geometry and the cap families

Write `phi=(1+sqrt(5))/2`. There are ten unit coordinate inradii, one `sqrt(2)` inradius, and eight golden inradii `phi`. Coordinates below divide each physical coordinate by its inradius. Thus the old box is `[-1,1]^19`, with spatial quadratic form

```text
diag(1 repeated 10 times, 2, 1+phi repeated 8 times).
```

Every old vertex has squared norm `N=20+8phi`. The retained core consists of the 38 vertical box walls and the spherical cap centered at zero with radius two. The half-vertex probes have normalized coordinates `(±1/2,...,±1/2)` and squared norm `N/4=5+2phi`.

In Lorentz coordinates `(x,D,U)` with square `x²-2DU`, a spherical cap of center `c` and radius `R` has normal `(c,(c²-R²)/2,1)`. At an ideal point of Euclidean coordinate `x`, its inequality is strict precisely when `|x-c|²-R²>0`. A cover by the closed cap balls is therefore necessary for finite volume. The theorem concerns finite coherent cap collections, with the core retained and every pair satisfying a Coxeter gate.

## The radius-two coordinate alphabet is complete for every label

For a radius-two cap, the two opposite physical gaps in a coordinate of inradius `b` sum to `2b`. Acuteness forces both gaps nonnegative. Their smaller gap is at most `b<2`, so it is either zero or `2cos(pi/m)`. Monotonicity in `m` restricts the smaller label to at most 3, 4, or 5 for `b=1`, `sqrt(2)`, or `phi`, respectively.

The only complementary gaps not already resolved by the frozen finite-label calculation are

```text
2sqrt(2)-1,       2phi-sqrt(2).
```

Both lie strictly between `91/50` and `46/25`. On the other hand,

```text
2cos(pi/7) < 2[1-(3/7)^2/2+(3/7)^4/24] < 91/50,
2cos(pi/8) = sqrt(2+sqrt(2)) > 46/25.
```

The first inequality uses `pi>3` and the strict fourth-order cosine upper bound. The checker verifies every remaining comparison by rational squaring and rational enclosures. Thus neither exceptional complement is an allowed finite angle; both are below two, so neither is parallel or dotted.

Consequently the complete physical center alphabet is

| Coordinate type | Allowed coordinates |
| --- | --- |
| Unit | `0, ±1` |
| Heavy | `0, ±sqrt(2)` |
| Golden | `0, ±(phi-1), ±phi` |

This removes the precursor's seed-label cutoff. Central and mutual radius-two cosines lie in `Q(phi)`; their only finite values are `0`, `1/2`, and `phi/2`. Indeed the degree of `2cos(pi/m)` is `EulerPhi(2m)/2` for `m>=3`; the degree-at-most-two possibilities are `m=3,4,5,6`, and the square-root-two and square-root-three cases are not in `Q(phi)`.

## Radius-two half-probe cutters and the cube bound

The replayed precursor classifies every radius-two cap that is compatible with the core and covers a fixed half-probe. There are exactly 88 recipes. Such a center has:

- one full golden coordinate, with the probe's sign;
- seven small golden coordinates `±(phi-1)`, with matching signs;
- a full heavy coordinate, with matching sign;
- matching full unit coordinates, with either one unit coordinate zero or none.

Every cutter has probe power `8-5phi<0`. A golden sign change raises the squared distance by at least two, so one cutter cannot cover two different golden sign patterns in a fixed unit/heavy fiber.

For probes differing in golden bit `j`, two selected cutters are compatible exactly when one of their full golden coordinates is `j` and the other is not, and their zero-unit indices are distinct and both present. Every such compatible pair is right. This exact rule is checked on all 88² pairs for one bit and transported by coordinate permutations and sign changes.

Now fix the eleven unit/heavy signs. Among the 256 golden sign patterns, let `A` be those covered by a coherent radius-two family, and select one cutter for each. Orient each induced cube edge toward the unique endpoint whose full golden coordinate matches that edge's direction. Each vertex accepts at most one edge. Thus, writing `v=|A|`, the induced edge count satisfies

```text
e(A) <= v.
```

Degree eight in the full cube also gives

```text
2e(A) = 8v - e(A,A^c) >= 8v-8(256-v),
e(A) >= 8v-1024.
```

Hence `7v<=1024`, so at most **146** probes per fiber are radius-two covered. There are 2,048 fibers, giving at most **299,008** radius-two-covered half-probes.

This is a global cube count, not a claim that every square is impossible. In fact a golden square is feasible. Around patterns `00,10,11,01`, choose `(full golden coordinate, zero unit coordinate)` as `(0,0),(1,1),(0,0),(1,1)`. The four edge pairs are right and the two diagonals have label three. The certificate checks this explicit four-cap configuration and its four strict probe cuts.

## Radius-four counting, and why all other vertex radii are irrelevant

For two signed vertices let

```text
H = h1 + 2h2 + (1+phi)h5,
```

where the Hamming differences count ten unit bits, one heavy bit and eight golden bits. A radius-four ball at a vertex contains a half-probe exactly when

```text
N/4 + 2H <= 16,
H <= 11/2-phi.
```

There is no boundary equality. The exact number of covered half-probes is **275**: 187 with no golden sign change, and 88 with one. The checker exhausts the 198 Hamming triples with binomial multiplicities.

The frozen compatibility bound allows at most **512** radius-four vertex caps in a coherent family. Briefly, at most four can have the same golden pattern by the weighted Plotkin bound; occupied golden patterns form an independent subset of the eight-cube, so at most 128 patterns occur. The proof uses all possible pair labels in `Q(phi)`, not a chosen bounded label list. Therefore radius-four caps cover at most **140,800** half-probes.

For a vertex-centered cap of arbitrary radius, the opposite unit wall has physical gap two. Its finite contact requires `2/R=cos(pi/m)`. Label three gives `R=4`; every label at least four gives `R²<=8`; a parallel or dotted contact gives `R<=2`. A right contact is impossible at positive gap two. Thus every seed-compatible vertex radius other than four has square at most eight. But every half-probe is at squared distance at least `5+2phi>8` from every vertex, so all those other balls strictly miss every half-probe.

## The mixed contradiction and volume scope

The central ball also misses every half-probe. Even ignoring overlaps, the entire mixed family covers at most

```text
299008 + 140800 = 439808 < 524288.
```

At least **84,480** half-probes remain outside every closed ball. Equivalently, full coverage would require the radius-four family to cover at least `2048*(256-146)=225280` probes, exceeding its 140,800 capacity.

Each missed probe is strictly inside the vertical box and has strictly positive power for every selected ball. For a finite collection these strict inequalities persist in a neighborhood and under a sufficiently small spacelike Lorentz perturbation. Thus there is an open ideal-boundary patch, and the resulting full-dimensional chamber has infinite volume. If a proposed collection has no full-dimensional interior, it already fails to be a finite-volume H20 polytope. The obstruction needs no assertion that every possible family member is mutually compatible: it applies to every finite coherent selection.

The next construction must therefore change at least one substantive premise: use a nonvertex radius other than two, change or remove the central cap, or change the retained box walls. Enlarging either of the two now-excluded dictionaries cannot close this core.
