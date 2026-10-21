# Independent exact audit of the 3,040 archived cap profiles

**PASS, with a census correction: the 3,040 profiles represent 3,034 distinct
cap mirrors.** Every archived normal has square two, admissible relations
with all sixteen seed walls, and the same strict timelike interior point.
The six duplicate pairs are exact repeated roots that the numerical scout
had separated slightly. No cap-cap compatibility, compactness, or finite
volume is concluded from this audit.

The [independent checker](../scripts/fresh_h9_split_leaf_3040_independent_20261013.py)
uses only the standard library and its own Fraction arithmetic. It imports
no producer or repository arithmetic implementation. Its
[result archive](fresh_h9_split_leaf_3040_independent_20261013.json)
records all counts and duplicate pairs. The audited
[exact normals](../data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json)
are anchored by the SHA-256 digest of the saved numerical profile source.

## Exact reconstruction

The checker independently rebuilds the sixteen-wall seed with
`a=phi+1/2`, `b=-5/22+phi/11`, and verifies every archived Gram entry.
An explicit sum/difference congruence decomposition gives inertia
`(9 positive,1 negative,6 zero)`. The all-one coefficient vector `S`
has square `-128a=-64-128phi` and pairs `-8a` with every seed wall.

Each profile has finite negative pairing magnitudes in `{0,1,phi}` and,
where present, one common unknown dotted magnitude `d`. Instead of using
the producer's symmetric coefficient formula, the independent checker
inverts a nonsingular ten-wall Gram minor by exact Gaussian elimination,
sets the complementary six coefficients to zero, and reconstructs the
same covector. This different coefficient choice differs only by the
Gram radical. It independently recovers the norm equation
`A*d^2+B*d+C=0` over `Q(phi)`.

For each archived radical, exact arithmetic checks its radicand, selected
root, full normal, direct multiplication by the sixteen-wall Gram, and
square two. Sign evaluation uses rational comparisons in `Q(phi)` and in
`Q(phi,sqrt(D))`, with the specified nonnegative square root. It also
handles `D=0` without treating a formal radical coefficient as a nonzero
real number.

The root-selection check uses exact rational intervals, replacing the
producer's floating comparison: every reconstructed pairing lies within
`1/1,000,000` of its saved decimal value. Finite profile entries are
recognized within `1/100,000,000` of their prescribed exact values.
These intervals establish correspondence with the numerical proposals;
they are not used to decide a geometric angle approximately.

The exact seed-pair census, counting all archived profiles, is:

| Relation | Count |
|---|---:|
| Right angle | 4,368 |
| pi/3 | 13,888 |
| pi/5 | 14,624 |
| Strictly dotted | 15,760 |
| Total | 48,640 |

There are no parallel cap-seed pairs. All 3,040 cap normals pair strictly
negatively with `S`, so their half-spaces have a common strict timelike
interior. The certificate has 89 literal radicands, which is not asserted
to be a count of distinct number fields. The branch counts are 192 without
a radical parameter, 2,812 positive branches, and 36 negative branches.

## Exact duplicate correction

The following zero-based profile pairs have identical exact mirrors; the
first index in each row is the canonical representative:

| Canonical index | Duplicate index |
|---:|---:|
| 1852 | 1858 |
| 1853 | 1859 |
| 1854 | 1860 |
| 1855 | 1861 |
| 1856 | 1862 |
| 1857 | 1863 |

Their repeated root is `d=phi^2`. The numerical scout had separated it
into values near `2.6180341863` and `2.6180337912`; exact reconstruction
identifies both with the same root. The archived normal coefficients and
radicands are already literally equal within each displayed pair.

The distinctness check also excludes coincidences beyond literal encoded
equality. Since `d>=2>phi`, a finite-versus-dotted position cannot be
confused. Different exact finite/dotted patterns therefore have different
covectors. Identical patterns give the same norm equation and radicand;
the checker compares their selected roots exactly. This leaves **3,034
distinct norm-two cap mirrors**. No claim of novelty relative to other
constructions, or completeness among all possible caps, follows.

## Replay and scope

```sh
python scripts/fresh_h9_split_leaf_3040_independent_20261013.py
python scripts/fresh_h9_split_leaf_3040_independent_20261013.py --verify
```

The audit rejects a changed normal and a forged branch indicator.
Its exact checks establish a valid input bank for subsequent selection.
Relations between different caps require separate certification, followed
by an actual finite-volume test for any selected chamber. The earlier
228-cap bank exclusion does not exclude this larger bank.
