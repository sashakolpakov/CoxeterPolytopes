# Independent review: integral A8 period and alpha2 interval worker

20 October 2026. **PASS for the conditional period theorem and source
coverage, before execution of the alpha2 worker.** No import, census,
scientific job, Lean or typesetting was performed by this reviewer.

Sources read in full:

* `results/compact_h9_sqrt3_a8_period_isometry_20261020.md`, SHA256
  `f4d80265de3c167c82e741af40d45cb62957d63ff020211e24b7616f4229137b`;
* `results/compact_h9_sqrt3_e8_alpha2_interval_20261020.md`, SHA256
  `c19b41db49c82c66e3a2a4b0bc2c14cac3a246af09071e2c0431a106a979edf1`;
* `scripts/check_compact_h9_sqrt3_e8_alpha2_edge_20261020.py`, SHA256
  `0e5712a0627d948a6d273f6cec1a9a0d7e858ad78125cfda532b68c35243b961`.

The earlier full integral splitting, both first-height walls, complete
eighteen-profile inventory and actual A7 generation/replay were separately
read in `review_compact_h9_sqrt3_e8_prefix_and_escape_finite57_20261020.md`.
Their mathematical boundaries remain explicit inputs here.

## Integral period and global implication

I independently recomputed the norm2 root r, its A8 contacts and the two
future corner vectors. Both have square9(1-s); V1's orthogonality follows
particularly from `14h+17=73+42s` and `lambda*a=16+10s`. The reversed
chain prescription therefore defines a unique K-isometry.

The original-coordinate expansion of a general lattice vector is correct.
Consequently integrality reduces to the displayed denominator9 glue
vector. Expanding the E8 inverse columns gives precisely its time
coefficient-21-12s and its eight listed integral root coefficients.
This proves A(M) is contained in M. The isometry determinant
is+1 or-1, so the integral adjugate proves surjectivity as well. Future
orientation follows from V0 mapping to V1. This is a full lattice
automorphism, not only an isometry of the index-nine corner sublattice.

I also rederived the trace calculation. The first seven diagonal entries
sum to-20-12s. Solving the two remaining image equations gives the v
coefficient69+40s and beta coefficient-44-25s, hence

    tr A=5+3s>9.

A finite-order real rank-nine matrix cannot have this trace. The map has
infinite order independently of the proposed first-wall premise.

The global argument is sound **if** the specified upper root is certified
as the first simple wall on the alpha2 edge. The previous A7 certificate
puts V0 and v in one full chamber; that additional premise puts V1 there
too, with exactly the displayed outward simple germs. A lattice
automorphism permutes the entire integral-reflector arrangement. Matching
these inward germs then identifies the whole image chamber with the same
chamber. It is not enough merely to match two equal-norm lattice vectors.

An infinite-order chamber stabilizer excludes compactness. The admissible
anisotropic form would make a cofinite integral reflection chamber compact,
so under the first-simple-wall premise M is nonreflective. The established
full scalar-facet transfer then excludes reflectivity of the original H9
host. Without that premise the map can move the chamber and gives no such
conclusion. No conclusion about other forms or compact H9 polytopes follows.

## Alpha2 interval: algebra and complete bounds

The upper root has norm2, is primitive through its contact-1, and is an
integral reflector. I checked its four cap products and the identity
`Z*^2-14lambda*a*^2=9`. The positive endpoint precedes the null point.
The eighteen profiles include the six mod-seven-impossible rows instead
of assuming their elimination. Every R, placement, norm and inverse-Cartan
contact matches the independently reviewed finite inventory.

For positive a,Z the collision ordering is a/Z. The upper ratio and the
norm equation give exactly `9a^2<=R*a*^2` and `9Z^2<=R*Z*^2`. The stated
rational bounds yield a<158,Z<500,|a'|<1,|Z'|<6. Taking physical/conjugate
half-sums gives the complete integer boxes in the source, including
`a0,a1>=0`. These nonnegative coefficients are consequences of the bounds,
not an additional root assumption. Divisibility by eta is exactly the
stored coefficient parity condition. The exact dictionary target uses
`1/(14lambda)=eta/28` with the correct signs.

Candidate reconstruction imposes h in O with h>=0, all eight old wall
products, the exact root square and reflection integrality. Primitivity
follows because a nonunit common factor would leave a represented unit
norm, already excluded for M. Cap products are computed and reported but
**never used to discard census candidates**. This preserves the difference
between a prefix facet and a full-chamber facet.

## Simultaneous walls, endpoint and execution boundary

For tied candidates the tested difference is
`a_left*right-a_right*left`. Its time coefficient is zero. Multiplying
its weight coefficients by C_E8 inverse correctly gives its simple-root
coefficients, with the released alpha2 coefficient zero by the collision
identity. Nonnegativity at the seven white roots is precisely the needed
inequality implication. The source claims a simple wall only when all
dominating choices give one positive projective ray; otherwise it records
UNKNOWN and exits2. The explicit-upper flag requires the returned root
itself, rather than an arbitrary simultaneous candidate.

The endpoint `Z*v-lambda*a*omega2` has square-lambda*R<0 and positive time.
Its old-wall, all-candidate and cap feasibility checks have the correct
directions. If the selected root is the explicit upper, this endpoint is
exactly V1 in the conditional period theorem. A mere PASS with a different
selected root would not establish that theorem's premise.

The worker pins its five written inputs and itself before/after, uses an
exclusive bounded JSON write and canonical replay, and has a45-second
alarm. It owns no shared lock; the root wrapper must enforce the global
lock, resource gates and memory cap. A replay uses the same arithmetic
implementation. No actual alpha2 execution result is asserted here.
