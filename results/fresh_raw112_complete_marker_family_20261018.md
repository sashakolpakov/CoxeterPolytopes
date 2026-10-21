# Complete raw112 extra-vertex family and the shared marker allowance

The complete196-class raw112 extension exports529 valid coupled
columns. Three violate the previous fixed-mass certificate. The least
componentwise-nondecreasing vertex-allowance repair has exact ratio

    14671331978990040 / 1224897727063181
    = 11.977597520869008…,

down from12.015522073308947…. This is a feasible dual point below12,
not a universal upper bound, a primal certificate or an exclusion of
compact arithmetic H13. It shows that this particular above12 witness
does not survive the specified monotone repair. Other face masses, or
a nonmonotone redistribution of vertex allowances, may still give a
ratio above12. No optimizer is invoked.

The scope is the current compact arithmetic H13/Q(sqrt5) model.
The [producer](../scripts/fresh_raw112_complete_marker_family_20261018.py),
[exact certificate](fresh_raw112_complete_marker_family_20261018.json)
and [standalone saved checker](../scripts/check_raw112_complete_marker_saved_20261018.py)
extend the frozen [raw108 joint family](fresh_raw108_complete_marker_family_20261018.md).
Every column is valid at arbitrary nonnegative face masses. Full
coverage here means all196 raw112 spectator classes and both extra
vertices, not all geometric families or optimal enlarged local cones.

## Actual vertices, actual simplex faces and source recovery

The frozen [four-core pilot](fresh_raw14_raw112_raw77_raw108_core_incidence_pilot_20261018.md)
proved the raw112 eight-wall core has physical inertia(6,1,1), conjugate
inertia(7,0,1), and trivial native automorphism group. Its three counted
quadrilaterals use six distinct finite supports. The complete28
two-wall-deletion tests give two additional positive definite supports:
omissions `{2,4}` and `{2,5}`. Spectators have rank7. Adding their whole
orthogonal factors gives a positive definite13-root frame at each
additional support.

These are actual finite vertices in any realizing Coxeter polytope,
including one with arbitrary additional facets. Let `h` be a strict
timelike interior vector and `G` the acute positive definite Gram
matrix of the13 actual support normals. With `b_i=h·n_i<0`, the
inverse-positivity property gives `c=G^{-1}b<=0`. Then

    x = h − sum(c_i*n_i)

is orthogonal to the support, has
`x²=h²−bᵀG^{-1}b<0`, and is in the same time cone as `h`.
For any other actual wall `m`, all `n_i·m<=0`, so
`x·m=h·m−sum(c_i*n_i·m)<0`. Thus no additional facet can truncate
the recovered finite vertex. Within this fixed core-plus-spectator
configuration, the complete inherited support test also rules out
additional supports formed by deleting a spectator root;
the remaining core has a nonpositive direction by interlacing.

The two marker types are

| Type | Omitted walls | Active finite support | Ordered roles |
|---|---|---|---|
| A112 | 2,4 | 2A1 + H4 | A1 roots(0,7); H4 path(1,3,6,5) |
| B112 | 2,5 | 3A1 + H3 | A1 roots(0,4,7); H3 path(1,3,6) |

Both are vertices of the same actual compact4-simplex face. Fix core
walls0,7 and the seven spectator walls. The five residual sides are
1,3,4,5,6. Their exact residual Gram has Lorentz signature and strictly
negative inverse, giving a compact simplex. Its five finite vertices
are actual by the preceding projection argument. The polytope's
intersection with the fixed support lies inside that simplex, while
convexity and the five actual vertices put its entire convex hull
inside the polytope. Hence this is the actual face, not merely a
compatible Gram diagram.

At each extra vertex the four incident simplex sides give six pairs
of ports, each belonging to a compact triangular2-face. An ordinary
quadrilateral star cannot occupy any such pair: at a simple finite
vertex that pair determines one actual2-face, already the forced
triangle. This is a compatibility exclusion only. No extra triangle
demand or second triangle payment is introduced.

Each marker records both omitted walls' complete contact covectors
against the actual13-root frame. All spectator coordinates are zero.
Norm2 normalization fixes each spacelike projection, and its positive
squared height with the common future sign fixes the remaining
timelike coefficient uniquely. The producer checks the native omitted
product with those same heights. Thus a marker recovers the two
actual missing normals, the full eight-wall aggregate, its supporting
walls and its three actual quadrilaterals. Their compact full stars
are determined by the recovered actual walls, rather than by abstract
congruence alone.

The connected nonfinite core and whole passive forest identify the
active decomposition. Its trivial native automorphism group and the
duplicate-free unordered covector keys remove any role/stabilizer
ambiguity. Each actual aggregate contributes **one A112 occurrence
and one B112 occurrence**, each at its own omitted support. It does
not contribute three of either merely because it contains three
counted quads. The two extra vertices are distinct. All11,012 native
role orientations enumerate possible marker placements; they are not
multipliers on the source demand. New native keys are also explicitly
distinct from ordinary stars and all prior marker keys at each target.

## Complete relative-role reduction

The complete spectator inventory has196 source classes,392 pointed
marker domains and383 target forests. A112 has ordered pairs of
distinct A1 roles and one H4 role; B112 has ordered triples of distinct
A1 roles and one H3 role. Every finite-component embedding is
enumerated, with count checked independently from component
automorphisms and multiplicities. The resulting11,012 native maps
have positive heights and the prescribed native omitted product.

For a pair of markers, only equalities among their active component
roles affect the projection. Orthogonal spectator components have
zero contacts. The complete relative domain therefore consists of

| Pair domain | Cases |
|---|---:|
| A112/A112 | 13 |
| B112/B112 | 67 |
| A112/B112 | 13 |
| New types against raw56/57/76 | 55 |
| New types against A108/B108 | 14 |
| Total | 162 |

All host forests have rank13 and are canonically sorted before looking
up a frozen marker frame. The source checks exact frame equality, not
just equality of forest names. The independent helper reconstructs
the full relative signatures from complete active embeddings.

Each wall pair gets its own native field and projection-determined
heights. The I2(10) inverse squareclass is included only when a
contacted block needs it. The exact gates test nonpositive product,
absence of the forbidden interval(-1,0), and the required integral
square in Q(sqrt5). Every surviving test remains compatible in the
relaxation; survival is not physical realization. Fourteen of the162
relative cases survive.

Write an A112 marker as `(p,q;H4)` and a B112 marker as `(p,r,q;H3)`.
The surviving rules are:

* Two A112 markers use the same ordered endpoint pair and different
  H4 components, or disjoint endpoint pairs.
* A B112 clique fixes its middle A1 role `r` and its H3 component;
  its ordered endpoint pairs are mutually disjoint.
* A112 and B112 may use exactly the same ordered endpoint pair with
  the middle role outside it, or disjoint sets of all A1 roles.
* A112 is compatible with raw56 or raw76 only on disjoint A1 roles.
  B112 is compatible with raw56 only on disjoint A1 roles and is
  incompatible with raw76. Both new types exclude raw57.
* Both new types are compatible with every A108 marker. A112 avoids
  the A1 role of B108; B112 requires that role to equal its middle
  role and requires the same H3 component.

These statements include all relative cases, not only positive-mass
classes. The rank13 hosts contain every equality pattern; embedding
them into an arbitrary spectator forest changes no contacted product.

## One common marker allowance

Let `a,h,H,t,s` count A1,H3,H4,I2(10),I2(5) components in a target
forest. Coordinates are ordered

    (raw56, raw57, raw76, A108, B108, A112, B112).

Unavailable marker families have zero coordinates in every formula
below. A108 contributes at most `2*H*s` and coexists with all the other
non-raw57 types. Raw57 is exclusive and has capacity `t`.

First consider a clique with neither B112 nor B108. If its raw76
count is `k>0`, those markers consume `k+1` A1 roles: one common first
role and `k` distinct second roles, with distinct H3 components. For
`k=0` they consume zero A1 roles. Let `g` be half the remaining A1
count, rounded down. A112 uses a matching of endpoint pairs, with at
most `H` markers on each pair. Raw56 uses disjoint matching pairs,
with at most `t` markers per pair. Thus every such clique lies
coordinatewise below a row

    ((g−m)*t, 0, k, 2*H*s, 0, m*H, 0),

where `0<=k<=min(a−1,h)`, `0<=m<=g`, and `g=floor(a/2)` for `k=0`
or `g=floor((a−k−1)/2)` for `k>0`. If A112 is unavailable, take `m=0`.

Next suppose B112 or B108 is present. Their common middle/center
role consumes one A1, and raw76 is excluded. If both are present,
the relative rule identifies the same center and H3 component. The
A112 and B112 endpoint pairs form one matching on the other A1
roles: pairs shared between the two types have the same orientation.
On each matching pair there can be at most `H` A112 markers and one
B112 marker. Raw56 uses disjoint remaining pairs. With
`g=floor((a−1)/2)`, every clique is bounded by a row

    ((g−m)*t, 0, 0, 2*H*s, 2*s, m*H, m),

with unavailable-family coordinates set to zero and `0<=m<=g`.
This branch includes B108 present with B112 absent. The first branch
retains the case where both are absent, so reserving a center is
optional. These are upper bounds even when one type uses fewer than
the available matching pairs: augmenting its count to the displayed
capacity only increases a nonnegative charge.

Together with the exclusive raw57 row, these formulas dominate every
clique. The source constructs a literal clique witness for every
displayed boundary row and checks every pair using the exact relative
rules. Such witnesses establish sharpness in the conservative graph,
not realizability in a Coxeter polytope. The independent checker
reconstructs every row and witness.

For nonnegative marker coefficients, it suffices to impose one
inequality per boundary row, with right-hand side1. All vertices of
this rational coefficient polytope are enumerated exactly, and its
coordinatewise nondominated vertices are retained. This provides one
allowance shared by all seven marker types; there are no independent
allowances that could charge the same target vertex repeatedly.

## Ordinary stars and exported face coefficients

The complete current ordinary domains are reconstructed chronologically
from the327-class/600-forest inputs. All origins, actual corner
multiplicities, first-owner assignments and old ordinary Q–Q graph
outcomes are unchanged, including every historical unknown. At targets
already present in an older marker packet, the entire ordinary native
record, group table, graph and history must agree exactly.

One marker per target/type is tested against every ordinary star,
using either the actual triangular-port exclusion or exact native
product gates. Complete frame automorphisms transport that anchor to
every marker. The source checks preservation of the full ordinary
native data and its owner, and checks the transported anchor's native
covectors and heights. It does not assume the historical Q–Q graph is
automorphism invariant. Old outcomes are simply preserved.

There are13,273 anchor tests and606 retained unknowns. Their exact
transport retains49,268 ordinary/marker pairs. The complete unreduced
domains would have1,146,756 ordinary/new-marker positions,
499,568 old/new-marker positions and919,466 new/new-marker positions.
The relative theorem and exact transports cover these domains without
repeating that pair census.

For each inherited nonnegative ordinary count-vector charge, set to
zero every coordinate whose orbit is compatible with any new or old
marker. Preserve all other coordinates. If a local configuration has
no marker, the old charge is at most1. If it has markers, every
ordinary star with a nonzero remaining coordinate is excluded, and
the joint marker charge is at most1 by the boundary theorem. Thus
their sum costs one target vertex. This argument uses no second
triangle charge and does not spend an old Q–Q budget twice.

Ordinary face coefficients use their full saved actual-corner
multiplicities. A112 and B112 each add their marker coefficient once
to their source aggregate face, by the recovery argument above. The
coefficient is not multiplied by the three quads or by the number of
native orientations. Duplicates are removed within this export only;
529 exports are not asserted to be529 new globally unmatched variables.

## Exact fixed-mass result and verification scope

The face masses remain those of the frozen4624-positive-face
certificate. Starting from the raw108 repaired vertex multipliers,
set each affected multiplier to the maximum of its old value and all
new column loads at that vertex. Exactly three targets increase:

| Target vertex | Increase |
|---|---:|
| 279 | 11744583/100000000000 |
| 495 | 3620121/2000000000 |
| 910 | 829263/625000000 |

The new normalizer is1224897727063181/1188000000000000. It exceeds
the target12 normalizer by2286728814011/1188000000000000. This is the
least monotone repair for these exports, not a minimization over all
enlarged local cones or unrestricted vertex multipliers.

All3,299 previous family columns and all14 full multi-owner costs are
explicitly checked at the repaired multipliers. The inherited full
matrix remains feasible by its established sign convention: every
nonzero vertex-budget coefficient is positive. Keeping face masses
fixed leaves every zero-budget left side unchanged; increasing vertex
allowances enlarges every other right side, including those with
signed face loads. No historical column or ownership is deleted.

The independent helper imports no producer. It implements ordered
Q(sqrt5) and multiquadratic arithmetic, independently reconstructs the
contacted-block scaled Gram inverse, checks all zero contacts on unused
orthogonal passive blocks, and verifies native products and heights.
It reconstructs the complete relative signatures and new native maps,
owner-preserving frame transports, boundary polytopes and literal
clique witnesses, saved ordinary charge arithmetic, source coefficients
and exact monotone repair. Historical catalogue completeness, ordinary
graph/vector completeness and owner chronology are inherited inputs;
the old global matrix is not reconstructed by this helper.

Generation session83347/chunkcc1a25, separate fresh replay
session69187/chunkb92335 and independent saved-checker
session49955/chunk5dea78 all exited0/PASS. The checker verified20,518
cross-product expressions in addition to the complete native-map,
transport, envelope and charge checks described above.
The early relative pilot and an unfrozen export bookkeeping failure
were superseded by this complete generation; no frozen input changed.
Root completed a full producer/report/helper mathematical review: PASS.
The read checked actual finite supports and compact simplex recovery,
one occurrence per extra vertex, the complete relative domain and
contacted I2(10) field choice, all seven-type matching branches including
the optional center and raw76 zero case, nonnegative one-allowance
charges, preservation of native owners and old Q–Q outcomes, all prior
and multi-owner costs, and the precise fixed-mass interpretation below12.
That review did not execute computations or replay historical inputs.
An independent full producer/report/standalone-helper mathematical read
by order16_fractional also passed, including rereads of the inherited
native-map, full-frame automorphism, current-domain transport and
contacted I2(10) projection interfaces. It checked the actual incidence
and source multiplicity, complete relative roles, seven-type boundary
and literal witnesses, rational coefficient vertices, owner-preserving
transport without Q–Q symmetry, ordinary-coordinate omission and all
prior costs. The helper's independent active/passive block, projection,
domain, transport, envelope and charge scope was also checked.
That reviewer did not execute pair tests, censuses, coefficient-polytope
or charge arithmetic, or a global-matrix replay.

Final integrity validation checks50 direct source bindings, five local
links and both Python ASTs. The packet is final/frozen. No global
optimization or new dimension conclusion follows from this freeze.
