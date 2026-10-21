# Complete raw108 extra-vertex marker family

This is a mass-independent coupled-constraint extension for the current
compact arithmetic hyperbolic 13-space model over Q(sqrt5). It covers all 196 raw108
spectator classes and both extra finite vertices of each three-quadrilateral
aggregate. It exports 415 valid one-vertex columns. The least monotone repair of
the preserved face masses has exact ratio

    14671331978990040 / 1221031586432741
      = 12.015522073308947 > 12.

The previous three-family ratio was12.018969384494238. Eight new columns
are strictly violated and eight vertex allowances rise. This is an exported
column repair, not an optimum of the enlarged cones or a global solve.
The above12 witness survives; no dimension or physical-realizability claim
follows. Every old outcome and all 14 multi-owner costs are preserved.

## Inputs and complete domain

The [producer](../scripts/fresh_raw108_complete_marker_family_20261018.py) and
[exact certificate](fresh_raw108_complete_marker_family_20261018.json) bind
the frozen [four-core pilot](fresh_raw14_raw112_raw77_raw108_core_incidence_pilot_20261018.md),
which proved the raw108 eight-wall geometry and enumerated its complete spectator
family. The current
[raw56/57/76 joint certificate](fresh_raw56_raw57_raw76_joint_marker_family_20261018.md)
supplies the prior marker constraints and fixed-mass allowance vector.
The ordinary domain remains the complete current327-class,600-forest domain;
new raw108 events are markers, not an assertion of a new ordinary-star census.

Raw108 is canonical core88, has trivial native automorphism group, and has
three quadrilateral partitions. Its six-normal finite supports obtained by
omitting two of the eight normals are exactly

    01,03,06,12,23,24,25,26.

The union of the three quadrilaterals uses six of these supports. The extra
supports omit24 and25. Denote their marker types by A=1084 and B=1085:

| Type | Omitted normals | Active finite forest | Role orientations |
|---|---|---|---|
| A | 2,4 | H4 + I2(5) | 2*h4*i25 |
| B | 2,5 | A1 + H3 + I2(5) | 2*a1*h3*i25 |

Here a1,h3,h4,i25 count whole finite components in the target forest. H3 and
H4 have trivial label-preserving diagram automorphisms; each I2(5) has two
orientations. The remaining spectator has rank 7. The complete 196 source
classes give 392 pointed target domains, 383 distinct target forests, and 1,676
native markers: 526 of type A and 1,150 of type B. Each source class has one
occurrence of each type, irrespective of the number of orientation maps.
These are complete family counts, not a positive-mass selection.

## Actual geometry and counting

The inherited raw108 native Gram has physical inertia(6,1,1) and conjugate
inertia(7,0,1). A positive-definite support together with its orthogonal rank 7
spectator consists of13 actual walls and yields an actual finite vertex in
any realizing acute Coxeter polytope. To see why arbitrary additional facets
cannot truncate it, choose a strict timelike interior center h, put
b_i=h.n_i<0 and c=G^{-1}b<=0 for the acute positive-definite Gram G, and set
x=h-sum(c_i n_i). Then x is future timelike, orthogonal to the support, and
x.m=h.m-sum(c_i(n_i.m))<=h.m<0 for every other wall m. The support therefore
meets at an actual simple finite vertex, not merely a formal Gram subset.

The fixed walls0,7 and the passive spectator define an actual compact
4-simplex with sides1,3,4,5,6; wall2 is the excluded wall. The saved residual
Gram has Lorentzian inertia and strictly negative inverse, and all five
candidate vertices are actual by the preceding argument. The polytope's
intersection with the fixed support lies in this simplex and contains its
five vertices; convexity gives equality. At either new vertex, every pair
among its four incident simplex sides therefore belongs to a compact
triangular2-face. The six resulting port pairs cannot be occupied by an
ordinary quadrilateral star at that vertex. This is a compatibility
exclusion and introduces no extra triangle charge.

For every complete active-component map, the producer reconstructs both
missing full-frame contact covectors, the positive timelike heights and the
prescribed omitted-wall product. A full positive-definite finite frame fixes
the spacelike projection of a normalized wall; its positive height and
common future sign uniquely determine the remaining timelike coordinate.
Thus a marker recovers both omitted actual walls and the full eight-wall
aggregate. The connected nonfinite core is distinguished from the passive
finite components. Its trivial native automorphism group prevents a further
core-role quotient, and duplicate-free unordered native covectors explicitly
check uniqueness in each complete map domain.

The source object is the full aggregate of three quadrilaterals. A marker
recovers that object, not three independently counted pointed quad records.
Each aggregate produces exactly one type-A and one type-B extra occurrence.
The two extra vertices are distinct supports and are counted separately,
even when two target descriptions have the same forest type. Accordingly,
each marker term has coefficient one times its selected charge weight on
the source aggregate face. Neither the source's three-quad demand nor the
number of orientations multiplies this coefficient.

## Complete relative-role reduction

The producer constructs finite rank 13 hosts large enough to contain every
relative role pattern, enumerates all whole-component embeddings, and
performs32 new exact pair tests:

| Pair type | Complete relative cases | Retained cases |
|---|---:|---:|
| A,A | 5 | 5 |
| B,B | 11 | 2 |
| A,B | 3 | 3 |
| A versus56,57,76 | 3 | 2 |
| B versus56,57,76 | 10 | 1 |

For I2(5), the relative choices are same orientation, reversed orientation,
or a different whole component. Equal/different H4 gives the A,A cases;
equal/different A1 and H3 gives the B,B cases. The identical marker is
removed. The A,B pair only shares a possible I2(5) component. Against old
families, only B can share A1 or H3: its A1 is the first, second, or neither
old ordered A1 role, and its H3 is the same or different where applicable.
These choices give exactly the listed domains. Orthogonal uncontacted
components affect none of the products. Complete native-frame transports
make this reduction independent of spectator size and labeling.

Every pair test uses a fresh immutable multiquadratic field for its two
external walls, strict positive heights, and the exact native projections.
A positive product, product in the forbidden(-1,0) gap, or failed necessary
arithmetic square condition excludes the pair. All other outcomes remain
unknown and compatible. No retained pair is asserted physically realizable.

The complete conservative marker graph has the following capacities.
All distinct A markers are mutually retained, so cA=2*h4*i25. A B clique
must have one common A1 role and one common H3 role; it can use both
orientations of every I2(5), so cB=2*i25. All A,B pairs are retained.
A is also compatible with all raw56 and raw76 markers. B is compatible
with a raw56 marker exactly when its A1 is distinct from both ordered
raw56 A1 roles; every B/raw76 pair is excluded. Raw57 is excluded by both
new types and remains excluded by old56/76.

Write a,h,t for the numbers of A1,H3,I2(10) components. In the coordinate
order(raw56,raw57,raw76,A,B), a complete upper boundary for marker-clique
counts is

    (floor(a/2)*t, 0, 0, cA, 0),
    (floor((a-k-1)/2)*t, 0, k, cA, 0), 1<=k<=min(a-1,h),
    (floor((a-1)/2)*t, 0, 0, cA, cB),
    (0,t,0,0,0).

Absent types give zero coordinates and remove their corresponding
branches. A raw76 clique consumes one common first A1 and k distinct
second A1 roles, leaving the first matching bound for raw56. A B clique
consumes its common A1 role and excludes raw76, leaving the second
matching bound. A can always use its full capacity alongside the other
non-raw57 markers. The raw57 branch is exclusive. These cases cover
all cliques, and every displayed row has a saved literal abstract clique
witness. They do not assert physical realizability.

The producer enumerates every vertex of the nonnegative coefficient
polytope given by these count inequalities, retaining its coordinatewise
maximal vertices. The independent checker separately reconstructs this
finite intersection calculation and checks all saved graph witnesses.
The resulting mixed coefficients use one allowance; the code does not
add independent raw108 and old-family allowances.

## Ordinary stars, owners, and coupled columns

All383 target frames are reconstructed. There are 6,687 ordinary stars in
these records. One anchor marker per pointed family requires 6,900 exact
ordinary tests instead of 60,784 ordinary/new-marker positions. The actual
full-frame automorphism from the anchor to every marker is saved. It
permutes all ordinary native keys and preserves every native field and
actual owner orbit. This does not assume that the historically tested old
Q-Q graph is invariant under that permutation; that graph is retained
literally. The transports yield 4,878 retained ordinary/marker pairs from 580
anchor unknowns. All are allowed. The 32 marker-relative tests cover 11,780
new/old-marker and 5,254 new/new-marker positions without an unreduced census.

For a saved nonnegative ordinary charge, zero every orbit compatible with
any old or new marker. This only weakens its inequality on the complete
old count-vector domain. With no marker present, the old ordinary inequality
bounds the charge by one. With a marker present, every positive-weight
ordinary star is excluded, and the mixed marker envelope
bounds the charge by one. Thus the sum is a single vertex-budget inequality,
not a sum of separately available allowances. Source coefficients are
reconstructed from complete actual-origin multiplicities, with the new
once-per-aggregate terms added as above. The 415 columns are deduplicated within this export, not against every
old global variable. They are valid for arbitrary nonnegative face masses; saved attaining ordinary
mixtures at earlier masses are not transported.

The exported repair retains all 2,884 prior marker columns and all 14 full
multi-owner costs. It fixes every face multiplier and sets each affected
vertex multiplier to the maximum of its prior value and its new column
loads. All full-model nonzero vertex-budget coefficients are positive.
Hence raising these multipliers enlarges every budgeted right-hand side,
including columns with signed face coefficients, while zero-budget left
sides remain unchanged. This is the inherited full-matrix monotonicity
argument; the present packet does not rerun that matrix or its optimizer.

The remaining exact normalizer gap to ratio12 is

    1579411816429 / 1188000000000000.

The retained above12 witness rules out a target12 primal certificate from
these exported constraints at this stage. No claim is made about the optimum
after further stronger constraints, complete enlarged cones, or a new global
solve. Raw112 is not included in this extension.

## Validation and scope

Corrected source-bound generation56256/chunkdba596 completed0/PASS.
The independent checker caught and corrected a draft-only host-labeling
mismatch: old saved marker frames use canonical sorted component order,
so every temporary relative-test host is now sorted and checked against
the exact saved frame before testing. The corrected cross-family outcomes
and mixed envelope above replace the draft exclusivity assumption. No
frozen predecessor changed. Separate fresh replay56955/chunkbca986 and
independent saved checker14072/chunkdebb6e both completed0/PASS. The latter
checked 32 relative cases, 1,676 maps, 6,900 anchor tests and 10,132 native
expressions, together with the 415 columns and full repair accounting.
Final integrity checks rehash all 47 direct source bindings, parse both
Python files and check the five local report links. These checks are
integrity-only and do not add a mathematical replay.

Root completed a full producer/report/standalone-helper mathematical
read: PASS. Its scope included the complete sorted32-case role domain,
full maps/heights/recovery and simplex incidence, conservative mixed
count boundary and clique witnesses, nonnegative coefficient vertices,
ordinary transports and owner invariance, omission union, once-per-
aggregate coefficients, and preservation of the2884 prior columns and14
multi-owner costs under monotone repair. Root did not execute the packet
or rerun the historical model.

The independent order16 peer completed a full producer/report/helper
mathematical read: PASS. It checked canonical sorted hosts, complete
relative signatures and component transports, independent field arithmetic
and passive zero-contact guards, mixed count bounds and graph witnesses,
ordinary-orbit zeroing with one allowance, source multiplicity and old
column/multi-owner preservation, together with actual finite-vertex/simplex
geometry and the fixed-mass scope. It did not execute pair tests, maps,
censuses, coefficient-polytope enumeration or historical matrix replay.

The standalone [saved checker](../scripts/check_raw108_complete_marker_saved_20261018.py)
uses no producer imports. It implements independent ordered Q(sqrt5) and
multiquadratic arithmetic, positive-definite support tests, contacted-block
inversion, complete component-role maps, all 32 relative gates, all 1,676
native marker recoveries, ordinary anchor gates and owner transports,
exact charge accounting and repair. Uncontacted passive blocks are replaced
by harmless diagonal blocks only after explicit zero-contact checks. The
historical catalogue, complete source-family enumeration, ordinary owner
chronology and old count-vector completeness remain frozen inputs. Actual
face existence and multiplicity are the mathematical argument above,
not conclusions of the arithmetic checker alone.
