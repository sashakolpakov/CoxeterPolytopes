# Depth17 suffices for the full high-domain corridor alternative

For every corner of the complete frozen 6,703-context high-overlap domain, the following
statement holds. If the actual whole finite A/B/D component has its two marked roots at
distance 8, and both marks have distance at least 17 from every endpoint, branch and
endpoint of a nonsimple bond, then the full actual supporting star has an opposite-black
path of length at most 7.

This is a **short-path alternative**. Only in the branch with no such path anywhere in
the full supporting star does it exclude the deep corridor. The absence of a path in the
saved old core alone is insufficient. There is no new pair weight, complete payment
system or dimension bound here.

The [source](../scripts/fresh_full_a39_corridor_high_domain_20261017.py) and [exact
certificate](fresh_full_a39_corridor_high_domain_20261017.json) cover all 26,812 corners
of the [inherited high domain](fresh_distance7_high_supporting_frontier_20261017.md),
including all four corners of the 251 formerly unresolved contexts and all earlier
structural outputs. No old principal matrix is required to be nonsingular or Lorentzian.
The proof instead uses a forced positive definite path in the proposed actual
completion.

## The forced path and every possible old embedding

Normalize all actual wall roots to squared norm 2. A feature of a finite A/B/D diagram
is a vertex of degree other than 2, or either endpoint of a bond whose label exceeds 3.
Depth at least 17 at each mark, together with marked distance 8, gives two outward
simple arms of length 17 and the simple path joining the marks. All vertices other than
the two marks are supporting whites. The two omitted black roots do not belong to this
actual corner.

Keep the first 15 outward whites on each arm. Together with the marks and the seven
whites between them, these form an induced A39 diagram. Number its vertices 0 through
38; the left and right marks are at 15 and 23. The unused two outward whites on each
side are retained as geometric inputs for the A17 contact estimates below, but are not
needed in the final matrix. A feature at outward distance exactly 17 is harmless: all
bonds inside the displayed A17 arm are simple.

The 342 old native cores are reconstructed from their frozen diagrams and positive
diagonal gauges. They have order at most 19, so any displayed corner has at most 17
roots. Their gauges are 1 or 2; dividing by the square roots of the gauges puts every
relevant product in Q(sqrt(2)). The source checks all four left/right choices in every
old context and independently reconstructs the displayed component partitions.

A displayed component containing a mark must be a simple path. An existing branch or
nonsimple bond remains a feature in the actual component and is at distance at most 16 <
17 from the mark. An existing cycle cannot occur inside a finite A/B/D tree. The checker
saves one such witness for each rejected nonsimple marked block.

If the two marks are in the same displayed path, that path must already give distance 8:
a tree cannot acquire a second, shorter route. Its embedding in A39 is then fixed once
the two marks are placed at 15 and 23. If the marks lie in different displayed paths,
orient each path from its outward endpoint toward the eventual joining segment. Let
t_L,t_R be the inward distances from the marks to the chosen endpoints. The connecting
segment must contain exactly

    k = 7 - t_L - t_R

additional whites, with k at least 1. The two displayed components are orthogonal, so no
zero-white connection is possible. Every orientation of both paths is enumerated, with
singleton orientations counted once. A connected displayed path has at most 17 roots and
already spans 8 edges between the marks, so its two outward tails together use at most 8
roots. For disconnected marked paths, their combined order is at most 17 and each is
nonempty; either path has at most 16 roots and outward extent at most 15. Thus every old
embedding fits inside the 15 retained outward whites on each side. The source asserts
these distinct bounds before checking every old pair product against the corresponding
A39 entries. This also covers paths joined through old spectator roots: those roots may
occupy positions not selected in the reduced old embedding. Their identities and extra
contacts are not assumed away. The Gram tests use a necessary subset of actual roots
only.

## All-label restrictions leave three white slots per omitted black

Write b for the omitted left black root and d for the omitted right black root. The
full-star no-short7 hypothesis implies b·d=0 and forbids any contact from b to the right
mark, or from d to the left mark. Their products with their own marks are fixed native
products and remain in the calculation.

For a norm-2 A17 arm, the inverse is positive and its diagonal at outward distance j is
j(18-j)/18. The omitted black and all these whites belong to an actual edge diagram,
which is strictly positive definite. Hence:

* At outward distance 1, only nonzero contacts of labels 3 or 4 are possible.
  A label at least 5 has squared magnitude greater than 5/2, and
  (17/18)(5/2)>2.
* At distance 2, only label 3 is possible, since (16/9)·2>2.
* At distances 3 through 15, even a unit contact gives projection greater
  than 2, so no contact is possible.

Additional nonpositive products cannot lower these projection bounds, because the
inverse is positive. Strictness excludes equality as well. The retained path stops at
outward distance 15. Contacts to the unused outer whites at distances 16 and 17 are not
claimed to vanish: those whites supply the A17 inverse bound and are omitted from the
final principal matrix. These statements have no finite label cutoff.

In the inner A7 segment, b contacting a white at distance j at least 2 from the left
mark gives a route to the right mark of length 9-j, at most 7. Only the first inner
white can therefore meet b in the no-short7 branch. Its A7 inverse diagonal is 7/8,
which excludes label 5 and every larger label by the same squared-magnitude comparison.
Contacts to the first two outward whites beyond the opposite mark also give short paths;
more distant contacts are already excluded by the other A17 arm. The argument for d is
the reflection of this argument.

Thus b can contact, apart from its own mark at 15, only positions 14,13,16, with allowed
nonzero labels respectively (3 or 4), 3, (3 or 4). Root d has the reflected slots
24,25,22. The source first checks every fixed old contact. A forbidden one supplies
either a recorded short path, an A17 projection obstruction, or a strict edge label
obstruction. Four orientations are rejected by the last reason; the remaining
fixed-contact exclusions use the former two reasons.

One omitted black cannot contact both consecutive outer whites. The positive A2 inverse
is (1/3)[[2,1],[1,2]], and two unit contacts already have projection square 2, contrary
to strict edge positivity. The remaining slot products are enumerated exactly as
0,1,sqrt(2), respecting every old fixed contact. They are an overestimate of realizable
completions, which suffices for an obstruction.

## Two necessary Gram tests

For a norm-2 simple path of order n, the inverse entries, using indices from 0, are

    J_n[i,j] = (min(i,j)+1)(n-max(i,j))/(n+1).

The source verifies this inverse formula by multiplication for n=2,7,15,17,39. Deleting
the two marked roots from A39 gives three independent white blocks A15+A7+A15. Let J_W
be their inverse. For the positive contact vectors x,y of b,d with the A39 path, strict
positive definiteness of each actual edge requires

    2 - x_W^T J_W x_W > 0,     2 - y_W^T J_W y_W > 0.

Every profile violating either strict inequality is impossible. Otherwise eliminate the
full positive definite A39 block. Because b·d=0, the remaining two-by-two Schur matrix
is

    [ 2 - x^T J_39 x       -x^T J_39 y     ]
    [   -x^T J_39 y      2 - y^T J_39 y    ].

In every surviving enumerated profile its first diagonal is strictly negative and its
determinant strictly positive. It therefore has two negative directions. By congruence
the full Gram matrix would have two negative directions, impossible for actual vectors
in a real Lorentz space. This argument requires no assumption on the old core's
signature. Adding further actual roots or spectator contacts cannot repair a principal
Gram matrix with two negative directions.

All signs and products use the [exact Q(sqrt(2))
helper](../scripts/fresh_quadratic_sqrt2_exact_20261017.py). No floating-point tolerance
or optimizer is used.

## Complete counts and scope

|Corner outcome|Corners|
|---|---:|
|Nonsimple marked block incompatible with depth 17|16,627|
|Connected marks at a fixed distance other than 8|2,187|
|No admissible distance-8 joining orientation|14|
|Complete corridor/contact calculation|7,984|
|Total|26,812|

The last row of cases has 13,644 complete path orientations. Fixed contacts exclude
8,337 orientations. The remaining 5,307 use 379 distinct contact templates, with every
template's multiplicity saved. There are 4,983 distinct joint contact profiles. With
orientation multiplicities restored, 57,582 profiles fail strict actual-edge positivity
and 50,137 give two negative directions. None remains unresolved within this
deep-corridor hypothesis.

This settles the proposed deep-pair alternative throughout the **high-overlap domain**,
including the formerly unresolved 251 contexts. It does not settle those contexts for
arbitrary pair shapes, and does not complete the finite-7 payment theorem. Low/disjoint
cases, the actual retained pair maps and all affected universal payment interfaces must
still be checked before changing a weight or importing a finite budget. No global
conclusion changes. No new finite-volume polytope is constructed.

## Conditional rank65 graph counts

A separate literal BFS count in whole finite A65/B65/D65 diagrams gives 23/22/21 pairs
at distance 8 whose two endpoints both have feature distance at least 17. The source
saves all eligible unordered pairs; the independent helper reconstructs the graphs and
all distances without producer imports. If an additional complete payment argument later
authorized halving exactly those pairs, the inherited same-component costs 486/486/487
would become 949/2,475,953/2, respectively. These are hypothetical arithmetic values,
not an authorized finite payment system or a finite-budget theorem. In particular, this
local high-domain alternative alone does not justify those changes.

## Verification and inherited scope

The uncompressed root pilot covered all 26,812 corners and found the same
orientation-weighted outcomes. The authoritative producer independently reconstructs the
old native matrices and all corner/path embeddings, compresses identical fixed-slot
domains, and saves every exact profile. Authoritative generation 93617 and separate
fresh saved replay 11768 both completed terminal exit0/PASS. No floating-point
tolerance, local optimizer or global solver is used in this packet.

The separate [saved-template
checker](../scripts/check_full_a39_corridor_saved_20261017.py) imports no project
modules. It computes generic rational Gauss–Jordan inverses for A39 and its 37 white
roots, checks the A17/A7 all-label inequalities directly, and decides Q(sqrt2) signs by
rational isolating intervals rather than the producer's squared-norm sign test. It
reconstructs every saved 379-slot domain and 4,983 joint profile, reconciles all 13,644
orientation multiplicities, and checks the optional rank65 graph-count table.
Independent saved checker 44786 completed terminal exit0/PASS; all 12 direct source
bindings were checked. It does not independently reconstruct old native cores,
fixed-contact exclusion witnesses, old-block embeddings, all-label geometry or the
inherited overlap census.

The [frozen A43 predecessor](fresh_full_a43_corridor_high_domain_20261017.md) remains
unchanged and is directly bound. The inherited all-label high-overlap census is a frozen
proof input. This successor reconstructs its native matrices, corner partitions and
complete new corridor/contact domain, but does not repeat the original overlap census.
Root and construction peer each completed a full producer, report and independent-helper
mathematical/source read. Both checked depth 17 feature persistence, the
connected/disconnected embedding bounds, unused outer whites, all-label slots, physical
gauges, strict edge positivity, the index-two contradiction, independent arithmetic
implementation and conditional rank65 counts. Neither review duplicated the 4,983
profile arithmetic or the inherited native/embedding/high-overlap census.

    python3 scripts/fresh_full_a39_corridor_high_domain_20261017.py
    python3 scripts/fresh_full_a39_corridor_high_domain_20261017.py --verify
    python3 scripts/check_full_a39_corridor_saved_20261017.py
