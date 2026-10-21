# Raw88/raw25 forces a new finite vertex, but its event has exact zero-cost covers

This is a completed local diagnostic in the compact arithmetic
\(\mathbb H^{13}\), \(K=\mathbb Q(\sqrt5)\) model. A particular pair of actual
raw88/raw25 prisms forces two finite vertices outside both prisms. Its common
intersection is a tetrahedron, and a complete finite-frame recovery at one of
the new vertices gives a capacity of twelve pointed events per target vertex.
Nevertheless, exact covers at both current face-mass choices avoid this event
altogether. The full ratio remains **11.99632388627268**, and the named248-Q
restriction remains **12.05165963034431**. There is no new dimension bound,
global optimization, pair exclusion, or change to a frozen owner outcome.

The [producer](../scripts/fresh_raw88_raw25_forced_vertex_frontier_20261018.py)
and [saved data](fresh_raw88_raw25_forced_vertex_frontier_20261018.json) bind the
[previous tetrahedron witness](fresh_shared_prism_tetrahedron_incidence_20261017.md),
the current literal-star machinery, and the inherited actual-vertex proof.
The [independent saved-data checker](../scripts/check_raw88_raw25_forced_vertex_saved_20261018.py)
checks the saved matrices, graph/event domain and rational covers without
importing project arithmetic.

## Literal geometry and actual faces

At source904, the full finite forest is
\(3A_1+A_2+2H_3+I_2(5)\). The two stars belong to current face45683 (raw88)
and face26832 (raw25), each with its unique historical provenance. Their
occupied-port complements are checked literally; they are not inferred from
an abstract graph isomorphism. Both are four-dimensional prisms with eight
finite corners.

After removing the common passive \(A_1+2H_3\), label the six source walls
\(0,\ldots,5\), and the three external walls \(6,7,8\). Raw88 uses walls6,8
and has fixed active support1,2; raw25 uses walls7,8 and has fixed support0,2.
The nonright products are the simple bonds
\(16,23,36,37,38,58\), the label5 bonds45,67, and the two dotted products
\(g_{06}=-\sqrt{5+2\sqrt5}\), \(g_{17}=-(1+\sqrt5)\).
Every unlisted product is zero. Scaling wall0 by
\(\sqrt{5+2\sqrt5}\) puts the complete Gram matrix inK. Its inertia is
\((6,1,2)\). Full future-side contact reconstruction at every transported
source verifies this same nine-wall configuration, including the three
external cross products and their strictly positive squared heights.

All84 six-element supports are examined exactly. Fourteen are positive
definite. Each prism has exactly eight such corners; four belong to both,
and two belong to neither. The two outside vertices have forests1325 and2252.
The former happens to have a previously used *forest type*; no old event or
capacity is identified with this new occurrence.

A positive-definite support here gives an actual finite vertex even with
arbitrary additional acute facets. Ifh is a strict future timelike interior
point andG is the selected positive-definite acute Gram, write
\(b_i=h\cdot n_i<0\), \(c=G^{-1}b\le0\), and
\(x=h-\sum c_i n_i\). Then \(x^2=h^2-b^TG^{-1}b<0\), and every other actual
wallm satisfies \(x\cdot m\le h\cdot m<0\). The full rank13 support therefore
cuts out exactly one actual finite vertex. This is the same strict-center
argument used in the frozen1325 packet; it is applied to actual walls, not
merely to a positive-definite abstract subdiagram.

The intersection of the two prisms fixes passive7 plus walls0,1,2, hence has
dimension3. Its four residual sides are3,4,5,8. Their Schur complement has
inertia(3,1), and its inverse is strictly negative. The four triple-side
supports are the actual vertices

| Forest | Active support |
|---|---|
|904|0,1,2,3,4,5|
|534|0,1,2,3,4,8|
|1022|0,1,2,3,5,8|
|690|0,1,2,4,5,8|

These sides bound a compact tetrahedron. The actual face is contained in that
simplex, and contains its four actual vertices; convexity gives equality.
Extra facets cannot truncate this tetrahedron. The same argument and the
saved eight finite corners identify each transported compact prism with the
same actual face. At a finite simple ambient vertex its side/support walls
are the incident walls, so a transport cannot substitute another outward
wall having only the same face trace.

Complete whole-component embeddings yield respectively24,48,12,72 literal
prism-pair records in these four forest frames. Every record verifies the
full native covectors, future heights, omitted product, face identity and
occupied-port support. All pair keys are distinct and their number equals
the number of embeddings. Conversely, actual outward-normal uniqueness and
the full-star identity recover the same prism pair at every tetrahedron
vertex. Thus an unlabelled actual pair is counted once at each of the four
distinct vertex types, with no orientation multiplier:
\[
X_{904}=X_{534}=X_{1022}=X_{690}=X.
\]

## Target-frame recovery and its limited arithmetic effect

The extra vertex2252 has active support2,4,5,6,7,8 and full forest
\(2A_1+3H_3+I_2(5)\). Select one of its twoA1 factors, one of its threeH3
factors, and an ordering of the I2(5) endpoints. This gives exactly
\(2\cdot3\cdot2=12\) markers. The markedH3 has its native label5 end4,
center5 and simple end8; its internal roles have no further symmetry.
The otherA1 and twoH3 components are passive and are not permuted as extra
markers.

For each marker, all three absent walls0,1,3 have complete contact covectors
in the full thirteen-wall target frame. Wall0 contacts only6, with squared
contact \(5+2\sqrt5\). Wall1 contacts6 by3 and7 with squared contact
\(6+2\sqrt5\). Wall3 contacts2,6,7,8 by3. Their squared timelike heights are,
respectively,
\[
5+3\sqrt5,\qquad11+5\sqrt5,\qquad(11+5\sqrt5)/4,
\]
all strictly positive. The three recovered walls are pairwise orthogonal.
All twelve complete maps, heights and cross-product checks are saved.

A full covector fixes the spacelike projection. Its positive height and the
common outward future-side convention fix the remaining timelike
coefficient uniquely. Thus each target marker determines at most one set of
three missing actual walls, one source904 vertex, and the two specified full
prism stars. Distinctness is explicit in the covector keys; the common
normal set and fixed supports determine the actual source pair. Since that
pair has exactly one target2252 support, this is an injection, giving
\(X\le12v_{2252}\). This is a local necessary incidence statement; no claim
is made that any prescribed target marker is actually realized.

The arithmetic obstruction is stronger than failure to violate that quota.
For both fixed face-mass choices, the saved local covers satisfy **X=0** at
all four source types. Even forbidding this entire event family leaves the
four local graph-cover costs within their existing vertex multipliers.
Reducing the factor12, or adding a second target allowance for these same
events, therefore cannot improve these fixed witnesses through this local
event subsystem.

## Exact cover certificate and scope

The four inherited conservative compatibility graphs are reconstructed with
all chronological owners and unknown outcomes unchanged. Their14,440
maximal cliques are exhausted, and all nonempty subsets give252 distinct
count/event vectors. Counting all subsets avoids any assumption that adding
a star preserves the event count. Every event value is0 or1. Coordinates
are actual corner-orbit occurrence counts, not the numbers of orientations.
The two demands are rebuilt separately from full and restricted face
multipliers.

Two small changed local LPs propose covers and charges. Exact rational
recovery checks every coordinate demand, cover mass, common event equality,
quota, charge inequality, and matching objective. There are58 positive
mixture terms across the two certificates. The minimum additional source
mass is exactly zero, with zero dual charges, zero event mass, and zero
increments. The current face and vertex multipliers are copied unchanged:
\[
R_{\rm full}=\frac{111428537994904130}{9288556982227791}<12,
\qquad
R_{\rm restricted}=\frac{18670750401222348}{1549226494433359}>12.
\]
The latter remains an obstruction to target12 in the represented relaxation.
All439,410 inherited model columns, earlier incidence columns and the two
owned unknown outcomes at1763 remain valid because neither multiplier vector
nor any owner outcome changes.

These certificates establish survival of the two saved witnesses under
this specific event subsystem. They do not construct a physical chamber or
prove simultaneous realizability of all local covers. They also do not
assert that all possible raw88/raw25 configurations, or all forced-vertex
constraints, have been exhausted.

## Verification

Final primary generation68248 and separate fresh saved replay43062 completed with
terminal exit0/PASS. The replay reconstructs the finite supports and source
maps, enumerates the same graphs, and checks saved rational certificates
without an optimizer. The independent stdlib-only saved checker97245 also
completed with terminal exit0/PASS: it independently checks all84 saved
support positivity tests, the tetrahedron residual/inverse, all168 saved
source/target frame height and cross-product identities, the14,440 maximal
cliques and252 subset/event vectors,504 charge inequalities,58 positive
cover terms and both unchanged witness ratios. It validates saved matrices
and contacts; it does not regenerate historical provenance, whole-component
orientation completeness or chronological pair ownership.

Root and a separate peer both completed full producer/source and mathematical
report reviews: actual PD-support vertices, compact tetrahedron/full-star
transports, selected-pattern target recovery, and exact zero-event cover
scope. Neither review reran native geometry, graph enumeration or the local
LPs. The peer additionally read the complete independent helper, including
its orderedK arithmetic, saved-frame identities, clique enumeration and
rational cover checks. A final schema guard explicitly requires all four
source mixtures in both producer and helper; the peer read that guard, and
all final runs above include it. The helper review did not execute its
arithmetic or regenerate canonical embeddings, provenance or owner history.
The root review did not include the helper.

All source/data/report files are final and frozen after these checks. No
additional geometric family, global solve, or optimistic-screen conclusion
is included in this packet.
