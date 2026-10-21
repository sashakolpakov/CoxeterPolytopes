# Scope re-audit of the selected disjoint-Lanner obstructions

Status: supporting hostile audit for the re-audited seven-selected-nerve
certificate.  The new aggregate freezes the symmetry/interlacing conclusions
below and supersedes both the original digest and the provisional two-nerve
demotion.

## Conclusion

The compact-Coxeter principal-block obstruction is unrestricted for each of
the selected nerves 30, 42, 43, 44, and 47.  The blanket statement that an
outside dotted or Lanner component can "mediate" between the selected
supports is not relevant to this obstruction: outside entries cannot alter
the inertia of a fixed principal Gram block, and they cannot repair a
non-elliptic principal subdiagram of an actual compact vertex diagram.

This conclusion is only about the displayed selected nerves and only about
compact hyperbolic Coxeter realizations.  It is not a whole-schedule theorem
for records 30 or 42, and it does not treat ideal vertices or finite volume.

## General certificate logic

Let `L` and `R` be disjoint minimal nonfaces.  In a compact Coxeter
realization their induced diagrams are Lanner diagrams.  Every pair in
`L x R` used below occurs together in an actual rank-nine vertex, so its
entry is a finite Coxeter label, not a dotted entry.  The census covers all
labelled Lanner placements, either literally or modulo a subgroup verified
on the complete minimal-nonface clutter, and every cross pair.

For each fixed placement pair, the census tests a single non-right cross
edge first at `m=3`.  Perron--Frobenius monotonicity for the nonnegative
Coxeter adjacency matrix makes this the extremal test.  The surviving
single edges are also tested at `m=4,5`; failure at `m=4` rules out every
higher finite label.  Every pair of individually surviving `m=3` edges is
tested; monotonicity then rules out every configuration with at least two
non-right cross edges.

With no cross edge, `G[L union R]` is the orthogonal sum of two Lanner
blocks and has negative index two.  Every allowed one-edge block is placed
in an exact weighted-graph class whose leading-principal-minor signs give
negative index two.  Cauchy interlacing therefore contradicts the one
negative direction of an ambient hyperbolic Gram matrix.  Labels on nodes
outside `L union R` cannot change this principal submatrix.

Likewise, when a trial fails an actual-vertex restriction, that restriction
is a principal subdiagram of the required positive-definite compact vertex
diagram.  Adding labels outside the restriction cannot make it elliptic.

For records 43, 44, and 47 one may fix any one audited right A5 support.  The
zero/one/at-least-two split covers every cross entry inside its ten-node
principal block.  Entries with an endpoint outside that block cannot alter
it; if another right support overlaps it, an edge whose endpoints both lie
in the fixed block is already among the audited cross entries.  Thus overlap
of the residual supports creates no missing case.

## Exact record audit

| record | decisive supports | verified normalization | cross census | terminal block |
|---|---|---|---|---|
| 30 | A4 `{3,4,5,6}` and A4 `{7,8,9,10}` | full `S4 x S4`, order 576; all 108 A4 placements are covered by 9 orbits on each side | `9 x 9 x 16 = 1296` label-three trials; zero survivors | always Lanner direct sum, negative index 2 |
| 42 | A4 `{2,3,4,5}` and A5 `{6,7,8,9,10}` | full `S4 x S5`, order 2880; 9 A4 and 5 A5 placement orbits cover 108 and 420 labelled placements | `9 x 5 x 20 = 900` label-three trials; two single-edge survivors; no label 4 or 5 survivor; their only residual double-edge trial fails | exact inertia `(7,2,0)` |
| 43, first residual nerve | A5 `{3,4,5,6,7}` and A5 `{1,2,8,9,10}` | full `S5` on the left support, order 120; its 5 orbits cover all 420 placements; all 420 right placements are literal | `5 x 420 x 25 = 52500` label-three trials; 300 single-edge survivors; no label 4 or 5 or double-edge survivor | two exact classes, each inertia `(8,2,0)` |
| 44 | A5 `{2,3,4,5,6}` and A5 `{1,7,8,9,10}` | same full-left-`S5` coverage | same 52,500/300/0 census | two exact classes, each inertia `(8,2,0)` |
| 47 | A5 `{3,4,5,6,7}` and A5 `{1,2,8,9,10}` | same full-left-`S5` coverage | same 52,500/300/0 census | two exact classes, each inertia `(8,2,0)` |

For order four, the 108 labelled Lanner placements split into nine `S4`
orbits, of sizes `6,6,12,12,12,12,12,12,24`.  For order five, the 420
placements split into five `S5` orbits, of sizes `60,60,60,120,120`.  The
canonical placement lists used by the diagnostics represent exactly these
orbits.

Record 43's second residual nerve has a separate, already standalone audit:
it uses the verified distinguished `S4 x S4` action, 20 placement orbits on
each support, all 25 cross pairs, and again ends with exact inertia
`(8,2,0)`.  Thus both residual record-43 nerves have compact principal-block
obstructions once the first nerve's full-`S5` justification above is made
explicit.

## Artifact status

The selected-nerve aggregate has been rebuilt with the symmetry orders and
principal-block scope lemma explicit.  The downstream whole-shadow JSON
certifiers for records 42, 43, 44, and 47 must be repointed to that new digest
and replayed before their old frozen hashes are cited again.
