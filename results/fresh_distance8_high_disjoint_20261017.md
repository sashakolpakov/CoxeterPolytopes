# All-label exclusion of disjoint high cores at distance greater than eight

**Complete local computer-assisted proof; primary generation and fresh replay
PASS.** Under the actual quadrilateral positivity conditions and absence of
an adjacent marked path of length at most8, its two minimal indefinite cores
cannot be disjoint if both have orders6 through10. This lemma does not
assert a dimension bound or control overlapping cores and their extensions.
It is supplementary to the new ideal-corner replacement route.

The [source](../scripts/fresh_distance8_high_disjoint_20261017.py) and
[certificate](fresh_distance8_high_disjoint_20261017.json) use exact Fraction
arithmetic, positive diagonal congruences and explicit nonnegative witnesses.
The [independent checker](fresh_distance8_high_disjoint_20261017_review.md)
reconstructs the native catalogue and arithmetic separately.

## Complete domain and inverse gates

Let W be the supporting white roots of an actual quadrilateral, with
opposite black pairs U,V. Normalize each root to square2. The required
conditions are W positive definite, W plus any one black positive definite,
W plus an adjacent black pair positive semidefinite, and the full star
having at most one negative direction. All off-diagonal entries are nonpositive.
The [frozen high-order catalogue](fresh_high_rank_quasilanner_census_20261017.md)
supplies every minimal indefinite diagram of orders6 through10.

Its384 marked contexts have maximum white-to-nearest-mark distance7.
Assuming no adjacent marked path of length at most8, a joining bond must
have white endpoints of distances summing to at least8, hence each distance
is at least1. Every external contact must pass both positive-definite
black-deletion inverse tests. Of2404 ports,754 have both normalized inverse
diagonals strictly below2 and can support a nonzero contact.

These ports yield12135 ordered candidate joins. For each, let x,y be the
inverse diagonals at its two endpoints in the separate white blocks.
A bond of squared magnitude b² requires b²xy<1 for white positive
definiteness. At a black deletion it also requires b² times each endpoint
inverse diagonal to be less than2. Nonnegative inverse entries make these
necessary conditions valid even when other contacts are present.

The simple-bond white test excludes3037 joins, each with a saved
nonnegative vector of square at most0. The other9098 joins admit label3.
Exactly178 also survive the necessary label4 tests. Every remaining label
at least4 is excluded by one of the two displayed inequalities. In those
178 cases, label at least5 has squared magnitude at least
(3+sqrt(5))/2>5/2, and the exact endpoint inverse test excludes it.
Parallel and dotted contacts have still larger magnitude. No numerical
upper bound on labels is imposed.

## Single joins and genuine corner obstructions

The resulting9276 single-join label cases have266 underlying labelled
endpoint types. For262 types, one endpoint has a zero diagonal in the
inverse of its full indefinite core: its deleted principal diagram is
affine, so the cofactor is zero. Minimal indefiniteness plus a finite
black deletion gives full signature(order−1,1,0). Eliminating this full
core across a single rank-one joining bond changes the other core by
zero. The joined matrix therefore has two negative directions, for
every bond magnitude. This argument also handles full cores with an odd
label4 cycle, where no rational square-root gauge exists.

Four remaining types have explicit rational congruences and exact inertia
checks. Altogether9256 single-join label cases have two negative directions.
The remaining20 cases have an adjacent corner with a saved nonnegative
vector of strictly negative square, contradicting positive semidefiniteness.
Strict negativity matters here: a null vector alone would not exclude an
ideal corner. The implementation promotes a zero-pivot witness to a
negative witness only when its product with another coordinate is negative;
the resulting quadratic value is checked exactly.

The earlier exploratory count22 included two joins already excluded by
white positivity. The final20 count is over the surviving single-join domain.

## Every multiple-cross-edge subset

There can be up to seven individually admissible cross edges for a fixed
ordered context pair. The source enumerates every subset of size at least2,
giving10990 subsets. For10987, already replacing every cross bond by a
simple bond produces a nonnegative nonpositive-square witness in an actual
edge principal matrix. Increasing any cross magnitude only decreases its
quadratic value, so that witness excludes every native label assignment.

The three remaining subsets have two edges each. Their individual label
gates leave six assignments in total. For L25 with black pair{0,1}, and
for black pair{1,4}, only(3,3) is possible; both have negative index2.
For black pair{1,2}, the assignments are(3,3),(3,4),(4,3),(4,4).
The first has negative index2; the other three have actual-edge witnesses.
Thus no multiple-edge assignment survives. With no cross edges the two
cores themselves give two negative directions. This completes the local
disjoint-core proof.

## Verification scope

Primary generation65276 and fresh saved replay6403 both terminated exit0
with `PASS_DISTANCE8_ALL_LABEL_HIGH_DISJOINT_EXCLUSION`. The replay recomputes
the full domain and compares it with the saved certificate. The separate
independent generation reconstructs26 native catalogue Grams,142 finite
deletions and49 ideal deletions, and checks all ports, label gates, edge
subsets, witnesses and inertia records without importing producer arithmetic.
Its final binding and replay are recorded in its linked report.

The catalogue classification and actual-face positivity are explicit
inherited mathematical inputs. This local result has no external refereeing
or proof-assistant formalization and makes no construction claim.
