# All remaining 89 raw4 classes: exact pair-only benefit barrier

This is a fixed-witness diagnostic for the current compact arithmetic
hyperbolic 13-dimensional model over Q(sqrt5). It leaves the authoritative
324-class, 600-forest domain and every old graph outcome unchanged. It proves
no new physical exclusion, packing optimum, realization, or dimension bound.

The complete remaining domain consists of 89 positive raw4 face classes:
18 have explicit current provenance and 71 have a unique saved context.
Their 534 actual corners meet 238 forest types, 183 already represented and 55
additional types. The current 5121 face masses are used throughout; historical
provenance masses are not reused.

The exact singleton-cover bound permits a ratio below 12 only for batches of
at least 42 classes. Even all 89 give a lower bound only slightly below 12:

| Domain | Largest permitted multiplier increase | Ratio lower bound |
|---|---:|---:|
| Best batch of 41 |126119770519/12375000000000|41113829770090560/3426149996940161 =12.000008699797|
| Best batch of 42 |126163919569/12375000000000|41113829770090560/3426161864204801 =11.999967135129|
| All 89 |14144447719463/1386000000000000|8222765954018112/685239135430273 =11.999848708079|

These are lower bounds on the repaired witness ratio, not attained packing
ratios. Passing below 12 in this table is only permission left by an upper
bound on possible improvement. It proves neither incompatibility of any new
pair nor a repair reaching 12.

## Actual corners and singleton extension

The five-node raw4 core has face obligation 3 and six actual corners. These
quantities are checked separately. Each corner retains a three-node finite
subdiagram; its two omitted walls have complete native contact covectors.
The producer reconstructs all six small templates, their finite component
types and automorphisms, exact spacelike projections, strictly positive
squared timelike heights, and the native omitted-wall product with a common
future-side sign. Adding the rank 10 orthogonal spectator changes none of
these identities.

For each of the 89 classes, the producer checks the current explicit
provenance or unique source context, core identity, spectator rank, and all
six actual forest occurrences. Whole-component embedding counts are obtained
from component multiplicities and template automorphisms. They bound the
number of possible new native stars; no complete embedding, native-overlap,
or new-orbit census is claimed.

Let alpha_f be the fixed face mass, m_v the current vertex multiplier, and
t_v the saved exact local covering optimum of forest v. Put t_v=0 where no
current graph exists, while retaining any current m_v. Let c_fv count actual
corners, and for a selected set S write

    e_v(S) = sum_(f in S) c_fv alpha_f.

The latest 600 covers are selected chronologically from the 290-class local
reconstruction, the 32-class historical extension, and the two new raw4
classes. Their exact masses obey t_v<=m_v. Every old-old graph outcome,
including every owned unknown, is retained.

An old cover extended by singleton mass for each added actual occurrence
covers the enlarged demands with total mass at most t_v+e_v(S). New native
stars may coincide with old stars or each other, and new occurrences may
merge into the same complete orbit. Such coincidences can only improve this
upper bound: the separate singleton contributions still cover every added
marginal. New-node pair exclusions cannot delete singleton configurations.
Consequently a feasible monotone repair satisfies

    m'_v <= max(m_v, t_v + e_v(S)),
    G(S) = sum_v max(0, t_v + e_v(S) - m_v).

Here the inequality means that a repair of at most that mass exists; it does
not constrain arbitrary unnecessarily enlarged choices of m'_v. With fixed
numerator N and baseline denominator M, this gives the surviving witness
ratio lower bound N/(M+G(S)). Every inherited signed, zero-budget and
multi-vertex column remains satisfied by this monotone repair.

The baseline is

    N =679792158897/55000000000,
    M =16961245013123269/16632000000000000,
    N/M =205569148850452800/16961245013123269.

The mass increase needed to reach 12 is
169517391081131/16632000000000000. For all 89, the bound exceeds that amount by
only 8639262097/665280000000000, approximately 0.0000129859038255.
A verified cover saving larger than this amount, measured after existing
multiplier slack is applied, would restore a strictly above 12 witness even
for the full 89-class pair-only extension.

## Complete cardinality calculation and workload

Connect two face classes when their actual corner domains share a forest.
The resulting face graph has 64 components, the largest containing 7 classes.
Components have disjoint forest supports, so G is additive across them.
The producer enumerates all 318 subsets within these components, saves each
exact gain, and convolves their cardinality maxima. The saved certificates
include every intermediate cardinality upper bound and an attaining subset.
Thus the 42 threshold concerns every subset of the 89 classes, rather than
only the mass-ranked prefix. The exact mass-ranked prefix first permits a
ratio below 12 at 45 classes. Its first four reproduce the frozen next-four
bound exactly.

The 42-class attaining subset is saved in the JSON. Its conservative workload
bound is 220350 new pair positions; the all 89 upper bound is 332054. These
numbers use whole-component embedding counts before orientation/native
identifications and can overcount substantially. They are not actual pair
domains or predicted clique counts. No pair census or covering LP is needed
for this barrier.

## Validation and limits

Generation 74952/chunk 09ab3f exited 0/PASS. Separate fresh producer replay 21231/chunk b08afb and independent saved
checker 85557/chunk 3d871c both exited 0/PASS. The checker verified 7,065
cardinality upper inequalities in addition to all 318 component subsets. The producer rebuilds the
six native templates and all 89 actual-corner interfaces; it does not rerun
old physical geometry, graph cliques, or covering LPs.

The parent completed a full producer, report, and saved-helper source and
mathematical review. It checked provenance and template interfaces, current
cover/slack transport, actual occurrence accounting, singleton feasibility
under orbit overlap, component factorization, complete cardinality
upper/attainment logic, and the stated workload and threshold scope. That
review did not repeat projections, embeddings, arithmetic, or the cardinality
calculation. A second peer review is not claimed.

The standalone checker uses only the Python standard library. It verifies
bound source hashes, current masses and provenance, actual-corner sums,
embedding-count bounds, latest saved cover references, independent
union-find components, all 318 subsets, every cardinality upper inequality
and attaining subset, all 89 prefix values, and exact ratio arithmetic.
Native template projections, complete inherited graph domains, and old
cover feasibility remain explicit theorem inputs. It does not independently
classify finite forests or reconstruct native embeddings.

This barrier applies only to fixed current face masses, preserved old-old
graph outcomes, and additional pair exclusions involving newly introduced
stars. It does not cover constraints excluding singletons, new cross-vertex
incidences, other face classes, or changed face masses. Current represented
counts remain 324 classes and 600 forests. No new optimization of the global
model is performed.

Files: [producer](../scripts/fresh_new5121_raw4_all89_cover_barrier_20261018.py),
[saved exact certificate](fresh_new5121_raw4_all89_cover_barrier_20261018.json),
[standard-library checker](../scripts/check_new5121_raw4_all89_barrier_saved_20261018.py).
Inputs include the [current 324-class covers](fresh_new5121_raw4_two_face_local_covers_20261018.md)
and [frozen next-four domain and ceiling](fresh_new5121_raw4_next_four_domain_20261018.md).
