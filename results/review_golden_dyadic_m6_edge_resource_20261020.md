# Outside-author pre-execution review: one actual M6 edge

**PASS for the written reduction and fixed checker source.** No program,
imported mathematical body, enumeration, scientific replay, Lean or TeX
job was executed in this review. Execution remains root's serial task.

|Reviewed input|SHA256|
|---|---|
|`results/golden_dyadic_m6_actual_short_edge_preparation_20261020.md`|`2afd545c3c3166cd07381949b972fb012280b18ed673a090025b6aac30537dd2`|
|`scripts/check_golden_dyadic_m6_actual_edge_20261020.py`|`9e48ab80a337f43925a1be7ee6761b02da5eaa0c60df6dfb649bc640ec5e6949`|
|Pure quadratic helper `check_arithmetic_dn_opposed_fork_models_20261019.py`|`185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9`|

I read both new sources in full, the imported field/inverse definitions
and their guarded entry point, the bounded runner/control code, and the
P6 lineage source. The four-edge/full-H3 seal was separately fully
reviewed at its unchanged pin d6b148954b09e259d0e0598cc4f754904ba1977b373a3f491196ef87f93f9501.
No correction is requested.

## Mathematical completeness of the finite interval

The unimodular M6 lattice has normalized primitive pool 1,2. At the
retained whole A5, a finite new simple root has at most one contact,
since two contacts make a cycle. The finite rank-six extensions are
exactly the two isolated rows, five A6/D6/E6 simple-three rows and two
B6 endpoint-four rows. The norm assignments and golden native products
give exactly the nine listed literal profiles; no H6 or mixed contact
is being suppressed. The ambient H3 seal removes transverse contacts.

For contact k, summing the six coordinates and computing their squares
gives precisely

`Z^2-6phi a^2=6q-k(6-k)`.

The full lattice glue is the integrality of the six reconstructed
coordinates `(y+1)^k,y^(6-k)`, with `y=(Z-k)/6` (also correct at k0).
Comparison with the explicit upper gives `a^2<=aH^2 R`; the conjugate
equation gives `a'^2<=phi R/6`. Positive Z and these bounds imply the
required collision comparison, explaining the checker's later assertion.

The rectangular coefficient box is complete. First, the maximum strip
has `a<59/2` and `|a'|<9/5`. Negative B is impossible. B at least 14
would give `a>14sqrt(5)-9/5>59/2`, using the displayed exact comparison
`5*140^2-313^2=31`. Then `A=a'+B/phi` gives `-1<=A<=9`.
The six R-values each test all 154 pairs, exactly 924 tests. These are
necessary interval bounds, not a freely chosen height cutoff.

The helper represents `a+b sqrt(5)` with Fractions. Its sign and
maximal-order integrality tests are exact. The square-root routine is
complete: rational inputs allow either rational or pure-sqrt(5) roots;
nonrational inputs use both rational norm signs and
`u^2=(x.a+/-sqrt(Nx))/2`. A root with nonzero irrational coefficient
has nonzero rational u, so the division cannot omit a solution. Both
signs of each resulting field root are returned; the later positive-Z
filter is the mathematically required orientation.

The code reconstructs every integral candidate and filters only by
already actual beta/u wall acuteness. It retains both expected vectors
h and `g=2h+beta1`; an unrecognized vector or an extra accepted profile
aborts. Their norm/profile identities and equal collision parameter are
checked explicitly, with no narrowing of the original domain or filters.
It does not assume the future wall w is actual to establish h. The
intended h has a unit spatial
coordinate, and its B6 corner has positive leading minors and determinant
one. The endpoint norm and all six incident products are checked exactly.

## Actual facets, reversal and the full M6 conclusion

The actual provenance of u is sound: in the intrinsic basis with
`e1=-c`, the old ambient a1 has coordinates `(phi-1;phi,0^5)`.
It is orthogonal to the whole H3. Its conjunction with those already
actual H3 walls is positive definite, so positive acute projection from
a strict point supplies its actual facet on the H3 face. Mere integral
root existence is not being substituted for this step.

An actual first wall exists before the finite integral upper h, even
though an upper mirror need not itself be simple. The two expected
accepted vectors have the same collision parameter. The extra g has
norm two and only A5 contact -1 at beta2, explaining its D6 profile.
Before h is known actual, its oriented integral arrangement halfspace
still contains the full chamber. Since beta1 is also a nonpositive
actual halfspace, equality in `g=2h+beta1` forces both independent
forms to vanish. This is codimension two, so g cannot be a facet.
Thus a successful corrected census leaves only h among possible first
simple walls. Once that result passes, h is actual. The determinant-one
endpoint and the full seal
give the full positive B6+H3 germ, not a subset of a hidden larger germ.

For the reversal basis `(u,beta1,...,beta6)`, elimination pivots are
`phi-1,-1,...,-1`; the imported nonpivot inverse is safe here. The
checker verifies the entire Gram identity, integral matrix, inverse via
`A^2=I`, future orientation, all prescribed images and trace one.
The orthogonal integral M6 splitting extends A by identity on the
positive H3 summand. It maps the actual positive vertex germ to the
actual germ in the same chamber and consequently preserves that full
chamber. Only then is `w=A(u)` declared an actual facet.

I checked the ten doubled-root matrix against the written P6_10 datum.
The ordinary labels follow from the checked native nine-root contacts.
The three reflected-pair identities give the two weights
`sqrt(2) phi^3` and `2phi^6-1=17+8sqrt(5)` and the remaining right pair,
so all dotted incidences also match. On the h-negative side the reflected
inequalities follow from the original ones by adding positive multiples
of h; hence the doubled region is the union of the two original halves.
The known compact double therefore proves compactness of this nine-wall
region without a classification-uniqueness input.

Once all nine roots are actual, every vertex of that compact region is
feasible in the full chamber by full-rank actual-wall positive projection.
Its vertex convex hull proves equality with the full M6 chamber. This
step is not a claim that any finite acute Gram prefix is automatically
complete. The inherited compact lineage source is
`p6_10_bugaenko_double_lineage_and_classification_audit.md`, SHA256
`651176b56e7d0d0a5e7596b1c889e187f0f5a9703b8336c456aa2e0b445d56c3`.
Its known-example compactness remains a written input; its unrelated
uniqueness claim is explicitly excluded.

## Control and output audit

The helper executes definitions only on import; its own main is guarded.
The new worker has no subprocess, optimizer, internal lock or unbounded
queue. It uses a 45-second child alarm under the existing parent-lock
59-second/1-GiB runner. Both write and verify recompute the reduction.
Generation exclusively creates the fixed certificate; verification bounds
its size and compares the entire reconstructed data. Output is capped at
128 KiB, below the runner's bound. The four direct proof/helper pins,
worker and runner are streamed-hashed before and after the scientific
body; unknown candidates and changed inputs fail closed.

The worker checks fixed arithmetic and matrices, not the published
compact-example theorem, finite-class classification, or the geometric
full-germ argument by executable means. Its certificate explicitly lists
those written inputs. No compact H9 example, global infinite period or
uniform dimension improvement follows. This review leaves all sources
unchanged and does not claim successful corrected generation or replay.

The initial root-owned bounded run failed closed at an overly strict
one-vector expected-count assertion; a bounded diagnostic identified h
and the simultaneous nonsimple g. No certificate was written. I then
fully reread the corrected sources at the pins above. Their mathematical
domain and candidate filters are unchanged, and the correction proves
redundancy instead of deleting a candidate. Successful corrected
generation and independent-process replay remain execution obligations
at this review's handoff.
