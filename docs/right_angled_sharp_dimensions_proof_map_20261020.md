# Sharp right-angled dimension bounds: proof and reproduction

For current reading, see the [current sharp right-angled manuscript](../arxiv/sharp_dimension_bounds_8_and_4/main.tex) and the
[submission guide](arxiv_submission_editions_20261028.md). The supporting
arguments and recorded verification scopes below remain distinct from those
current submission builds. Standalone draft typesetting routes have been retired.

20 October2026. One focused manuscript proves both bounds:

|Class|Published ceiling|New internally reviewed ceiling|Largest published example|
|---|---:|---:|---:|
|Finite-volume right-angled, all vertex types|12|8|8|
|Entirely ideal finite-volume right-angled|6|4|4|
|Compact right-angled|4|4, unchanged published result|4|

The new bounds are sharp. They do not assume arithmeticity, a restricted
facet spectrum, a combinatorial enumeration or an optimizer result.
Written mathematical review is distinct from external refereeing and
from a full formalization of hyperbolic geometry.

## One primary paper

The supporting argument
and PDF
share the finite-point projection lemma. Its complete root and two
independent manuscript reviews are recorded separately. The general
proof has a complete [authoritative note](../results/right_angled_dimension8_descending_height_proof_20261020.md)
with full [root](../results/review_right_angled_dimension8_descent_root_20261020.md)
and two independent reviews
([first](../results/review_right_angled_dimension8_descent_finite57_20261020.md),
[second](../results/review_right_angled_dimension8_descending_height_batch_guard_20261020.md)).

The earlier four-page ideal cusp-height manuscript is supporting
history. The completed12-page mixed10 manuscript and older ideal5/6
necessary-condition papers are preserved predecessors. They are not
additional current proposed papers, and none is needed by the new
dimension-eight proof. The separate ideal4 volume/structure manuscript
keeps its existing scope.

## External inputs

|Input|Where used|Exact source and scope|
|---|---|---|
|At least27 facets for every finite-volume right-angled six-polytope|General8|Dufour2010 Lemma1 and Proposition4. The exact universal statement is explicitly given in [Alexandrov2023 Section4](https://arxiv.org/html/2303.09533v1#S4); Section6 derives27 from nu5>=26 and F6>=1+nu5. The original Dufour full text was not retrieved in this pass. This is a facet minimum, not merely a27-facet construction.|
|At most24 vectors in a60-degree spherical code in R4|Ideal4 only|[Musin2008 Main theorem, printed p4](https://annals.math.princeton.edu/wp-content/uploads/annals-v168-n1-p01.pdf). No kissing-code uniqueness is assumed.|
|At least24 facets for ideal right-angled4, with metric24-cell equality|Ideal4 only|[Kolpakov2012 Theorem4, printed p1945](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf).|
|Actual finite-volume right-angled examples through dimension8|Sharpness only|[Potyagailo--Vinberg2005 Section3 and Proposition3.1, pp66–67](https://ems.press/content/serial-article-files/43032). Its P8 has240 facets and2160 cusps. This paper is not being cited for the later universal27-facet minimum.|

Standard foundations are all-wall acuteness, finite-volume/right-angle
inheritance by actual faces of dimension at least three, and inheritance
of entire ideality. The all-wall sign is explicitly recalled in
Vinberg1985, printed p42, Section2.1. No descent reaches a polygon.

## Internal argument

At a finite relative-interior point q, put t(m)=-q.m for each actual
outward unit wall normal. Choose n maximizing t=t(n) and project to

    q'=(q+t*n)/sqrt(1+t²).

All-wall acuteness proves strict feasibility at every other wall.
Orthogonality makes the neighboring normals already intrinsic unit
normals. Their heights are at most t/sqrt(1+t²). Repeating gives
A_j<=A0/sqrt(1+j A0²)<1/sqrt j.

For a hypothetical nine-polytope, three actual projections give a
six-face with all heights below1/sqrt3. Its radial wall directions
are distinct vectors in S5 with products strictly below1/4. The
paper proves the cubic identity

    F(x)=(x+1/2)²(x-1/4)=x³+(3/4)x²-1/16

as an explicit sum of tensor squares. It yields N<=27; at equality
all off-diagonal products would be-1/2. Four such vectors would
have squared sum-2, so N<=26. This contradicts the published27
facet minimum. Higher dimensions are excluded by nine-face inheritance.

For a hypothetical entirely ideal five-polytope, one projection gives
an ideal four-facet with all heights strictly below1. The radial
code has at most24 points by Musin; Kolpakov gives equality and the
regular metric24-cell. Its opposite normal pairs exclude a point
whose heights are all strictly below1. Five-face inheritance gives
the higher-dimensional tail.

## Exact script and conditional Lean boundaries

The [SOS worker](../scripts/check_right_angled_sharp_dimension_sos_20261020.py)
works in a literal polynomial ring over Fraction in twelve independent
variables. It expands all36 rank-two and216 rank-three tensor entries,
checks their products and explicit unit-sphere ideal remainders, the
scalar factorization, positive coefficients, ratio27 and strict endpoint.
This is universal algebra, not sampling or code enumeration. Its
generation and replay records must be read separately from the proof.

The [24-cell worker](../scripts/check_ideal_right_angled_dimension4_model_20261020.py)
checks the displayed normal model,276 Gram pairs,12 opposite pairs,
24 displayed null vertices/cubical links and actual facet witnesses.
It does not claim complete vertex enumeration or arbitrary-polytope
classification. Neither worker proves real norm positivity, actual
face geometry, Musin's theorem or the geometric facet minima.

The [Std-only Lean module](../formal/CoxeterProofs/RightAngledSharpDimensions.lean)
checks exact cubic coefficient lists, the concrete doubled Gram
quadratic form-4, N²<=27N implying N<=27 and the strict cap26,
composition of three actual projection steps, and both dimension
implications. Only actual facet counts are natural numbers. Real
height bounds are abstract geometric predicates; no real height or
norm is replaced by an integer variable. Published inputs, actual
geometry, the universal real tensor argument and real Gram positivity
remain separately typed explicit premises. Universal8 is independent
of Musin/Kolpakov and the entire-ideal4 theorem.

## Reproduce

From the repository root, the catalogue entry is `right_angled_sharp`:

```sh
python3 -B scripts/verify_manuscript.py right_angled_sharp --mode core --pdf
```

This runs the two exact workers, one capped Lean module with its
12-theorem printed-axiom audit, and a temporary PDF compilation,
sequentially. The normal shared lock, one-thread limits,40GiB disk
and6GiB available-memory start gates remain in force. A resource
deferral means no verification pass is claimed for the deferred job.
Temporary build files and inspected figure previews are removed;
proof sources, compact certificates, final PDF and ledgers are kept.

## Actual final verification

The [complete four-job core](../results/right_angled_sharp_core_verification_20261020.json)
passes both exact scripts, the new conditional Lean module and a fresh
temporary PDF build. The SOS certificate was generated in a separate
[bounded batch](../results/right_angled_sharp_sos_generation_20261020.json)
and independently recomputed by the core replay. Its generator used
under56MB process-tree RSS. Both scripts use exact arithmetic.

The five-page final PDF
has zero warnings, with source-bound [build record](../results/right_angled_sharp_dimension_bounds_build_20261018.json).
Its sole figure was directly viewed and its disposable preview removed;
the [visual record](../results/right_angled_sharp_dimensions_figure_visual_review_20261020.json)
retains its hash and scope. The core's temporary recompilation preserves
that published PDF. PDF timestamps are not a byte-for-byte reproducibility
claim; the exact certificate replay and source hashes are deterministic.

The [conditional module](../formal/CoxeterProofs/RightAngledSharpDimensions.lean)
compiled with one thread and a2048MiB cap. Its separate
[12-theorem printed-axiom audit](../results/right_angled_sharp_lean_axiom_audit_20261020.json)
passes. Coefficient/Gram value theorems and the projection-chain
composition are axiom-free; the other printed dependencies contain
only propext and Quot.sound, with no sorryAx. Explicit geometric and
published hypotheses are still hypotheses. No complete geometric
theorem is labelled end-to-end Lean-verified.

The catalogue and assumption register as used by the core are retained
as byte-verified compressed snapshots. Subsequent catalogue status
updates do not silently rebind those historical hashes. The earlier
ideal-cusp module's resource-deferred core record remains preserved;
it is not used as verification of the new module or result.

The [immutable package record](../results/right_angled_sharp_dimensions_20261020_freeze.json)
binds the final sources, PDF, exact certificate, actual ledgers, reviews,
proof map and core-time input snapshots. Copies of the Lean compile and
audit records live under `results/`; future ordinary replays may refresh
the runner's mutable `.proof-batches` records without changing that evidence.
