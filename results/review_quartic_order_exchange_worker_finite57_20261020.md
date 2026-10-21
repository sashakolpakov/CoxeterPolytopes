# Source review of the fixed quartic order-exchange checker

20 October 2026. PASS by full written-source inspection; no import or
execution of the checker/helper, and no independent AST or proof job.

Reviewed `scripts/check_quartic_order_exchange_period_20261020.py`,
SHA256 `a30ea0e6b2b8e5e71f94ab3fd9dc5acab58a320a0a624c42929cceb9b08a9340`.
I read the complete original source (SHA prefix `92435ac0`), then checked the sole final
delta: all collected source paths are rehashed before the result returns.
I also read the complete imported helper
`scripts/fresh_sqrt2sqrt5_lanner_catalogue_20261016_review.py`, SHA256
`122cfeea5190517d9203a61d933fbe964485652be04c5fdca157e492eca2c811`.

**Disclosure.** I authored the first/second-edge and order-exchange paper
sources, with mathematical contributions explicitly communicated by root
and resource_cleanup. Root authored this verifier. This is an independent
implementation review, not another outside-author mathematical review of
my own paper. The separate batch_guard paper review is pinned by the code.

## Arithmetic and import boundary

The helper's K class uses Fraction coefficients in the internal basis
`(1,sqrt5,sqrt2,sqrt10)`. Its opposite quadratic-tower multiplication,
three-conjugate inverse and two-stage exact sign test are correct for the
operations used here. The conversion to the full integral basis is
`(A-B,2B,C-D,2D)`, and both `integral()` and certificate encoding use
precisely that conversion. Thus an integral entry is not confused with
an entry integral only in the smaller square-root basis. The phi
constructor and four embedding signs are correctly ordered.

The helper performs only class/function definitions and a small fixed
constant initialization on import. Its catalogue, support enumeration,
legacy certificate reads and run() are behind uncalled functions and
the main guard. The new worker calls none of them, including root(),
field-degree routines or rank/elimination. It uses no floats or sampling.

## Matrix and output checks

Both five-column frames exactly reproduce the written A/B lists in
the paper bound to SHA prefix `043f1f1f`. The row/column transposes are
consistent. The two
positive Gram blocks and the negative time block are correct. The
displayed block inverse has the correct signs and denominators.
Gram adjoints provide both frame inverses; integral entries and BOTH
inverse identities are asserted before the period is formed. The
period/inverse, Lorentz Gram identity and every frame image are then
checked independently. These checks establish invertibility without
an unsafe pivot assumption on an indefinite matrix.

The four positive simple images include the correct Q expressions
`-F` and `C*E+U*F`. The future-point images and their positive physical
time coordinates are checked; their negative norms follow from the
asserted common Gram matrix and positive physical a. The active trace
is computed from the actual period diagonal, compared with both its
quartic polynomial and its K-cubic expression, and proved greater
than5 by the helper's exact sign test. The extra five fixed coordinates
give the recorded N9 trace. Both forbidden B3 projections are formed
as full vectors and checked nonintegral in O, rather than inspected
through a single guessed coefficient.

Pinned papers, reviews, helper and wrapper are checked initially; the
worker adds its own source digest, then rechecks EVERY collected path
after arithmetic. Fresh writing uses exclusive creation, and replay
recomputes before comparing the entire saved JSON. No mode modifies an
existing certificate. Resource gates and the parent lock remain the
explicit responsibility of the pinned bounded wrapper, not this worker.

## Precise scope

This verifies fixed five-coordinate algebra only. It does not rerun the
first/second-edge finite strips, prove arithmetic class number or norm
pool completeness, classify finite groups, prove either full face
identity, identify actual full vertex germs, or prove the arithmetic
lattice theorem. All those are named written inputs, not encoded as
unexplained computation. The geometry-to-chamber and infinite-index
arguments remain in the paper and its outside-author review. There is
no new global dimension bound, numerical root search or end-to-end
formalization in this checker.
