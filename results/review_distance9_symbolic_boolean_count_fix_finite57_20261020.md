# Distance-nine replay: exact Boolean-count correction

Source-only review, 20 October 2026. I authored the batched replay and
independently compared this correction with the immutable producer; this is
not an outside-author review of that worker's original design.

Reviewed current worker:
`scripts/check_selective_distance9_batches_20261020.py`,
SHA256 `126d5fe3588cfbb56009bff1f748dbab2154582e3df814feff50b32003cf023b`.
Its byte-preserved predecessor decodes to SHA256
`d9ab46e4404ed8dc02b0166f220cb0e63bc62eb42beaf68a5a7823fa65f84745`
from
`manuscripts/source_archives/check_selective_distance9_batches_20261020_before_symbolic_boolean_count_20261020.py.gz`.
I checked the exact diff: its sole change is line227,

```python
sum(x*y >= 1 for _, _, x, y in rows)
# becomes
sum(1 for _, _, x, y in rows if x*y >= 1)
```

**Equivalence conclusion:** the correction preserves the complete candidate
domain and every classification predicate. It replaces invalid arithmetic
on SymPy BooleanAtoms with addition of Python integer ones.

## Why this is the producer's same count

The immutable producer
`scripts/fresh_selective_distance9_dimension108_20261017.py`
has unchanged SHA256
`f5cd7d0181cf9d5f96b10b6d4d6e7654caa3369b18ffa7473a16033f3124c81a`.
Its `disjoint` function computes the same white inverse
diagonals by the same positive diagonal congruence and exact rational
inverse. For each ordered eligible port pair it branches on `x*y >= 1`
and creates `white_rejection` exactly in that branch.
Its line196 then counts the Python membership predicates
`sum('white_rejection' in r for r in rows)`.

The batched worker's `join_domain` computes those identical
\(x,y\) but needs only the reconstructed predicate to count rejected joins.
The new filter uses that same truth evaluation and adds a Python `1`
for each true predicate. Thus it is precisely the producer's rejection
count, including a boundary product \(xy=1\). The checks of 4442 candidates,
1103 rejections, the remaining grouped joins and 85 inertia types are
unchanged. No candidate, ordering, matching size, port or context is removed.
An undecidable symbolic predicate would still fail closed in the filter,
as it would in the producer's `if`; no approximate Boolean
conversion or tolerance is introduced.

The recorded domains failure is a runtime type error in counting, not a
failed mathematical comparison. The runtime ledger reports the predecessor's
metadata pilot succeeded and the domains pilot stopped at this sum.
Neither that metadata pilot nor this source review establishes complete
replay coverage.

## Exhaustive reduction/type check of the mutable worker

I read every `sum` and scalar reduction in the predecessor,
then confirmed no other line changed:

- Lines240,248--249,270--272 count distance predicates. These distances are
  Python integers: `old.P.distances` resolves to the BFS in
  `fresh_selective_distance13_dimension146_20261017.py`,
  initialized with Python zero and incremented by Python one. Its predicates
  are built-in booleans, so these existing sums are valid.
- The combinations total at273 sums Python `math.comb`
  integers. Counts at348 are Python Counter values. The metadata increment
  at323 compares a JSON string with a string and adds a built-in Boolean;
  its other counts and maxima are JSON integers.
- Other SymPy comparisons occur in `if`, `assert`,
  `all`, `any`, or numerical `max/min`
  contexts, where decided exact predicates are truth-tested rather than
  added. Rationalization asserts every matrix entry is rational, and the
  inherited inertia and corner helpers operate on those exact entries.
  I found no second arithmetic accumulation of a SymPy BooleanAtom.
- The facet polynomial explicitly converts its exact integer coefficients
  using `int`, matching the producer. Profile arithmetic uses
  Python integers and `fractions.Fraction`. Inertia signs are
  converted by the producer's explicit `int(bool(...))`.
  JSON-facing matrices and inverse entries remain strings; no new SymPy
  scalar is sent to the row equality serializer.

I found no additional concrete scalar-type blocker requiring a change.
There is therefore no reason to broaden the minimal patch to the unchanged
enumeration, exact comparisons or frozen helpers.

## Scope and companion bindings

Current fixture pin read:
`scripts/check_distance9_array_scanner_fixtures_20261020.py`,
SHA256 `1f0bf6d952d536a0eb7316c088d8dceb23c52687c8e7fdc2536547b5371cf09b`.
Current coverage-checker pin read:
`scripts/check_partitioned_replay_coverage_20261020.py`,
SHA256 `41f3c82b86ad05aa0877a34cf4718f92d5253e72fe669c083dd161a22c4342c6`.
Companion fixture and coverage reruns belong to root's bounded serial
execution. This review ran no mathematical worker, import, fixture,
compilation or solver. It does not certify the remaining phase indices,
their resource performance, or the geometric theorem behind the frozen
packet. All eight complete phase domains and their existing freshness,
input-pin and output-equality obligations remain in force.
