# Hostile audit: Raw Stage 0 to PB-Atlas Stage A

**Verdict (2026-08-10): PASS after one provenance/scope wording repair.**
For every accepted index-two Stage-0 candidate, the PB-Atlas boundary weight
is exactly the physical scale recovered by adjunction.  The canonical class,
boundary indices, and weights in all four calibration payloads are derived
without reading their frozen PB counterparts.  The Gram and wall arrays have
shared fixture provenance and were therefore removed from the claim of
independent inverse evidence; their equality remains a useful round-trip
regression check.

Audited theorem:
[`raw_to_pb_atlas_stage0_a_bridge.md`](raw_to_pb_atlas_stage0_a_bridge.md).

## Pinned audited artifacts

| artifact | SHA-256 |
|---|---|
| `results/raw_to_pb_atlas_stage0_a_bridge.md` | `9aa8dca6d1272e171b77b207bbf6c0330a4e538e068ae131b02ea4efdcf781b2` |
| `scripts/certify_raw_to_pb_atlas_stage0_a_bridge.py` | `c5be3a7bc45d686b05e571335929e224b4406bf48fe469a5d04b64d6a02777e1` |
| `data/enumerations/raw_to_pb_atlas_stage0_a_bridge.json` | `ec75b3d30389644fb811e9be2905858b1746566e104382594babf23d57e5448a` |
| `tests/test_raw_to_pb_atlas_stage0_a_bridge.py` | `d1e8d6ccbbb30c743dc9d43706c3668bd521ef4c253b7787e9dc7c9bf39dcc10` |
| `tests/test_raw_to_pb_atlas_stage0_a_bridge_hostile_audit.py` | `b723e0bfc3c9e8475404a57875b48e7f15b9287691eb007250e4a582fde2af7c` |
| `data/fixtures/raw_adjunction_log_stage0_inputs.json` | `c197624eb2986e44c2a991f3433cfe81d5b16ca00593282b481844721dcc4113` |
| `data/fixtures/pb_atlas_numerical_stage_a_inputs.json` | `5b31cb4bf3ef53032f1ce1ff837c2223df40ca6d1cac01a82e42d6ce6aaebd75` |
| `scripts/pb_atlas_contraction_core.py` | `3787c63a1beabd0b8737da1c66293280b6848fc2761af6abd8ce1a4e57d3107d` |

## 1. Why the boundary weight is the physical scale

Let (r_i) be a primitive supplied wall and let Raw Stage 0 recover the
physical smooth-rational class (B_i=d_i r_i).  On an accepted index-two
branch, every inferred ramification index is two and the already-checked log
equation is

\[
       K+\frac12\sum_{i\in I}B_i=0.
\]

Substitution gives

\[
       \sum_{i\in I}d_i r_i=-2K.                 \tag{1}
\]

The PB-Atlas input does not store (B_i) as a separate vector.  It stores
the primitive wall (r_i) and a positive `boundary_weight` (q_i), with the
semantics that (q_i r_i) is the physical anti-bicanonical boundary class.
Consequently (q_i=d_i), not merely because those numbers happen to solve
the aggregate equation, but because they are two encodings of the same
physical class.  Equation (1) is the independent parser check on that
fieldwise identification.

The hostile suite recomputes every (d_i r_i), checks that its square is
(-4), sums the physical boundary directly, and obtains (-2K) in all four
rows.  Increasing any single derived weight by one makes the independent
PB-Atlas parser reject the packet at its weighted-boundary equation.  The
exact results are

| packet | inferred physical boundary | PB weight |
|---|---|---|
| `H10/30` | (B=2r), (r^2=-1) | `[2]` |
| `H11/26` | two primitive (B_i^2=-4) | `[1,1]` |
| `H14/60` | five primitive (B_i^2=-4) | five `1`s |
| `H17/48` | eight primitive (B_i^2=-4) | eight `1`s |

This is exactly the primitive/twice-primitive distinction required by the
half-lattice theorem.

## 2. What is and is not circular in the four equalities

The frozen raw fixtures record, in their provenance map, that their Gram and
wall arrays were projected from
`data/fixtures/pb_atlas_numerical_stage_a_inputs.json`.  Therefore equality
of these two pass-through fields after conversion cannot be independent
evidence that Raw Stage 0 reconstructed anything.  The original report's
unqualified `byte-for-structure` phrase was repaired to state this
explicitly.

The remaining three fields are noncircular:

1. `solve_raw_input` receives only the exact raw object with keys `schema`,
   `picard_gram`, and `wall_vectors`.
2. It derives (K), every scale (d_i), and the nonweak boundary index set
   before the bridge is called.
3. `candidate_to_pb_input` receives only that raw object and the accepted
   candidate.  Its function body has no PB fixture or calibration alias
   input.
4. Only after conversion does `build_certificate` use `CASE_TO_ALIAS` to
   select a frozen native payload for comparison.

The raw fixture's separate `expected.canonical_class` is likewise never
passed to the solver.  The hostile test replaces both modules' PB-fixture
globals by forbidden objects and replaces the bridge loader by a function
that raises on access; Raw Stage 0 and all four conversions still return the
same derived payloads.  In the opposite direction, it corrupts the native
canonical class, boundary indices, and boundary weights.  Construction is
unchanged and the final structural equality becomes false.  Thus the
expected object is an oracle used only after construction, not an ingredient
of construction.

A second hostile path bypasses `candidate_to_pb_input` entirely.  It
assembles each six-field payload independently from the raw arrays and the
Stage-0 certificate, passes it through `parse_physical_branch`, and obtains
the same native structure.  The four equalities therefore decompose as

```text
Gram, walls                copied; shared-provenance regression
K, boundary indices, q_i   reconstructed; noncircular inverse evidence
schema                     fixed interface constant
```

That is the exact evidentiary strength—neither less nor more.

## 3. Routing and generality

The converter itself is row-agnostic and works candidate by candidate.
`CASE_TO_ALIAS` occurs only in the calibration harness.  The four Coble
fixtures happen to have one accepted index-two candidate each, so each emits
one PB-Atlas packet.  A general raw input with several accepted candidates
would emit a finite packet list; the repaired theorem no longer implies an
unstated uniqueness or selection rule.

The toric `H6/9` accepted candidate has the weak `(-1)/(-2)` alphabet and an
empty inferred boundary, so it is routed to the direct weak theorem.  The
maximal (\operatorname{Bl}_p(\mathbb F_4)) triangle has no accepted Stage-0
candidate because (K+B/2\ne0), so it never reaches PB-Atlas.  Neither
control is converted by relabelling or by a maximality assumption.

## 4. Exact boundary of the PASS

This audit proves a closed **numerical data interface**.  It does not prove
that any wall is an effective irreducible curve, that a rational surface or
incidence stratum exists, that the supplied walls exhaust the Mori cone, or
that the eventual construction returns the same polytope.  In particular,
the four exact calibration equalities are not a dimension-bound existence
theorem.  Those obligations remain in PB-Atlas Stage B and the geometric
return theorem.

## Replay

```bash
/Users/sasha/.venvs/global/bin/python \
  scripts/certify_raw_to_pb_atlas_stage0_a_bridge.py --check

/Users/sasha/.venvs/global/bin/python -m pytest -q \
  tests/test_raw_to_pb_atlas_stage0_a_bridge.py \
  tests/test_raw_to_pb_atlas_stage0_a_bridge_hostile_audit.py
```

The bridge archive rebuilds exactly.  The combined bridge suite has 11
passing tests: five owner tests and six hostile provenance/algebra tests.
