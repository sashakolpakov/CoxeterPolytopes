# Compact arithmetic rank21: exact charge upper40

The new same-hub angle pools admit an exact rational weight vector whose
maximum charge over every allowed finite rank21 forest is **40**. Every
enlarged native quadrilateral inequality is satisfied. This is a checked
full-model statement; promotion to a geometric dimension theorem additionally
uses the inherited complete contact domains and face-payment composition.
Independent replay and the final composition review are recorded separately.

## The mathematical change

For each of `A1`, `I2(5)` and `H3`, let `k` be its actual number of whole
components in a vertex forest and `r` its rank. If `k>=2`, give every angle
whose two marked roots lie in distinct copies of that family the additional
weight

\[
  \frac{u_H(a,b,h,d)}{r^2\binom{k}{2}}.
\]

Here `(a,b,h,d)` is the full tuple of numbers of `A1`, `I2(5)`, `H3` and
`A4` components. The total new charge at the vertex is exactly one copy of
`u_H(a,b,h,d)` for each eligible family. All these weights are nonnegative.
The earlier D5 and H4 face-dependent transfers remain in force. The full
formulation and its geometric coverage argument are in
`compact21_self_hub_full_model_formulation_20261018.md`.

The implementation retains every literal native face before deduplicating
identical enlarged inequalities. It has1979 native faces,102917 spectator
cases and53495 distinct enlarged face rows. The finite-forest recurrence
uses all68 allowed irreducible types,729 exact hub tuples and243 capped
groups. New pool costs enter each complete forest budget once, and never
the residual-component Bellman items.

## Exact reconstruction

The single guarded IPM solve is retained, including its raw hexadecimal
primal and numerical dual, in
`fresh_compact21_self_hub_epigraph_minimum_20261018.json`.
Its numerical objective was40.00000000000242. Upward rounding of every
positive value to the1e-10 grid produced a genuinely checked full-model
upper40+3e-10, attained by `D21`; that rounded vector does not prove40.

The solver-free reconstruction changes exactly five entries:

| Global column | Rounded proposal | Reconstructed value |
|---|---:|---:|
|19152|1/10000000000|0|
|19170|1/10000000000|0|
|19188|1/10000000000|0|
|95219|360000000001/10000000000|36|
|95228|390000000001/10000000000|39|

The three raw zero residues are between5.9e-13 and7.6e-13. The two pool
values differ numerically from36 and39 by about2.42e-12. Every other
positive raw entry is greater than1e-4. All other rounded bonuses and the
transfer `14030552/21687803` are unchanged. There are30596 positive old
bonuses and591 positive self-hub pools after reconstruction.

These numerical observations only select a rational candidate. Acceptance
uses exact arithmetic throughout: each enlarged face row is evaluated over
the rationals, and the complete finite-forest maximum is recomputed from
the rational component charges. Numerical auxiliaries and numerical dual
feasibility are not used.

An earlier diagnostic removing only the three zero residues failed: all
faces passed, but `16 A1 + H3 + I2(6)` had charge40+1e-10 because the pool
value36 was still rounded upward. Its compact failure ledger is retained.
The successful reconstruction includes both nearly integral pool values.

## Checked outcome and scope

The successful generation was session29967, exit0, terminal chunk565ddb.
It took3.776s and at most189382656bytes of process-tree RSS, with one
thread, a59s wall cap and1GiB RSS cap. Its start resource gate measured
195.71GiB free disk and6.32GiB available RAM. Temporary candidate files
were removed after atomic publication.

The resulting certificate is
`fresh_compact21_self_hub_exact40_reconstruction_20261018.json`.
It verifies every enlarged face row and gives exact maximum **40**, attained
by `D21`, with zero self-hub addition. No second optimizer was run.

The independent checker subsequently passed in root session59683,
chunk7766fb. It imports none of the new model or producer helpers and
reconstructs the native angle and forest calculations separately. It
checked1979 faces,102917 spectator cases,729 exact hub tuples and850
forest offers, with minimum payment `60007431727/60000000000>1` and
exact maximum40 at D21. Its ledger is
`fresh_compact21_self_hub_exact40_independent_check_20261018.json`.

This value is also the exact optimum of this weight model: the `D21`
ordinary forest has baseline40, receives no exceptional transfer, and all
additional bonuses are nonnegative. In its tree there are20 adjacent
unordered root pairs and20 pairs at graph distance two.

For geometric nonexistence, the inherited strict rank21 h-vector inequality
requires average vertex charge greater than40 once the local face-payment
conditions hold. Therefore an upper bound of40 suffices; a strict upper
bound below40 is unnecessary. The field reduction to `Q(sqrt(5))`, complete
finite/dotted contact domains, the D5/H4 payment lemmas and this global
composition remain explicit dependencies. This file does not replace their
proofs or claim Lean verification. The final geometric composition and
fresh replay are separate records. The producer's fresh `--verify` replay
passed in session9122, chunk5e5c28. Root's streamed final composition
passed in chunk64777c, checking all13 dimension21 field branches and24
higher-dimensional branches,107 recursive manifests and228 distinct files.
Its certificate is `compact_arithmetic20_composition_20261018.json`; the
full root review is `review_compact_arithmetic20_composition_root_20261018.md`.
Together these discharge the internal compact arithmetic ceiling20
composition. The noncompact and combined arithmetic ceilings remain21.
The JSON's null dimension field records its original generation scope;
the separately bound composition certificate supplies the dimension result.

Reproduce the exact check without invoking an optimizer:

```
python3 -B scripts/fresh_compact21_self_hub_exact40_reconstruction_20261018.py --verify
```

The replay checks source/data SHA256 pins, reconstructs the vector from the
retained raw proposal, evaluates all face rows and recomputes the complete
forest maximum. The independent checker uses the same certificate via its
`--candidate` argument and reconstructs the pool and forest calculations
without importing the new model engine.
