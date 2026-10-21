# Static hostile audit: `JK-H12` five-vertex fork, `B=3`

Status (2026-08-10): **SUPERSEDED BY TERMINAL PASS**.  This document records
the read-only launch-time audit of the then-active run.  The run subsequently
terminated green, and every item in the terminal checklist was discharged by
the [hostile terminal acceptance certificate](../data/enumerations/jacobian_kummer_h12_five_vertex_fork_b3_terminal_acceptance.json).
The completed theorem is recorded in the
[`B=3` freeze report](jacobian_kummer_h12_five_vertex_fork_b3_complete_census.md).

## Frozen code and inputs

The audited top-level builder is
`scripts/certify_jacobian_kummer_h12_five_vertex_fork_streaming_census.py`
with SHA-256

```text
64dd51e4410a379723669eee4002f16029c4cab296ba8693e4483f1eb7f6627e
```

Its current imported exact builders have SHA-256 values

```text
d7c5c02803911f8d205a63c59027d8d9584cce998de5efa0a7acb134d56d00e4  P5 orbit/boundary builder
bbb99d59726e0b6bcf0d6b3748823ca9bc8f1efe5ecbe5c29a16f4a5db6a53fb  star builder
474c796521edeb0f7460906d3779cb0e540d8cbebef0f97176fbe812ec1ba78c  color-refinement builder
```

The selected frozen archives hash to

```text
2306ebe0048b576a5ba1aa145d5af500b6efdb13f16ffacb8174c7984d1d5d10  all-shape tree denominator
32943e8531bdaa1874af41d45635ea792c9337694adc4f20dfc152aac3d7713a  complete P5 B<=2 census
b535bbe5d1639358e5a336ce606a4d7eae8dedd0b55835f9ee057321ac394299  P5 ABBB summary/rows
04c9ff621b6ff47e690902e872204790982040cc725548ed3abd808670219e22  P5 BABB rows
21fbd0946b910360cdd4096392e4e629907f971be85359d5b87e312a78d4f309  complete star census
```

The independently archived fork pattern breakdown used below hashes to
`99294000d296e82d072a886c298265082acc20d2f912d2291e267c74bf3e7be7`.

The high-`B` archive hashes agree with the constants enforced by the active
builder.  The builder also rechecks the lower-`B` P5 archive, star archive,
tree denominator, stock Gram/action gate, and A/B preservation before doing
the census.  `python3 -m py_compile` passed for this exact top-level snapshot.
These hashes validate the source snapshot; because the emitted format does
not yet record builder hashes, the terminal handoff should retain the
launch-time hash evidence rather than infer process provenance from the later
workspace state.

## Losslessness and exact arithmetic

The `B=3` fork denominator is exactly

```text
BB|AB  36,247
BA|BB  18,055
AB|BB  18,456
total   72,758
```

For each pattern, the program independently requires equality of (i) the raw
proper-colouring count, (ii) the sum of explicit orbit sizes, and (iii) the
Burnside orbit count.  The three `B=3` patterns therefore form a lossless
partition of the selected fork layer.  Every resulting row must have
`752+9*3=779` facets.

The corrected completed same-facet prior iterator has the exact disjoint
source arithmetic

```text
P5 ABBB   36,247
P5 BABB   36,648
star B=3   6,117
total      79,012
```

No other completed row can have the same facet number: all three shapes use
the exact formula `752+9B`.  Thus restricting the cross-shape comparison to
these `79,012` rows is exact, not a heuristic sample.  If the terminal run
passes, the already certified `237,068`-type cross-shape union increases by
exactly `72,758` to `309,826`; that aggregate must not be claimed before the
terminal archive is checked.

## Collision and refinement logic

The nonisomorphism logic is sound.

1. SHA-256 is used only as a negative filter.
2. Every equal current degree digest is replayed to literal degree-profile
   bytes and partitioned by those bytes.
3. Every exact literal degree collision is replayed to the canonical
   color-refinement byte string.  Equal refinement bytes cause an `UNKNOWN`
   abort; different bytes are a valid nonisomorphism certificate even when
   refinement is nondiscrete.
4. Against prior shapes, every equal archived digest triggers literal replay
   on both sides.  Every literal match then triggers exact refinement on both
   sides, again aborting on equality.

The boundary gate is likewise exact: the five transforms
`1,Ra,Rb,Rc,Ra Rd` are checked to have precisely the four marked internal
adjacencies; those primitive hyperplanes are deleted; all surviving walls are
oriented by a strict summed interior.  Each output then passes primitive
integral reflection, timelike-interior, rank-13 Smith-`1` core, complete
ordinary-Coxeter, and signed-int64 bound checks.  The row transcript is
canonical JSONL in deterministic gzip and records independent uncompressed,
compressed, gate-stream, and degree-stream hashes.

## Concrete defects found

### 1. Prior colour/digest association is not self-authenticating

This is a real provenance/losslessness weakness, but it does not invalidate a
run made from the frozen snapshot above.  The prior iterator reads only the
archived degree-digest sequence, regenerates the P5/star colour sequence with
the current imported builders, and attaches the two sequences by positional
`zip`.  The archives themselves already store each row's colours, but those
colours are not compared.  The active output also does not pin the imported
P5 orbit-builder hash.  A later ordering-only builder edit could therefore
pair a valid archived digest with the wrong colour; a wrong noncandidate could
then be skipped before literal replay.

For this active snapshot, the relevant P5 builder is exactly the frozen
`d7c5...00e4` version, so the present ordering is the archived ordering.  The
defect is reproducibility under a future code change, not a detected mismatch
in the active run.

Minimal post-run repair: read `(archived colours, digest)` together from every
prior archive and never reconstruct their association by position.  At a
minimum, pin the imported-builder hashes and assert every regenerated colour
tuple equals the archived tuple before yielding it.

### 2. One summary counter can undercount

`color_refined_discrete_current_count` is incremented during within-fork
collision refinement only.  A current row first refined solely because of a
cross-shape literal match is added to `refined_current_indices`, but its
discrete flag is not added to that scalar.  This cannot create a false
nonisomorphism result; it can only make the diagnostic count too small.

Minimal post-run repair: after both collision loops, recompute the scalar from
the final rows, e.g. count rows with a non-null refinement and true discrete
flag.

### 3. Builder provenance is absent from the emitted archive

`pinned_inputs` records data archives, but not the top-level builder or its
imported exact builders.  Add their hashes after this run; doing so now would
change the active archive target.

Minor nonsemantic hygiene: one image-pair assignment and one theorem mapping
key are duplicated, and `io` is unused.

## Terminal acceptance checklist

Do not promote this audit beyond **PASS-PENDING-TERMINAL** unless the completed
summary and row archive satisfy all of the following:

- `selected_B_range == [3,3]` and the selected patterns are exactly
  `BB|AB`, `BA|BB`, `AB|BB`;
- row/type count `72,758`, orbit-size sum equal to the raw denominator, and
  facet histogram exactly `{779: 72758}`;
- prior same-facet count `79,012`, with an iterator-only audit recovering the
  source split `36,247 + 36,648 + 6,117`;
- every digest candidate is literally replayed, every literal collision is
  refinement-separated, and unresolved cross-shape pair count is zero;
- `all_checks_pass` is true, row count is `72,758`, deterministic transcript
  hashes are present and independently rehashed, and peak RSS remains within
  the stated 1 GiB ceiling.

Recommended lightweight regression test after termination: exercise only the
prior iterator (no boundary rebuild), compare each yielded colour/digest pair
directly with its archived row, and assert the three source counts, total,
`B=3`, source-local sequential indices, and uniqueness of `(source,index)`.
