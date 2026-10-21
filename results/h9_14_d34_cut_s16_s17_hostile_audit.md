# Hostile audit: the variable `S16/S17` D34-cut theorem

Status: **independent raw-denominator and finite-volume audit PASS**
(2026-08-08).  The K3 overclosure, Torelli, and complete-nef results are
hash-pinned but explicitly **transported**, not presented as a second proof.

## Verdict

No discrepancy was found in the new variable two-support theorem.  A
clean-room program importing neither the production denominator, the old D34
corpus/compiler, a T2 loader, nor the `D16+/D17+` theorem synthesizer gives:

| stage | `S16` | `S17` | total |
|---|---:|---:|---:|
| support automorphism order | 48 | 8 | — |
| labelled eight-`m4` assignments | 12,870 | 24,310 | **37,180** |
| support-automorphism orbits | 376 | 3,234 | **3,610** |
| cut-normalizable labelled assignments | 3,270 | 1,519 | **4,789** |
| cut-normalizable diagram orbits | 136 | 238 | **374** |
| scaled-dual oriented Gram lanes | 272 | 476 | **748** |
| Lorentzian corank-four diagram orbits | 1 | 1 | **2** |

The orbit counts agree by two routes: explicit group-action partition and
Burnside averaging.  The graph groups themselves come from a custom exhaustive
adjacency-preserving bijection search; no graph-isomorphism package or archived
automorphism list is used.

## The cut criterion is exhaustive

For every one of the 37,180 raw masks, the audit independently solves the xor
system

```text
m3: b_u xor b_v = 0,
m4: b_u xor b_v = 1.
```

It separately enumerates all vertex cuts with vertex 14 fixed outside.  The two
sets agree exactly.  Thus a primitive `{2,4}` norm assignment exists if and only
if the eight `m4` edges equal `delta(X)`.  Connectedness makes `X` unique up to
complement, and the two complements are exactly the scaled-dual orientations.

This is both a proof from the coboundary equations and a full computational
check on the stated denominator—not a pattern inferred from the survivors.

## All 748 exact inertia decisions

The audit constructs every integral Gram directly from the cut bits and uses
exact rational symmetric congruence, including singular two-by-two pivots.
There are no floating eigenvalue decisions.

| support | exact oriented-lane histogram |
|---|---|
| `S16` | `216 rank14/(2,12,0)`, `46 rank13/(2,11,1)`, `8 rank12/(2,10,2)`, **`2 rank10/(1,9,4)`** |
| `S17` | `326 rank14/(2,12,0)`, `102 rank14/(3,11,0)`, `36 rank13/(2,11,1)`, `4 rank13/(3,10,1)`, `6 rank12/(2,10,2)`, **`2 rank10/(1,9,4)`** |

The full lane records are committed through deterministic inventory hashes in
the JSON certificate.  Each pair of Lorentzian lanes is one scaled-dual diagram
orbit, and its mask is fixed by the entire support group.

The oriented selector has a useful sharper denominator:

| selector stage | `S16` | `S17` |
|---|---:|---:|
| `short4/long10` lanes before inertia | 25 | 40 |
| `short4/long10` Lorentzian lanes | **1** | **1** |

Thus `D16+` and `D17+` each select exactly one marked Lorentzian Gram lane; the
opposite scaled-dual lane has `short10/long4`.

## Independent cone and Vinberg replay

For each selected survivor, a lexicographic nonsingular principal root-span
basis is found from the generated Gram.  Every extreme ray is then regenerated
from all rank-nine active-wall subsets.  The program verifies a strict positive
interior, no hyperideal rays, and a positive-square relative-interior witness
for every one of the fourteen facets.

Those active-set endpoints agree exactly with the independently invoked generic
Vinberg endpoint/full-neighbourhood engine:

| survivor | exact rays | finite | ideal | elliptic rank-8 edges | endpoint types `(0,2)/(1,1)/(2,0)` | bad edges | bad full cusp neighbourhoods |
|---|---:|---:|---:|---:|---|---:|---:|
| `S16` | 110 | 90 | 20 | 683 | `97/362/224` | 0 | 0 |
| `S17` | 127 | 112 | 15 | 707 | `53/300/354` | 0 | 0 |

Both are therefore finite-volume, noncompact Coxeter polytopes in `H9`.  This
does not confuse the existence of finitely many rays with finite volume: every
elliptic edge has exactly two finite/ideal endpoints and every connected affine
component has a valid full zero-neighbourhood.

## Public-matrix lineage and synthesis checks

Only after completing the raw census, the audit reads separately hash-pinned
public matrix arrays.

- The generated oriented `S16` Gram equals the public T2 Gram entry-for-entry.
- The generated `S17` Coxeter relation matrix equals the public T1 relation
  matrix label-for-label.  The public factory Gram is a different root scaling,
  so entrywise equality is neither expected nor claimed.
- The T2 fixed Gram, rank-18 Gram, involution, and ordered 24-root hashes are
  recomputed directly from the public arrays and equal the synthesis manifest.
- The frozen T2 bridge matrix has determinant `+1`, integral inverse, and its
  published string-matrix digest replays exactly.
- The two T1 rank-18 hash manifests in the obstruction archive agree with one
  another and with the synthesis constants.

## Exact independence boundary

Independently proved here:

- both raw support graphs and automorphism groups;
- all 37,180 masks and 3,610 raw orbits;
- the cut-if-and-only-if normalization criterion and all 374 compatible orbits;
- all 748 exact Grams/inertias and the two oriented selectors;
- both survivor facet Grams, irredundant cones, exact endpoint sets, and full
  Vinberg finite-volume decisions.

Independently matrix/hash-checked from frozen arrays:

- T2 facet/package hashes and bridge-matrix unimodularity;
- T1 labelled Coxeter identity;
- all separately pinned synthesis-source file hashes.

Transported, not re-proved:

- the `16/5` stable-overclosure inventories and uniqueness of the weak K3
  packages;
- the remaining full T2 marked-bridge equalities;
- period/Torelli existence and the complete T2 invariant-nef equality `C`;
- the rank-18 arithmetic proving that the five T1 hidden normals lift to
  exchanged `(-2)`-root pairs.

The last limitation is unavoidable from the selected public input: the T1
obstruction JSON exposes two agreeing rank-18 hash manifests but not the full
rank-18 matrices.  The audit therefore records the `-6,+6` witness metadata and
pins the source, but does not pretend to have recomputed those hidden-root
pairings.  Consequently this is an independent proof of the variable
surface-free arithmetic and `V`, plus a hostile synthesis cross-check—not a
second self-contained proof of all K3 obligations `U/R/C`.

## Replay

```bash
PYTHONPATH=.:src python3 scripts/audit_h9_14_d34_cut_s16_s17_hostile.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_h9_14_d34_cut_s16_s17_hostile_audit.py
python3 -m py_compile scripts/audit_h9_14_d34_cut_s16_s17_hostile.py
```

The focused suite has **6 passing tests**.  Frozen SHA-256 values:

- script: `a62dc8a100f60fe6f95f5722fb87c2df0df3c79b8438db86ff98c3a9d0f38efe`;
- JSON: `b3c26174158f34c48b176b735896a623db372b455bdf19e5bed36f5a25f8cf3e`;
- tests: `e9aa6a8734e0a77ef78effca0f7a8522759a5942c859b84b8e94c8c71edfe490`.

Artifacts:

- [`audit_h9_14_d34_cut_s16_s17_hostile.py`](../scripts/audit_h9_14_d34_cut_s16_s17_hostile.py)
- [`h9_14_d34_cut_s16_s17_hostile_audit.json`](../data/enumerations/h9_14_d34_cut_s16_s17_hostile_audit.json)
- [`test_h9_14_d34_cut_s16_s17_hostile_audit.py`](../tests/test_h9_14_d34_cut_s16_s17_hostile_audit.py)
