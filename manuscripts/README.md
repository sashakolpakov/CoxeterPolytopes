# Six papers: reading and verification

These are the six TeX papers prepared for audit. All six are **written,
compiled and internally reviewed**, with matching PDFs at their stated
input boundaries. This is a manuscript status, not a claim of complete
computational replay, external refereeing or full geometric Lean verification.
The arithmetic paper now includes two necessary field proofs previously
cited only from research notes; its completed 21-page PDF replaces the
preserved 17-page predecessor. The other 27 TeX files are preserved predecessors,
supporting notes, input fragments and generated atlases.

| Paper and TeX | Result in this edition | PDF | Proof and reproduction map | Verification ID |
| --- | --- | --- | --- | --- |
| [General finite volume](finite_volume57.tex) | Dimension at most 57 | [48 pages](build/finite_volume57.pdf): 12 main + 36 essential proof appendices | [Corner payments, forests and cusps](../docs/finite_volume57_lemma_map.md) | `finite57` |
| [Compact arithmetic](compact_arithmetic_dimension_bound.tex) | Dimension at most 20; field proofs completed; parked refinements excluded | [21 pages](build/compact_arithmetic_dimension_bound.pdf) | [Field reductions and exact angle rules](../docs/compact_arithmetic20_proof_map.md) | `compact20` |
| [General compact](compact_farthest_vertex_bounds.tex) | Dimension at most 23 | [36 pages](build/compact_farthest_vertex_bounds.pdf) | [Farthest vertex and component energies](../docs/compact_farthest_vertex_bounds_proof_map_20261020.md) | `compact_farthest_bounds` |
| [Sharp right-angled dimensions](right_angled_sharp_dimension_bounds_publication.tex) | At most 8 for finite volume; at most 4 if entirely ideal; both sharp | [5 pages](build/right_angled_sharp_dimension_bounds_publication.pdf) | [Projection, heights and spherical codes](../docs/right_angled_sharp_dimensions_publication_proof_map_20261021.md) | `right_angled_sharp` |
| [Ideal right-angled four-volume and structure](ideal_ra4_volume_and_structure_publication.tex) | Published minimum; second volume and its unique double; restricted 24-cell chains | [17 pages](build/ideal_ra4_volume_and_structure_publication.pdf) | [Volume, equality and decomposition](../docs/ideal_ra4_volume_and_structure_proof_map.md) | `ideal_ra4_volume_and_structure` |
| [Compact right-angled four-volume constraints](compact_ra4_volume_constraints.tex) | Facet capacity, neighboring excess and recognition results | [14 pages](build/compact_ra4_volume_constraints.pdf) | [Capacity and excess](../docs/compact_ra4_volume_constraints_reproduce_20261019.md) | `compact_ra4` |

The priority reading is the 57 paper and the sharp right-angled paper.
Both now have fresh adversarial mathematical reviews: the
[finite57 composition](../results/review_finite57_geometric_composition_new_audit_20261021.md),
[sharpRA geometry](../results/review_sharp_ra_geometric_foundations_adversarial_20261021.md)
and [sharpRA tensor/equality argument](../results/review_sharp_ra_tensor_ideal_equality_adversarial_20261021.md).
No blocking gap was found at their explicit finite-data/published-input
boundaries. These are internal mathematical audits, distinct from full
producer replay, geometric formalization and external refereeing.
The remaining papers also have complete written reviews. The fresh
[ideal RA4 audit](../results/review_ideal_ra4_volume_structure_fresh_adversarial_20261021.md)
checks the universal equality reduction and actual cutting/gluing proof;
the fresh [compact RA4 audit](../results/review_compact_ra4_volume_constraints_fresh_adversarial_20261021.md)
checks its complete capacity, metric and graph cases. The
[completed-audit checkpoint](../results/six_paper_completed_audits_and_restricted_continuations_20261021.json)
binds these reviews and the unchanged six editions. Written audit
completion does not turn partial reproduction or prepared Lean into
executed verification.
The former establishes one corner-payment rule before separating finite
vertex forests from ideal cusp links. Its appendices contain necessary
proofs and algorithm specifications, rather than raw execution logs.
The latter shares one projection lemma between two short branches:
a direct cubic tensor argument for dimension eight, and the
Musin/Kolpakov inputs for entirely ideal dimension four.

The compact arithmetic and general compact PDFs retain the 20 and 23
theorems requested for this finishing phase. Later reviewed written
research now gives 17 and 19, respectively. The
[complete new composition](../results/compact_dimension19_and_arithmetic_dimension17_global_composition_20261021.md)
and its full internal source/composition reviews do not turn those stronger
statements into theorems of the older PDFs.
Any integration belongs in these same two papers, with one coherent
argument rather than chapters for successive increments. The
[paper architecture](../docs/paper_architecture.md) explains the division
by method. Compact right-angled minimum and second-volume questions
remain open; the sixth paper does not settle them.

## What the audit editions assume

The manuscripts distinguish published results, internal geometric
arguments and exact computed propositions. In particular, the arithmetic
20 paper identifies its internal companion dimension/field reductions;
the finite-volume paper identifies its inherited finite tables. Their
proof maps give the precise statements and proof homes.

Every proof-relevant algorithm is to be read as mathematics: its domain,
steps, termination, soundness/completeness and exact output are separate
obligations. The shared
[computational proof contract](../docs/computational_proof_contract.md)
sets out this convention. A hash match or a successful checker run does
not establish geometric coverage by itself.

The parked field-specific arithmetic refinements are excluded from the
current universal paper. Both
[eight-field](appendices/compact_arithmetic_eight_field19.tex) and
[K7](appendices/compact_arithmetic_k7_field19.tex) appendices, the
[combined predecessor](compact_arithmetic20.tex), its
[20-page PDF](build/compact_arithmetic20.pdf) and verification evidence
are preserved. These are not additional papers for this audit.

## Reproduce the declared checks

Run from the repository root. Inspect a plan first; it performs no proof
computation. For the two priorities:

```sh
python3 -B scripts/verify_manuscript.py finite57 --mode full --pdf --plan
python3 -B scripts/verify_manuscript.py right_angled_sharp --mode core --pdf --plan
```

For sharpRA, remove `--plan` to execute its core checks sequentially.
The finite57 full command is an inventory: its two eager loader commands
still need the reviewed bounded replacement fully executed and integrated.
Its
[reproduction guide](../docs/finite_volume57_reproduce.md) gives the
smaller runnable core and the completed bounded parts. The IDs for the
other papers are in the table. Use `--ledger results/NEW_NAME.json` to
save a new receipt; existing receipts are never overwritten. The
[verification guide](../docs/manuscript_verification.md) explains
`scripts`, `logic`, `core` and `full`, including inherited input scopes.
The compact RA4 selection currently checks typesetting only.

Finite57 has completed bounded replacement replays for D13 (19 parts)
and distance nine (2,196 units in 56 batches), plus its other recorded
checks. There is no combined 123-job full-manuscript runtime PASS;
some transitive finite-data inputs remain assumed. The new
[bounded low8 replay](../results/finite57_audit_edition_low8_bounded_replay_20261021.json)
passes all2,552 cusp profiles. Matching bounded records now cover111/113
selected script commands. The two eager commands still await replacement
in the selection. Their reviewed25-unit bounded replacement now has three
actual PASS records (parser and the first two B13 pin groups). The
[successor ledger](../results/finite57_loader_replay_quadratic_b16_20261022.json)
preserves the earlier resource refusal and records the later successful
second pin group; it is not a complete replay. These
execution boundaries do not alter the stated theorem.
The sharp right-angled
package has exact algebra/model checks and selected conditional Lean
verification. Historical core runs retain their original source bindings;
editorial successors have separate source-bound builds and delta reviews.
The proof maps record that distinction.

Selected Lean modules check internal implications with explicit external,
script and geometric premises. They do not formalize every geometric
argument. The [assumption register](../docs/manuscript_assumptions.md)
and per-kernel maps identify those premises; printed axiom lists alone
are not an assumption inventory.

All scientific, Lean and TeX work shares one exclusive lock. Batches
require 40 GiB free disk and 6 GiB available RAM; script/TeX jobs retain
their 59-second/1-GiB caps and Lean its 2,048-MiB cap. Temporary build
products are removed after each batch. Below 20 GiB free disk, cleanup
starts before further jobs. Proof sources, certificates, PDFs and frozen
evidence are retained.

The [audit-edition receipt](../results/six_papers_audit_edition_20261021.json)
binds the preserved six source/PDF pairs before the arithmetic field-proof
completion. The arithmetic predecessor is preserved in its
[translation record](../results/compact_arithmetic20_internal_field_translation_preservation_20261021.json);
the other five pairs are unchanged. The
[current delivery record](../results/six_paper_final_tex_delivery_20261021.json)
binds all six matching pairs and the separate verification scopes.
These are internal audit editions, not
externally refereed or end-to-end Lean-verified publications.
