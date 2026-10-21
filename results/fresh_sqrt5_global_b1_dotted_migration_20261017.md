# Exact comparison before removing a redundant affine-subdiagram scan

The stopped dotted attachment search has6,537 processed cores and8,149 marked
states. Every processed core has been freshly evaluated by the successor
search. All14,848 accepted attachments agree, together with the full contact
domains, unchanged arithmetic tests and retained child maps.

The [conjugate radical lemma](fresh_conjugate_radical_lemma_20261017.md)
shows why the former exponential search over affine principal subsets is
redundant: conjugate positivity transports a principal kernel into the full
kernel. The existing exact full-matrix tests already reject any resulting
nonnegative relation, as well as the other incompatible signatures.

The [comparison checker](../scripts/fresh_sqrt5_global_b1_dotted_migration_20261017.py)
separately checks each old affine witness. Every one remains rejected by the
successor, with2,655 changed rejection explanations:

| New rejection reason | Count |
|---|---:|
| Two physical negative directions |1,585|
| Conjugate matrix fails the required positivity |419|
| Cyclic product fails the coefficient-field test |651|

All other individual trial records agree exactly. Only after the full
comparison passed was the stopped database copied to a separate successor
database and its processed records replaced by the newly checked versions.
The predecessor is immutable and will not be resumed. The
[saved comparison](fresh_sqrt5_global_b1_dotted_migration_20261017.json)
binds that stopped database, its output and the exact sources. The successor
database is deliberately outside this frozen comparison because its queue
continues to grow.

This is equivalence of a completed **partial prefix**, not exhaustion of the
attachment search or a new face-payment theorem. The successor retains the
same377 initial parameter seeds, full finite-corner contact domains, marked
isomorphisms and pruning assumptions. Singular stopping configurations are
still recorded and expanded; no face obligation is removed by this change.

Generation and a fresh saved replay both passed. The replay recomputes every
processed predecessor decision again without touching the growing successor:

```sh
python3 scripts/fresh_sqrt5_global_b1_dotted_migration_20261017.py --verify
```

The active successor is
[the separate fast queue](../scripts/fresh_sqrt5_global_b1_dotted_fast_20261017.py).
Its eventual exhaustion and complete saved replay remain required before
calling the attachment census complete. Any weighting proof also needs the
actual face payments and all relevant field branches.
