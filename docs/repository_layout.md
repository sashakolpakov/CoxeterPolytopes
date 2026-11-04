# Repository layout

The publication layout follows the separation used in `PseudoArithmetic`:
manuscripts, local upload packages, current release metadata, audits and
research history have distinct homes.

| Location | Contents |
|---|---|
| `arxiv/` | The three paper families and their [submission selection](../arxiv/SUBMISSION.md) |
| `arxiv/hyperbolic_coxeter_dimension_bounds/` | Concise 57/17/15 manuscript, matching PDF, bibliography and exact proof-addendum map |
| `manuscripts/addenda/dimension_bounds_20261030/` | Byte-preserved full dimension proof: source, bibliography and 551-page PDF |
| `dist/` | Latest local upload ZIP and checksum; generated and ignored by Git |
| `release/manifest.json` | Current manuscript and package index |
| `release/dimension_bounds/` | Concise-edition preparation and typesetting receipts |
| `repro/` | Reproduction instructions, pinned requirements, input and archive manifests |
| `formal/` | Lean sources, project pins, certificates and retained verification inputs |
| `results/`, `data/` | Scientific results, certificates, exact inputs and execution evidence at their recorded paths |
| `audit/` | Repository organization checks, priority evidence and typesetting admission records |
| `notes/` | Guide to research notes; local exploratory scripts under ignored `notes/experiments/` |
| `docs/history/` | Ignored local session history and byte-verified archives of old publication administration |
| `.publication/` | Protected prepared bulk archives, primary reference inputs, current delivery pointers and live publication work |

## Historical records and fixed input paths

The full proof addenda, six accepted TeX editions and exact checker inputs
remain at their established paths. This retains source hashes and the
meaning of the historical execution records. `GOAL.md` and
`RESEARCH_PLAN.md` are ignored local history but are also named by retained
checkpoint scripts; their original paths are preserved.

Old publication staging lists, commit-message drafts, workflow watches and
assembly records have been removed from the live staging area. Their exact
bytes are retained in:

```text
docs/history/publication_admin_20261104.tar.gz
docs/history/publication_admin_20261104_manifest.json
```

The manifest records every original path, length and SHA-256 hash. Every
archived member was read back and compared before its loose copy was
removed. The [cleanup receipt](../audit/repository_cleanup_20261104.json)
records the operation and protected manuscript hashes. Historical paths
inside that archive can be restored when needed; the archive is local
administration and is not part of the public proof package.

The cancelled long expansion and superseded presentation bookkeeping are
also retained in the ignored `docs/history/cancelled_dimension_expansion_20261104.tar.gz`,
with a per-file manifest and verified byte restoration. They are not current
submission sources. The original complete proof edition remains available
in `manuscripts/addenda/dimension_bounds_20261030/`.

## Generated products

Lean objects, TeX intermediates, Python caches, temporary sites and local
upload archives remain untracked. Final PDFs, scientific certificates,
verification stdout, failures and proof sources are retained. The cleanup
found no remaining Lean or TeX compiler products in this checkout.

The separate 936-file/33-archive research payload was not altered or
uploaded. The [public-data guide](../repro/PUBLIC_DATA.md) records its
availability independently of source packaging.
