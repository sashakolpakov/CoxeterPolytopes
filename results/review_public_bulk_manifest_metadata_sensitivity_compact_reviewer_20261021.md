# Bulk publication inventory: metadata and source review

21 October 2026. This review supplies read-only evidence following the automatic rejection of the archive upload. It does not approve an upload or certify that arbitrary payload bytes are non-sensitive. I read the complete small inventory, packaging and inspection sources, and selected provenance sources. I did not read/hash the 17 GB of original payloads, open the bulk archives, execute/import a worker, use the network, or alter the Git index.

## Bound inventory and sources

* `repro/public-data-manifest.json`: SHA-256 `57dfd7fee51da317b7482dda3401b330e84c7b36fbc8bb12d8e0cf0110ff0e6d`, 213,687 bytes.
* `tools/prepare_public_data.py`: `f3d064395928ea2fb2e29883fcd3f5dfb6a161ab0490ec3c3700906a2760b6eb`, 11,692 bytes.
* `tools/audit_public_data_contents.py`: `c5fbdc9a5316b16d6e7aedd98aeb7ced6c0dac6d058a9409267e9da18d2ac4c7`, 6,868 bytes. Its execution/results are outside this review.
* `results/public_data_reference_exclusion_20261021.json`: `860f24cac4ae9ace910c4c6642041fa7fea1854fd158a86d9262b93f996eaa60`, 1,519 bytes.

The current manifest has **936 unique canonical relative paths**, 17,665,895,736 original bytes, and 33 archives totaling 2,177,183,089 compressed bytes. Every listed length is an integer at least 1 MiB and every original digest has the expected SHA-256 syntax. These are inventory checks, not fresh original/archive hash checks.

The exact path partition is:

| Path class | Files | Original bytes | Name/source provenance |
| --- | ---: | ---: | --- |
| `data/candidates/` | 4 | 25,346,330 | Silver K3 Vinberg root-search diagnostics; corresponding diagnostic/certifier sources exist. |
| `data/certificates/` | 1 | 7,483,982 | Dimension58 arithmetic export; its exporter explicitly constructs integer/rational tables and source hashes. |
| `data/enumerations/`, including three subdirectories | 213 | 2,614,475,096 | Mathematical search/census data: Coxeter roots, walls, corners, finite diagrams, surface/lattice chamber computations, SAT instances/proofs, numerical arrays and caches. |
| `data/known_examples/` | 1 | 3,641,778 | General 15-nodal quartic chamber stock. The producing source identifies published Dolgachev–Shimada/Shimada lattice data as its input; this is externally derived mathematical data, not a new privacy or redistribution determination. |
| `results/` directly | 717 | 15,014,948,550 | 716 `fresh_` research outputs plus `fixed_point_symmetry_and_contact_replacements_2026_10_14_audit.json`. The latter's source emits relative source/hash/link provenance for mathematical records, not session transcripts. |

The seven extension classes exhaust those 936 paths: **887 JSON, 14 JSON.gz, 11 JSONL.gz, 13 NPZ, 9 SQLite, one CNF.zst, one DRAT.zst**. The 13 arrays are named H20–H22 geometry banks; the nine databases are named square-root-five root/contact/triangle searches. The SAT pair is the H9 record30 dual exhaustion. Such names identify a mathematical purpose; they do not establish the contents of every cell/string. As an additional filename-level provenance check, 645 of the 717 result paths have an exactly matching Python source basename under `scripts/`; the other 72 include parameterized batch/summary outputs and the separately named audit producer. I do not assert a complete producer execution audit from that check.

## Specific exception and historical material

The original **937-file** manifest `b4373eef8bc086e51fc8d42b98e789bbe9342f4a9fbf85f93288610c4cadcccf` contained the third-party literature PDF `data/references/fresh_all_ideal_seven_simplex_literature_20261017_chen_labbe.pdf` (1,619,220 bytes). It is absent from the current inventory. The current metadata equals the old metadata on every retained original; no original was added, and only archive006's record changed. The old manifest is byte-recovered from `results/source_archive/public_data_before_content_review_20261021/public-data-manifest.json.gz` (compressed SHA `3559d2c43cca166e13d50fee5aab0d7ecf0592f27744b65f1d2e786dcabb3e7a`, 57,165 bytes; decoded 213,980 bytes, exact old digest). I verified this small gzip round trip; the receipt's large archive round-trip assertion is reported evidence, not my execution.

Historical mathematical work remains intentionally included: one `_scratch_fermat...` enumeration, databases with `_running` names, initial/stopped/timeout records, and historical corner/frontier outputs. Those names are not evidence of success or of live jobs. None of the current paths is under `docs/history`, `data/external`, `source_archive(s)`, `.git`, `.aws`, `.ssh`, `.codex`, `.agents`, `.publication`, or `.proof-batches`; none has an AGENTS/session/transcript/credential/private-key-named component. This rules out those *named classes*, not private strings embedded in scientific records.

## Packaging and inspection boundary

The packer's actual discovery is restricted to files at least 1 MiB under `data/` and `results/`; it excludes caches, leaf symlinks, administrative filenames, external downloads and now `data/references/*.pdf`. It does not package the repository, home directories or administrative history wholesale. It constructs tar headers from relative paths with fixed mode/time, rather than copying filesystem owner names and timestamps. One-MiB streaming reads, original/archive digests and round-trip comparisons establish identity and coverage, not privacy or mathematical correctness. The revised manifest, rather than regrouping the old preparation state, is the release identity.

The inspection source has no network calls, subprocess launch, mathematical-worker imports or mutation of original data. It reads original bytes, decoded gzip/Zstandard/NPZ members and read-only immutable SQLite schema, then writes a local report under `.publication/`. Payload reads use 1-MiB chunks; NumPy headers are capped at 64 KiB and each inspected stream at 64 GiB. ZIP directory/report/schema metadata are not independently fixed-size, so the publication watchdog remains necessary; this is not a strict allocator-memory proof. JSON receives a prefix check, not full grammatical validation.

Its token/key/credential patterns are heuristic: for example, the bound source does not have a PGP-private-key pattern and cannot recognize every encoding, password or personal datum. A negative marker report cannot establish universal non-sensitivity. The local report also retains decoded prefixes, array headers and database schema; it should remain local until reviewed before any public sharing. No matched credential values are intentionally printed in progress/findings.

**Result:** the complete revised inventory has mathematical/research names and the concrete literature-PDF exception has been removed. No additional administrative/private *path class* was identified. Payload sensitivity, all externally derived data permissions, actual scanner findings and upload authorization remain separate questions. The automatic rejection has not been bypassed by this metadata review.
