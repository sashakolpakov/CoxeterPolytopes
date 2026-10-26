# Reproduction and verification

Start with the **[detailed reproduction guide](https://sashakolpakov.github.io/CoxeterPolytopes/reproducibility.html)**
and [mathematical algorithm contracts](https://sashakolpakov.github.io/CoxeterPolytopes/algorithms.html).
Plans launch no scientific child:

```sh
python3 -B scripts/verify_manuscript.py --list
python3 -B scripts/verify_manuscript.py finite57 --mode full --pdf --plan
python3 -B scripts/fetch_public_data.py --plan
```

| Recorded catalogue edition | ID |
|---|---|
| Finite volume 57 | `finite57` |
| Compact 23 | `compact_farthest_bounds` |
| Compact arithmetic 20 | `compact20` |
| Sharp right-angled 8/ideal 4 | `right_angled_sharp` |
| Ideal RA4 volume and structure | `ideal_ra4_volume_and_structure` |
| Compact RA4 constraints | `compact_ra4` |

These IDs select their recorded editions. Later compact 17/arithmetic 15 and
illustrated revisions have separate source/build records. Finite57's current
64-shard/259-unit conjunction passed at explicit premises; no combined fresh
146-job replay or complete transitive regeneration is claimed. Its
[evidence guide](../docs/finite_volume57_reproduce.md) preserves exact scopes.

Bulk assets (936 files in 33 archives) remain **unuploaded**. See
[availability](PUBLIC_DATA.md), the [manifest](public-data-manifest.json) and
[third-party retrieval](EXTERNAL_INPUTS.md). Byte restoration is not mathematical
verification. Pinned packages are in [requirements-proof](requirements-proof-20261018.txt);
Lean configuration and [separate compilation/audit instructions](../formal/README.md)
are preserved.

All scientific, Lean and TeX work runs serially under the shared lock and
monitoring. Admit at 2% free disk and 15% available RAM; stop below 1% disk
or 6% RAM, and start cleanup below 2% disk. Use exact byte ratios against
total filesystem capacity and physical RAM; record totals and the policy hash.
Scientific/TeX batches use 59 seconds, 1 GiB and one thread; Lean uses one
module, 2048 MiB and a separate audit. Retain sources, certificates, PDFs,
failures and receipts; remove only owned unused build products.
Do not rerun successful unchanged jobs for a documentation update.

Build Sphinx with warnings as errors; it launches no mathematical checks.
The [site-build instructions](https://sashakolpakov.github.io/CoxeterPolytopes/reproducibility.html#building-only-this-site)
and GitHub Pages workflow describe deployment.
