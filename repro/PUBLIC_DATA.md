# Bulk research inputs

The Git repository contains the code, manuscripts, six current PDFs, proof
notes, verification records and smaller certificates. Larger scientific
files are distributed as assets of the
[research-data release](https://github.com/sashakolpakov/CoxeterPolytopes/releases/tag/research-data-2026-10-21-r1).
The versioned `public-data-manifest.json` records their original paths,
lengths and SHA-256 hashes, and the corresponding archive hashes.

**Availability:** no bulk assets have been uploaded. The revised manifest
describes 936 research files in 33 prepared archives; publication of that
payload is awaiting explicit authorization. The source repository, six
paper PDFs and documentation site are public. The restore commands below
require the release assets and do not yet establish a successful public
restore.

Inspect the local plan first. This command performs no network requests or
mathematical checks:

```sh
python3 -B scripts/fetch_public_data.py --plan
```

Restore all missing large inputs sequentially:

```sh
python3 -B scripts/fetch_public_data.py --restore
```

Or restore one exact input named in a proof's plan:

```sh
python3 -B scripts/fetch_public_data.py --restore \
  --path results/fresh_ambient57_payment_frontier_20261018.json
```

Use repeated `--path` options for several inputs. Selection is by exact
path, not by a shell glob. An archive may contain additional files, but
only the selected paths are restored. Already downloaded assets can be
used with `--from-directory /path/to/assets`.

The helper requires Python 3.11 or later and uses the standard library.
It checks the compressed archive, permits only manifest-listed regular
files, checks every restored file, and refuses to overwrite a different
existing file. It rejects path traversal and symlink destinations.
Downloads and restoration run one archive at a time with 1-MiB buffers;
downloaded archives are removed after use unless `--keep-cache` is given.
Plan enough disk for the missing original data and the largest selected
archive; restoration also maintains a 20-GiB free-space floor at admission.
The full original bulk collection is about 17 GiB.

Archive restoration proves byte identity only. It does not establish the
correctness, completeness or geometric interpretation of an enumeration.
The [verification guide](../docs/manuscript_verification.md) and
[finite57 reproduction guide](../docs/finite_volume57_reproduce.md) retain
their separate mathematical and execution boundaries. In particular, the
full finite57 selection remains a **plan**, not a completed safe replay.

Administrative session history, local build caches and temporary outputs
are outside the public dataset. They remain preserved in the local working
archive; they are not inputs to a new mathematical theorem.
Third-party downloads under `data/external/` and the Chen–Labbé reference PDF
also remain local. The [external-input guide](EXTERNAL_INPUTS.md) and their
per-source README files record provenance and retrieval references; those
downloads are not bundled in the release assets. The six primary papers'
selected verification plans do not reference those excluded downloads.

Release availability and checksums must be verified before a public restore
is claimed successful. Preparing this guide alone does not upload data.
