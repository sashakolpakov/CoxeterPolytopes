#!/usr/bin/env python3
"""One-time, source-bound removal of an upstream PDF from unuploaded assets.

Preserves the old manifest/archive locally, changes only archive006, and
retains every other scientific input. No upload or original-file mutation.
"""
import gzip
import hashlib
import json
import os
from pathlib import Path

import prepare_public_data as packing

ROOT = Path(__file__).resolve().parents[1]
PREDECESSOR = "b4373eef8bc086e51fc8d42b98e789bbe9342f4a9fbf85f93288610c4cadcccf"
PDF = "data/references/fresh_all_ideal_seven_simplex_literature_20261017_chen_labbe.pdf"
PDF_HASH = "1d7a1fd1f20421f9cc7fdeb9b3deacdc05f6ce814e37b7849cd377f72d7f5e99"


def main():
    if os.environ.get("COXETER_PUBLICATION_IO_GUARD") != "1":
        raise SystemExit("Use the publication I/O watchdog")
    source = ROOT / "repro/public-data-manifest.json"
    before = source.read_bytes()
    if hashlib.sha256(before).hexdigest() != PREDECESSOR:
        raise ValueError("Expected exact unuploaded predecessor manifest")
    manifest = json.loads(before)
    matches = [(i, p, r) for i, p in enumerate(manifest["parts"])
               for r in p["files"] if r["path"] == PDF]
    if len(matches) != 1:
        raise ValueError("Unexpected reference inventory")
    index, old, row = matches[0]
    if old["name"] != "research-data-006.tar.gz" or row["bytes"] != 1619220 or row["sha256"] != PDF_HASH:
        raise ValueError("Reference identity differs")
    original_archive = ROOT / ".publication/data" / old["name"]
    if packing.digest_file(original_archive) != old["sha256"]:
        raise ValueError("Old archive differs")
    preserved = ROOT / "results/source_archive/public_data_before_content_review_20261021"
    preserved.mkdir(parents=True, exist_ok=True)
    saved = preserved / "public-data-manifest.json.gz"
    with saved.open("xb") as target:
        target.write(gzip.compress(before, compresslevel=9, mtime=0))
    if gzip.decompress(saved.read_bytes()) != before:
        raise ValueError("Manifest preservation differs")
    work = ROOT / ".publication/revised-data"
    work.mkdir(exist_ok=True)
    successor = packing.pack([r for r in old["files"] if r["path"] != PDF], work / old["name"])
    if successor["files"] != [r for r in old["files"] if r["path"] != PDF]:
        raise ValueError("Retained original identities differ")
    history = ROOT / ".publication/pre-content-review"
    history.mkdir(exist_ok=True)
    if (history / old["name"]).exists():
        raise ValueError("Old archive already preserved")
    original_archive.rename(history / old["name"])
    (work / old["name"]).rename(original_archive)
    manifest["parts"][index] = {key: successor[key] for key in ("name", "bytes", "sha256", "files")}
    manifest["original_bytes"] -= row["bytes"]
    manifest["release_tag"] = "research-data-2026-10-21-r1"
    after = (json.dumps(manifest, indent=2) + "\n").encode()
    source.write_bytes(after)
    receipt = {"status": "REFERENCE_PDF_EXCLUDED_BEFORE_ANY_ASSET_UPLOAD",
               "old_manifest_sha256": PREDECESSOR,
               "new_manifest_sha256": hashlib.sha256(after).hexdigest(),
               "preserved_manifest": str(saved.relative_to(ROOT)),
               "preserved_manifest_gzip_sha256": hashlib.sha256(saved.read_bytes()).hexdigest(),
               "excluded_reference": row,
               "unchanged_archive_count": len(manifest["parts"]) - 1,
               "changed_archive_before": {k: old[k] for k in ("name", "bytes", "sha256")},
               "changed_archive_after": {k: successor[k] for k in ("name", "bytes", "sha256")},
               "remaining_files": sum(len(p["files"]) for p in manifest["parts"]),
               "original_bytes": manifest["original_bytes"],
               "compressed_bytes": sum(p["bytes"] for p in manifest["parts"]),
               "local_original_pdf_preserved": packing.digest_file(ROOT / PDF) == PDF_HASH,
               "archive006_full_roundtrip_verified": successor["roundtrip_verified"],
               "no_network_or_original_scientific_mutation": True,
               "fresh_packing_note": "The original preparation state is historical; asset identities for this revision are fixed by the revised public manifest, not by rediscovering group boundaries."}
    (ROOT / "results/public_data_reference_exclusion_20261021.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main()
