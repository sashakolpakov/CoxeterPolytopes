#!/usr/bin/env python3
"""Prepare only the two reviewed manuscript editions; never upload or build."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAGING = ROOT / ".publication/arxiv/manuscript_release_20261028"
REVIEW = ROOT / ".publication/arxiv/presentation_review_20261028"
REPOSITORY = "sashakolpakov/CoxeterPolytopes"
TAG = "arxiv-papers-2026-10-28"
URL = f"https://github.com/{REPOSITORY}/releases/tag/{TAG}"
MANIFEST = "results/arxiv_publication_release_manifest_20261028.json"
PAPERS = (
    ("right_angled_sharp", "Sharp dimension bounds for right-angled hyperbolic polytopes",
     "sharp_dimension_bounds_8_and_4", "layout_final", "right_angled_sharp_dimension_bounds",
     "672b148a0733f8d4a9597c13a7e3b5807646b9ed55df6b9cba2d9ad147005747",
     "a4bca1c0acd85842792ee7057a7e3d217ce265ed89fa6d276195595a565d2daa", 7, 2,
     "results/arxiv_sharp_submission_build_20261028.json"),
    ("ideal_ra4_volume_and_structure", "Volumes and octahedral chains of ideal right-angled hyperbolic four-polytopes",
     "ideal_four_dimensional_volume_and_structure", "historical_context", "ideal_ra4_volume_and_structure",
     "ddc3431d435803c490b620f8c8652f32d383de79553950e2adc0bad7dd9385a3",
     "368d2fb3b84b39df1880ba1b25a29318e9775772e813db52cd08f1a821c42be7", 19, 5,
     "results/arxiv_ideal_submission_build_20261028.json"),
)


def checked(path: Path) -> Path:
    path.relative_to(ROOT)
    for candidate in (path, *path.parents):
        if candidate.is_symlink():
            raise ValueError(f"Refusing symlink: {candidate}")
        if candidate == ROOT:
            break
    return path


def read(path: Path, cap: int = 4 * 1024**2) -> bytes:
    path = checked(path)
    if not path.is_file() or not 0 < path.stat().st_size < cap:
        raise ValueError(f"Unexpected file size or type: {path}")
    return path.read_bytes()


def sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def save(path: Path, payload: bytes) -> None:
    checked(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if read(path) != payload:
            raise ValueError(f"Refusing to replace a different existing file: {path}")
    else:
        with path.open("xb") as output:
            output.write(payload)
    if read(path) != payload:
        raise RuntimeError(f"Copy readback failed: {path}")


def identity(path: Path) -> dict:
    payload = read(path)
    return {"path": str(path.relative_to(ROOT)), "bytes": len(payload), "sha256": sha(payload)}


def encoded(value: dict) -> bytes:
    return (json.dumps(value, indent=2) + "\n").encode()


def main() -> None:
    if os.environ.get("COXETER_PUBLICATION_IO_GUARD") != "1":
        raise SystemExit("Run through tools/run_publication_io.py")
    records, assets = [], []
    for suite, title, folder, tag, basename, source_sha, pdf_sha, pages, figures, export in PAPERS:
        report_path = REVIEW / f"arxiv_{folder}_20261028_{tag}_build.json"
        report_bytes = read(report_path, 128 * 1024)
        report = json.loads(report_bytes)
        source = ROOT / "arxiv" / folder / "main.tex"
        pdf = ROOT / "manuscripts/build" / f"arxiv_{folder}_20261028_{tag}.pdf"
        source_bytes, pdf_bytes = read(source, 256 * 1024), read(pdf)
        if (report.get("status") != "PASS_ISOLATED_TEX_BUILD_ONLY"
                or report.get("warnings") != []
                or report.get("undefined_references_or_citations") is not False
                or report.get("repository_inputs_read_by_tex") != []
                or report.get("source") != str(source.relative_to(ROOT))
                or report.get("pdf") != str(pdf.relative_to(ROOT))
                or report.get("pages") != pages
                or report.get("source_figure_environments") != figures):
            raise ValueError(f"The accepted source/PDF/build record differs: {suite}")
        if not (sha(source_bytes) == source_sha == report["source_sha256"]
                and sha(pdf_bytes) == pdf_sha == report["pdf_sha256"]):
            raise ValueError(f"Source or PDF hash differs: {suite}")
        edition = ROOT / "manuscripts/revisions" / f"{basename}_arxiv_20261028.tex"
        exported = ROOT / export
        save(edition, source_bytes)
        save(exported, report_bytes)
        record = {"id": suite, "title": title, "source": identity(edition),
                  "pdf": identity(pdf), "build_record": identity(exported),
                  "pages": pages, "inline_figures": figures,
                  "author": "Alexander Kolpakov", "ai_acknowledgment_included": True}
        records.append(record)
        for extension, payload in (("tex", source_bytes), ("pdf", pdf_bytes)):
            asset = STAGING / f"{basename}_arxiv_20261028.{extension}"
            save(asset, payload)
            assets.append({"name": asset.name, "bytes": len(payload), "sha256": sha(payload)})
    manifest = {
        "schema_version": 1, "repository": REPOSITORY, "release_tag": TAG,
        "release_url": URL, "edition_date": "2026-10-28",
        "documentation_url": "https://sashakolpakov.github.io/CoxeterPolytopes/",
        "companion_snapshot": "04c9bab7fec573b7f315f070395780632acbfdf5",
        "papers": records,
        "verification_scope": {
            "source_review_and_isolated_pdf_builds": True,
            "new_mathematical_check_executions": 0, "new_lean_executions": 0,
            "end_to_end_geometric_formalization_claimed": False,
            "selected_small_check_inputs": "Tracked in the public repository; these two script suites need no unuploaded bulk archives.",
            "guarded_replay_platform": "The pinned repository launcher uses macOS resource probes.",
            "catalogue_pdf_scope": "The catalogue --pdf job builds its preserved publication edition, not these submission revisions.",
        },
        "six_frozen_editions_preserved": True,
        "arxiv_submission_performed": False, "submission_zips_created": False,
        "bulk_assets_included": 0,
        "bulk_availability": "The separate 936-file/33-archive research-data collection remains unuploaded and requires specific authorization.",
    }
    payload = encoded(manifest)
    save(ROOT / MANIFEST, payload)
    manifest_asset = STAGING / Path(MANIFEST).name
    save(manifest_asset, payload)
    assets.append({"name": manifest_asset.name, "bytes": len(payload), "sha256": sha(payload)})
    sums = "".join(f"{row['sha256']}  {row['name']}\n" for row in sorted(assets, key=lambda row: row["name"])).encode()
    save(STAGING / "SHA256SUMS", sums)
    assets.append({"name": "SHA256SUMS", "bytes": len(sums), "sha256": sha(sums)})
    save(STAGING / "release_assets.json", encoded({"repository": REPOSITORY, "tag": TAG, "assets": assets}))
    notes = """# Right-angled arXiv manuscript editions — 28 October 2026

This manuscript-only release supplies two self-contained TeX sources, their reading PDFs, a source/build manifest and SHA-256 checksums:

* **Sharp dimension bounds for right-angled hyperbolic polytopes:** 7 pages, 2 inline figures; finite-volume dimension 8 and entirely ideal dimension 4.
* **Volumes and octahedral chains of ideal right-angled hyperbolic four-polytopes:** 19 pages, 5 inline figures; the second ideal volume and unique reflected double, with a restricted octahedral-chain theorem.

Both editions name Alexander Kolpakov and include the requested AI-use acknowledgment. They are prepared for arXiv submission; this release does not claim an arXiv submission has occurred. No submission ZIPs are supplied.

The ideal-volume introduction situates the result after the Potyagailo–Vinberg compact facet question, Guo–Zheng's compact facet/volume minima (arXiv:2610.33053), and Kolpakov's ideal facet/volume minima.

Both sources were reviewed and compiled alone in fresh temporary directories with zero final warnings. This release updates exposition, figures, citations and reproduction mapping. It does not claim new mathematical-check or Lean execution, external refereeing, or an end-to-end formalization of the geometry. Published geometric inputs and finite-model interfaces remain explicit.

Companion suites are `right_angled_sharp` and `ideal_ra4_volume_and_structure`, pinned in the papers to public revision `04c9bab7`. Their small selected inputs are tracked; the guarded launcher currently uses macOS resource probes. The six frozen paper editions are preserved.

Documentation: https://sashakolpakov.github.io/CoxeterPolytopes/

The separate 936-file/33-archive bulk research-data payload remains unuploaded. Its existing data-release tags and manifests are unchanged; this manuscript release supplies no bulk archives and makes no public bulk-restoration claim.
"""
    save(STAGING / "release_notes.md", notes.encode())
    print(json.dumps({"status": "PREPARED_MANUSCRIPT_ONLY_RELEASE_NOT_UPLOADED",
                      "manifest": MANIFEST, "asset_directory": str(STAGING.relative_to(ROOT)),
                      "release_tag": TAG, "assets": assets}, indent=2))


if __name__ == "__main__":
    main()
