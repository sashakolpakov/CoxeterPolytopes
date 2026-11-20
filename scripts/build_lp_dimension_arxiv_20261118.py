#!/usr/bin/env python3
"""Build and package the LP article through the existing serial TeX worker.

Only publication typesetting and three-file source packaging are performed.
The caller supplies the existing 59-second/1-GiB/one-thread runner and lock.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import signal
import tempfile
import zipfile

import build_arxiv_submission_20261028 as publication

ROOT = Path(__file__).resolve().parents[1]
PAPER = "linear_programming_dimension_bounds"
SOURCE = ROOT / "arxiv" / PAPER
RELEASE = ROOT / "release" / PAPER
DIST = ROOT / "dist"
MEMBERS = ("main.tex", "references.bib", "main.bbl")
PACKAGE = DIST / "CoxeterPolytopes-linear-programming-dimension-bounds-arxiv.zip"


def bound(path: Path) -> dict:
    publication.checked(path)
    return {"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
            "sha256": publication.digest(path)}


def replace_pair(replacements: list[tuple[Path, Path]]) -> None:
    """Replace the small package/checksum pair, restoring old bytes on error."""
    previous = {target: target.read_bytes() if target.exists() else None
                for _, target in replacements}
    try:
        for staged, target in replacements:
            os.replace(staged, target)
    except BaseException:
        for target, data in previous.items():
            if data is None:
                target.unlink(missing_ok=True)
            else:
                with tempfile.NamedTemporaryFile(prefix=".lp-arxiv-rollback-", dir=DIST,
                                                 delete=False) as stream:
                    restore = Path(stream.name)
                    stream.write(data)
                try:
                    os.replace(restore, target)
                finally:
                    restore.unlink(missing_ok=True)
        raise


def package_sources(receipt: Path, build_receipt: Path, tag: str) -> None:
    report = {"status": "FAILED_LOCAL_SOURCE_PACKAGING", "typesetting_only": True,
              "arxiv_submission_or_network_upload_performed": False,
              "build_receipt": bound(build_receipt), "tag": tag,
              "worker": bound(Path(__file__)), "members": []}
    staged_paths = []
    try:
        build = json.loads(build_receipt.read_text())
        if (build["status"] != "PASS_ISOLATED_TEX_BUILD_ONLY" or build["warnings"]
                or build["bibtex_warnings"] or build["undefined_references_or_citations"]):
            raise RuntimeError("A clean isolated TeX/BibTeX build is required before packaging")
        if build["font_subtypes"] != ["/Type1"]:
            raise RuntimeError("The submission PDF must use Type1 fonts")
        bindings = [bound(SOURCE / name) for name in MEMBERS]
        expected_hashes = (build["source_sha256"], build["bibliography_sha256"],
                           build["compiled_bibliography_sha256"])
        if tuple(row["sha256"] for row in bindings) != expected_hashes:
            raise RuntimeError("Current source bytes do not match the successful build")
        if publication.digest(SOURCE / "main.pdf") != build["pdf_sha256"]:
            raise RuntimeError("Current PDF does not match the successful build")
        with tempfile.NamedTemporaryFile(prefix=".lp-arxiv-upload-", suffix=".tmp",
                                         dir=DIST, delete=False) as stream:
            staged = Path(stream.name)
        staged_paths.append(staged)
        with zipfile.ZipFile(staged, "w", compression=zipfile.ZIP_DEFLATED,
                             compresslevel=9) as archive:
            for name in MEMBERS:
                info = zipfile.ZipInfo(name, date_time=(2026, 11, 18, 0, 0, 0))
                info.create_system = 3
                info.external_attr = 0o100644 << 16
                info.compress_type = zipfile.ZIP_DEFLATED
                archive.writestr(info, (SOURCE / name).read_bytes(), compresslevel=9)
        if not 0 < staged.stat().st_size <= 1024**2:
            raise RuntimeError("Unexpected source ZIP size")
        with zipfile.ZipFile(staged) as archive:
            if tuple(archive.namelist()) != MEMBERS or archive.testzip() is not None:
                raise RuntimeError("Unexpected archive members or CRC error")
            for name in MEMBERS:
                if archive.read(name) != (SOURCE / name).read_bytes():
                    raise RuntimeError("Source ZIP byte round trip failed")
        checksum = PACKAGE.with_name(PACKAGE.name + ".sha256")
        with tempfile.NamedTemporaryFile(mode="w", prefix=".lp-arxiv-checksum-", suffix=".tmp",
                                         dir=DIST, delete=False) as stream:
            staged_checksum = Path(stream.name)
            stream.write(publication.digest(staged) + "  " + PACKAGE.name + "\n")
        staged_paths.append(staged_checksum)
        publication.checked(PACKAGE)
        publication.checked(checksum)
        if [bound(SOURCE / name) for name in MEMBERS] != bindings:
            raise RuntimeError("Sources changed during package construction")
        replace_pair([(staged, PACKAGE), (staged_checksum, checksum)])
        report.update(status="PASS_LOCAL_THREE_FILE_SOURCE_PACKAGE", members=bindings,
                      package=bound(PACKAGE), checksum=bound(checksum),
                      byte_round_trip_verified=True, crc_check_passed=True,
                      deterministic_metadata=True, pages=build["pages"],
                      pdf=bound(SOURCE / "main.pdf"),
                      excluded="PDF, README, TeX intermediates, scientific data and Lean products")
    except BaseException as error:
        report["error"] = str(error)
        raise
    finally:
        for path in staged_paths:
            path.unlink(missing_ok=True)
        publication.write_report(receipt, report)
    print(json.dumps({"status": report["status"], "pages": report["pages"],
                      "package": report["package"]}))


def main() -> None:
    if "COXETER_CHECK_OWNER" not in os.environ:
        raise SystemExit("Use scripts/run_bounded_research_check_20261018.py --parent-lock")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--render-pages", help="Only render the PDF bound by this tag's build receipt")
    args = parser.parse_args()
    if not re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,31}", args.tag):
        raise ValueError("Use a fresh safe ASCII tag")

    def interrupted(signum, frame):
        raise RuntimeError(f"Interrupted bounded LP publication job: signal {signum}")

    for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP, signal.SIGALRM):
        signal.signal(signum, interrupted)
    signal.setitimer(signal.ITIMER_REAL, 53)
    publication.checked(RELEASE).mkdir(parents=True, exist_ok=True)
    publication.checked(DIST).mkdir(exist_ok=True)
    receipt = RELEASE / f"build_{args.tag}.json"
    package_receipt = RELEASE / f"package_{args.tag}.json"
    publication.CURRENT_PDFS[PAPER] = SOURCE / "main.pdf"
    if args.render_pages:
        publication.REVIEW = ROOT / ".publication/arxiv" / PAPER
        publication.checked(publication.REVIEW).mkdir(parents=True, exist_ok=True)
        publication.render(PAPER + "_" + args.tag, SOURCE / "main.pdf", receipt,
                           args.render_pages)
        return
    if receipt.exists() or package_receipt.exists():
        raise ValueError("Build/package receipt exists; use a fresh tag")
    publication.build(PAPER, PAPER + "_" + args.tag, SOURCE / "main.pdf", receipt,
                      replace_current=True, source_directory=SOURCE)
    package_sources(package_receipt, receipt, args.tag)


if __name__ == "__main__":
    main()
