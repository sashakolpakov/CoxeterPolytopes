#!/usr/bin/env python3
"""Typeset or preview the complete dimension-bound arXiv source, serially.

This reuses the isolated submission builder, with bounded limits appropriate
to the three complete case appendices. It does not run mathematical checks.
The existing 59-second/1-GiB/single-thread research wrapper is unchanged.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import signal

from pypdf import PdfReader
import build_arxiv_submission_20261028 as submission

PAPER = "hyperbolic_coxeter_dimension_bounds"
ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / ".publication/arxiv/dimension_bounds_work"
SOURCE_LIMIT = 4 * 1024**2
PDF_LIMIT = 32 * 1024**2
PAGE_LIMIT = 800


def describe_pdf(path: Path) -> dict:
    if not 0 < path.stat().st_size <= PDF_LIMIT:
        raise ValueError("Large submission PDF exceeds 32 MiB")
    fonts, part_pages = set(), {}
    with path.open("rb") as stream:
        reader = PdfReader(stream)
        count = len(reader.pages)
        if not 1 <= count <= PAGE_LIMIT:
            raise ValueError("Unexpected combined-manuscript page count")
        for page in reader.pages:
            resources = page.get("/Resources", {})
            resources = resources.get_object() if hasattr(resources, "get_object") else resources
            font_map = resources.get("/Font", {})
            font_map = font_map.get_object() if hasattr(font_map, "get_object") else font_map
            for font in font_map.values():
                fonts.add(str(font.get_object().get("/Subtype", "unknown")))
        for name, destination in reader.named_destinations.items():
            if name.startswith("part."):
                part_pages[name] = reader.get_destination_page_number(destination) + 1
    return {"pages": count, "font_subtypes": sorted(fonts),
            "part_start_pages": part_pages,
            "pdf_inspection": "All page resources inspected; full-page text extraction omitted."}


def main() -> None:
    if "COXETER_CHECK_OWNER" not in os.environ:
        raise SystemExit("Run through scripts/run_proof_batch_20261024.py check")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--render-pages")
    parser.add_argument("--replace-current", action="store_true")
    args = parser.parse_args()
    if not re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,31}", args.tag):
        raise ValueError("Unsafe tag")

    def interrupted(signum, frame):
        raise RuntimeError(f"Interrupted bounded submission I/O: signal {signum}")
    for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP, signal.SIGALRM):
        signal.signal(signum, interrupted)
    submission.SOURCE_LIMIT = SOURCE_LIMIT
    submission.PDF_LIMIT = PDF_LIMIT
    submission.REVIEW = WORK
    submission.pdf_description = describe_pdf
    submission.checked(WORK).mkdir(parents=True, exist_ok=True)
    stem = f"arxiv_{PAPER}_20261028_{args.tag}"
    pdf = submission.checked(ROOT / "manuscripts/build" / (stem + ".pdf"))
    if args.replace_current:
        pdf = submission.checked(submission.CURRENT_PDFS[PAPER])
    report = submission.checked(WORK / (stem + "_build.json"))
    if args.render_pages:
        submission.render(stem, pdf, report, args.render_pages)
        return
    profile = {"status": "TYPESETTING_CONFIGURATION_ONLY",
               "worker": str(Path(__file__).relative_to(ROOT)),
               "worker_sha256": submission.digest(Path(__file__)),
               "shared_builder_sha256": submission.digest(Path(submission.__file__)),
               "source_limit_bytes": SOURCE_LIMIT, "pdf_limit_bytes": PDF_LIMIT,
               "page_limit": PAGE_LIMIT, "per_pass_timeout_seconds": 18,
               "max_tex_passes": 3, "bibtex_passes": 1, "scientific_checks_run": False,
               "resource_policy": "Unchanged existing exclusive 59-second/1-GiB/one-thread wrapper"}
    submission.write_report(WORK / (stem + "_profile.json"), profile)
    submission.build(PAPER, stem, pdf, report, replace_current=args.replace_current)
    print(json.dumps(profile, indent=2))


if __name__ == "__main__":
    main()
