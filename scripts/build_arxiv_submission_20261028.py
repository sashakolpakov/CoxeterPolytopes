#!/usr/bin/env python3
"""Build or preview one submission copy inside the existing bounded wrapper.

This is typesetting and visual-inspection I/O, never mathematical verification.
The public run_proof_batch check entry supplies the exclusive lock, 59-second
watchdog, process-tree RSS cap, resource gates and single-thread environment.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import tempfile

from pypdf import PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject

ROOT = Path(__file__).resolve().parents[1]
REVIEW = ROOT / ".publication/arxiv/presentation_review_20261028"
PAPERS = (
    "sharp_dimension_bounds_8_and_4",
    "ideal_four_dimensional_volume_and_structure",
)
SOURCE_LIMIT = 256 * 1024
PDF_LIMIT = 8 * 1024**2


def checked(path: Path) -> Path:
    path.relative_to(ROOT)
    for ancestor in (path, *path.parents):
        if ancestor.is_symlink():
            raise ValueError(f"Refusing symlink: {ancestor}")
        if ancestor == ROOT:
            break
    return path


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_report(path: Path, report: dict) -> None:
    with checked(path).open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")


def pdf_description(path: Path) -> dict:
    if not 0 < path.stat().st_size <= PDF_LIMIT:
        raise ValueError("PDF exceeds the bounded inspection size")
    figure_pages, fonts = [], set()
    with path.open("rb") as stream:
        reader = PdfReader(stream)
        if not 1 <= len(reader.pages) <= 128:
            raise ValueError("Unexpected submission page count")
        for number, page in enumerate(reader.pages, 1):
            captions = sorted(set(map(int, re.findall(
                r"Figure\s*(\d+)\s*:", page.extract_text() or ""))))
            if captions:
                figure_pages.append({"page": number, "figures": captions})
            resources = page.get("/Resources", {})
            resources = resources.get_object() if hasattr(resources, "get_object") else resources
            font_map = resources.get("/Font", {})
            font_map = font_map.get_object() if hasattr(font_map, "get_object") else font_map
            for font in font_map.values():
                fonts.add(str(font.get_object().get("/Subtype", "unknown")))
        return {"pages": len(reader.pages), "figure_pages": figure_pages,
                "font_subtypes": sorted(fonts)}


def build(paper: str, stem: str, pdf: Path, report_path: Path) -> None:
    source = checked(ROOT / "arxiv" / paper / "main.tex")
    if not source.is_file():
        raise ValueError(f"Missing submission source: {source}")
    with source.open("rb") as stream:
        payload = stream.read(SOURCE_LIMIT + 1)
    if not payload or len(payload) > SOURCE_LIMIT:
        raise ValueError("Submission source exceeds 256 KiB")
    original = hashlib.sha256(payload).hexdigest()
    if pdf.exists() or report_path.exists():
        raise ValueError("Build output already exists; choose a fresh --tag")
    engine = shutil.which("pdflatex")
    if engine is None:
        raise ValueError("pdflatex is required")
    report = {
        "status": "INCOMPLETE_TYPESETTING_ONLY", "source": str(source.relative_to(ROOT)),
        "source_sha256": original, "source_bytes": len(payload),
        "pdf": str(pdf.relative_to(ROOT)), "passes": [],
        "builder_sha256": digest(Path(__file__)),
        "mathematical_verification_claimed": False,
        "isolated_input": "Only the exact main.tex bytes are copied to a fresh temporary directory.",
        "temporary_auxiliaries_removed": False,
    }
    try:
        with tempfile.TemporaryDirectory(prefix="coxeter-arxiv-", dir="/private/tmp") as temporary:
            folder = Path(temporary)
            (folder / "main.tex").write_bytes(payload)
            environment = dict(os.environ)
            for name in ("TEXINPUTS", "BIBINPUTS", "BSTINPUTS", "TEXMFOUTPUT"):
                environment.pop(name, None)
            command = [engine, "-no-shell-escape", "-interaction=nonstopmode",
                       "-halt-on-error", "-file-line-error", "-recorder", "main.tex"]
            for iteration in range(3):
                # Inherit the wrapper's owned process group; no detached child.
                result = subprocess.run(command, cwd=folder, env=environment,
                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                        timeout=18)
                report["passes"].append({"pass": iteration + 1,
                                         "exit_code": result.returncode,
                                         "stdout_bytes": len(result.stdout)})
                if result.returncode:
                    raise RuntimeError(result.stdout[-7000:].decode(errors="replace"))
                log_path = folder / "main.log"
                if log_path.stat().st_size > 2 * 1024**2:
                    raise RuntimeError("Unexpectedly large TeX log")
                log = log_path.read_text(errors="replace")
                needs_rerun = bool(re.search(
                    r"Rerun to get|Label\(s\) may have changed|rerunfilecheck Warning", log))
                if iteration >= 1 and not needs_rerun:
                    break
            if needs_rerun:
                raise RuntimeError("Cross-references still request a rerun after three passes")
            report["warnings"] = [line.strip() for line in log.splitlines()
                                  if re.search(r"Warning|Overfull|Underfull", line)]
            report["undefined_references_or_citations"] = bool(re.search(
                r"(?:Reference|Citation).*undefined|undefined references|undefined citations", log))
            if report["undefined_references_or_citations"]:
                raise RuntimeError("Undefined references or citations in final TeX pass")
            fls = (folder / "main.fls").read_text(errors="replace")
            unexpected = []
            for line in fls.splitlines():
                if line.startswith("INPUT "):
                    item = Path(line[6:])
                    resolved = (item if item.is_absolute() else folder / item).resolve()
                    if resolved.is_relative_to(ROOT):
                        unexpected.append(str(resolved.relative_to(ROOT)))
            report["repository_inputs_read_by_tex"] = sorted(set(unexpected))
            if unexpected:
                raise RuntimeError("TeX read repository files outside the isolated submission")
            built = folder / "main.pdf"
            report.update(pdf_description(built))
            report["pdf_sha256"] = digest(built)
            report["pdf_bytes"] = built.stat().st_size
            report["source_figure_environments"] = len(re.findall(
                rb"\\begin\{figure\*?\}", payload))
            if digest(source) != original:
                raise RuntimeError("Submission source changed during compilation")
            with checked(pdf).open("xb") as output, built.open("rb") as incoming:
                shutil.copyfileobj(incoming, output, length=1024 * 1024)
        report["temporary_auxiliaries_removed"] = True
        report["status"] = "PASS_ISOLATED_TEX_BUILD_ONLY"
        write_report(report_path, report)
        print(json.dumps(report, indent=2))
    except BaseException as error:
        report["status"] = "FAILED_ISOLATED_TEX_BUILD"
        report["error"] = str(error)
        report["temporary_auxiliaries_removed"] = True
        if not report_path.exists():
            write_report(report_path, report)
        print(json.dumps(report, indent=2))
        raise


def render(stem: str, pdf: Path, build_report: Path, pages_text: str) -> None:
    if not re.fullmatch(r"\d+(?:,\d+){0,5}", pages_text):
        raise ValueError("Choose one to six one-based comma-separated pages")
    pages = [int(value) for value in pages_text.split(",")]
    if len(set(pages)) != len(pages):
        raise ValueError("Duplicate preview pages")
    saved = json.loads(checked(build_report).read_text())
    if saved.get("status") != "PASS_ISOLATED_TEX_BUILD_ONLY" or digest(pdf) != saved["pdf_sha256"]:
        raise ValueError("Existing PDF does not match its successful build report")
    if not 0 < pdf.stat().st_size <= PDF_LIMIT:
        raise ValueError("PDF exceeds the bounded rendering size")
    if not Path("/usr/bin/sips").is_file():
        raise ValueError("This optional preview mode requires macOS /usr/bin/sips")
    output = Path(tempfile.mkdtemp(prefix=stem + "_preview_", dir=REVIEW))
    rows = []
    with pdf.open("rb") as source, tempfile.TemporaryDirectory(
            prefix="coxeter-arxiv-preview-", dir="/private/tmp") as temporary:
        reader = PdfReader(source)
        if any(not 1 <= page <= len(reader.pages) for page in pages):
            raise ValueError("Preview page outside the PDF")
        for number in pages:
            writer = PdfWriter()
            writer.add_page(reader.pages[number - 1])
            canvas = writer.pages[0]
            content = canvas.get_contents()
            stream = DecodedStreamObject()
            width, height = float(canvas.mediabox.width), float(canvas.mediabox.height)
            stream.set_data(f"q 1 1 1 rg 0 0 {width} {height} re f Q\n".encode()
                            + (content.get_data() if content is not None else b""))
            canvas.replace_contents(stream)
            single = Path(temporary) / f"page{number}.pdf"
            writer.write(single)
            target = output / f"page{number}.png"
            converted = subprocess.run(
                ["/usr/bin/sips", "-s", "format", "png", "-Z", "1400",
                 str(single), "--out", str(target)], capture_output=True, timeout=8)
            if converted.returncode or not 0 < target.stat().st_size <= 2 * 1024**2:
                raise RuntimeError("Page preview failed or exceeded 2 MiB")
            rows.append({"page": number, "png": str(target.relative_to(ROOT)),
                         "bytes": target.stat().st_size, "sha256": digest(target)})
    report = {"status": "RENDERED_FOR_VISUAL_REVIEW_NOT_PROOF",
              "pdf": str(pdf.relative_to(ROOT)), "pdf_sha256": saved["pdf_sha256"],
              "source_sha256": saved["source_sha256"], "pages": rows,
              "temporary_page_pdfs_removed": True}
    write_report(output / "render.json", report)
    print(json.dumps(report, indent=2))


def main() -> None:
    if "COXETER_CHECK_OWNER" not in os.environ:
        raise SystemExit("Run with scripts/run_proof_batch_20261024.py check -- this script")
    def interrupted(signum, frame):
        raise RuntimeError(f"Interrupted bounded submission I/O: signal {signum}")
    for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP, signal.SIGALRM):
        signal.signal(signum, interrupted)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paper", choices=PAPERS)
    parser.add_argument("--tag", default="", help="Fresh safe ASCII version tag for retries")
    parser.add_argument("--render-pages", help="Render existing PDF only; one-based comma-separated pages")
    args = parser.parse_args()
    if args.tag and not re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,31}", args.tag):
        raise ValueError("Unsafe --tag")
    stem = "arxiv_" + args.paper + "_20261028" + ("_" + args.tag if args.tag else "")
    pdf = checked(ROOT / "manuscripts/build" / (stem + ".pdf"))
    report = checked(REVIEW / (stem + "_build.json"))
    checked(REVIEW).mkdir(parents=True, exist_ok=True)
    checked(pdf.parent).mkdir(exist_ok=True)
    if args.render_pages:
        render(stem, pdf, report, args.render_pages)
    else:
        build(args.paper, stem, pdf, report)


if __name__ == "__main__":
    main()
