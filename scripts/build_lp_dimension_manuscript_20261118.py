#!/usr/bin/env python3
"""Typeset only the new LP companion, inside the existing bounded runner.

The caller supplies --parent-lock to run_bounded_research_check_20261018.py.
This owner runs no mathematical checker. Temporary TeX products are removed;
the final PDF is written beside its source only after successful compilation.
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
import time

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "manuscripts/linear_programming_dimension_bounds"
INPUTS = ("main.tex", "finite_volume.tex", "energy_framework.tex", "compact_cases.tex")
STARTED = time.monotonic()


def require(ok, message):
    if not ok:
        raise ValueError(message)
    if time.monotonic() - STARTED >= 53:
        raise TimeoutError("Typesetting entry allowance exhausted")


def bound(path):
    with path.open("rb") as stream:
        sha = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
            "sha256": sha}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--replace", action="store_true")
    parser.add_argument("--preview-pages", help="One to six one-based pages for visual review")
    args = parser.parse_args()
    report_path = (ROOT / args.report).resolve()
    require(report_path.parent == ROOT / "results" and not report_path.exists(),
            "Choose a fresh results report")
    require(FOLDER.resolve() == FOLDER, "Aliased manuscript directory")
    target = FOLDER / "main.pdf"
    require(target.resolve() == target and (args.replace or not target.exists()),
            "Existing PDF needs the explicit replacement flag")
    engine = shutil.which("pdflatex")
    require(engine is not None, "pdflatex is unavailable")
    bindings, payloads = {}, {}
    for name in INPUTS:
        path = FOLDER / name
        require(path.resolve() == path and 0 < path.stat().st_size <= 512 * 1024,
                "Missing, aliased or oversized manuscript input")
        bindings[name] = bound(path)
        payloads[name] = path.read_bytes()
        require(hashlib.sha256(payloads[name]).hexdigest() == bindings[name]["sha256"],
                "Source changed during read")
    previous = bound(target) if target.exists() else None
    environment = dict(os.environ)
    for key in ("TEXINPUTS", "BIBINPUTS", "BSTINPUTS", "TEXMFOUTPUT"):
        environment.pop(key, None)
    passes, previews = [], []
    with tempfile.TemporaryDirectory(prefix="coxeter-lp-paper-", dir="/private/tmp") as temp:
        folder = Path(temp)
        for name, payload in payloads.items():
            (folder / name).write_bytes(payload)
        command = [engine, "-no-shell-escape", "-interaction=nonstopmode",
                   "-halt-on-error", "-file-line-error", "-recorder", "main.tex"]
        for iteration in range(1, 5):
            remaining = 53 - (time.monotonic() - STARTED)
            require(remaining > 5, "Insufficient remaining typesetting allowance")
            result = subprocess.run(command, cwd=folder, env=environment,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    timeout=min(12, remaining - 3))
            passes.append({"pass": iteration, "exit_code": result.returncode,
                           "stdout_bytes": len(result.stdout)})
            if result.returncode:
                raise RuntimeError(result.stdout[-9000:].decode(errors="replace"))
            log = (folder / "main.log").read_text(errors="replace")
            require(len(log) <= 2 * 1024**2, "Typesetting log exceeds cap")
            rerun = bool(re.search(r"Rerun to get|Label\(s\) may have changed|rerunfilecheck Warning|Rerun to get cross-references|Table widths have changed", log))
            if iteration >= 3 and not rerun:
                break
        require(not rerun, "Cross-references or table widths still require a rerun")
        require(not re.search(r"(?:Reference|Citation).*undefined|undefined references|undefined citations", log),
                "Undefined reference or citation in final pass")
        warnings = [line.strip() for line in log.splitlines()
                    if re.search(r"Warning|Overfull|Underfull", line)]
        pdf = folder / "main.pdf"
        require(pdf.is_file() and 0 < pdf.stat().st_size <= 16 * 1024**2,
                "Missing or oversized PDF")
        from pypdf import PdfReader
        reader = PdfReader(pdf)
        pages = len(reader.pages)
        require(1 <= pages <= 160, "Unexpected manuscript page count")
        text = "\n".join(page.extract_text() or "" for page in reader.pages)
        require("Linear programming certificates" in text
                and "377" in text and "arithmetic" in text.lower(),
                "PDF text does not identify the expected manuscript")
        fonts = set()
        for page in reader.pages:
            resources = page.get("/Resources", {})
            if hasattr(resources, "get_object"):
                resources = resources.get_object()
            font_map = resources.get("/Font", {})
            if hasattr(font_map, "get_object"):
                font_map = font_map.get_object()
            for font in font_map.values():
                fonts.add(str(font.get_object().get("/Subtype", "unknown")))
        if args.preview_pages:
            require(re.fullmatch(r"\d+(?:,\d+){0,5}", args.preview_pages),
                    "Choose one to six preview page numbers")
            selected = [int(item) for item in args.preview_pages.split(",")]
            require(len(selected) == len(set(selected))
                    and all(1 <= page <= pages for page in selected), "Invalid preview pages")
            from pypdf import PdfWriter
            from pypdf.generic import DecodedStreamObject
            preview_folder = ROOT / ".publication" / (report_path.stem + "_previews")
            preview_folder.mkdir(exist_ok=False, parents=True)
            for number in selected:
                writer = PdfWriter()
                writer.add_page(reader.pages[number - 1])
                canvas = writer.pages[0]
                content = canvas.get_contents()
                width, height = float(canvas.mediabox.width), float(canvas.mediabox.height)
                white = DecodedStreamObject()
                white.set_data(f"q 1 1 1 rg 0 0 {width} {height} re f Q\n".encode()
                               + (content.get_data() if content is not None else b""))
                canvas.replace_contents(white)
                single = folder / f"preview{number}.pdf"
                writer.write(single)
                target_png = preview_folder / f"page{number}.png"
                rendered = subprocess.run(
                    ["/usr/bin/sips", "-s", "format", "png", "-Z", "1400",
                     str(single), "--out", str(target_png)],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=8)
                require(rendered.returncode == 0 and target_png.is_file()
                        and 0 < target_png.stat().st_size <= 2 * 1024**2,
                        "PDF preview failed or exceeds cap")
                previews.append({"page": number, **bound(target_png)})
        for name, binding in bindings.items():
            require(bound(FOLDER / name) == binding, "Source changed during build")
        data = pdf.read_bytes()
    partial = target.with_suffix(".pdf.partial")
    require(not partial.exists(), "Owned PDF partial already exists")
    try:
        with partial.open("xb") as stream:
            stream.write(data)
        os.replace(partial, target)
    finally:
        partial.unlink(missing_ok=True)
    for name, binding in bindings.items():
        require(bound(FOLDER / name) == binding, "Source changed before final build seal")
    report = {
        "schema": "lp-dimension-manuscript-typesetting-v1",
        "status": "PASS_TYPESETTING_ONLY", "sources": bindings,
        "builder": bound(Path(__file__).resolve()), "pdf": bound(target),
        "previous_pdf": previous, "passes": passes, "pages": pages,
        "font_subtypes": sorted(fonts), "warnings": warnings,
        "previews": previews,
        "undefined_references_or_citations": False,
        "temporary_auxiliaries_removed": True,
        "mathematical_checks_or_solvers_run": False,
        "elapsed_seconds": round(time.monotonic() - STARTED, 3),
    }
    with report_path.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps({key: report[key] for key in (
        "status", "pages", "pdf", "warnings", "elapsed_seconds")}))


if __name__ == "__main__":
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError("Typesetting deadline")))
    signal.setitimer(signal.ITIMER_REAL, 53)
    main()
