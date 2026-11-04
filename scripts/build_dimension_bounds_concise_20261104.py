#!/usr/bin/env python3
"""Typeset the concise edition inside the existing exclusive bounded runner.

This worker performs publication typesetting only. It uses the established
submission PDF inspector and previewer, and runs no mathematical checkers.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import tempfile
import zipfile

import build_arxiv_submission_20261028 as publication

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "arxiv/hyperbolic_coxeter_dimension_bounds"
OUTPUT = ROOT / "release/dimension_bounds"
DIST = ROOT / "dist"


def main() -> None:
    if "COXETER_CHECK_OWNER" not in os.environ:
        raise SystemExit("Use scripts/run_proof_batch_20261024.py check -- this worker")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--render-pages")
    args = parser.parse_args()
    if not re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,31}", args.tag):
        raise ValueError("Use a fresh safe ASCII tag")

    def interrupted(signum, frame):
        raise RuntimeError(f"Interrupted typesetting: signal {signum}")

    for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP, signal.SIGALRM):
        signal.signal(signum, interrupted)

    publication.checked(OUTPUT).mkdir(parents=True, exist_ok=True)
    publication.checked(DIST).mkdir(exist_ok=True)
    pdf = SOURCE / "main.pdf"
    receipt = OUTPUT / f"build_{args.tag}.json"
    if args.render_pages:
        publication.REVIEW = OUTPUT
        publication.render("concise_" + args.tag, pdf, receipt, args.render_pages)
        return
    if receipt.exists():
        raise ValueError("Build receipt exists; choose a fresh tag")

    source, bibliography = SOURCE / "main.tex", SOURCE / "references.bib"
    inputs = {}
    for path in (source, bibliography):
        publication.checked(path)
        if not 0 < path.stat().st_size <= 256 * 1024:
            raise ValueError(f"Unexpected source size: {path}")
        inputs[str(path.relative_to(ROOT))] = publication.digest(path)
    engine, bibtex = shutil.which("pdflatex"), shutil.which("bibtex")
    if not engine or not bibtex:
        raise RuntimeError("pdflatex and bibtex are required")
    report = {
        "status": "INCOMPLETE_TYPESETTING_ONLY",
        "inputs_sha256": inputs,
        "source": str(source.relative_to(ROOT)),
        "source_sha256": publication.digest(source),
        "pdf": str(pdf.relative_to(ROOT)),
        "builder_sha256": publication.digest(Path(__file__)),
        "shared_inspector_sha256": publication.digest(Path(publication.__file__)),
        "passes": [],
        "mathematical_verification_claimed": False,
        "full_proof_addenda_snapshot": "36eeaddf646eff25d4990c9ab730a20ada21f80f",
        "temporary_auxiliaries_removed": False,
    }
    try:
        with tempfile.TemporaryDirectory(prefix="coxeter-concise-", dir="/private/tmp") as work:
            folder = Path(work)
            shutil.copyfile(source, folder / "main.tex")
            shutil.copyfile(bibliography, folder / "references.bib")
            environment = dict(os.environ)
            for name in ("TEXINPUTS", "BIBINPUTS", "BSTINPUTS", "TEXMFOUTPUT"):
                environment.pop(name, None)
            tex = [engine, "-no-shell-escape", "-interaction=nonstopmode",
                   "-halt-on-error", "-file-line-error", "-recorder", "main.tex"]
            pipeline = [("pdflatex", tex), ("bibtex", [bibtex, "main"]),
                        ("pdflatex", tex), ("pdflatex", tex)]
            for number, (program, command) in enumerate(pipeline, 1):
                result = subprocess.run(command, cwd=folder, env=environment,
                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                        timeout=5 if program == "bibtex" else 18)
                report["passes"].append({"pass": number, "program": program,
                                         "exit_code": result.returncode})
                if result.returncode:
                    raise RuntimeError(result.stdout[-7000:].decode(errors="replace"))
            log = (folder / "main.log").read_text(errors="replace")
            biblog = (folder / "main.blg").read_text(errors="replace")
            report["warnings"] = [line.strip() for line in log.splitlines()
                                  if re.search(r"Warning|Overfull|Underfull", line)]
            report["bibtex_warnings"] = [line.strip() for line in biblog.splitlines()
                                         if "Warning" in line]
            report["undefined_references_or_citations"] = bool(re.search(
                r"(?:Reference|Citation).*undefined|undefined references|undefined citations", log))
            if report["undefined_references_or_citations"] or report["bibtex_warnings"]:
                raise RuntimeError("Undefined reference or bibliography warning")
            if re.search(r"Rerun to get|Label\(s\) may have changed|rerunfilecheck Warning", log):
                raise RuntimeError("References require another pass")
            unexpected = []
            for line in (folder / "main.fls").read_text().splitlines():
                if line.startswith("INPUT "):
                    path = Path(line[6:])
                    resolved = (path if path.is_absolute() else folder / path).resolve()
                    if resolved.is_relative_to(ROOT):
                        unexpected.append(str(resolved.relative_to(ROOT)))
            if unexpected:
                raise RuntimeError(f"Unexpected repository TeX inputs: {unexpected}")
            built = folder / "main.pdf"
            report.update(publication.pdf_description(built))
            report["pdf_sha256"] = publication.digest(built)
            report["pdf_bytes"] = built.stat().st_size
            report["source_figure_environments"] = 0
            for name, original in inputs.items():
                if publication.digest(ROOT / name) != original:
                    raise RuntimeError("A source changed while typesetting")
            for incoming, target in [(folder / "main.bbl", SOURCE / "main.bbl"),
                                     (built, pdf)]:
                publication.checked(target)
                with tempfile.NamedTemporaryFile(prefix=".concise-", dir=target.parent,
                                                 delete=False) as staged:
                    stage = Path(staged.name)
                    with incoming.open("rb") as stream:
                        shutil.copyfileobj(stream, staged, length=1024 * 1024)
                try:
                    os.replace(stage, target)
                finally:
                    stage.unlink(missing_ok=True)
            report["compiled_bibliography_sha256"] = publication.digest(SOURCE / "main.bbl")
            package = DIST / "CoxeterPolytopes-dimension-bounds-arxiv.zip"
            names = ["main.tex", "references.bib", "main.bbl"]
            with tempfile.NamedTemporaryFile(prefix=".concise-upload-", dir=DIST,
                                             delete=False) as temporary_zip:
                stage = Path(temporary_zip.name)
            try:
                with zipfile.ZipFile(stage, "w", compression=zipfile.ZIP_DEFLATED) as archive:
                    for name in names:
                        archive.write(SOURCE / name, arcname=name)
                with zipfile.ZipFile(stage) as archive:
                    for name in names:
                        if archive.read(name) != (SOURCE / name).read_bytes():
                            raise RuntimeError("Submission ZIP byte round trip failed")
                os.replace(stage, package)
            finally:
                stage.unlink(missing_ok=True)
            package.with_name(package.name + ".sha256").write_text(
                publication.digest(package) + "  " + package.name + "\n")
            report["submission_zip"] = {"path": str(package.relative_to(ROOT)),
                                         "sha256": publication.digest(package),
                                         "entries": names,
                                         "byte_round_trip_verified": True}
        report["temporary_auxiliaries_removed"] = True
        report["status"] = "PASS_ISOLATED_TEX_BUILD_ONLY"
    except BaseException as error:
        report["status"] = "FAILED_ISOLATED_TEX_BUILD"
        report["error"] = str(error)
        report["temporary_auxiliaries_removed"] = True
        publication.write_report(receipt, report)
        print(json.dumps(report, indent=2))
        raise
    publication.write_report(receipt, report)
    manifest_path = ROOT / "release/manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text())
        for paper in manifest["papers"]:
            if paper["id"] == "dimension_bounds":
                paper["submission_status"] = "TYPESET_PACKAGE_PREPARED"
                paper["submission_source"] = {
                    "path": str(source.relative_to(ROOT)),
                    "bytes": source.stat().st_size,
                    "sha256": report["source_sha256"]}
                paper["bibliography"] = {
                    "path": str(bibliography.relative_to(ROOT)),
                    "bytes": bibliography.stat().st_size,
                    "sha256": publication.digest(bibliography)}
                catalogue = SOURCE / "ADDENDA.md"
                paper["addendum_map"] = {
                    "path": str(catalogue.relative_to(ROOT)),
                    "bytes": catalogue.stat().st_size,
                    "sha256": publication.digest(catalogue)}
                paper["pdf"] = {"path": str(pdf.relative_to(ROOT)),
                                 "bytes": report["pdf_bytes"],
                                 "sha256": report["pdf_sha256"]}
                paper["pages"] = report["pages"]
                paper["current_pdf_role"] = "Concise consolidated article; complete proofs and computations in the GitHub supplement"
                paper.pop("future_pdf", None)
                paper.pop("commit_and_push_held_pending_expanded_proof_review", None)
                paper.pop("extra_length_allowed_to_retain_essential_proofs", None)
                paper["local_upload_package"] = {
                    "path": str(package.relative_to(ROOT)),
                    "bytes": package.stat().st_size,
                    "sha256": publication.digest(package)}
                paper["package_entries"] = names
                paper["compiled_bibliography"] = {
                    "path": str((SOURCE / "main.bbl").relative_to(ROOT)),
                    "sha256": report["compiled_bibliography_sha256"]}
                paper["typesetting_receipt"] = str(receipt.relative_to(ROOT))
                break
        manifest["publication_state"] = "Concise submission package prepared; arXiv submission and bulk upload have not been performed."
        with tempfile.NamedTemporaryFile(mode="w", prefix=".manifest-",
                                         dir=manifest_path.parent, delete=False) as output:
            staged_manifest = Path(output.name)
            json.dump(manifest, output, indent=2)
            output.write("\n")
        try:
            os.replace(staged_manifest, manifest_path)
        finally:
            staged_manifest.unlink(missing_ok=True)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
