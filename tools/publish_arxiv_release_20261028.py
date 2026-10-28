#!/usr/bin/env python3
"""Publish only the fixed, reviewed manuscript assets after an audited push."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
STAGING = ROOT / ".publication/arxiv/manuscript_release_20261028"
REPOSITORY = "sashakolpakov/CoxeterPolytopes"
TAG = "arxiv-papers-2026-10-28"
ASSETS = {
    "right_angled_sharp_dimension_bounds_arxiv_20261028.tex": "672b148a0733f8d4a9597c13a7e3b5807646b9ed55df6b9cba2d9ad147005747",
    "right_angled_sharp_dimension_bounds_arxiv_20261028.pdf": "a4bca1c0acd85842792ee7057a7e3d217ce265ed89fa6d276195595a565d2daa",
    "ideal_ra4_volume_and_structure_arxiv_20261028.tex": "ddc3431d435803c490b620f8c8652f32d383de79553950e2adc0bad7dd9385a3",
    "ideal_ra4_volume_and_structure_arxiv_20261028.pdf": "368d2fb3b84b39df1880ba1b25a29318e9775772e813db52cd08f1a821c42be7",
    "arxiv_publication_release_manifest_20261028.json": None,
    "SHA256SUMS": None,
}
DATA_TAGS = {
    "research-data-2026-10-21-r1": "b8d2cb73a395417e7cb1fb7b4a56a2040cc4a595",
    "research-data-2026-10-21": "32b68355cdb13ca8629f3fa8b98d424dd70472c7",
}
PUBLIC_PAPERS = {
    "right_angled_sharp": (
        "manuscripts/revisions/right_angled_sharp_dimension_bounds_arxiv_20261028.tex",
        "manuscripts/build/arxiv_sharp_dimension_bounds_8_and_4_20261028_layout_final.pdf",
        "results/arxiv_sharp_submission_build_20261028.json"),
    "ideal_ra4_volume_and_structure": (
        "manuscripts/revisions/ideal_ra4_volume_and_structure_arxiv_20261028.tex",
        "manuscripts/build/arxiv_ideal_four_dimensional_volume_and_structure_20261028_historical_context.pdf",
        "results/arxiv_ideal_submission_build_20261028.json"),
}


def checked(path: Path) -> Path:
    path.relative_to(ROOT)
    for candidate in (path, *path.parents):
        if candidate.is_symlink():
            raise ValueError(f"Refusing symlink: {candidate}")
        if candidate == ROOT:
            break
    return path


def payload(path: Path, cap: int = 4 * 1024**2) -> bytes:
    path = checked(path)
    if not path.is_file() or not 0 < path.stat().st_size < cap:
        raise ValueError(f"Unexpected publication file: {path}")
    return path.read_bytes()


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def run(arguments: list[str]) -> bytes:
    return subprocess.check_output(arguments, cwd=ROOT, timeout=90)


def git(*arguments: str) -> bytes:
    return run(["git", "-c", "core.packedGitLimit=16m", "-c", "core.packedGitWindowSize=1m",
                "-c", "core.deltaBaseCacheLimit=8m", *arguments])


def api(endpoint: str, missing: bool = False):
    result = subprocess.run(["gh", "api", endpoint], cwd=ROOT, capture_output=True, timeout=90)
    if result.returncode:
        if missing and b"HTTP 404" in result.stderr:
            return None
        raise RuntimeError("GitHub API failed: " + result.stderr[-1200:].decode(errors="replace"))
    return json.loads(result.stdout)


def release_assets(release: dict) -> list:
    return api(f"repos/{REPOSITORY}/releases/{release['id']}/assets?per_page=100")


def data_snapshot() -> dict:
    snapshot = {}
    for tag, expected in DATA_TAGS.items():
        release = api(f"repos/{REPOSITORY}/releases/tags/{tag}")
        actual = api(f"repos/{REPOSITORY}/commits/{tag}")["sha"]
        if actual != expected or release_assets(release):
            raise ValueError("Existing bulk-data tag or empty-asset state differs: " + tag)
        snapshot[tag] = {"release_id": release["id"], "commit": actual, "assets": 0}
    return snapshot


def matches(asset: dict, expected: dict) -> bool:
    return (asset.get("name") == expected["name"] and asset.get("size") == expected["bytes"]
            and asset.get("digest") == "sha256:" + expected["sha256"]
            and asset.get("state") == "uploaded")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--publish", action="store_true", required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    if os.environ.get("COXETER_PUBLICATION_IO_GUARD") != "1":
        raise SystemExit("Run through tools/run_publication_io.py")
    if not re.fullmatch(r"[0-9a-f]{40}", args.commit):
        raise ValueError("Supply the exact audited 40-digit source commit")
    inventory = json.loads(payload(STAGING / "release_assets.json", 128 * 1024))
    if inventory.get("repository") != REPOSITORY or inventory.get("tag") != TAG:
        raise ValueError("Unexpected manuscript publication destination")
    rows = inventory["assets"]
    if len(rows) != 6 or {row["name"] for row in rows} != set(ASSETS):
        raise ValueError("Only the exact six manuscript assets may be uploaded")
    total = 0
    for row in rows:
        data = payload(STAGING / row["name"])
        total += len(data)
        if len(data) != row["bytes"] or sha(data) != row["sha256"]:
            raise ValueError("Local asset differs from prepared inventory: " + row["name"])
        if ASSETS[row["name"]] and sha(data) != ASSETS[row["name"]]:
            raise ValueError("Asset differs from the accepted paper edition: " + row["name"])
    if total >= 4 * 1024**2:
        raise ValueError("Manuscript release exceeds the fixed 4-MiB aggregate bound")
    manifest = json.loads(payload(STAGING / "arxiv_publication_release_manifest_20261028.json", 128 * 1024))
    if manifest.get("repository") != REPOSITORY or manifest.get("release_tag") != TAG or manifest.get("bulk_assets_included") != 0:
        raise ValueError("Unexpected public manuscript manifest")
    papers = manifest.get("papers", [])
    if len(papers) != 2 or {paper.get("id") for paper in papers} != set(PUBLIC_PAPERS):
        raise ValueError("The manifest must describe exactly the two accepted papers")
    for paper in papers:
        actual_paths = tuple(paper[key]["path"] for key in ("source", "pdf", "build_record"))
        if actual_paths != PUBLIC_PAPERS[paper["id"]]:
            raise ValueError("Unexpected public manuscript or build-record path")
        stem = Path(paper["source"]["path"]).stem
        for key, extension in (("source", ".tex"), ("pdf", ".pdf")):
            name = stem + extension
            row = next(item for item in rows if item["name"] == name)
            if (paper[key]["sha256"] != row["sha256"]
                    or paper[key]["sha256"] != ASSETS[name]
                    or paper[key]["bytes"] != row["bytes"]):
                raise ValueError("Public manifest and uploaded paper asset disagree: " + name)
    sums = "".join(f"{row['sha256']}  {row['name']}\n" for row in sorted(rows, key=lambda row: row["name"]) if row["name"] != "SHA256SUMS").encode()
    if payload(STAGING / "SHA256SUMS") != sums:
        raise ValueError("SHA256SUMS differs from the five manuscript/manifest assets")
    if git("rev-parse", "HEAD").decode().strip() != args.commit:
        raise ValueError("Working HEAD is not the exact audited commit")
    for paper in manifest["papers"]:
        for key in ("source", "pdf", "build_record"):
            item = paper[key]
            content = git("show", args.commit + ":" + item["path"])
            if len(content) != item["bytes"] or sha(content) != item["sha256"]:
                raise ValueError("A public paper input does not match the committed source")
    committed_manifest = git("show", args.commit + ":results/arxiv_publication_release_manifest_20261028.json")
    if committed_manifest != payload(STAGING / "arxiv_publication_release_manifest_20261028.json"):
        raise ValueError("The release manifest is not the exact committed manifest")
    if api(f"repos/{REPOSITORY}/commits/{args.commit}")["sha"] != args.commit:
        raise ValueError("The audited source commit is not available on GitHub")
    before = data_snapshot()
    tag_ref = api(f"repos/{REPOSITORY}/git/ref/tags/{TAG}", missing=True)
    tag_commit = api(f"repos/{REPOSITORY}/commits/{TAG}") if tag_ref is not None else None
    if tag_commit is not None and tag_commit["sha"] != args.commit:
        raise ValueError("Existing manuscript tag points to another commit; refusing to retag")
    release = api(f"repos/{REPOSITORY}/releases/tags/{TAG}", missing=True)
    if release:
        if tag_commit is None or release["draft"] or release["prerelease"]:
            raise ValueError("Existing manuscript release is not a published matching-tag edition")
        expected = {row["name"]: row for row in rows}
        existing = release_assets(release)
        if len(existing) != len({asset["name"] for asset in existing}):
            raise ValueError("Duplicate existing release asset names")
        if any(asset["name"] not in expected or not matches(asset, expected[asset["name"]]) for asset in existing):
            raise ValueError("Existing assets differ; no asset will be overwritten or deleted")
    notes = checked(STAGING / "release_notes.md")
    payload(notes, 128 * 1024)
    if release is None:
        run(["gh", "release", "create", TAG, "--repo", REPOSITORY, "--target", args.commit,
             "--title", "Right-angled arXiv manuscript editions — 28 October 2026",
             "--notes-file", str(notes), "--latest"])
    else:
        run(["gh", "release", "edit", TAG, "--repo", REPOSITORY,
             "--title", "Right-angled arXiv manuscript editions — 28 October 2026",
             "--notes-file", str(notes), "--draft=false", "--latest"])
    release = api(f"repos/{REPOSITORY}/releases/tags/{TAG}")
    verified = []
    for row in rows:
        current = [asset for asset in release_assets(release) if asset["name"] == row["name"]]
        if not current:
            run(["gh", "release", "upload", TAG, str(STAGING / row["name"]), "--repo", REPOSITORY])
            current = [asset for asset in release_assets(release) if asset["name"] == row["name"]]
        if len(current) != 1 or not matches(current[0], row):
            raise ValueError("Remote asset digest/size/state differs: " + row["name"])
        verified.append({"name": row["name"], "bytes": row["bytes"], "sha256": row["sha256"],
                         "id": current[0]["id"], "url": current[0]["browser_download_url"]})
        print(json.dumps({"verified_manuscript_asset": row["name"], "completed": len(verified)}), flush=True)
    if api(f"repos/{REPOSITORY}/commits/{TAG}")["sha"] != args.commit:
        raise ValueError("Published manuscript tag does not resolve to the audited commit")
    final_assets = release_assets(release)
    if len(final_assets) != 6 or {asset["name"] for asset in final_assets} != set(ASSETS):
        raise ValueError("Final release has an unexpected asset inventory")
    after = data_snapshot()
    if after != before:
        raise ValueError("Existing bulk-data releases changed during manuscript publication")
    receipt = {"status": "PASS_MANUSCRIPT_ONLY_REMOTE_RELEASE_DIGESTS", "repository": REPOSITORY,
               "tag": TAG, "commit": args.commit, "url": release["html_url"], "assets": verified,
               "bulk_data_releases_unchanged": after, "bulk_assets_uploaded": 0}
    output = checked(STAGING / "remote_release_receipt.json")
    if output.exists() and payload(output) != (json.dumps(receipt, indent=2) + "\n").encode():
        raise ValueError("Refusing to replace a different existing remote release receipt")
    if not output.exists():
        with output.open("xb") as stream:
            stream.write((json.dumps(receipt, indent=2) + "\n").encode())
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main()
