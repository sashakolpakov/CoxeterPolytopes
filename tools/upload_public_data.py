#!/usr/bin/env python3
"""Serial release-data upload with local hashes and GitHub digest checks.

Publication I/O only. The release must already exist at the intended commit.
Existing assets are verified and reused, never overwritten or deleted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def gh_json(endpoint):
    return json.loads(subprocess.check_output(["gh", "api", endpoint]))


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024**2), b""):
            result.update(block)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--plan", action="store_true")
    mode.add_argument("--upload", action="store_true")
    args = parser.parse_args()
    manifest_path = ROOT / "repro/public-data-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    repo, tag = manifest["repository"], manifest["release_tag"]
    if repo != "sashakolpakov/CoxeterPolytopes" or tag != "research-data-2026-10-21":
        raise ValueError("Unexpected publication destination")
    parts = manifest["parts"]
    if not 0 < len(parts) < 100 or len({p["name"] for p in parts}) != len(parts):
        raise ValueError("Invalid asset inventory")
    if args.plan:
        print(json.dumps({"repository": repo, "tag": tag, "archives": len(parts),
                          "bytes": sum(p["bytes"] for p in parts), "network": False}))
        return
    if os.environ.get("COXETER_PUBLICATION_IO_GUARD") != "1":
        raise ValueError("Use tools/run_publication_io.py for monitored upload")
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                     text=True).strip()
    release = gh_json(f"repos/{repo}/releases/tags/{tag}")
    target = gh_json(f"repos/{repo}/commits/{tag}")["sha"]
    if target != commit or release["draft"]:
        raise ValueError("Release must be published at the current source commit")
    endpoint = f"repos/{repo}/releases/{release['id']}/assets?per_page=100"
    report = {"publication_io_only": True, "repository": repo, "tag": tag,
              "commit": commit, "manifest_sha256": digest(manifest_path), "assets": []}
    report_path = ROOT / ".publication/upload_report.json"
    for part in parts:
        name = part["name"]
        if Path(name).name != name or not name.endswith(".tar.gz"):
            raise ValueError("Unsafe archive filename")
        local = ROOT / ".publication/data" / name
        if local.is_symlink() or local.stat().st_size != part["bytes"] or digest(local) != part["sha256"]:
            raise ValueError(f"Local asset changed: {name}")
        assets = gh_json(endpoint)
        if len(assets) >= 100:
            raise ValueError("Asset listing requires pagination; refusing incomplete checks")
        matches = [a for a in assets if a["name"] == name]
        if len(matches) > 1:
            raise ValueError(f"Duplicate remote asset: {name}")
        if not matches:
            subprocess.run(["gh", "release", "upload", tag, str(local),
                            "--repo", repo], check=True)
            matches = [a for a in gh_json(endpoint) if a["name"] == name]
        if len(matches) != 1:
            raise ValueError(f"Uploaded asset missing: {name}")
        asset = matches[0]
        if asset["size"] != part["bytes"] or asset.get("digest") != "sha256:" + part["sha256"] or asset["state"] != "uploaded":
            raise ValueError(f"Remote digest/size/state differs: {name}")
        report["assets"].append({"name": name, "bytes": asset["size"],
                                  "digest": asset["digest"], "id": asset["id"],
                                  "url": asset["browser_download_url"]})
        temp = report_path.with_suffix(".tmp")
        temp.write_text(json.dumps(report, indent=2) + "\n")
        temp.replace(report_path)
        print(json.dumps({"verified_asset": name, "completed": len(report["assets"]),
                          "total": len(parts)}), flush=True)
    print(json.dumps({"status": "PASS_REMOTE_ASSET_DIGESTS", "assets": len(parts),
                      "report": str(report_path.relative_to(ROOT))}), flush=True)


if __name__ == "__main__":
    main()
