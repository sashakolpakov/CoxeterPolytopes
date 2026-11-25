#!/usr/bin/env python3
"""Execute one reviewed Hexagon API step under the publication-I/O watchdog.

No timers or background workers are installed. The default command only
reports missing declarations. Only the three files in a ready queue entry can
be uploaded. API responses are kept in an ignored, owner-only local receipt;
the credential is never an argument, a payload field or a recorded value.
"""
from __future__ import annotations

import argparse
import contextlib
import datetime as dt
import difflib
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
import unicodedata
import uuid

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ORIGIN = "https://hexagonmath.org"
KEY_PATH = ROOT / "hexagon_key.local"
RECEIPT_PATH = HERE / "execution.local"
PREVIEW_MARKER = HERE / "unresolved-preview.local"
COMMIT_MARKER = HERE / "unresolved-commit.local"
QUEUE_CONFIGURATION = HERE / "scheduled-execution.local"
LAUNCH_LABEL = "org.coxeterpolytopes.hexagon-submissions-20261123"
MAX_RESPONSE = 2 * 1024 * 1024
# Hexagon rejected the generated .bbl as an unsupported source-asset type.
# Keep it in the canonical release checks, but send only the TeX and BibTeX inputs.
HEXAGON_SOURCE_FILES = {"main.tex", "references.bib"}
PAPERS = (
    "right_angled_sharp",
    "ideal_ra4_volume_and_structure",
    "linear_programming_dimension_bounds",
)
ACTIONS = (
    "prepare", "create-draft", "draft-status", "create-upload", "put-part",
    "complete-upload", "source-check", "preview", "preview-status", "commit",
    "resolve-preview", "preview-artifacts", "inspect-preview", "run-queue",
    "inspect-queue", "apply-five-per-day", "refresh-statuses", "inspect-service", "completion-report",
    "inspect-update-api",
    "update-lp-vinberg",
)
DAILY_COMMIT_LIMIT = 5
DAILY_LIMIT_AUTHORIZATION = (
    'User instruction, 23 November 2026: the side hexagon agent may proceed '
    'with the remaining submissions since now we can submit up to 5 per day.'
)


class Stop(Exception):
    """A local check failed; the message never includes a credential."""


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode("utf-8")


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_json(path: Path) -> object:
    return json.loads(path.read_text(encoding="utf-8"))


def payload_snapshot(queue: dict) -> dict:
    fields = ("id", "metadata_template", "upload_files", "reference_pdf", "pages",
              "source_and_pdf_hashes_match_release_manifest", "matching_build_receipt")
    return {"source_manifest": queue["source_manifest"],
            "active_upload_names": sorted(HEXAGON_SOURCE_FILES),
            "papers": [{key: next(row for row in queue["entries"] if row["id"] == paper)[key]
                        for key in fields} for paper in PAPERS]}


def preparation(paper: str) -> tuple[dict, dict, dict[str, bytes], list[str]]:
    queue = load_json(HERE / "queue.json")
    entry = next(row for row in queue["entries"] if row["id"] == paper)
    manifest_descriptor = queue["source_manifest"]
    if manifest_descriptor["path"] != "release/manifest.json":
        raise Stop("Release manifest path differs from the reviewed path")
    manifest_bytes = (ROOT / manifest_descriptor["path"]).read_bytes()
    if (len(manifest_bytes) != manifest_descriptor["bytes"] or
            digest(manifest_bytes) != manifest_descriptor["sha256"]):
        raise Stop("Release manifest differs from the reviewed snapshot")
    release_paper = next(row for row in json.loads(manifest_bytes)["papers"]
                         if row["id"] == paper)
    descriptor_keys = ("path", "bytes", "sha256")
    released_sources = {row["path"]: row for row in
        (release_paper["submission_source"], release_paper["bibliography"],
         release_paper["compiled_bibliography"])}
    for row in entry["upload_files"]:
        if {key: row[key] for key in descriptor_keys} != released_sources.get(row["path"]):
            raise Stop("Upload descriptor differs from the release manifest")
    reference_pdf = {key: entry["reference_pdf"][key] for key in descriptor_keys}
    if reference_pdf != release_paper["pdf"]:
        raise Stop("Reference PDF differs from the release manifest")
    pdf_bytes = (ROOT / reference_pdf["path"]).read_bytes()
    if len(pdf_bytes) != reference_pdf["bytes"] or digest(pdf_bytes) != reference_pdf["sha256"]:
        raise Stop("Reference PDF bytes differ from the reviewed release")
    if entry.get("source_and_pdf_hashes_match_release_manifest") is not True:
        raise Stop("Selected paper has no matching source/PDF release record")
    metadata = load_json(HERE / entry["metadata_template"])
    files = {}
    expected = {"main.tex", "main.bbl", "references.bib"}
    if {row["relativePath"] for row in entry["upload_files"]} != expected:
        raise Stop("The reviewed three-file upload manifest has changed")
    for row in entry["upload_files"]:
        path = ROOT / row["path"]
        if path.is_symlink() or not path.resolve().is_relative_to(ROOT / "arxiv"):
            raise Stop("An upload source is outside the approved source directory")
        if path.name != row["relativePath"] or path.suffix == ".local":
            raise Stop("An upload source has an unexpected name")
        if path.stat().st_size > 4 * 1024 * 1024:
            raise Stop("An upload source exceeds the documented text-file limit")
        data = path.read_bytes()
        if len(data) != row["bytes"] or digest(data) != row["sha256"]:
            raise Stop("An upload source differs from the reviewed manifest")
        files[row["relativePath"]] = data
    missing = []
    if metadata.get("license") != "CC BY-SA 4.0":
        missing.append("user-selected license")
    if metadata.get("authorshipMode") != "mixed":
        missing.append("user-selected mixed text origin")
    if metadata.get("humanUnderstanding") != "some":
        missing.append("user-selected human understanding: some")
    if metadata.get("humanUnderstanding") == "some" and not metadata.get("provenance"):
        missing.append("scope of human understanding")
    for name in ("correctnessConfirmed", "communityInterestConfirmed",
                 "citationsConfirmed", "rightsConfirmed"):
        if metadata.get(name) is not True:
            missing.append(name)
    return entry, metadata, files, missing


def credential() -> str:
    info = KEY_PATH.lstat()
    if not stat.S_ISREG(info.st_mode) or info.st_mode & 0o077:
        raise Stop("Credential must be a regular owner-only file")
    if info.st_uid != os.getuid() or info.st_size > 4096:
        raise Stop("Credential ownership or size is unexpected")
    text = KEY_PATH.read_text(encoding="utf-8").strip()
    # Accept plain token, a Bearer header value, or one named assignment.
    # Never evaluate shell syntax or load arbitrary dotenv/environment files.
    match = re.fullmatch(r"(?:HEXAGON_API_KEY|HEXAGON_TOKEN)\s*=\s*(.+)", text)
    if match:
        text = match.group(1).strip()
        if len(text) >= 2 and text[0] == text[-1] and text[0] in "\"'":
            text = text[1:-1]
    if text.startswith("Bearer "):
        text = text[7:]
    if not re.fullmatch(r"[A-Za-z0-9._~+/=-]{16,4096}", text):
        raise Stop("Credential format is unsupported; contents were not displayed")
    return text


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def safe_response(value: object, token: str) -> object:
    if isinstance(value, dict):
        return {key.replace(token, "[redacted]"): ("[redacted]" if re.search(
            r"authorization|credential|secret|token|cookie", key, re.I)
            else safe_response(item, token)) for key, item in value.items()}
    if isinstance(value, list):
        return [safe_response(item, token) for item in value]
    if isinstance(value, str):
        return value.replace(token, "[redacted]")
    return value


def request(method: str, path: str, token: str, body: bytes | None,
            content_type: str, idempotency: str | None,
            artifact_path: Path | None = None) -> tuple[int, object]:
    if not re.fullmatch(r"/api/v1/submissions(?:/[A-Za-z0-9_./-]+)?", path):
        raise Stop("Request path is outside the fixed submission API")
    headers = {"Authorization": "Bearer " + token, "Accept": "application/json",
               "User-Agent": "CoxeterPolytopes-Hexagon-publication-client/1"}
    if body is not None:
        headers["Content-Type"] = content_type
    if idempotency:
        headers["Idempotency-Key"] = idempotency
    req = urllib.request.Request(ORIGIN + path, data=body, headers=headers,
                                 method=method)
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}), NoRedirect())
    try:
        try:
            response = opener.open(req, timeout=45)
        except urllib.error.HTTPError as error:
            response = error
        with response:
            raw = response.read(MAX_RESPONSE + 1)
            code = response.code
        if len(raw) > MAX_RESPONSE:
            raise Stop("Response exceeded the bounded receipt size")
        if artifact_path is not None and code == 200:
            if artifact_path.parent != HERE or artifact_path.suffix != ".local":
                raise Stop("Artifact output must remain in the private task directory")
            if token.encode("utf-8") in raw:
                raise Stop("Server response contained credential material; no artifact saved")
            if path.endswith("/pdf") and not raw.startswith(b"%PDF-"):
                raise Stop("Preview PDF response has no PDF header")
            if artifact_path.exists():
                if artifact_path.is_symlink() or artifact_path.read_bytes() != raw:
                    raise Stop("Existing private artifact differs from this response")
            else:
                fd = os.open(artifact_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
                with os.fdopen(fd, "wb") as artifact:
                    artifact.write(raw)
                    artifact.flush()
                    os.fsync(artifact.fileno())
            return code, {"artifact_path": str(artifact_path), "bytes": len(raw),
                          "sha256": digest(raw)}
        try:
            value = json.loads(raw) if raw else {}
        except (ValueError, UnicodeError):
            value = {"non_json_response": True, "response_bytes": len(raw)}
        return code, safe_response(value, token)
    except (urllib.error.URLError, TimeoutError, OSError):
        # A failed response does not prove that a mutation failed. No retry.
        raise Stop("Network result uncertain; inspect draft state before retrying") from None


def component(value: str | None, name: str) -> str:
    if not value or not re.fullmatch(r"[A-Za-z0-9_-]{1,128}", value):
        raise Stop("A valid " + name + " is required")
    return value


def explicit_utc_timestamp(value: str) -> dt.datetime:
    try:
        result = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    except (AttributeError, TypeError, ValueError):
        raise Stop("Timestamp is not a valid explicit UTC value") from None
    if result.tzinfo is None or result.utcoffset() != dt.timedelta(0):
        raise Stop("Timestamp has no explicit UTC timezone")
    return result


def bind_draft(draft: str, paper: str, metadata: dict) -> None:
    records = load_json(RECEIPT_PATH) if RECEIPT_PATH.exists() else []
    for row in records:
        response = row.get("response") or {}
        if (row.get("action") == "create-draft" and
                row.get("http_status") in {200, 201} and
                row.get("paper") == paper and
                row.get("metadata_sha256") == digest(canonical(metadata)) and
                row.get("payload_sha256") == digest(canonical(metadata)) and
                isinstance(response, dict) and response.get("draftId") == draft):
            return
    raise Stop("Draft is not bound to this paper and its current metadata")


def bind_verified_metadata(draft: str, paper: str, metadata: dict, entry: dict) -> None:
    records = load_json(RECEIPT_PATH) if RECEIPT_PATH.exists() else []
    expected_files = digest(canonical(entry["upload_files"]))
    for row in reversed(records):
        response = row.get("response") or {}
        if (row.get("action") != "draft-status" or row.get("http_status") != 200
                or row.get("paper") != paper or not isinstance(response, dict)
                or response.get("id") != draft):
            continue
        if (row.get("files_sha256") != expected_files or
                row.get("metadata_sha256") != digest(canonical(metadata)) or
                response.get("status") != "uploading"):
            raise Stop("Draft verification does not match the current paper/files")
        remote = json.loads(response["metadataJson"])
        for key, value in metadata.items():
            actual = response.get("arxivId") if key == "arxivId" else remote.get(key)
            if actual != value:
                raise Stop("Verified remote metadata differs from the selected paper")
        return
    raise Stop("Read and verify this draft's metadata before uploading")


def upload_descriptor(entry: dict, filename: str) -> dict:
    row = next(row for row in entry["upload_files"] if row["relativePath"] == filename)
    return {"kind": row["kind"], "filename": filename, "relativePath": filename,
            "contentType": row["contentType"], "sizeBytes": row["bytes"]}


def bind_upload(draft: str, upload: str, paper: str, filename: str,
                metadata: dict, entry: dict) -> None:
    records = load_json(RECEIPT_PATH)
    expected_payload = digest(canonical(upload_descriptor(entry, filename)))
    for row in reversed(records):
        response = row.get("response") or {}
        if (row.get("action") == "create-upload" and row.get("http_status") in {200, 201}
                and row.get("paper") == paper and row.get("manifest_file") == filename
                and row.get("api_path") == "/api/v1/submissions/" + draft + "/uploads"
                and row.get("metadata_sha256") == digest(canonical(metadata))
                and row.get("files_sha256") == digest(canonical(entry["upload_files"]))
                and row.get("payload_sha256") == expected_payload
                and isinstance(response, dict) and response.get("uploadId") == upload
                and isinstance(response.get("partSizeBytes"), int)
                and response["partSizeBytes"] >= upload_descriptor(entry, filename)["sizeBytes"]):
            return
    raise Stop("Upload identifier is not bound to this paper, draft and file")


def require_uploaded_part(draft: str, upload: str, paper: str, filename: str,
                          metadata: dict, entry: dict, data: bytes) -> None:
    for row in reversed(load_json(RECEIPT_PATH)):
        if (row.get("action") == "put-part" and row.get("http_status") == 200
                and row.get("paper") == paper and row.get("manifest_file") == filename
                and row.get("api_path") == "/api/v1/submissions/" + draft +
                    "/uploads/" + upload + "/parts/1"
                and row.get("metadata_sha256") == digest(canonical(metadata))
                and row.get("files_sha256") == digest(canonical(entry["upload_files"]))
                and row.get("payload_sha256") == digest(data)):
            return
    raise Stop("The exact source bytes have not completed their upload")


def require_completed_files(draft: str, paper: str, metadata: dict, entry: dict) -> None:
    expected = HEXAGON_SOURCE_FILES
    complete = set()
    for row in load_json(RECEIPT_PATH):
        if (row.get("action") != "complete-upload" or row.get("http_status") != 200
                or row.get("paper") != paper or row.get("manifest_file") not in expected
                or row.get("metadata_sha256") != digest(canonical(metadata))
                or row.get("files_sha256") != digest(canonical(entry["upload_files"]))):
            continue
        prefix = "/api/v1/submissions/" + draft + "/uploads/"
        path = row.get("api_path", "")
        if not path.startswith(prefix) or not path.endswith("/complete"):
            continue
        upload = component(path[len(prefix):-len("/complete")], "upload identifier")
        bind_upload(draft, upload, paper, row["manifest_file"], metadata, entry)
        complete.add(row["manifest_file"])
    if complete != expected:
        raise Stop("Both accepted source files must be completed before preview")


def validated_source_digest(draft: str, paper: str, metadata: dict, entry: dict) -> str:
    bind_verified_metadata(draft, paper, metadata, entry)
    records = load_json(RECEIPT_PATH)
    responses = {}
    paths = {"draft-status": "/api/v1/submissions/" + draft,
             "source-check": "/api/v1/submissions/" + draft + "/preview/check",
             "preview-status": "/api/v1/submissions/" + draft + "/preview"}
    for row in records:
        action = row.get("action")
        if (action in paths and row.get("http_status") == 200 and row.get("paper") == paper
                and row.get("api_path") == paths[action]
                and row.get("metadata_sha256") == digest(canonical(metadata))
                and row.get("files_sha256") == digest(canonical(entry["upload_files"]))):
            responses[action] = row["response"]
    if set(responses) != set(paths):
        raise Stop("Read draft, source-check and preview status before compiling")
    uploads = responses["draft-status"].get("uploads", [])
    if len(uploads) != len(HEXAGON_SOURCE_FILES):
        raise Stop("Remote draft has a different number of source files")
    found = set()
    for upload in uploads:
        filename = upload.get("relativePath")
        if filename not in HEXAGON_SOURCE_FILES or filename in found:
            raise Stop("Remote draft has an unexpected or duplicate source path")
        found.add(filename)
        expected = upload_descriptor(entry, filename)
        if any(upload.get(key) != expected[key] for key in
               ("filename", "relativePath", "kind", "contentType")):
            raise Stop("Remote source descriptor differs from the reviewed file")
        if (upload.get("status") != "complete" or
                upload.get("sizeBytes") != expected["sizeBytes"] or
                upload.get("expectedSizeBytes") != expected["sizeBytes"]):
            raise Stop("Remote source is incomplete or has a different byte count")
        bind_upload(draft, component(upload.get("id"), "upload identifier"),
                    paper, filename, metadata, entry)
        source = next(row for row in entry["upload_files"] if row["relativePath"] == filename)
        require_uploaded_part(draft, upload["id"], paper, filename,
                              metadata, entry, (ROOT / source["path"]).read_bytes())
        complete_path = "/api/v1/submissions/" + draft + "/uploads/" + upload["id"] + "/complete"
        if not any(row.get("action") == "complete-upload" and row.get("http_status") == 200
                   and row.get("paper") == paper and row.get("manifest_file") == filename
                   and row.get("api_path") == complete_path
                   and row.get("metadata_sha256") == digest(canonical(metadata))
                   and row.get("files_sha256") == digest(canonical(entry["upload_files"]))
                   and row.get("response") == {"uploadId": upload["id"],
                       "status": "complete", "sizeBytes": source["bytes"]}
                   for row in records):
            raise Stop("Current remote upload lacks its own successful completion receipt")
        # The live API currently returns sha256:null. The byte binding then
        # rests on successful exact-body PUT receipts and the immutable upload
        # identity; do not report a server-provided content hash in that case.
        if upload.get("sha256") not in {None, source["sha256"]}:
            raise Stop("Remote source hash differs from the reviewed source")
    check, preview = responses["source-check"], responses["preview-status"]
    tree = check.get("digest")
    if (not isinstance(tree, str) or not re.fullmatch(r"[a-f0-9]{64}", tree)
            or check.get("ready") is not True or check.get("treeChecked") is not True
            or check.get("truncated") is not False or check.get("warnings") != []
            or preview.get("draftId") != draft or preview.get("digest") != tree
            or preview.get("ready") is not True or preview.get("engine") != "pdflatex"
            or preview.get("engineRefusal") is not None):
        raise Stop("Remote static checks or compiler-source digest do not match")
    return tree


def _append_receipt_locked(receipt: dict) -> None:
    if RECEIPT_PATH.is_symlink():
        raise Stop("Local receipt must not be a symlink")
    records = load_json(RECEIPT_PATH) if RECEIPT_PATH.exists() else []
    records.append(receipt)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=HERE,
                prefix="receipt-", suffix=".local", delete=False) as stream:
            temporary = Path(stream.name)
            os.fchmod(stream.fileno(), 0o600)
            json.dump(records, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, RECEIPT_PATH)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def save_receipt(receipt: dict) -> None:
    fd = os.open(HERE / "receipts.lock.local", os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "a") as lock:
        os.fchmod(lock.fileno(), 0o600)
        fcntl.flock(lock, fcntl.LOCK_EX)
        _append_receipt_locked(receipt)


def write_private_json_once(path: Path, value: dict) -> None:
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


def exchange(args, entry: dict, metadata: dict, token: str, action: str,
             method: str, path: str, body: bytes | None = None,
             content_type: str = "application/json", idem: str | None = None,
             artifact_path: Path | None = None):
    base = {"request_id": str(uuid.uuid4()), "paper": args.paper,
            "action": action, "method": method, "manifest_file": args.file,
            "api_path": path, "metadata_sha256": digest(canonical(metadata)),
            "files_sha256": digest(canonical(entry["upload_files"])),
            "payload_sha256": digest(body) if body is not None else None}
    save_receipt(dict(base, utc_time=dt.datetime.now(dt.timezone.utc).isoformat(),
                      http_status=None, state="REQUEST_PENDING_OR_RESULT_UNKNOWN"))
    code, result = request(method, path, token, body, content_type, idem, artifact_path)
    save_receipt(dict(base, utc_time=dt.datetime.now(dt.timezone.utc).isoformat(),
                      http_status=code, state="HTTP_RESPONSE_RECEIVED", response=result))
    print(json.dumps({"paper": args.paper, "action": action, "http_status": code,
                      "receipt_path": str(RECEIPT_PATH)}), flush=True)
    return code, result


def preview_with_lock(args, entry: dict, metadata: dict) -> int:
    draft = component(args.draft_id, "draft identifier")
    base = "/api/v1/submissions/" + draft
    inherited_lock = getattr(args, "scientific_lock", None)
    lock_context = (contextlib.nullcontext(inherited_lock) if inherited_lock is not None
                    else (ROOT / "formal/.proof-batches/exclusive.lock").open("a"))
    with lock_context as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise Stop("Scientific/typesetting lane is occupied; preview deferred") from None
        if PREVIEW_MARKER.exists():
            raise Stop("A remote preview is unresolved; inspect its recorded state before any new preview")
        token = credential()
        for action, suffix in (("draft-status", ""), ("source-check", "/preview/check"),
                               ("preview-status", "/preview")):
            code, _ = exchange(args, entry, metadata, token, action, "GET", base + suffix)
            if code != 200:
                raise Stop("Live source validation failed; no preview requested")
        tree = validated_source_digest(draft, args.paper, metadata, entry)
        descriptor = {"paper": args.paper, "draft_id": draft, "source_digest": tree,
                      "metadata_sha256": digest(canonical(metadata)),
                      "files_sha256": digest(canonical(entry["upload_files"])),
                      "state": "REMOTE_PREVIEW_RESULT_UNRESOLVED",
                      "created_utc": dt.datetime.now(dt.timezone.utc).isoformat()}
        fd = os.open(PREVIEW_MARKER, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
        with os.fdopen(fd, "w", encoding="utf-8") as marker:
            json.dump(descriptor, marker, indent=2)
            marker.write("\n")
            marker.flush()
            os.fsync(marker.fileno())
        code, _ = exchange(args, entry, metadata, token, "preview", "POST", base + "/preview")
        if code not in {200, 202}:
            raise Stop("Preview was not accepted; inspect the actual response")
        deadline = time.monotonic() + 240
        while time.monotonic() < deadline:
            time.sleep(5)
            code, value = exchange(args, entry, metadata, token, "preview-status", "GET", base + "/preview")
            if code != 200 or not isinstance(value, dict):
                raise Stop("Preview state is uncertain; coordinator must check before releasing lane")
            if value.get("draftId") != draft or value.get("digest") != tree:
                raise Stop("Preview source digest changed; keep publication on hold")
            preview = value.get("preview") or {}
            state = preview.get("status") if isinstance(preview, dict) else None
            print(json.dumps({"preview_status": state, "source_digest": tree}), flush=True)
            if state == "succeeded":
                # The first actual response established this terminal schema.
                # The resolver verifies the exact job, digest, completion,
                # PDF/log availability and gate before clearing the marker.
                resolve_preview(args, entry, metadata)
                return 0
            if not PREVIEW_MARKER.exists():
                print(json.dumps({"preview_terminal_confirmed_by_coordinator": True}), flush=True)
                return 0
        raise Stop("Preview result remains uncertain; coordinator must check before releasing lane")


def resolve_preview(args, entry: dict, metadata: dict) -> int:
    """Record the observed successful-terminal schema without reading a key."""
    marker_bytes = PREVIEW_MARKER.read_bytes()
    marker = json.loads(marker_bytes)
    draft = component(args.draft_id, "draft identifier")
    tree = validated_source_digest(draft, args.paper, metadata, entry)
    if (marker.get("paper") != args.paper or marker.get("draft_id") != draft or
            marker.get("source_digest") != tree or
            marker.get("metadata_sha256") != digest(canonical(metadata)) or
            marker.get("files_sha256") != digest(canonical(entry["upload_files"]))):
        raise Stop("Unresolved marker does not match this paper and source set")
    records = load_json(RECEIPT_PATH)
    matching = [row for row in records if row.get("paper") == args.paper
                and row.get("api_path") == "/api/v1/submissions/" + draft + "/preview"
                and row.get("metadata_sha256") == digest(canonical(metadata))
                and row.get("files_sha256") == digest(canonical(entry["upload_files"]))]
    submitted = [row for row in matching if row.get("action") == "preview"
                 and row.get("http_status") in {200, 202}]
    observed = [row for row in matching if row.get("action") == "preview-status"
                and row.get("http_status") == 200]
    if not submitted or not observed:
        raise Stop("Missing actual preview request or terminal observation")
    value = observed[-1]["response"]
    preview, gate = value.get("preview") or {}, value.get("gate") or {}
    observed_id = component(preview.get("id"), "observed preview identifier")
    submitted_id = component(submitted[-1]["response"].get("preview", {}).get("id"),
                             "submitted preview identifier")
    if (value.get("draftId") != draft or value.get("digest") != tree or
            observed_id != submitted_id or
            preview.get("digest") != tree or preview.get("current") is not True or
            preview.get("status") != "succeeded" or not preview.get("completedAt") or
            preview.get("errorCode") is not None or preview.get("hasLog") is not True or
            preview.get("hasPdf") is not True or not isinstance(preview.get("pdfSizeBytes"), int) or
            preview["pdfSizeBytes"] <= 0 or gate.get("state") != "previewed" or
            gate.get("allowed") is not True or gate.get("warn") is not False):
        raise Stop("Actual preview evidence does not establish successful termination")
    proof = {"paper": args.paper, "draft_id": draft, "preview_id": preview["id"],
             "source_digest": tree, "preview_status": "succeeded",
             "pdf_bytes": preview["pdfSizeBytes"], "completed_at": preview["completedAt"],
             "status_receipt_sha256": digest(canonical(observed[-1])),
             "status_receipt": observed[-1], "bibliography_reviewed": False}
    proof_path = HERE / (args.paper + ".preview-resolution.local")
    fd = os.open(proof_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8") as stream:
        json.dump(proof, stream, indent=2)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    if PREVIEW_MARKER.read_bytes() != marker_bytes:
        raise Stop("Preview marker changed during terminal confirmation")
    PREVIEW_MARKER.unlink()
    print(json.dumps({"paper": args.paper, "preview_terminal": "succeeded",
                      "resolution_path": str(proof_path)}), flush=True)
    return 0


def preview_artifacts(args, entry: dict, metadata: dict) -> int:
    draft = component(args.draft_id, "draft identifier")
    bind_draft(draft, args.paper, metadata)
    proof = load_json(HERE / (args.paper + ".preview-resolution.local"))
    tree = validated_source_digest(draft, args.paper, metadata, entry)
    if proof.get("draft_id") != draft or proof.get("paper") != args.paper or proof.get("source_digest") != tree:
        raise Stop("Terminal preview proof differs from this paper/source set")
    preview_id = component(proof.get("preview_id"), "resolved preview identifier")
    base = "/api/v1/submissions/" + draft + "/preview"
    token = credential()
    def check_current() -> None:
        code, value = exchange(args, entry, metadata, token, "preview-status", "GET", base)
        if code != 200 or not isinstance(value, dict) or not isinstance(value.get("preview"), dict):
            raise Stop("Preview status did not return the expected object")
        item = value["preview"]
        if (value.get("draftId") != draft or value.get("digest") != tree
                or item.get("id") != preview_id or item.get("digest") != tree
                or item.get("status") != "succeeded" or item.get("current") is not True
                or item.get("hasLog") is not True or item.get("hasPdf") is not True
                or item.get("pdfSizeBytes") != proof.get("pdf_bytes")):
            raise Stop("Successful preview changed before or during artifact retrieval")
    check_current()
    artifacts = {}
    for kind in ("log", "pdf"):
        target = HERE / (args.paper + "." + preview_id + ".preview." + kind + ".local")
        code, descriptor = exchange(args, entry, metadata, token, "download-preview-" + kind,
                                     "GET", base + "/" + kind, artifact_path=target)
        if code != 200:
            raise Stop("Preview artifact was not available; inspect the actual response")
        if kind == "pdf" and descriptor.get("bytes") != proof["pdf_bytes"]:
            raise Stop("Preview PDF byte count differs from the terminal proof")
        artifacts[kind] = descriptor
    check_current()
    result = {"paper": args.paper, "draft_id": draft, "preview_id": preview_id,
              "source_digest": tree, "artifacts": artifacts, "bibliography_reviewed": False}
    target = HERE / (args.paper + ".preview-artifacts.local")
    fd = os.open(target, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    print(json.dumps({"paper": args.paper, "preview_artifacts": str(target)}), flush=True)
    return 0


def inspect_preview(args, entry: dict, metadata: dict) -> int:
    """Read the two PDFs and existing log; no key, network or TeX execution."""
    from pypdf import PdfReader
    artifact_record = load_json(HERE / (args.paper + ".preview-artifacts.local"))
    proof = load_json(HERE / (args.paper + ".preview-resolution.local"))
    draft = component(args.draft_id, "draft identifier")
    bind_draft(draft, args.paper, metadata)
    preview_id = component(proof.get("preview_id"), "preview identifier")
    tree = validated_source_digest(draft, args.paper, metadata, entry)
    recorded = proof.get("status_receipt") or {}
    if (artifact_record.get("paper") != args.paper or artifact_record.get("draft_id") != draft
            or artifact_record.get("preview_id") != preview_id
            or artifact_record.get("source_digest") != tree or proof.get("source_digest") != tree
            or recorded.get("metadata_sha256") != digest(canonical(metadata))
            or recorded.get("files_sha256") != digest(canonical(entry["upload_files"]))):
        raise Stop("Preview artifact record does not match the terminal proof")
    paths = {}
    for kind in ("pdf", "log"):
        row = artifact_record["artifacts"][kind]
        path = HERE / (args.paper + "." + preview_id + ".preview." + kind + ".local")
        if row.get("artifact_path") != str(path) or path.is_symlink():
            raise Stop("Unexpected preview artifact path")
        raw = path.read_bytes()
        if len(raw) != row["bytes"] or digest(raw) != row["sha256"]:
            raise Stop("Preview artifact differs from its download receipt")
        paths[kind] = path
    original = PdfReader(ROOT / entry["reference_pdf"]["path"])
    remote = PdfReader(paths["pdf"])
    texts = ["\n".join(page.extract_text() or "" for page in reader.pages)
             for reader in (original, remote)]
    def normalize(text: str) -> str:
        return re.sub(r"\s+", "", unicodedata.normalize("NFKC", text))
    bibliographies = [text.rsplit("References", 1)[1] if "References" in text else None
                      for text in texts]
    date_equivalent = False
    title_dates = None
    source_descriptor = next(row for row in entry["upload_files"] if row["relativePath"] == "main.tex")
    source_text = (ROOT / source_descriptor["path"]).read_text(encoding="utf-8")
    if re.search(r"(?m)^[ \t]*\\date[ \t]*\{[ \t]*\\today[ \t]*\}[ \t]*(?:%[^\n]*)?$", source_text):
        lines = [text.splitlines() for text in texts]
        candidates = [[(i, line) for i, line in enumerate(rows[:12])
                       if re.fullmatch(r"[A-Z][a-z]+ \d{1,2}, \d{4}", line.strip())]
                      for rows in lines]
        if all(len(items) == 1 for items in candidates):
            positions = [items[0][0] for items in candidates]
            dates = [dt.datetime.strptime(items[0][1].strip(), "%B %d, %Y").date()
                     for items in candidates]
            completed_day = dt.datetime.strptime(proof["completed_at"], "%Y-%m-%d %H:%M:%S").date()
            if positions[0] == positions[1] and dates[1] == completed_day:
                for rows, position in zip(lines, positions):
                    rows[position] = "<TITLE_DATE_FROM_TODAY>"
                date_equivalent = normalize("\n".join(lines[0])) == normalize("\n".join(lines[1]))
                title_dates = [day.isoformat() for day in dates]
    log = paths["log"].read_text(encoding="utf-8", errors="replace")
    final_pass = log.rsplit("This is pdfTeX", 1)[-1]
    result = {"paper": args.paper, "draft_id": draft, "preview_id": preview_id,
              "source_digest": proof["source_digest"],
              "metadata_sha256": digest(canonical(metadata)),
              "files_sha256": digest(canonical(entry["upload_files"])),
              "pdf_sha256": artifact_record["artifacts"]["pdf"]["sha256"],
              "reference_pdf_sha256": entry["reference_pdf"]["sha256"],
              "reference_pages": len(original.pages), "preview_pages": len(remote.pages),
              "page_count_matches": len(original.pages) == len(remote.pages) == entry["pages"],
              "normalized_full_text_equal": normalize(texts[0]) == normalize(texts[1]),
              "normalized_text_equal_except_explicit_today": date_equivalent,
              "reference_and_preview_title_dates": title_dates,
              "normalized_bibliography_equal": all(text is not None for text in bibliographies)
                  and normalize(bibliographies[0]) == normalize(bibliographies[1]),
              "bibtex_executed": "This is BibTeX" in log,
              "final_pass_has_undefined_references": "undefined" in final_pass.lower(),
              "visual_layout_review_performed": False,
              "comparison_scope": "Extracted PDF text, page count and bibliography"}
    result["bibliography_reviewed"] = (result["normalized_bibliography_equal"]
        and result["bibtex_executed"] and not result["final_pass_has_undefined_references"])
    result["content_review_passed"] = (result["page_count_matches"]
        and (result["normalized_full_text_equal"] or date_equivalent) and result["bibliography_reviewed"])
    if not result["normalized_full_text_equal"]:
        result["text_diff"] = list(difflib.unified_diff(texts[0].splitlines(), texts[1].splitlines(),
                                                       fromfile="reference", tofile="preview"))[:250]
    target = HERE / (args.paper + ".preview-content-review.local")
    fd = os.open(target, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    print(json.dumps({key: value for key, value in result.items()
                      if key not in {"text_diff", "draft_id", "preview_id"}}), flush=True)
    return 0


def commit_publication(args, entry: dict, metadata: dict) -> int:
    """One authorized commit, bound to the inspected current preview; no retry."""
    draft = component(args.draft_id, "draft identifier")
    bind_draft(draft, args.paper, metadata)
    if getattr(args, "operation_lock", None) is None:
        raise Stop("Commit requires the shared private operation lock")
    if PREVIEW_MARKER.exists() or COMMIT_MARKER.exists():
        raise Stop("An earlier preview or commit is unresolved; inspect its actual state")
    records = load_json(RECEIPT_PATH)
    if any(row.get("action") == "commit" and row.get("paper") == args.paper for row in records):
        raise Stop("This paper already has a commit attempt; read its state instead of retrying")
    today = dt.datetime.now(dt.timezone.utc).date().isoformat()
    used = sum(row.get("action") == "commit" and row.get("http_status") == 202
               and row.get("utc_time", "").startswith(today) for row in records)
    if used >= DAILY_COMMIT_LIMIT:
        raise Stop("The user-authorized daily deposit allowance is exhausted")
    review_path = HERE / (args.paper + ".preview-content-review.local")
    review_bytes = review_path.read_bytes()
    review = json.loads(review_bytes)
    artifact_record = load_json(HERE / (args.paper + ".preview-artifacts.local"))
    proof = load_json(HERE / (args.paper + ".preview-resolution.local"))
    preview_id = component(review.get("preview_id"), "reviewed preview identifier")
    meta_sha, files_sha = digest(canonical(metadata)), digest(canonical(entry["upload_files"]))
    if (review.get("paper") != args.paper or review.get("draft_id") != draft
            or review.get("metadata_sha256") != meta_sha or review.get("files_sha256") != files_sha
            or review.get("reference_pdf_sha256") != entry["reference_pdf"]["sha256"]
            or review.get("content_review_passed") is not True
            or review.get("page_count_matches") is not True
            or review.get("bibliography_reviewed") is not True
            or review.get("normalized_bibliography_equal") is not True
            or review.get("bibtex_executed") is not True
            or review.get("final_pass_has_undefined_references") is not False
            or not (review.get("normalized_full_text_equal") is True or
                    review.get("normalized_text_equal_except_explicit_today") is True)):
        raise Stop("Current paper lacks the required successful PDF content comparison")
    for record in (artifact_record, proof):
        if (record.get("paper") != args.paper or record.get("draft_id") != draft
                or record.get("preview_id") != preview_id
                or record.get("source_digest") != review.get("source_digest")):
            raise Stop("PDF review, artifact and terminal proof are not the same preview")
    for kind in ("pdf", "log"):
        descriptor = artifact_record["artifacts"][kind]
        expected_path = HERE / (args.paper + "." + preview_id + ".preview." + kind + ".local")
        if descriptor.get("artifact_path") != str(expected_path) or expected_path.is_symlink():
            raise Stop("Reviewed preview artifact path is unexpected")
        data = expected_path.read_bytes()
        if len(data) != descriptor["bytes"] or digest(data) != descriptor["sha256"]:
            raise Stop("Reviewed preview artifact has changed")
    if review.get("pdf_sha256") != artifact_record["artifacts"]["pdf"]["sha256"]:
        raise Stop("Content review is not bound to the downloaded PDF")
    base = "/api/v1/submissions/" + draft
    token = credential()
    current = None
    for action, suffix in (("draft-status", ""), ("source-check", "/preview/check"),
                           ("preview-status", "/preview")):
        code, current = exchange(args, entry, metadata, token, action, "GET", base + suffix)
        if code != 200:
            raise Stop("Pre-commit live validation failed; no commit requested")
    tree = validated_source_digest(draft, args.paper, metadata, entry)
    if not isinstance(current, dict) or not isinstance(current.get("preview"), dict):
        raise Stop("Current preview status is missing")
    preview, gate = current["preview"], current.get("gate") or {}
    if (tree != review.get("source_digest") or preview.get("id") != preview_id
            or preview.get("digest") != tree or preview.get("current") is not True
            or preview.get("status") != "succeeded" or preview.get("errorCode") is not None
            or preview.get("completedAt") != proof.get("completed_at")
            or preview.get("hasPdf") is not True or preview.get("hasLog") is not True
            or preview.get("pdfSizeBytes") != artifact_record["artifacts"]["pdf"]["bytes"]
            or gate.get("state") != "previewed" or gate.get("allowed") is not True
            or gate.get("warn") is not False):
        raise Stop("Live preview differs from the inspected successful preview")
    evidence = {"paper": args.paper, "draft_id": draft, "preview_id": preview_id,
                "source_digest": tree, "metadata_sha256": meta_sha, "files_sha256": files_sha,
                "content_review_sha256": digest(review_bytes),
                "requested_utc": dt.datetime.now(dt.timezone.utc).isoformat()}
    write_private_json_once(COMMIT_MARKER, dict(evidence, state="COMMIT_RESULT_UNRESOLVED"))
    code, value = exchange(args, entry, metadata, token, "commit", "POST", base + "/commit")
    if code != 202:
        raise Stop("Commit result requires reconciliation; do not retry automatically")
    receipt = dict(evidence, http_status=202, state="ACCEPTED_FOR_ASYNCHRONOUS_PROCESSING",
                   response=value, accepted_utc=dt.datetime.now(dt.timezone.utc).isoformat())
    target = HERE / (args.paper + ".submission-receipt.local")
    write_private_json_once(target, receipt)
    COMMIT_MARKER.unlink()
    print(json.dumps({"paper": args.paper, "submission_accepted": True,
                      "publication_status": "processing", "receipt_path": str(target)}), flush=True)
    return 0


def run_queue(args) -> int:
    """Submit the remaining authorized papers within the updated daily limit."""
    if os.environ.get("COXETER_PUBLICATION_IO_GUARD") != "1":
        raise Stop("Scheduled execution requires the publication watchdog")
    config = load_json(QUEUE_CONFIGURATION)
    if (config.get("client_sha256") != digest(Path(__file__).read_bytes()) or
            config.get("watchdog_sha256") != digest((ROOT / "tools/run_publication_io.py").read_bytes()) or
            config.get("label") != LAUNCH_LABEL or config.get("papers") != list(PAPERS)):
        raise Stop("Scheduled source/monitor configuration differs from the reviewed version")
    if (config.get("max_accepted_submissions_per_utc_day") != DAILY_COMMIT_LIMIT or
            config.get("daily_limit_authorization") != DAILY_LIMIT_AUTHORIZATION):
        raise Stop("Updated daily allowance has not been bound to the user's authorization")
    snapshot_bytes = (HERE / "approved_payloads.json").read_bytes()
    if (digest(snapshot_bytes) != config.get("approved_payloads_sha256") or
            payload_snapshot(load_json(HERE / "queue.json")) != json.loads(snapshot_bytes) or
            digest((ROOT / "release/manifest.json").read_bytes()) != config.get("release_manifest_sha256")):
        raise Stop("Scheduled payload descriptors or release manifest have changed")
    if (HERE / "queue-hold.local").exists():
        print(json.dumps({"queue_held": "A previous scheduled step needs receipt reconciliation"}), flush=True)
        return 0
    for paper in PAPERS:
        path = HERE / (paper + ".metadata.template.json")
        if (config.get("metadata_sha256") or {}).get(paper) != digest(path.read_bytes()):
            raise Stop("Scheduled metadata differs from the user's reviewed declarations")
    # The first completed submission anchors both the authorized order and UTC days.
    accepted = {}
    entries = {row["id"]: row for row in load_json(HERE / "queue.json")["entries"]}
    history = load_json(RECEIPT_PATH)
    for paper in PAPERS:
        receipt_path = HERE / (paper + ".submission-receipt.local")
        if receipt_path.exists():
            receipt_bytes = receipt_path.read_bytes()
            receipt = json.loads(receipt_bytes)
            metadata = load_json(HERE / (paper + ".metadata.template.json"))
            meta_sha = digest(canonical(metadata))
            files_sha = digest(canonical(entries[paper]["upload_files"]))
            reply = receipt.get("response") or {}
            server = reply.get("receipt") or {}
            if (receipt.get("paper") != paper or receipt.get("http_status") != 202
                    or receipt.get("state") != "ACCEPTED_FOR_ASYNCHRONOUS_PROCESSING"
                    or receipt.get("metadata_sha256") != meta_sha or receipt.get("files_sha256") != files_sha
                    or reply.get("draftId") != receipt.get("draft_id")
                    or server.get("draftId") != receipt.get("draft_id")
                    or server.get("title") != metadata.get("title")
                    or not isinstance(server.get("identifier"), str)
                    or not re.fullmatch(r"\d{4}\.\d{5}", server["identifier"])
                    or server.get("versionId") != server["identifier"] + "v1"
                    or reply.get("versionId") != server.get("versionId")
                    or reply.get("identifier") != "hexagon:" + server["identifier"]):
                raise Stop("A prior submission receipt is inconsistent")
            if not any(row.get("action") == "commit" and row.get("http_status") == 202
                       and row.get("paper") == paper and row.get("method") == "POST"
                       and row.get("api_path") == "/api/v1/submissions/" + receipt["draft_id"] + "/commit"
                       and row.get("metadata_sha256") == meta_sha and row.get("files_sha256") == files_sha
                       and row.get("response") == reply for row in history):
                raise Stop("Accepted receipt has no matching successful commit execution")
            if paper == PAPERS[0]:
                first_descriptor = config.get("first_accepted_receipt") or {}
                if (first_descriptor.get("filename") != receipt_path.name
                        or first_descriptor.get("bytes") != len(receipt_bytes)
                        or first_descriptor.get("sha256") != digest(receipt_bytes)
                        or server.get("versionId") != config.get("first_accepted_submission")
                        or dt.datetime.fromisoformat(receipt["accepted_utc"]).date().isoformat()
                            != config.get("first_accepted_utc_day")):
                    raise Stop("The first accepted submission differs from the schedule anchor")
            accepted[paper] = receipt
    if PAPERS[0] not in accepted:
        raise Stop("The first accepted submission must anchor this scheduled queue")
    if list(accepted) != list(PAPERS[:len(accepted)]):
        raise Stop("Accepted submissions do not follow the reviewed queue order")
    accepted_days = [dt.datetime.fromisoformat(row["accepted_utc"]).date() for row in accepted.values()]
    if any(right < left for left, right in zip(accepted_days, accepted_days[1:])):
        raise Stop("Accepted receipts violate the chronological queue order")
    if any(accepted_days.count(day) > DAILY_COMMIT_LIMIT for day in set(accepted_days)):
        raise Stop("Accepted receipts exceed the user-authorized daily allowance")
    if len(accepted) == len(PAPERS):
        result = subprocess.run(["/bin/launchctl", "disable", f"gui/{os.getuid()}/{LAUNCH_LABEL}"],
                                capture_output=True, timeout=10)
        if result.returncode != 0:
            raise Stop("Submissions are accepted, but disabling the completed schedule failed")
        print(json.dumps({"queue_complete": True, "schedule_disabled": True}), flush=True)
        return 0
    first_day = dt.datetime.fromisoformat(accepted[PAPERS[0]]["accepted_utc"]).date()
    today = dt.datetime.now(dt.timezone.utc).date()
    paper = next(paper for paper in PAPERS if paper not in accepted)
    not_before = first_day
    used_today = sum(dt.datetime.fromisoformat(row["accepted_utc"]).date() == today
                     for row in accepted.values())
    if today < not_before or used_today >= DAILY_COMMIT_LIMIT:
        print(json.dumps({"queue_deferred": "UTC eligibility/daily allowance", "paper": paper,
                          "not_before_utc_date": not_before.isoformat()}), flush=True)
        return 0
    with (ROOT / "formal/.proof-batches/exclusive.lock").open("a") as science_lock:
        try:
            fcntl.flock(science_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print(json.dumps({"queue_deferred": "scientific/typesetting lane occupied", "paper": paper}), flush=True)
            return 0
        current = argparse.Namespace(paper=paper, action="prepare", file=None, draft_id=None,
                                     upload_id=None, operation_lock=args.operation_lock,
                                     scientific_lock=science_lock)
        entry, metadata, _, missing = preparation(paper)
        if missing:
            raise Stop("Scheduled paper has incomplete declarations")
        def perform(action: str, filename: str | None = None, upload: str | None = None) -> None:
            current.action, current.file, current.upload_id = action, filename, upload
            rows = load_json(RECEIPT_PATH)
            answered = {row.get("request_id") for row in rows
                        if row.get("state") == "HTTP_RESPONSE_RECEIVED"}
            if action == "create-draft":
                unknown = [row for row in rows if row.get("paper") == paper and row.get("action") == action
                           and row.get("state") == "REQUEST_PENDING_OR_RESULT_UNKNOWN"
                           and row.get("request_id") not in answered]
                if any(dt.datetime.now(dt.timezone.utc) - dt.datetime.fromisoformat(row["utc_time"])
                       >= dt.timedelta(hours=24) for row in unknown):
                    raise Stop("An unknown draft creation is outside its documented idempotency window")
            if action in {"create-upload", "put-part", "complete-upload", "preview", "commit"} and any(row.get("paper") == paper and
                    row.get("action") == action and row.get("manifest_file") == filename and
                    row.get("state") == "REQUEST_PENDING_OR_RESULT_UNKNOWN" and
                    row.get("request_id") not in answered for row in rows):
                raise Stop("An earlier request has no known result; read its state before retrying")
            try:
                if execute_args(current) != 0:
                    raise Stop("Scheduled API step failed; inspect its receipt before resuming")
            except Stop as error:
                if action in {"create-draft", "draft-status", "source-check", "preview-status"}:
                    # GETs can be repeated; draft creation reuses the same documented idempotency key.
                    raise
                hold = {"paper": paper, "action": action, "manifest_file": filename,
                        "reason": str(error), "utc_time": dt.datetime.now(dt.timezone.utc).isoformat()}
                write_private_json_once(HERE / "queue-hold.local", hold)
                raise
        def successful(action: str, filename: str | None = None, upload: str | None = None):
            rows = load_json(RECEIPT_PATH)
            base = "/api/v1/submissions/" + (current.draft_id or "")
            if action == "create-draft":
                expected_path, expected_body = "/api/v1/submissions", digest(canonical(metadata))
            elif action == "create-upload":
                expected_path, expected_body = base + "/uploads", digest(canonical(upload_descriptor(entry, filename)))
            elif action == "put-part":
                expected_path = base + "/uploads/" + component(upload, "resume upload identifier") + "/parts/1"
                expected_body = next(row["sha256"] for row in entry["upload_files"] if row["relativePath"] == filename)
            elif action == "complete-upload":
                expected_path = base + "/uploads/" + component(upload, "resume upload identifier") + "/complete"
                expected_body = None
            else:
                raise Stop("Unexpected resumable queue action")
            return [row for row in rows if row.get("action") == action
                    and row.get("paper") == paper and row.get("http_status") in {200, 201, 202}
                    and row.get("metadata_sha256") == digest(canonical(metadata))
                    and row.get("files_sha256") == digest(canonical(entry["upload_files"]))
                    and row.get("api_path") == expected_path and row.get("payload_sha256") == expected_body
                    and (filename is None or row.get("manifest_file") == filename)]
        created = successful("create-draft")
        if not created:
            perform("create-draft")
            created = successful("create-draft")
        current.draft_id = component(created[-1]["response"].get("draftId"), "queued draft identifier")
        expires = created[-1]["response"].get("expiresAt")
        if not isinstance(expires, str) or dt.datetime.now(dt.timezone.utc) >= explicit_utc_timestamp(expires):
            raise Stop("Stored private draft has expired; reconcile it before creating another")
        perform("draft-status")
        for filename in sorted(HEXAGON_SOURCE_FILES):
            registered = successful("create-upload", filename)
            if not registered:
                perform("create-upload", filename)
                registered = successful("create-upload", filename)
            upload = component(registered[-1]["response"].get("uploadId"), "queued upload identifier")
            if not successful("put-part", filename, upload):
                perform("put-part", filename, upload)
            if not successful("complete-upload", filename, upload):
                perform("complete-upload", filename, upload)
        perform("draft-status")
        perform("source-check")
        perform("preview-status")
        if not (HERE / (paper + ".preview-resolution.local")).exists():
            perform("preview")
        if not (HERE / (paper + ".preview-artifacts.local")).exists():
            perform("preview-artifacts")
        if not (HERE / (paper + ".preview-content-review.local")).exists():
            perform("inspect-preview")
        review = load_json(HERE / (paper + ".preview-content-review.local"))
        if review.get("content_review_passed") is not True:
            raise Stop("Scheduled PDF content comparison needs review; no commit requested")
        perform("commit")
        print(json.dumps({"scheduled_submission_accepted": paper}), flush=True)
    return run_queue(args)  # Reconcile receipts before the next paper or schedule disable.


def inspect_queue() -> int:
    """Print only submission identifiers/statuses and nonsecret queue state."""
    if os.environ.get("COXETER_PUBLICATION_IO_GUARD") != "1":
        raise Stop("Queue inspection requires the publication watchdog")
    records = load_json(RECEIPT_PATH) if RECEIPT_PATH.exists() else []
    entries = []
    for paper in PAPERS:
        accepted_path = HERE / (paper + ".submission-receipt.local")
        accepted = load_json(accepted_path) if accepted_path.exists() else {}
        reply = accepted.get("response") or {}
        server = reply.get("receipt") or {}
        completed = [row for row in records if row.get("paper") == paper
                     and row.get("state") == "HTTP_RESPONSE_RECEIVED"]
        statuses = [row for row in completed if row.get("action") == "draft-status"
                    and row.get("http_status") == 200]
        last_status = statuses[-1].get("response", {}) if statuses else {}
        entries.append({"paper": paper, "accepted_receipt_present": accepted_path.exists(),
            "accepted_utc": accepted.get("accepted_utc"),
            "version_id": server.get("versionId"),
            "last_known_status": last_status.get("status"),
            "draft_created": any(row.get("action") == "create-draft"
                and row.get("http_status") in {200, 201} for row in completed),
            "commit_attempts": sum(row.get("action") == "commit" for row in completed),
            "preview_review_present": (HERE / (paper + ".preview-content-review.local")).exists()})
    answered = {row.get("request_id") for row in records if row.get("state") == "HTTP_RESPONSE_RECEIVED"}
    unknown = [{"paper": row.get("paper"), "action": row.get("action")}
               for row in records if row.get("state") == "REQUEST_PENDING_OR_RESULT_UNKNOWN"
               and row.get("request_id") not in answered]
    print(json.dumps({"queue_entries": entries, "unanswered_requests": unknown,
        "holds": {name: (HERE/name).exists() for name in
                  ("unresolved-preview.local", "unresolved-commit.local", "queue-hold.local")}}), flush=True)
    return 0


def apply_five_per_day() -> int:
    """Bind the existing approved queue to the new explicit user instruction."""
    if os.environ.get("COXETER_PUBLICATION_IO_GUARD") != "1":
        raise Stop("Policy update requires the publication watchdog")
    config_bytes = QUEUE_CONFIGURATION.read_bytes()
    config = json.loads(config_bytes)
    if (config.get("label") != LAUNCH_LABEL or config.get("papers") != list(PAPERS)
            or config.get("watchdog_sha256") != digest((ROOT/"tools/run_publication_io.py").read_bytes())
            or config.get("approved_payloads_sha256") != digest((HERE/"approved_payloads.json").read_bytes())
            or payload_snapshot(load_json(HERE/"queue.json")) != load_json(HERE/"approved_payloads.json")):
        raise Stop("Existing approved payload or monitor configuration changed")
    backup = HERE/"scheduled-execution.before_five_per_day.local"
    if not backup.exists():
        fd = os.open(backup, os.O_WRONLY|os.O_CREAT|os.O_EXCL|os.O_NOFOLLOW, 0o600)
        with os.fdopen(fd, "wb") as stream:
            stream.write(config_bytes)
            stream.flush()
            os.fsync(stream.fileno())
        if backup.read_bytes() != config_bytes:
            raise Stop("Prior schedule configuration backup failed byte verification")
    for paper in PAPERS:
        if preparation(paper)[3]:
            raise Stop("Missing declarations in an approved paper")
    config.update(client_sha256=digest(Path(__file__).read_bytes()),
        max_accepted_submissions_per_utc_day=DAILY_COMMIT_LIMIT,
        daily_limit_authorization=DAILY_LIMIT_AUTHORIZATION,
        daily_limit_updated_utc=dt.datetime.now(dt.timezone.utc).isoformat())
    fd = os.open(QUEUE_CONFIGURATION, os.O_WRONLY | os.O_TRUNC | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8") as stream:
        json.dump(config, stream, indent=2)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    print(json.dumps({"daily_limit": DAILY_COMMIT_LIMIT, "payloads_unchanged": True,
        "authorization": DAILY_LIMIT_AUTHORIZATION, "client_sha256": config["client_sha256"]}), flush=True)
    return 0


def refresh_statuses() -> int:
    """Refresh actual accepted submissions without displaying private IDs."""
    if os.environ.get("COXETER_PUBLICATION_IO_GUARD") != "1":
        raise Stop("Status refresh requires the publication watchdog")
    for paper in PAPERS:
        path = HERE/(paper+".submission-receipt.local")
        if not path.exists():
            continue
        receipt = load_json(path)
        args = argparse.Namespace(paper=paper, action="draft-status", file=None,
            draft_id=receipt["draft_id"], upload_id=None)
        if execute_args(args) != 0:
            raise Stop("Accepted submission status refresh failed")
    return inspect_queue()


def inspect_service() -> int:
    """Fetch the current official public schema, without an account credential."""
    if os.environ.get("COXETER_PUBLICATION_IO_GUARD") != "1":
        raise Stop("Service inspection requires the publication watchdog")
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}), NoRedirect())
    request_ = urllib.request.Request(ORIGIN+"/api/v1/openapi.json", headers={
        "Accept": "application/json", "User-Agent": "CoxeterPolytopes-publication-client/1"})
    with opener.open(request_, timeout=30) as response:
        raw = response.read(MAX_RESPONSE+1)
        if response.code != 200 or len(raw) > MAX_RESPONSE:
            raise Stop("Current public API schema could not be read")
    schema = json.loads(raw)
    paths = schema.get("paths") or {}
    quota_paths = [path for path in paths if re.search(r"account|quota|limit|profile", path, re.I)]
    target = HERE/"current-api-schema-five-per-day.local"
    if not target.exists():
        fd = os.open(target, os.O_WRONLY|os.O_CREAT|os.O_EXCL|os.O_NOFOLLOW, 0o600)
        with os.fdopen(fd, "wb") as stream:
            stream.write(raw)
    print(json.dumps({"api_version": schema.get("info", {}).get("version"),
        "schema_bytes": len(raw), "schema_sha256": digest(raw),
        "account_quota_paths": quota_paths,
        "daily_allowance_source": DAILY_LIMIT_AUTHORIZATION,
        "server_enforcement_retained": True}), flush=True)
    return 0


def completion_report() -> int:
    """Write a public summary of actual receipts, excluding private draft data."""
    refresh_statuses()
    history = load_json(RECEIPT_PATH)
    papers = []
    for paper in PAPERS:
        entry, metadata, _, missing = preparation(paper)
        receipt_path = HERE/(paper+".submission-receipt.local")
        if missing or not receipt_path.exists():
            raise Stop("The queue is not complete")
        receipt = load_json(receipt_path)
        server = (receipt.get("response") or {}).get("receipt") or {}
        if receipt.get("http_status") != 202:
            raise Stop("A submission lacks its actual accepted receipt")
        states = [row for row in history if row.get("paper") == paper
                  and row.get("action") == "draft-status" and row.get("http_status") == 200]
        latest = states[-1]
        status = (latest.get("response") or {}).get("status")
        review = load_json(HERE/(paper+".preview-content-review.local"))
        papers.append({"paper": paper, "version_id": server.get("versionId"),
            "identifier": server.get("identifier"), "commit_http_status": 202,
            "accepted_utc": receipt.get("accepted_utc"), "current_status": status,
            "status_checked_utc": latest.get("utc_time"),
            "publication_confirmed": status == "published",
            "public_record_url": ORIGIN+"/abs/"+server["identifier"],
            "accepted_receipt_sha256": digest(receipt_path.read_bytes()),
            "source_hashes_match_approved_release": True,
            "preview": {key: review.get(key) for key in (
                "reference_pages", "preview_pages", "page_count_matches",
                "normalized_full_text_equal", "normalized_text_equal_except_explicit_today",
                "normalized_bibliography_equal", "bibtex_executed",
                "final_pass_has_undefined_references", "content_review_passed")}})
    disabled = subprocess.run(["/bin/launchctl", "print-disabled", f"gui/{os.getuid()}"],
                              capture_output=True, text=True, timeout=10)
    is_disabled = (disabled.returncode == 0 and
        re.search(r'"'+re.escape(LAUNCH_LABEL)+r'"\s*=>\s*(?:disabled|true)\b', disabled.stdout) is not None)
    if not is_disabled:
        raise Stop("Completed schedule is not confirmed disabled")
    holds = {name: (HERE/name).exists() for name in
             ("unresolved-preview.local", "unresolved-commit.local", "queue-hold.local")}
    if any(holds.values()):
        raise Stop("A remote operation remains unresolved")
    report = {"schema": "hexagon-authorized-five-per-day-completion-v1",
        "status": "ALL_THREE_APPROVED_PAPERS_ACCEPTED_SCHEDULE_DISABLED",
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "daily_limit_authorization": DAILY_LIMIT_AUTHORIZATION,
        "local_max_accepted_per_utc_day": DAILY_COMMIT_LIMIT,
        "account_quota_endpoint_available": False,
        "account_allowance_source": "User's current account report; server checks remain active",
        "official_api_schema_sha256": digest((HERE/"current-api-schema-five-per-day.local").read_bytes()),
        "client_sha256": digest(Path(__file__).read_bytes()),
        "submitted_client_sha256": digest((HERE/"client.submitted_five_per_day.local").read_bytes()),
        "prior_client_sha256": digest((HERE/"client.before_five_per_day.local").read_bytes()),
        "prior_configuration_bytes_preserved": True,
        "approved_payloads_unchanged_sha256": digest((HERE/"approved_payloads.json").read_bytes()),
        "papers": papers, "schedule_disabled": is_disabled,
        "unresolved_markers": holds,
        "original_supporting_bbl_rejection_preserved": any(row.get("http_status") == 400
            and row.get("manifest_file") == "main.bbl" for row in history),
        "held_geometric_paper_submitted": False,
        "held_geometric_reason": "Source/bibliography differ from recorded release/build package",
        "credential_contents_recorded": False, "bulk_assets_uploaded": False,
        "scientific_or_local_TeX_jobs_launched": False,
        "remote_TeX_preview_executed": True}
    target = HERE/"completion_five_per_day_20261123.json"
    if target.exists():
        raise Stop("Public completion report already exists")
    target.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"completion_report": str(target), "papers": [
        {key: row[key] for key in ("paper", "version_id", "current_status", "publication_confirmed")}
        for row in papers], "schedule_disabled": True}), flush=True)
    return 0


def inspect_update_api() -> int:
    """Read current official revision/correction schemas without credentials."""
    if os.environ.get("COXETER_PUBLICATION_IO_GUARD") != "1":
        raise Stop("API inspection requires the publication watchdog")
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}), NoRedirect())
    request_ = urllib.request.Request(ORIGIN+"/api/v1/openapi.json", headers={
        "Accept": "application/json", "User-Agent": "CoxeterPolytopes-publication-client/1"})
    with opener.open(request_, timeout=30) as response:
        raw = response.read(MAX_RESPONSE+1)
        if response.code != 200 or len(raw) > MAX_RESPONSE:
            raise Stop("Current public API schema could not be read")
    schema = json.loads(raw)
    paths = schema.get("paths") or {}
    selected = {path: value for path, value in paths.items()
        if path in {"/api/v1/submissions", "/api/v1/submissions/{draftId}"}
        or re.search(r"/submissions/.*/(?:revise|corrections|lifecycle)$", path)
        or path == "/api/v1/submissions/{draftId}/uploads/{uploadId}"}
    target = HERE/"lp-remark-current-api-schema.local"
    if target.exists():
        if target.read_bytes() != raw:
            raise Stop("Previously inspected revision schema has changed")
    else:
        fd = os.open(target, os.O_WRONLY|os.O_CREAT|os.O_EXCL|os.O_NOFOLLOW, 0o600)
        with os.fdopen(fd, "wb") as stream:
            stream.write(raw)
    print(json.dumps({"api_version": schema.get("info", {}).get("version"),
        "schema_sha256": digest(raw), "paths": selected,
        "schemas": schema.get("components", {}).get("schemas", {})}), flush=True)
    return 0


def execute_args(args) -> int:
    entry, metadata, files, missing = preparation(args.paper)
    if args.action == "prepare":
        print(json.dumps({"paper": args.paper, "source_hashes_match": True,
                          "missing_declarations": missing, "network_requests": 0}))
        return 0
    if os.environ.get("COXETER_PUBLICATION_IO_GUARD") != "1":
        raise Stop("Network actions require tools/run_publication_io.py")
    if missing:
        raise Stop("Required personal declarations remain pending")
    if args.action == "resolve-preview":
        return resolve_preview(args, entry, metadata)
    if args.action == "preview-artifacts":
        return preview_artifacts(args, entry, metadata)
    if args.action == "inspect-preview":
        return inspect_preview(args, entry, metadata)
    if args.action == "commit":
        return commit_publication(args, entry, metadata)
    # Every action loads reviewed source bytes before the credential. No ZIP,
    # arbitrary path, bulk archive or local settings file can enter a payload.
    method, body, content_type, idem = "GET", None, "application/json", None
    if args.action == "create-draft":
        path = "/api/v1/submissions"
        method, body = "POST", canonical(metadata)
        idem = "coxeter-" + digest(canonical({"metadata": metadata,
                                              "files": entry["upload_files"]}))
    else:
        draft = component(args.draft_id, "draft identifier")
        bind_draft(draft, args.paper, metadata)
        path = "/api/v1/submissions/" + draft
        if args.action == "create-upload":
            bind_verified_metadata(draft, args.paper, metadata, entry)
            if args.file is None:
                raise Stop("A manifest file is required")
            path += "/uploads"
            method, body = "POST", canonical(upload_descriptor(entry, args.file))
        elif args.action in {"put-part", "complete-upload"}:
            bind_verified_metadata(draft, args.paper, metadata, entry)
            if args.file is None:
                raise Stop("A manifest file is required")
            upload = component(args.upload_id, "upload identifier")
            bind_upload(draft, upload, args.paper, args.file, metadata, entry)
            path += "/uploads/" + upload
            if args.action == "put-part":
                if args.file is None:
                    raise Stop("A manifest file is required")
                path += "/parts/1"
                method, body, content_type = "PUT", files[args.file], "application/octet-stream"
            else:
                require_uploaded_part(draft, upload, args.paper, args.file,
                                      metadata, entry, files[args.file])
                path += "/complete"
                method = "POST"
        elif args.action == "source-check":
            require_completed_files(draft, args.paper, metadata, entry)
            path += "/preview/check"
        elif args.action in {"preview", "preview-status"}:
            require_completed_files(draft, args.paper, metadata, entry)
            path += "/preview"
            if args.action == "preview":
                method = "POST"
    if args.action == "preview":
        return preview_with_lock(args, entry, metadata)
    token = credential()
    code, _ = exchange(args, entry, metadata, token, args.action, method, path,
                       body, content_type, idem)
    return 0 if 200 <= code < 300 else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper", choices=PAPERS)
    parser.add_argument("--action", choices=ACTIONS, default="prepare")
    parser.add_argument("--draft-id")
    parser.add_argument("--upload-id")
    parser.add_argument("--file", choices=tuple(sorted(HEXAGON_SOURCE_FILES)))
    parser.add_argument("--update-step", choices=("prepare", "status", "run", "finish"), default="prepare")
    args = parser.parse_args()
    if args.action not in {"run-queue", "inspect-queue", "apply-five-per-day", "refresh-statuses", "inspect-service", "completion-report", "inspect-update-api", "update-lp-vinberg"} and args.paper is None:
        parser.error("--paper is required except for run-queue")
    if args.action == "inspect-queue":
        return inspect_queue()
    if args.action == "refresh-statuses":
        return refresh_statuses()
    if args.action == "inspect-service":
        return inspect_service()
    if args.action == "completion-report":
        return completion_report()
    if args.action == "inspect-update-api":
        return inspect_update_api()
    if args.action == "update-lp-vinberg":
        import importlib.util
        path = ROOT/"release/hexagon_lp_vinberg_update_20261123/workflow.py"
        specification = importlib.util.spec_from_file_location("hexagon_lp_vinberg_update", path)
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
        return module.run(sys.modules[__name__], args.update_step)
    mutating = {"create-draft", "create-upload", "put-part", "complete-upload", "preview", "commit", "run-queue", "apply-five-per-day"}
    if args.action not in mutating:
        return execute_args(args)
    fd = os.open(HERE / "submission.lock.local", os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "a") as lock:
        os.fchmod(lock.fileno(), 0o600)
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise Stop("Another mutating submission operation is active") from None
        if PREVIEW_MARKER.exists() or COMMIT_MARKER.exists():
            raise Stop("A remote operation is unresolved; read its state before any mutation")
        args.operation_lock = lock
        if args.action == "apply-five-per-day":
            return apply_five_per_day()
        if args.action == "run-queue":
            return run_queue(args)
        return execute_args(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Stop as error:
        print(json.dumps({"stopped": str(error)}))
        raise SystemExit(1)
    except Exception:
        # Do not print tracebacks or arbitrary server/credential-derived text.
        print(json.dumps({"stopped": "Unexpected client error; no credentials printed"}))
        raise SystemExit(1)
