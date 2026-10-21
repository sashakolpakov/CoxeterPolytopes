#!/usr/bin/env python3
"""Small publication-I/O regressions; root must run under its I/O admission.

This fixture uses temporary synthetic bytes, never repository research data.
It is not a mathematical checker. No network is used. Importing this file does
not import either utility; main() does that only when explicitly executed.
"""
from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha(data):
    return hashlib.sha256(data).hexdigest()


class PublicDataIO(unittest.TestCase):
    packer = fetcher = None

    def setUp(self):
        self.area = tempfile.TemporaryDirectory(prefix="public-data-io-")
        self.addCleanup(self.area.cleanup)
        self.root = Path(self.area.name)
        for directory in ("data", "results", "repro"):
            (self.root / directory).mkdir()
        self.addCleanup(mock.patch.stopall)
        mock.patch.object(self.packer, "ROOT", self.root).start()

    def file(self, name="data/example.bin", data=b"synthetic fixture\n"):
        target = self.root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        return {"path": name, "bytes": len(data), "sha256": sha(data)}

    def archive(self, rows, *, member_type=None):
        target = self.root / "part.tar.gz"
        with tarfile.open(target, "w:gz") as output:
            for row in rows:
                info = tarfile.TarInfo(row["path"])
                if member_type is not None:
                    info.type, info.linkname = member_type, "../outside"
                    output.addfile(info)
                else:
                    payload = row.pop("payload", b"fixture")
                    info.size = len(payload)
                    output.addfile(info, io.BytesIO(payload))
        return target, {"name": target.name, "bytes": target.stat().st_size,
                        "sha256": self.fetcher.digest(target), "files": rows}

    def call_main(self, module, args):
        with mock.patch.object(sys, "argv", [module.__file__, *args]), \
                contextlib.redirect_stdout(io.StringIO()):
            return module.main()

    def test_pack_restore_and_originals_unchanged(self):
        payload = b"original data\0" * 1000
        row = self.file(data=payload)
        archive = self.root / "roundtrip.tar.gz"
        part = self.packer.pack([row], archive)
        destination = self.root / "restored"
        destination.mkdir()
        self.assertTrue(part["roundtrip_verified"])
        self.assertEqual(part["files"][0]["sha256"], sha(payload))
        self.assertEqual(self.fetcher.restore_part(
            archive, part, destination, {row["path"]}), 1)
        self.assertEqual((destination / row["path"]).read_bytes(), payload)
        self.assertEqual((self.root / row["path"]).read_bytes(), payload)
        self.assertEqual(self.fetcher.restore_part(
            archive, part, destination, {row["path"]}), 0)

    def test_local_directory_restore_keeps_asset(self):
        row = self.file()
        asset = self.root / "asset.tar.gz"
        part = self.packer.pack([row], asset)
        manifest = self.root / "manifest.json"
        manifest.write_text(json.dumps({"schema_version": 1,
            "repository": "sashakolpakov/CoxeterPolytopes",
            "release_tag": "fixture", "original_bytes": row["bytes"],
            "scope": "Synthetic publication-I/O fixture only", "parts": [part]}))
        destination = self.root / "local-restore"
        with mock.patch.object(self.fetcher.urllib.request, "urlopen",
                               side_effect=AssertionError("network forbidden")):
            self.call_main(self.fetcher, ["--restore", "--manifest", str(manifest),
                "--destination", str(destination), "--from-directory", str(self.root)])
        self.assertTrue(asset.is_file())
        self.assertEqual((destination / row["path"]).read_bytes(),
                         (self.root / row["path"]).read_bytes())

    def test_existing_mismatch_is_preserved(self):
        row = self.file(data=b"right")
        target = self.root / row["path"]
        target.write_bytes(b"wrong")
        with self.assertRaises(ValueError):
            self.fetcher.existing_matches(target, row)
        self.assertEqual(target.read_bytes(), b"wrong")

    def test_traversal_and_symlink_destinations(self):
        for name in ("../outside", "/tmp/outside", "data/../../outside"):
            with self.assertRaises(ValueError):
                self.fetcher.safe_target(self.root, name)
        (self.root / "data/link").symlink_to(self.root / "results", target_is_directory=True)
        with self.assertRaises(ValueError):
            self.fetcher.safe_target(self.root, "data/link/file")

    def test_link_archive_member_is_rejected(self):
        row = {"path": "data/link", "bytes": 0, "sha256": sha(b"")}
        archive, part = self.archive([row], member_type=tarfile.SYMTYPE)
        with self.assertRaises(ValueError):
            self.fetcher.restore_part(archive, part, self.root, {row["path"]})
        self.assertFalse((self.root / row["path"]).is_symlink())

    def test_wrong_restored_digest_leaves_no_output_or_partial(self):
        row = {"path": "data/bad", "bytes": 7, "sha256": sha(b"different")}
        archive, part = self.archive([row])
        with self.assertRaises(ValueError):
            self.fetcher.restore_part(archive, part, self.root, {row["path"]})
        self.assertFalse((self.root / "data/bad").exists())
        self.assertFalse((self.root / "data/bad.public-data.partial").exists())

    def test_preexisting_dangling_download_is_not_deleted(self):
        target = self.root / "asset.tar.gz"
        partial = target.with_suffix(target.suffix + ".download")
        partial.symlink_to(self.root / "absent")
        with mock.patch.object(self.fetcher.urllib.request, "urlopen",
                               side_effect=OSError("fixture refusal")):
            with self.assertRaises((OSError, ValueError)):
                self.fetcher.download("https://invalid.invalid/", target,
                    {"name": target.name, "bytes": 0, "sha256": sha(b"")})
        self.assertTrue(partial.is_symlink())

    def test_restore_open_race_does_not_delete_unowned_partial(self):
        row = {"path": "data/race", "bytes": 7, "sha256": sha(b"fixture")}
        archive, part = self.archive([row])
        partial = self.root / "data/race.public-data.partial"
        original_open = Path.open

        def competing_open(path, *args, **kwargs):
            if path == partial and args and args[0] == "xb":
                with original_open(path, "wb") as output:
                    output.write(b"other writer")
                raise FileExistsError("fixture concurrent creator")
            return original_open(path, *args, **kwargs)

        with mock.patch.object(Path, "open", competing_open):
            with self.assertRaises(FileExistsError):
                self.fetcher.restore_part(archive, part, self.root, {row["path"]})
        self.assertEqual(partial.read_bytes(), b"other writer")

    def seed_incomplete_pack(self):
        rows = [self.file("data/a.bin", b"a" * self.packer.MIB),
                self.file("data/b.bin", b"b" * self.packer.MIB)]
        # Exercise main's state/resume contract using two SMALL synthetic parts.
        # The grouping formula itself is checked separately below.
        patch = mock.patch.object(self.packer, "groups", side_effect=
            lambda actual: [[r] for r in actual])
        patch.start()
        self.addCleanup(patch.stop)
        self.call_main(self.packer, ["--max-parts", "1"])
        self.assertFalse((self.root / "repro/public-data-manifest.json").exists())
        return rows

    def test_same_length_source_change_cannot_earn_resume_credit(self):
        rows = self.seed_incomplete_pack()
        (self.root / rows[0]["path"]).write_bytes(b"z" * self.packer.MIB)
        with self.assertRaises((RuntimeError, ValueError, OSError)):
            self.call_main(self.packer, ["--max-parts", "1"])
        self.assertFalse((self.root / "repro/public-data-manifest.json").exists())

    def test_missing_completed_archive_cannot_earn_resume_credit(self):
        self.seed_incomplete_pack()
        (self.root / ".publication/data/research-data-001.tar.gz").unlink()
        with self.assertRaises((RuntimeError, ValueError, OSError)):
            self.call_main(self.packer, ["--max-parts", "1"])
        self.assertFalse((self.root / "repro/public-data-manifest.json").exists())

    def test_earlier_source_changed_during_last_part_blocks_manifest(self):
        rows = self.seed_incomplete_pack()
        original_pack = self.packer.pack

        def mutate_after_last_part(group, target):
            part = original_pack(group, target)
            if target.name == "research-data-002.tar.gz":
                (self.root / rows[0]["path"]).write_bytes(b"z" * self.packer.MIB)
            return part

        with mock.patch.object(self.packer, "pack", mutate_after_last_part):
            with self.assertRaises((RuntimeError, ValueError, OSError)):
                self.call_main(self.packer, ["--max-parts", "1"])
        self.assertFalse((self.root / "repro/public-data-manifest.json").exists())

    def test_whole_file_grouping_never_splits_large_member(self):
        rows = [{"path": "a", "bytes": 600 * self.packer.MIB},
                {"path": "b", "bytes": self.packer.MIB}]
        self.assertEqual(self.packer.groups(rows), [[rows[0]], [rows[1]]])


def main():
    PublicDataIO.packer = load("publication_packer_fixture", ROOT / "tools/prepare_public_data.py")
    PublicDataIO.fetcher = load("publication_fetcher_fixture", ROOT / "scripts/fetch_public_data.py")
    unittest.main(argv=[sys.argv[0]])


if __name__ == "__main__":
    main()
