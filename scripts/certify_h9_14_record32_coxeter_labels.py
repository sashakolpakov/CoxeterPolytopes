#!/usr/bin/env python3
"""Build or verify the record-32 unrestricted compact-Coxeter no-go."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record32_coxeter_labels import (
    audit_record32_coxeter_labels,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record32_coxeter_labels.json"
MODULE = (
    ROOT
    / "src/coxeter_polytopes/combinatorics/h9_14_record32_coxeter_labels.py"
)
DEPENDENCIES = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_disjoint_lanner_gate.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record32.py",
)
POLYTOPALITY = ROOT / "data/enumerations/h9_14_record32_polytopality.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "15bf62a37b5f65fdf7118e80be347482d30b0e47b84eeab2dd831f089b5b4412"
)


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def certificate_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
    ).hexdigest()


def generate() -> dict[str, object]:
    audit = audit_record32_coxeter_labels()
    document: dict[str, object] = {
        "schema_version": 1,
        "object": "unrestricted compact-Coxeter obstruction for H9/14 record 32",
        "source": {
            "audit_module": {
                "path": "src/coxeter_polytopes/combinatorics/h9_14_record32_coxeter_labels.py",
                "sha256": _file_sha256(MODULE),
            },
            "imported_dependency_sha256": {
                str(path.relative_to(ROOT)): _file_sha256(path)
                for path in DEPENDENCIES
            },
            "ordinary_polytopality_certificate": {
                "path": "data/enumerations/h9_14_record32_polytopality.json",
                "sha256": _file_sha256(POLYTOPALITY),
            },
        },
        "theorem": audit.to_record(),
        "proof_scope": {
            "compact_coxeter_nerve_excluded": True,
            "ordinary_polytopality_retained": True,
            "dotted_values_enumerated": False,
            "why_dotted_values_are_irrelevant": (
                "The decisive negative-index-two principal block contains no "
                "dotted pair, and the alternative cross-edge obstruction is "
                "already a finite actual-vertex obstruction."
            ),
        },
        "nonclaims": [
            "No compact H9 Coxeter polytope is constructed.",
            "No finite-volume noncompact chamber with a changed ideal-face structure is classified.",
            "The theorem is scoped to the explicit record-32 nerve.",
        ],
    }
    document["certificate_sha256"] = certificate_digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    source = document.get("source")
    if not isinstance(source, Mapping):
        raise ValueError("record-32 label certificate lacks source metadata")
    module = source.get("audit_module")
    dependencies = source.get("imported_dependency_sha256")
    polytope = source.get("ordinary_polytopality_certificate")
    if not isinstance(module, Mapping) or module.get("sha256") != _file_sha256(MODULE):
        raise ValueError("record-32 Coxeter-label module digest changed")
    expected_dependencies = {
        str(path.relative_to(ROOT)): _file_sha256(path) for path in DEPENDENCIES
    }
    if dependencies != expected_dependencies:
        raise ValueError("record-32 Coxeter-label dependency digest changed")
    if (
        not isinstance(polytope, Mapping)
        or polytope.get("sha256") != _file_sha256(POLYTOPALITY)
    ):
        raise ValueError("record-32 polytopality input digest changed")
    if document.get("certificate_sha256") != certificate_digest(document):
        raise ValueError("record-32 Coxeter-label certificate digest changed")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and document.get("certificate_sha256") != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-32 Coxeter-label certificate digest is not frozen")
    rebuilt = generate()
    if document != rebuilt:
        raise ValueError("record-32 Coxeter-label certificate does not replay")
    theorem = document.get("theorem")
    if not isinstance(theorem, Mapping) or theorem.get("verified") is not True:
        raise ValueError("record-32 Coxeter-label theorem is not verified")
    no_go = theorem.get("unrestricted_compact_coxeter_no_go")
    if not isinstance(no_go, Mapping) or no_go.get("proved") is not True:
        raise ValueError("record-32 unrestricted no-go disappeared")
    return {
        "a4_compatible_states": theorem["overlapping_order_four_gate"][
            "compatible_join_count"
        ],
        "a4_passing_states": theorem["overlapping_order_four_gate"][
            "passing_state_count"
        ],
        "sparse_states": theorem["sparse_locally_elliptic_state_count"],
        "cross_tests": theorem["cross_support_gate"]["test_count"],
        "cross_survivors": theorem["cross_support_gate"]["passing_count"],
        "unrestricted_compact_coxeter_no_go": True,
        "certificate_sha256": document["certificate_sha256"],
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    if args.verify:
        document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    else:
        document = generate()
    if args.write:
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(verify(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
