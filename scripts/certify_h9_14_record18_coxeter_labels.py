#!/usr/bin/env python3
"""Freeze the selected record-18 compact-Coxeter no-go certificate."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record18_coxeter_labels import (
    audit_record18_compact_coxeter_labels,
)
from scripts.independent_h9_14_record18_coxeter_audit import (
    EXPECTED_SUMMARY_SHA256,
    replay as independent_replay,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record18_coxeter_labels.json"
SELECTION_CERTIFICATE = ROOT / "data/enumerations/h9_14_record18_selection.json"
SELECTION_FILE_SHA256 = (
    "3791c7f0a20bc6bfae365a59c54ad27ae3f585aa52e7febf5fb9b5976f37e808"
)
SELECTION_CERTIFICATE_SHA256 = (
    "27193c530804237c8505b7ce7400e7551a6b152d5a4391f2dc624371d8f54406"
)
PROOF_SOURCES = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record18_coxeter_labels.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record18.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_disjoint_lanner_gate.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_selected_nerve_coxeter_no_go.py",
    ROOT / "scripts/independent_h9_14_record18_coxeter_audit.py",
)
EXPECTED_PROOF_SOURCE_SHA256 = {
    "src/coxeter_polytopes/combinatorics/h9_14_record18_coxeter_labels.py": (
        "ed0fafb5c58a46e07657bc1896cac9bbbbea675a0541b6ebefee1bf7e8d63739"
    ),
    "src/coxeter_polytopes/combinatorics/h9_14_record18.py": (
        "c96d21616004da447a43c6c803a1b27583a2c997454eba64ccafc2eeee512e15"
    ),
    "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py": (
        "53d7e82a66d3cbb3b091109a6b79ffbdab271d15c9d098f1862a56355f126c43"
    ),
    "src/coxeter_polytopes/combinatorics/h9_14_disjoint_lanner_gate.py": (
        "3a1b58fe874f84b4a3bacf5e23654b3c4832246f0ae3e0c9b688c84fa3da68ce"
    ),
    "src/coxeter_polytopes/combinatorics/h9_14_selected_nerve_coxeter_no_go.py": (
        "c87e4239ef69beabf97856c64d876e6e5369dbc47c83f89cfc490c01db7507ec"
    ),
    "scripts/independent_h9_14_record18_coxeter_audit.py": (
        "eb66db5b487a5073fdb8866f5b73b2639acb97859e517585becf0dee0ca56b06"
    ),
}
EXPECTED_CERTIFICATE_SHA256 = (
    "0b168250c951ed814f0678b01122f67de8da8943a3d37bd4f5f1d881aed6c88a"
)


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def certificate_digest(document: Mapping[str, object]) -> str:
    normalized = dict(document)
    normalized.pop("certificate_sha256", None)
    return sha256(
        json.dumps(normalized, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def proof_source_records() -> list[dict[str, str]]:
    return [
        {"path": str(path.relative_to(ROOT)), "sha256": file_sha256(path)}
        for path in PROOF_SOURCES
    ]


def build() -> dict[str, object]:
    if file_sha256(SELECTION_CERTIFICATE) != SELECTION_FILE_SHA256:
        raise AssertionError("record-18 selected-nerve certificate file changed")
    selection = json.loads(SELECTION_CERTIFICATE.read_text(encoding="utf-8"))
    if (
        selection.get("certificate_sha256") != SELECTION_CERTIFICATE_SHA256
        or selection.get("status") != "abstract_audited_nerve_candidate"
    ):
        raise AssertionError("record-18 selected nerve source changed")
    production = audit_record18_compact_coxeter_labels()
    if not production.unrestricted_compact_coxeter_no_go:
        raise AssertionError("record-18 compact no-go disappeared")
    independent = independent_replay()
    if independent["summary_sha256"] != EXPECTED_SUMMARY_SHA256:
        raise AssertionError("record-18 compact hostile replay changed")

    document: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record18_selected_nerve_compact_coxeter_no_go",
        "status": "certified_selected_nerve_compact_coxeter_no_go",
        "claim": (
            "The explicit abstract record-18 nerve candidate is not the "
            "nerve of a compact hyperbolic Coxeter 9-polytope."
        ),
        "scope": (
            "One selected record-18 incidence completion only; other "
            "completions of the graph shadow remain outside the theorem."
        ),
        "nonclaims": [
            "No other completion of record 18 is excluded.",
            "No PL-sphericity or ordinary convex-polytopality theorem is added.",
            "No finite-volume noncompact Coxeter realization is excluded.",
            "No Vinberg parabolic-cusp finite-volume audit is performed.",
            "No global Coxeter labeling or Lorentzian Gram realization is supplied.",
        ],
        "selected_nerve_source": {
            "path": str(SELECTION_CERTIFICATE.relative_to(ROOT)),
            "file_sha256": SELECTION_FILE_SHA256,
            "certificate_sha256": SELECTION_CERTIFICATE_SHA256,
            "status": selection["status"],
        },
        "proof_sources": proof_source_records(),
        "production_audit": production.to_record(),
        "independent_hostile_replay": independent,
        "proof_logic": {
            "compact_vertex_condition": (
                "Every actual rank-nine vertex diagram must be elliptic."
            ),
            "minimal_nonface_condition": (
                "Each displayed minimal five-face supports a labelled "
                "order-five Lanner diagram in a compact realization."
            ),
            "monotonicity": (
                "For normalized Gram I-A, increasing a finite cross label "
                "or adding cross edges increases A entrywise, so a local "
                "nonellipticity failure at m=3 cannot be repaired."
            ),
            "inertia_obstruction": (
                "Zero cross edges give an orthogonal sum of two Lanner "
                "blocks; every surviving one-edge block has exact negative "
                "index two. Principal interlacing contradicts ambient "
                "compact H9 Gram inertia (9,1,4)."
            ),
        },
        "vinberg_scope": {
            "compact_case_excluded": True,
            "finite_volume_noncompact_checked": False,
            "elliptic_edge_endpoint_audit_run": False,
            "maximal_parabolic_cusp_audit_run": False,
        },
    }
    document = json.loads(json.dumps(document, sort_keys=True))
    document["certificate_sha256"] = certificate_digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record18_selected_nerve_compact_coxeter_no_go":
        raise ValueError("wrong record-18 compact-gate certificate kind")
    for path in PROOF_SOURCES:
        relative = str(path.relative_to(ROOT))
        if file_sha256(path) != EXPECTED_PROOF_SOURCE_SHA256[relative]:
            raise ValueError(f"record-18 compact proof source changed: {relative}")
    rebuilt = build()
    if document != rebuilt:
        raise ValueError("record-18 compact certificate did not reproduce")
    digest = certificate_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("record-18 compact certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("record-18 compact certificate is not frozen")
    return {
        "record_index": 18,
        "selected_nerve_only": True,
        "unrestricted_compact_coxeter_no_go": True,
        "finite_volume_noncompact_checked": False,
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    if args.verify:
        document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    else:
        document = build()
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(verify(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
