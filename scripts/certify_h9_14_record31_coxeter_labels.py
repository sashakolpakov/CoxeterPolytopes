#!/usr/bin/env python3
"""Freeze the selected record-31 compact-Coxeter no-go certificate."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record31_coxeter_labels import (
    audit_record31_compact_coxeter_labels,
)
from scripts.independent_h9_14_record31_coxeter_audit import (
    EXPECTED_SUMMARY_SHA256,
    replay as independent_replay,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record31_coxeter_labels.json"
SELECTION_CERTIFICATE = ROOT / "data/enumerations/h9_14_record31_polytopality.json"
SELECTION_FILE_SHA256 = (
    "c7ac4053b47d8ae7dd863de620297684a7f75227e85c53aef564fe7b040ad83b"
)
SELECTION_CERTIFICATE_SHA256 = (
    "12b01aec7f3de0bb4e71040291753e8362787b9c7063704877bd24693884a799"
)
PROOF_SOURCES = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record31_coxeter_labels.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record31.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_disjoint_lanner_gate.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_selected_nerve_coxeter_no_go.py",
    ROOT / "scripts/independent_h9_14_record31_coxeter_audit.py",
)
EXPECTED_PROOF_SOURCE_SHA256 = {
    "src/coxeter_polytopes/combinatorics/h9_14_record31_coxeter_labels.py": (
        "22ab79f6e8b7d0e15372299a5f56dc918cae572b9a2030effd3aaf092f98ee0e"
    ),
    "src/coxeter_polytopes/combinatorics/h9_14_record31.py": (
        "ad30ec2b567c865165d805f1f2375211cf03aec766b574e8fe2fcf6cb987f890"
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
    "scripts/independent_h9_14_record31_coxeter_audit.py": (
        "cb910ac7c3dab214d49738a520a0b39e7b4a3cdfdf9a56e4a94f2c3c8a5ece24"
    ),
}
EXPECTED_CERTIFICATE_SHA256 = (
    "92aacb8b7d4ab3ca9faa0eefeba972d1736e6e919376bd77409af64424857b97"
)


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def certificate_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def proof_source_records() -> list[dict[str, str]]:
    return [
        {"path": str(path.relative_to(ROOT)), "sha256": file_sha256(path)}
        for path in PROOF_SOURCES
    ]


def build() -> dict[str, object]:
    if file_sha256(SELECTION_CERTIFICATE) != SELECTION_FILE_SHA256:
        raise AssertionError("record-31 selected-nerve file changed")
    selection = json.loads(SELECTION_CERTIFICATE.read_text(encoding="utf-8"))
    if (
        selection.get("certificate_sha256") != SELECTION_CERTIFICATE_SHA256
        or selection.get("theorem_status")
        != "certified_explicit_ordinary_polytopal"
    ):
        raise AssertionError("record-31 selected nerve source changed")

    production = audit_record31_compact_coxeter_labels()
    if not production.unrestricted_compact_coxeter_no_go:
        raise AssertionError("record-31 compact no-go disappeared")
    independent = independent_replay()
    if independent["summary_sha256"] != EXPECTED_SUMMARY_SHA256:
        raise AssertionError("record-31 hostile replay changed")

    document: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record31_selected_nerve_compact_coxeter_no_go",
        "status": "certified_selected_nerve_compact_coxeter_no_go",
        "claim": (
            "The explicit ordinary-polytopal record-31 nerve is not the "
            "nerve of a compact hyperbolic Coxeter 9-polytope."
        ),
        "scope": (
            "One selected record-31 incidence completion only; other "
            "completions of the graph shadow remain outside the theorem."
        ),
        "nonclaims": [
            "No other completion of graph-shadow record 31 is excluded.",
            "No finite-volume noncompact Coxeter realization is excluded.",
            "No Vinberg parabolic-cusp finite-volume audit is performed.",
            "No global Coxeter labeling or Lorentzian Gram realization is supplied.",
        ],
        "selected_nerve_source": {
            "path": str(SELECTION_CERTIFICATE.relative_to(ROOT)),
            "file_sha256": SELECTION_FILE_SHA256,
            "certificate_sha256": SELECTION_CERTIFICATE_SHA256,
            "status": selection["theorem_status"],
        },
        "proof_sources": proof_source_records(),
        "production_audit": production.to_record(),
        "independent_hostile_replay": independent,
        "proof_logic": {
            "compact_vertex_condition": (
                "Every actual rank-nine vertex diagram must be elliptic."
            ),
            "minimal_nonface_condition": (
                "All nine displayed missing five-faces must simultaneously "
                "support labelled order-five Lanner diagrams."
            ),
            "symmetry_coverage": (
                "Four independent clone transpositions preserve the complete "
                "selected nerve and reduce each 420-state Lanner support to "
                "110 representatives; their actions on the two supports "
                "factor, so both placements can be normalized simultaneously."
            ),
            "cross_edge_exhaustion": (
                "No cross edge gives two orthogonal negative directions. "
                "At least two non-right cross edges fail an actual-vertex "
                "test, as do all labels m=4,5 in the one-edge case."
            ),
            "one_edge_residue": (
                "The 1,932 locally surviving m=3 representative trials split "
                "into four rooted block classes. The two large classes, of "
                "total size 1,860, have exact principal-block inertia (8,2,0)."
            ),
            "simultaneous_lanner_csp": (
                "The remaining 72 representative trials admit zero mutually "
                "consistent choices from the complete 420-placement tables "
                "on the other seven missing five-faces."
            ),
            "ambient_inertia": (
                "A compact fourteen-wall H9 chamber needs Gram inertia "
                "(9,1,4); a principal block with negative index two is "
                "impossible by Cauchy interlacing."
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
    if document.get("kind") != "h9_14_record31_selected_nerve_compact_coxeter_no_go":
        raise ValueError("wrong record-31 compact-gate certificate kind")
    for path in PROOF_SOURCES:
        relative = str(path.relative_to(ROOT))
        if file_sha256(path) != EXPECTED_PROOF_SOURCE_SHA256[relative]:
            raise ValueError(f"record-31 compact proof source changed: {relative}")
    rebuilt = build()
    if document != rebuilt:
        raise ValueError("record-31 compact certificate did not reproduce")
    digest = certificate_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("record-31 compact certificate digest mismatch")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-31 compact certificate is not frozen")
    return {
        "record_index": 31,
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
