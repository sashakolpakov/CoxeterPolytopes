#!/usr/bin/env python3
"""Freeze the selected record-14 compact-Coxeter no-go certificate."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record14_coxeter_labels import (
    audit_record14_compact_coxeter_labels,
)
from scripts.independent_h9_14_record14_coxeter_audit import (
    EXPECTED_SUMMARY_SHA256,
    replay as independent_replay,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record14_coxeter_labels.json"
POLYTOPALITY_CERTIFICATE = (
    ROOT / "data/enumerations/h9_14_record14_polytopality.json"
)
POLYTOPALITY_FILE_SHA256 = (
    "01ac2e8ca83a72f40dd813ee0051081911bab203be01966902882eb61d751dc1"
)
POLYTOPALITY_CERTIFICATE_SHA256 = (
    "0fba773d1449780048625c0243eb4740efc0ccd23d87b266e5f3ca06deb7edb1"
)
PROOF_SOURCES = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record14_coxeter_labels.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record14.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_disjoint_lanner_gate.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_selected_nerve_coxeter_no_go.py",
    ROOT / "scripts/independent_h9_14_record14_coxeter_audit.py",
)
EXPECTED_PROOF_SOURCE_SHA256 = {
    "src/coxeter_polytopes/combinatorics/h9_14_record14_coxeter_labels.py": (
        "9e47fd8f4f16f65d7596adf43562e3842ca092b97b474ac98ccec63f924b9771"
    ),
    "src/coxeter_polytopes/combinatorics/h9_14_record14.py": (
        "9c83d807df43f364455164e7a120c2600558fb965400b2aa61d166de2a7c1891"
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
    "scripts/independent_h9_14_record14_coxeter_audit.py": (
        "672eba2fc3091c46af3001fa72656d1dffe99eb3426a88949dab6dbbebc2cc9a"
    ),
}
EXPECTED_CERTIFICATE_SHA256 = (
    "c31d43853470a12c6b88fa95e01bbb7b8078ce314a00b9342cc1fae2c1ad8f4e"
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


def _proof_sources() -> list[dict[str, str]]:
    return [
        {"path": str(path.relative_to(ROOT)), "sha256": file_sha256(path)}
        for path in PROOF_SOURCES
    ]


def build() -> dict[str, object]:
    if file_sha256(POLYTOPALITY_CERTIFICATE) != POLYTOPALITY_FILE_SHA256:
        raise AssertionError("record-14 polytopality certificate file changed")
    polytopality = json.loads(POLYTOPALITY_CERTIFICATE.read_text(encoding="utf-8"))
    if (
        polytopality.get("certificate_sha256")
        != POLYTOPALITY_CERTIFICATE_SHA256
        or polytopality.get("status")
        != "certified_ordinary_polytopal_completion"
    ):
        raise AssertionError("record-14 selected nerve source changed")

    production = audit_record14_compact_coxeter_labels()
    if not production.unrestricted_compact_coxeter_no_go:
        raise AssertionError("record-14 compact no-go disappeared")
    independent = independent_replay()
    if independent["summary_sha256"] != EXPECTED_SUMMARY_SHA256:
        raise AssertionError("record-14 compact hostile replay changed")

    document: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record14_selected_nerve_compact_coxeter_no_go",
        "status": "certified_selected_nerve_compact_coxeter_no_go",
        "claim": (
            "The explicitly Gale-realized record-14 completion is not the "
            "nerve of a compact hyperbolic Coxeter 9-polytope."
        ),
        "scope": (
            "One selected record-14 incidence completion only; other "
            "completions of the graph shadow remain outside the theorem."
        ),
        "nonclaims": [
            "No other completion of record 14 is excluded.",
            "No finite-volume noncompact Coxeter realization is excluded.",
            "No Vinberg parabolic-cusp finite-volume audit is performed.",
            "No global Coxeter labeling or Lorentzian Gram realization is supplied.",
        ],
        "selected_nerve_source": {
            "path": str(POLYTOPALITY_CERTIFICATE.relative_to(ROOT)),
            "file_sha256": POLYTOPALITY_FILE_SHA256,
            "certificate_sha256": POLYTOPALITY_CERTIFICATE_SHA256,
        },
        "proof_sources": _proof_sources(),
        "production_audit": production.to_record(),
        "independent_hostile_replay": independent,
        "proof_logic": {
            "compact_vertex_condition": (
                "Every actual rank-nine vertex diagram must be elliptic."
            ),
            "minimal_nonface_condition": (
                "Each of the two displayed minimal five-faces supports a "
                "labelled order-five Lanner diagram."
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
    if document.get("kind") != "h9_14_record14_selected_nerve_compact_coxeter_no_go":
        raise ValueError("wrong record-14 compact-gate certificate kind")
    for path in PROOF_SOURCES:
        relative = str(path.relative_to(ROOT))
        if file_sha256(path) != EXPECTED_PROOF_SOURCE_SHA256[relative]:
            raise ValueError(f"record-14 compact proof source changed: {relative}")
    rebuilt = build()
    if document != rebuilt:
        raise ValueError("record-14 compact certificate did not reproduce")
    digest = certificate_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("record-14 compact certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("record-14 compact certificate is not frozen")
    return {
        "record_index": 14,
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
