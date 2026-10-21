#!/usr/bin/env python3
"""Freeze and verify the explicit record-30 a=0 completion and label no-go."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

from scripts.independent_h9_14_record30_a0_candidate_audit import audit


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/independent_h9_14_record30_a0_candidate_audit.py"
ARCHIVE = ROOT / "data/enumerations/h9_14_record30_a0_candidate.json"
EXPECTED_SOURCE_SHA256 = "ba1c71006c2768802ceb5108f15c6297ce2457b0de44f19fb3c1094633c7ecc2"
EXPECTED_CERTIFICATE_SHA256 = "e40b22cc094901627c31186405faf6bd238dfd0633440677e788d02b5003b5fd"


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if sha256(SOURCE.read_bytes()).hexdigest() != EXPECTED_SOURCE_SHA256:
        raise ValueError("record-30 a0 independent source digest changed")
    if document.get("kind") != "independent_h9_14_record30_a0_candidate_audit":
        raise ValueError("wrong record-30 a0 candidate kind")
    if document.get("certificate_sha256") != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("record-30 a0 candidate certificate digest changed")
    critical = dict(document)
    digest = critical.pop("certificate_sha256", None)
    rebuilt = sha256(
        json.dumps(critical, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if digest != rebuilt:
        raise ValueError("record-30 a0 candidate internal digest is invalid")
    if (
        document.get("minimal_transversal_size_histogram") != {"5": 132}
        or document.get("face_vector")
        != [14, 88, 328, 806, 1364, 1596, 1248, 594, 132]
        or document.get("ridge_degree_histogram") != {"2": 594}
        or document.get("dual_component_count") != 1
        or document.get("abstract_mod2_homology_sphere_gate_passes") is not True
    ):
        raise ValueError("record-30 a0 abstract-sphere audit changed")
    multiwedge = document.get("multiwedge")
    if not isinstance(multiwedge, Mapping) or (
        multiwedge.get("exact_multiwedge_identity") is not True
        or multiwedge.get("base_f_vector") != [10, 42, 88, 90, 36]
        or multiwedge.get("convex_polytopality_proved") is not False
    ):
        raise ValueError("record-30 a0 multiwedge audit changed")
    coxeter = document.get("coxeter_label_obstruction")
    if not isinstance(coxeter, Mapping) or (
        coxeter.get("terminal_state_count") != 0
        or coxeter.get("unrestricted_compact_coxeter_no_go_for_this_completion")
        is not True
        or [row["locally_elliptic_state_count"] for row in coxeter["production_order"]["stages"]]
        != [420, 96, 6, 0]
        or [row["locally_elliptic_state_count"] for row in coxeter["independent_reverse_order"]["stages"]]
        != [420, 96, 0]
    ):
        raise ValueError("record-30 a0 Coxeter-label obstruction changed")
    if (
        document.get("convex_polytopality_proved") is not False
        or document.get("finite_volume_noncompact_chamber_excluded") is not False
        or document.get("compact_only") is not True
    ):
        raise ValueError("record-30 a0 candidate scope was over-promoted")
    return {
        "facet_count": 132,
        "ridge_count": 594,
        "base_multiwedge_facet_count": multiwedge["base_facet_count"],
        "compact_label_terminal_states": 0,
        "convex_polytopality_proved": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }


def fresh_record() -> dict[str, object]:
    result = audit()
    verify(result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args()
    if args.refresh:
        document = fresh_record()
        ARCHIVE.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    else:
        document = json.loads(ARCHIVE.read_text())
        verify(document)
        if args.recompute and fresh_record() != document:
            raise ValueError("fresh record-30 a0 replay differs from archive")
    print(json.dumps(verify(document), sort_keys=True))


if __name__ == "__main__":
    main()
