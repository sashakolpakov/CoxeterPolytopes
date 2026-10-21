#!/usr/bin/env python3
"""Freeze the exact compact-Coxeter label obstruction for H9/14 record 35."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

from coxeter_polytopes.combinatorics.h9_14_record35_coxeter_labels import (
    audit_record35_coxeter_labels,
)


ROOT = Path(__file__).resolve().parents[1]
MODULE = (
    ROOT
    / "src/coxeter_polytopes/combinatorics/h9_14_record35_coxeter_labels.py"
)
INCIDENCE_MODULE = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record35.py"
)
OUTPUT = ROOT / "data/enumerations/h9_14_record35_coxeter_labels.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "fcb5b0374c1ee68b32193491dcb1d905a7ac39ef61e7c56a7bbfe6bf6b56d23c"
)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_sha256(value: object) -> str:
    normalized = json.loads(json.dumps(value))
    return sha256(
        json.dumps(normalized, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _certificate_digest(manifest: Mapping[str, object]) -> str:
    critical = dict(manifest)
    critical.pop("certificate_sha256", None)
    return _canonical_sha256(critical)


def build() -> dict[str, object]:
    audit = audit_record35_coxeter_labels()
    if not audit.verified:
        raise AssertionError("record-35 label audit failed")
    census = audit.to_record()
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record35_compact_coxeter_label_obstruction",
        "claim": (
            "The ordinary-polytopal H9/14 record 35 cannot be the nerve of a "
            "compact hyperbolic Coxeter polytope.  Exhausting the 420 labelled "
            "order-five Lanner placements on each of three common-core missing "
            "faces leaves no elliptic completion at the rank-nine vertices. "
            "The all-right shadow obstruction is monotone under every possible "
            "additional non-right finite Coxeter label."
        ),
        "theorem": {
            "record_index": 35,
            "dimension": 9,
            "facet_count": 14,
            "compact_coxeter_realization_exists": False,
            "proof_scope": (
                "unrestricted in the finite Coxeter labels on meeting pairs "
                "outside the five forced Lanner supports"
            ),
            "not_claimed": [
                "a nonexistence theorem for finite-volume noncompact chambers",
                "a global Gram-matrix theorem for another nerve",
                "a new Coxeter polytope",
            ],
        },
        "proof": {
            "compact_inputs": [
                "each order-five minimal nonface supports a Lanner diagram",
                "each rank-nine vertex diagram is elliptic",
                "the classical order-five Lanner table has 420 labelled placements",
            ],
            "exact_reduction": (
                "For the supports {0,5,7,9,10}, {1,5,7,9,10}, and "
                "{4,5,7,9,10}, join the placement tables on their common "
                "four-node core.  Twelve rank-five restrictions leave 96 "
                "states.  Under S3 x S4 these form orbits of sizes 24 and 72; "
                "rank-six restrictions exhibit affine D4+A1 and affine D5."
            ),
            "monotone_shadow_lemma": (
                "For normalized Coxeter Gram G=I-A, A is symmetric and "
                "entrywise nonnegative.  Ellipticity is lambda_max(A)<1. "
                "Changing an external m=2 entry to any m>=3 increases A "
                "entrywise, hence cannot decrease lambda_max by the "
                "Perron-Frobenius monotonicity theorem.  Therefore no extra "
                "non-right labels can repair the certified non-elliptic "
                "principal subdiagram."
            ),
        },
        "census": census,
        "sources": {
            "label_gate_module": {
                "path": str(MODULE.relative_to(ROOT)),
                "sha256": _file_sha256(MODULE),
            },
            "record35_incidence_module": {
                "path": str(INCIDENCE_MODULE.relative_to(ROOT)),
                "sha256": _file_sha256(INCIDENCE_MODULE),
            },
        },
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != "h9_14_record35_compact_coxeter_label_obstruction":
        raise ValueError("unexpected record-35 Coxeter-label certificate kind")
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("record-35 Coxeter-label certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-35 Coxeter-label certificate digest mismatch")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-35 Coxeter-label certificate is not frozen")
    cluster = rebuilt["census"]["three_support_common_core_cluster"]
    return {
        "record_index": rebuilt["theorem"]["record_index"],
        "joined_common_core_states": cluster["join_counts"][-1],
        "rank_five_survivors": cluster["rank_five_restrictions"][
            "passing_states"
        ],
        "compact_coxeter_realization_exists": rebuilt["theorem"][
            "compact_coxeter_realization_exists"
        ],
        "certificate_sha256": digest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
        print(json.dumps(verify(manifest), indent=2, sort_keys=True))
        return
    manifest = build()
    if args.write:
        OUTPUT.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(
        json.dumps(
            {
                "output": str(OUTPUT.relative_to(ROOT)),
                "written": args.write,
                "certificate_sha256": manifest["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
