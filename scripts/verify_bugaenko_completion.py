"""Verify Bugaenko's search transcript and compact chamber together."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.verify import (
    CertificateError,
    VinbergTranscriptError,
    verify_bugaenko_vinberg_completion,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rerun the exact Vinberg enumeration and compactness verifier."
    )
    parser.add_argument("transcript", type=Path)
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()
    with args.transcript.open(encoding="utf-8") as handle:
        transcript = json.load(handle)
    with args.certificate.open(encoding="utf-8") as handle:
        certificate = json.load(handle)
    try:
        result = verify_bugaenko_vinberg_completion(transcript, certificate)
    except (CertificateError, VinbergTranscriptError) as exc:
        raise SystemExit(f"REJECTED: {exc}") from exc
    print(
        json.dumps(
            {
                "status": "verified_complete_compact_reflection_chamber",
                "dimension": result.chamber.dimension,
                "facet_count": result.chamber.facet_count,
                "finite_vertex_count": len(result.chamber.finite_vertices),
                "cusp_count": len(result.chamber.cusps),
                "vinberg_height_count": len(result.search.levels),
                "positive_height_root_counts": [
                    len(level.roots) for level in result.search.levels
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
