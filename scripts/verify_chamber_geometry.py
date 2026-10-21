"""Verify the exact geometry portion of a chamber certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.verify import CertificateError, verify_chamber_geometry


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Independently recompute Gram data, Coxeter relations, and all "
            "finite/ideal/hyperideal extreme rays. Vinberg wall completeness "
            "is a separate certificate stage."
        )
    )
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()
    with args.certificate.open(encoding="utf-8") as handle:
        certificate = json.load(handle)
    try:
        result = verify_chamber_geometry(certificate)
    except CertificateError as exc:
        raise SystemExit(f"REJECTED: {exc}") from exc
    print(
        json.dumps(
            {
                "status": "geometry_verified_completeness_not_checked",
                "dimension": result.dimension,
                "facet_count": result.facet_count,
                "finite_vertices": [list(face) for face in result.finite_vertices],
                "cusps": [list(face) for face in result.cusps],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
