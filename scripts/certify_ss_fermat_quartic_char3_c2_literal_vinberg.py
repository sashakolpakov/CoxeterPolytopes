"""Recompute the literal Vinberg audits for the two cyclic Fermat folds."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.enumerate_jacobian_kummer_h10_class0007_components import (
    enumerate_components,
)
from scripts.verify_jacobian_kummer_h10_class0007_vinberg import verify


ROOT = Path(__file__).resolve().parents[1]


def paths(dimension: int) -> tuple[Path, Path, Path]:
    stem = ROOT / "data/enumerations" / f"ss_fermat_quartic_char3_c2_h{dimension}"
    return (
        Path(f"{stem}_coxeter_screen.json"),
        Path(f"{stem}_component_audit.json"),
        Path(f"{stem}_vinberg_certificate.json"),
    )


def normalized_components(payload: dict[str, object], dimension: int) -> dict[str, object]:
    payload["schema"] = f"ss_fermat_quartic_char3_c2_h{dimension}_component_audit_v1"
    return payload


def normalized_vinberg(payload: dict[str, object], dimension: int) -> dict[str, object]:
    payload["schema"] = (
        f"ss_fermat_quartic_char3_c2_h{dimension}_literal_vinberg_volume_certificate_v1"
    )
    payload["elliptic_edge_rank"] = dimension - 1
    payload["elliptic_edge_count"] = payload.pop("elliptic_rank_9_edge_count")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dimension", type=int, choices=(10, 11), required=True)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    screen, component_path, vinberg_path = paths(args.dimension)
    components = normalized_components(enumerate_components(screen), args.dimension)
    component_text = json.dumps(components, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if component_path.read_text(encoding="utf-8") != component_text:
            raise SystemExit(f"component archive differs: {component_path}")
    else:
        component_path.write_text(component_text, encoding="utf-8")

    # ``verify`` consumes the archived component path.  In write mode it was
    # just regenerated; in replay mode it was proved byte-for-byte identical.
    vinberg = normalized_vinberg(
        verify(screen_path=screen, component_path=component_path), args.dimension
    )
    vinberg_text = json.dumps(vinberg, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if vinberg_path.read_text(encoding="utf-8") != vinberg_text:
            raise SystemExit(f"Vinberg archive differs: {vinberg_path}")
    else:
        vinberg_path.write_text(vinberg_text, encoding="utf-8")
    print(
        json.dumps(
            {
                "dimension": args.dimension,
                "facets": vinberg["facet_count"],
                "edges": vinberg["elliptic_edge_count"],
                "finite_vertices": vinberg["finite_vertex_count"],
                "ideal_vertices": vinberg["ideal_vertex_count"],
                "decision": vinberg["vinberg_decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
