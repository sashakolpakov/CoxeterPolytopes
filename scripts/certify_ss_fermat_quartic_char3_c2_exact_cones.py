"""GMP exact cone audit for the cyclic Fermat-quartic Coxeter folds."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import sympy as sp

from scripts.certify_ss_fermat_quartic_char3_v4_exact_cones import (
    _fraction_text,
    analyze_rays,
    parse_cdd_rays,
    run_cdd,
)


ROOT = Path(__file__).resolve().parents[1]
COMBINED_SCREEN = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_survivor_screen.json"
)


def paths(dimension: int) -> tuple[Path, Path, Path, Path]:
    stem = ROOT / "data/enumerations" / f"ss_fermat_quartic_char3_c2_h{dimension}"
    return (
        Path(f"{stem}_coxeter_screen.json"),
        Path(f"{stem}_exact_cone.ine"),
        Path(f"{stem}_exact_cone_rays.jsonl.gz"),
        Path(f"{stem}_exact_cone_certificate.json"),
    )


def rebuild_cone(dimension: int) -> dict[str, object]:
    screen_path, _, _, _ = paths(dimension)
    screen = json.loads(screen_path.read_text(encoding="utf-8"))
    combined = json.loads(COMBINED_SCREEN.read_text(encoding="utf-8"))
    matches = [
        cast
        for cast in combined["casts"]
        if int(cast["hyperbolic_dimension"]) == dimension
    ]
    if len(matches) != 1:
        raise ArithmeticError(f"combined screen has no unique H{dimension} cast")
    cast = matches[0]
    if not (
        cast["exact_complete"]
        and cast["unresolved_count"] == 0
        and cast["coxeter_pair_gate"]["all_pairs_coxeter"]
        and screen["status"]["all_facet_pairs_are_coxeter"]
    ):
        raise ArithmeticError("cyclic cast has not passed the complete Coxeter gate")

    positive = sp.Matrix(
        [[sp.Rational(value) for value in row] for row in screen["positive_wall_gram"]]
    )
    raw = -positive
    _, basis_positions = positive.rref()
    if len(basis_positions) != dimension + 1:
        raise ArithmeticError("facet Gram has wrong Lorentzian rank")
    inequality = raw[:, list(basis_positions)]
    basis_gram = raw.extract(basis_positions, basis_positions)
    if inequality.rank() != dimension + 1:
        raise ArithmeticError("cyclic chamber cone is not pointed")

    witness_slacks = sp.Matrix(
        [sp.Rational(item["ample_pairing"]) for item in cast["strict_facet_witnesses"]]
    )
    ample_coordinates = basis_gram.inv() * witness_slacks[list(basis_positions), :]
    if inequality * ample_coordinates != witness_slacks:
        raise ArithmeticError("archived ample slacks do not reconstruct in the facet basis")
    if any(value <= 0 for value in witness_slacks):
        raise ArithmeticError("ample point is not strictly feasible")

    lines = [
        "H-representation",
        "begin",
        f"{inequality.rows} {inequality.cols + 1} rational",
    ]
    lines.extend(
        "0 "
        + " ".join(
            _fraction_text(inequality[row, column])
            for column in range(inequality.cols)
        )
        for row in range(inequality.rows)
    )
    lines.append("end")
    return {
        "dimension": dimension,
        "facet_indices": tuple(range(positive.rows)),
        "basis_positions": tuple(basis_positions),
        "basis_gram": basis_gram,
        "inequality_matrix": inequality,
        "inequality_text": "\n".join(lines) + "\n",
        "minimum_ample_slack": min(witness_slacks),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dimension", type=int, choices=(10, 11), required=True)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument(
        "--cdd-executable", default=shutil.which("cddexec_gmp") or "cddexec_gmp"
    )
    args = parser.parse_args()
    cone = rebuild_cone(args.dimension)
    _, inequality_path, ray_path, certificate_path = paths(args.dimension)
    inequality_path.write_text(str(cone["inequality_text"]), encoding="utf-8")
    print(
        f"prepared H{args.dimension}: {cone['inequality_matrix'].rows} facets -> "
        f"{inequality_path}",
        flush=True,
    )
    if args.prepare_only:
        return
    cdd_text = run_cdd(
        str(cone["inequality_text"]), args.cdd_executable, args.timeout
    )
    rays = parse_cdd_rays(cdd_text, args.dimension + 1)
    summary, archive = analyze_rays(
        cone,
        rays,
        schema_prefix="ss_fermat_quartic_char3_c2",
    )
    exact = summary["exact_double_description"]
    literal_path = (
        ROOT
        / "data/enumerations"
        / f"ss_fermat_quartic_char3_c2_h{args.dimension}_vinberg_certificate.json"
    )
    literal = json.loads(literal_path.read_text(encoding="utf-8"))
    literal_match = (
        exact["finite_vertex_count"] == literal["finite_vertex_count"]
        and exact["ideal_vertex_count"] == literal["ideal_vertex_count"]
    )
    summary = {
        "schema": f"ss_fermat_quartic_char3_c2_h{args.dimension}_exact_cone_certificate_v1",
        "coxeter_screen": str(paths(args.dimension)[0].relative_to(ROOT)),
        "literal_vinberg_certificate": str(literal_path.relative_to(ROOT)),
        "method": "cddlib GMP rational double description with exact integer replay",
        "dimension": args.dimension,
        "facet_count": int(cone["inequality_matrix"].rows),
        "coordinate_basis_facet_positions": [
            index + 1 for index in cone["basis_positions"]
        ],
        "inequality_scale": exact["inequality_integral_scale"],
        "basis_gram_scale": exact["basis_gram_integral_scale"],
        "maximum_primitive_coordinate": exact[
            "primitive_ray_coordinate_absolute_maximum"
        ],
        "ray_count": exact["extreme_ray_count"],
        "finite_ray_count": exact["finite_vertex_count"],
        "ideal_ray_count": exact["ideal_vertex_count"],
        "hyperideal_ray_count": exact["hyperideal_vertex_count"],
        "active_facet_count_histogram": exact["active_facet_count_histogram"],
        "literal_vinberg_endpoint_counts_match": literal_match,
        "finite_volume_certified": (
            exact["hyperideal_vertex_count"] == 0 and literal_match
        ),
    }
    ray_path.write_bytes(archive)
    certificate_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "finite_volume_certified": summary["finite_volume_certified"],
                "finite_rays": summary["finite_ray_count"],
                "ideal_rays": summary["ideal_ray_count"],
                "hyperideal_rays": summary["hyperideal_ray_count"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
