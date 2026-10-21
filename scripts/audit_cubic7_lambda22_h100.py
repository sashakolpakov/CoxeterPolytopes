#!/usr/bin/env python3
"""Freeze the exact ``<-theta> + Lambda22`` Vinberg prefix through H=100."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import groupby
import json
from pathlib import Path
import shutil
import subprocess

from sympy import Matrix, kronecker_product

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_lambda22 import (
    GRAM_FILE_SHA256,
    GRAM_PATH,
    GRAM_STREAM_SHA256,
    LAMBDA22_UNIT_SEED,
    PARI_QFMINIM_CERTIFICATES,
    SOURCE_HTML_SHA256,
    SOURCE_URL,
    TRACE_GRAM,
    Lambda22GramType,
    lambda22_basis_audit,
    lambda22_control_levels,
    lambda22_gram,
    lambda22_gram_file_sha256,
    lambda22_gram_stream_sha256,
    lambda22_level_audit,
    lambda22_norm_gate,
    lambda22_pairing_radical_basis,
    lambda22_pairing_radical_gram,
    lambda22_rank3_gram_types,
    lambda22_smith_invariants,
)


HEIGHT_BOUND = 100
DEFAULT_OUTPUT = Path("data/enumerations/cubic7_lambda22_h100.json")
EXPECTED_CONTROL_COUNTS = {1: 5, 2: 16, 4: 5, 3: 1, 6: 2, 12: 1}


def _value(value: Cubic7) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(place) for place in range(3)],
    }


def _gram_type(value: Lambda22GramType) -> dict[str, object]:
    return {
        "tuple_U_V_W_uv_uw_vw_det": [
            value.u_square,
            value.v_square,
            value.w_square,
            value.u_dot_v,
            value.u_dot_w,
            value.v_dot_w,
            value.determinant,
        ],
        "rank": value.rank,
    }


def _gp_matrix(matrix: tuple[tuple[int, ...], ...]) -> str:
    return "[" + ";".join(
        ",".join(str(value) for value in row) for row in matrix
    ) + "]"


def _pari_qfminim(
    gp: str, matrix: tuple[tuple[int, ...], ...]
) -> dict[str, object]:
    program = (
        f"r=qfminim({_gp_matrix(matrix)});"
        'print(r[1]);print(r[2]);print(matsize(r[3])[2]);print(r[3]);\n'
    )
    result = subprocess.run(
        [gp, "-q", "-s", "2000000000"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=300,
    )
    lines = result.stdout.splitlines(keepends=True)
    if len(lines) != 4:
        raise AssertionError("PARI qfminim emitted an unexpected record")
    vector_stream = lines[3].encode("utf-8")
    return {
        "kissing_number": int(lines[0]),
        "minimum": int(lines[1]),
        "representatives_mod_sign": int(lines[2]),
        "vector_stream_bytes": len(vector_stream),
        "vector_stream_sha256": sha256(vector_stream).hexdigest(),
    }


def _pari_replay() -> dict[str, object]:
    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP is required to regenerate the Lambda22 artifact")
    version_result = subprocess.run(
        [gp, "--version"], text=True, capture_output=True, check=True, timeout=30
    )
    version_lines = (version_result.stdout + version_result.stderr).splitlines()
    version = version_lines[0] if version_lines else "PARI/GP (version output unavailable)"
    trace = Matrix(TRACE_GRAM)
    rows = []
    for modulus in (1, 2, 3, 6):
        ordinary_gram = lambda22_pairing_radical_gram(modulus)
        trace_matrix = kronecker_product(trace, Matrix(ordinary_gram))
        trace_gram = tuple(
            tuple(int(trace_matrix[i, j]) for j in range(trace_matrix.cols))
            for i in range(trace_matrix.rows)
        )
        ordinary = _pari_qfminim(gp, ordinary_gram)
        trace_audit = _pari_qfminim(gp, trace_gram)
        expected = PARI_QFMINIM_CERTIFICATES[modulus]
        if ordinary != {
            "kissing_number": expected["ordinary_kissing_number"],
            "minimum": expected["ordinary_minimum"],
            "representatives_mod_sign": expected[
                "ordinary_representatives_mod_sign"
            ],
            "vector_stream_bytes": expected["ordinary_vector_stream_bytes"],
            "vector_stream_sha256": expected["ordinary_vector_stream_sha256"],
        }:
            raise AssertionError("an ordinary PARI Lambda22 certificate changed")
        if trace_audit != {
            "kissing_number": expected["trace_kissing_number"],
            "minimum": expected["trace_minimum"],
            "representatives_mod_sign": expected[
                "trace_representatives_mod_sign"
            ],
            "vector_stream_bytes": expected["trace_vector_stream_bytes"],
            "vector_stream_sha256": expected["trace_vector_stream_sha256"],
        }:
            raise AssertionError("a trace PARI Lambda22 certificate changed")
        basis = Matrix(lambda22_pairing_radical_basis(modulus))
        rows.append(
            {
                "pairing_modulus": modulus,
                "definition": f"N_{modulus}={{y in Lambda22 : (y,Lambda22) subset {modulus}Z}}",
                "index": abs(int(basis.det())),
                "basis_determinant": int(basis.det()),
                "ordinary_rank": 22,
                "trace_rank": 66,
                "ordinary": ordinary,
                "trace_over_OK": trace_audit,
            }
        )
    return {
        "engine": version,
        "command_shape": "qfminim(M) with signed count r[1], minimum r[2], and representative matrix r[3]",
        "coordinate_order_for_trace": "three 22-coordinate blocks u,v,w; Gram=TraceGram tensor ordinaryGram",
        "shell_orbit_assumption": False,
        "radical_lattices": rows,
    }


def build_audit() -> dict[str, object]:
    if not lambda22_basis_audit():
        raise AssertionError("the vendored Lambda22 Gram audit failed")
    if lambda22_smith_invariants() != (1,) * 20 + (2, 6):
        raise AssertionError("the Lambda22 Smith invariants changed")
    seed = LAMBDA22_UNIT_SEED
    if not seed.admissible:
        raise AssertionError("a=theta lost its distinguished Lorentzian signs")
    levels = lambda22_control_levels(HEIGHT_BOUND)
    norm_counts = Counter(int(level.norm.a) for level in levels)
    if len(levels) != 30 or dict(norm_counts) != EXPECTED_CONTROL_COUNTS:
        raise AssertionError("the Lambda22 H100 control census changed")

    audits = tuple(lambda22_level_audit(level) for level in levels)
    if any(audit.gated_gram_types for audit in audits):
        raise AssertionError("the Lambda22 H100 prefix acquired a possible wall")
    height_zero = []
    for norm in seed.canonical_mirror_norms:
        gate = lambda22_norm_gate(norm)
        types = lambda22_rank3_gram_types(
            norm, gate.coefficient_lattice_minimum, gate.pairing_modulus
        )
        if types:
            raise AssertionError("Lambda22 unexpectedly acquired a height-zero wall")
        height_zero.append(
            {
                "norm": _value(norm),
                "pairing_modulus": gate.pairing_modulus,
                "target_trace": int(norm.trace()),
                "coefficient_trace_lattice_minimum": PARI_QFMINIM_CERTIFICATES[
                    gate.pairing_modulus
                ]["trace_minimum"],
                "necessary_gram_type_count": 0,
            }
        )

    level_rows = []
    for audit in audits:
        level = audit.level
        gate = lambda22_norm_gate(level.norm)
        level_rows.append(
            {
                "height": _value(level.height),
                "norm": _value(level.norm),
                "control": _value(level.control),
                "spatial_square_target": _value(level.target),
                "target_trace": int(level.target.trace()),
                "gate": {
                    "pairing_modulus": gate.pairing_modulus,
                    "coefficient_lattice_index": gate.coefficient_lattice_index,
                    "coefficient_lattice_minimum": gate.coefficient_lattice_minimum,
                },
                "ambient_minimum_four_gram_types": [
                    _gram_type(value) for value in audit.ambient_gram_types
                ],
                "pairing_congruence_gram_types": [
                    _gram_type(value) for value in audit.congruence_gram_types
                ],
                "fully_gated_gram_types": [
                    _gram_type(value) for value in audit.gated_gram_types
                ],
                "candidate_root_count": 0,
            }
        )

    height_packets = []
    for height, iterator in groupby(audits, key=lambda audit: audit.level.height):
        group = tuple(iterator)
        height_packets.append(
            {
                "height": _value(height),
                "control_count": len(group),
                "norms": [_value(audit.level.norm) for audit in group],
                "fully_gated_gram_type_count": sum(
                    len(audit.gated_gram_types) for audit in group
                ),
                "accepted_wall_count": 0,
            }
        )

    rank_counts = Counter(
        gram_type.rank
        for audit in audits
        for gram_type in audit.ambient_gram_types
    )
    congruence_rank_counts = Counter(
        gram_type.rank
        for audit in audits
        for gram_type in audit.congruence_gram_types
    )
    pari = _pari_replay()
    last = levels[-1]
    return {
        "schema": "coxeter.cubic7.lambda22_h100.v1",
        "field": "Q(theta), theta^3+theta^2-2theta-1=0",
        "hyperbolic_dimension": 22,
        "height_bound": HEIGHT_BOUND,
        "spatial_lattice": {
            "name": "Lambda22 laminated lattice",
            "rank": 22,
            "rootless": True,
            "catalogue_source": SOURCE_URL,
            "retrieved_utc_date": "2026-08-12",
            "source_html_sha256": SOURCE_HTML_SHA256,
            "vendored_gram_path": str(GRAM_PATH.relative_to(Path.cwd())),
            "vendored_gram_file_sha256": lambda22_gram_file_sha256(),
            "normalized_full_gram_stream_sha256": lambda22_gram_stream_sha256(),
            "catalogue_and_exact_matrix_invariants": {
                "determinant": 12,
                "minimum": 4,
                "kissing_number": 49_896,
                "smith_invariants": [1] * 20 + [2, 6],
                "even": True,
                "integral": True,
            },
            "source_pins_match": (
                lambda22_gram_file_sha256() == GRAM_FILE_SHA256
                and lambda22_gram_stream_sha256() == GRAM_STREAM_SHA256
            ),
        },
        "pari_qfminim_replay": pari,
        "seed": {
            "seed_id": seed.seed_id,
            "timelike_coefficient": _value(seed.timelike_coefficient),
            "support": seed.support,
            "admissible_at_distinguished_place": True,
            "canonical_mirror_norms": [
                _value(norm) for norm in seed.canonical_mirror_norms
            ],
            "complete_local_norm_proof": seed.local_scope,
            "norm_gates": [
                {
                    "norm": _value(norm),
                    "pairing_modulus": lambda22_norm_gate(norm).pairing_modulus,
                    "coefficient_lattice_minimum": lambda22_norm_gate(
                        norm
                    ).coefficient_lattice_minimum,
                    "coefficient_lattice_index": lambda22_norm_gate(
                        norm
                    ).coefficient_lattice_index,
                }
                for norm in seed.canonical_mirror_norms
            ],
        },
        "height_zero": {
            "norm_class_audits": height_zero,
            "reflecting_wall_count": 0,
            "statement": (
                "For every complete norm class, the exact trace target is below "
                "the PARI-certified trace minimum of the required pairing-radical lattice."
            ),
        },
        "positive_prefix": {
            "control_count": len(levels),
            "distinct_height_count": len({level.height for level in levels}),
            "control_counts_by_norm": [
                {
                    "norm": _value(norm),
                    "count": norm_counts[int(norm.a)],
                }
                for norm in seed.canonical_mirror_norms
            ],
            "control_levels": level_rows,
            "height_packets": height_packets,
            "ambient_gram_type_count": sum(rank_counts.values()),
            "ambient_gram_type_rank_histogram": {
                str(rank): rank_counts[rank] for rank in sorted(rank_counts)
            },
            "pairing_congruence_gram_type_count": sum(
                congruence_rank_counts.values()
            ),
            "pairing_congruence_rank_histogram": {
                str(rank): congruence_rank_counts[rank]
                for rank in sorted(congruence_rank_counts)
            },
            "fully_gated_gram_type_count": 0,
            "accepted_positive_wall_count": 0,
        },
        "frontier": {
            "all_controls_at_height_at_most_100_exhausted": True,
            "last_control_height": _value(last.height),
            "last_control_norm": _value(last.norm),
            "last_control": _value(last.control),
            "last_spatial_square_target": _value(last.target),
        },
        "cone_status": {
            "height_zero_wall_count": 0,
            "positive_wall_count_through_H100": 0,
            "polyhedral_chamber_closed": False,
            "finite_volume_test_applicable": False,
            "compactness_test_applicable": False,
            "hyperideal_extreme_ray_test_applicable": False,
            "reason": (
                "No reflecting wall occurs in the certified prefix, so there is no "
                "positive-dimensional wall cone on which to run an extreme-ray closure test."
            ),
        },
        "scope": (
            "Exact necessary-and-realized Vinberg prefix only through H<=100 for "
            "<-theta> + Lambda22 at the distinguished embedding. It proves no "
            "reflecting wall in this bounded prefix. It makes no finite-volume, "
            "compactness, reflectivity, or global nonreflectivity conclusion."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
