"""Independent GMP cddlib cone-ray audit for a Fermat-char-3 volume input."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import argparse
import json
import math
from pathlib import Path
import subprocess

import numpy as np
import sympy as sp


def fraction_text(value: Fraction | sp.Rational | int) -> str:
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def primitive_ray(tokens: list[str]) -> tuple[int, ...]:
    values = tuple(Fraction(token) for token in tokens)
    denominator = math.lcm(*(value.denominator for value in values))
    integers = tuple(
        value.numerator * (denominator // value.denominator) for value in values
    )
    divisor = math.gcd(*(abs(value) for value in integers))
    if not divisor:
        raise ArithmeticError("cddlib emitted the zero ray")
    # Preserve cddlib's feasible orientation; unlike an unoriented projective
    # line, a ray and its negative are not interchangeable in a pointed cone.
    return tuple(value // divisor for value in integers)


def integral_matrix(matrix: sp.Matrix) -> tuple[np.ndarray, int]:
    scale = math.lcm(*(int(sp.Rational(value).q) for value in matrix))
    return (
        np.asarray(
            [[int(scale * sp.Rational(value)) for value in row] for row in matrix.tolist()],
            dtype=np.int64,
        ),
        scale,
    )


def audit(screen_path: Path, vinberg_path: Path | None) -> dict[str, object]:
    screen = json.loads(screen_path.read_text(encoding="utf-8"))
    positive = sp.Matrix(
        [[sp.Rational(value) for value in row] for row in screen["positive_wall_gram"]]
    )
    raw = -positive
    dimension = int(screen["group_class"]["hyperbolic_dimension"])
    _, pivots = positive.rref()
    basis = tuple(pivots)
    if len(basis) != dimension + 1:
        raise ArithmeticError("facet Gram has wrong Lorentzian rank")
    inequalities = raw[:, list(basis)]
    if inequalities.rank() != dimension + 1:
        raise ArithmeticError("cone inequalities have lineality")
    basis_gram = raw.extract(basis, basis)

    lines = [
        "H-representation",
        "begin",
        f"{inequalities.rows} {inequalities.cols + 1} rational",
    ]
    lines.extend(
        "0 " + " ".join(fraction_text(value) for value in inequalities.row(row))
        for row in range(inequalities.rows)
    )
    lines.append("end")
    completed = subprocess.run(
        ["cddexec_gmp", "--rep"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        timeout=1800,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(completed.stdout[-2000:] + completed.stderr[-2000:])
    output = completed.stdout.splitlines()
    begin = output.index("begin")
    count_text, width_text, arithmetic = output[begin + 1].split()
    count, width = int(count_text), int(width_text)
    if width != dimension + 2 or arithmetic != "rational":
        raise ArithmeticError("unexpected cddlib output header")
    rows = output[begin + 2 : begin + 2 + count]
    if len(rows) != count or output[begin + 2 + count].strip() != "end":
        raise ArithmeticError("truncated cddlib output")
    rays = []
    for row in rows:
        tokens = row.split()
        if len(tokens) != width or Fraction(tokens[0]) != 0:
            raise ArithmeticError("homogeneous cone emitted a non-ray row")
        rays.append(primitive_ray(tokens[1:]))
    if len(set(rays)) != len(rays):
        raise ArithmeticError("cddlib emitted duplicate projective rays")

    ray_matrix = np.asarray(rays, dtype=np.int64)
    inequality_integer, inequality_scale = integral_matrix(inequalities)
    basis_gram_integer, basis_gram_scale = integral_matrix(basis_gram)
    coordinate_bound = int(np.max(np.abs(ray_matrix), initial=0))
    arithmetic_bound = max(
        inequalities.cols * coordinate_bound * int(np.max(np.abs(inequality_integer))),
        inequalities.cols**2
        * coordinate_bound**2
        * int(np.max(np.abs(basis_gram_integer))),
    )
    if arithmetic_bound >= np.iinfo(np.int64).max:
        raise OverflowError("int64 exact replay bound exceeded")
    norms = np.einsum("bi,ij,bj->b", ray_matrix, basis_gram_integer, ray_matrix)
    active_histogram: Counter[int] = Counter()
    finite_count = ideal_count = hyperideal_count = 0
    for start in range(0, len(rays), 1024):
        slacks = ray_matrix[start : start + 1024] @ inequality_integer.T
        if np.any(slacks < 0):
            raise ArithmeticError("cddlib ray violates an exact cone inequality")
        for local, row in enumerate(slacks):
            active_histogram[int(np.count_nonzero(row == 0))] += 1
            norm = int(norms[start + local])
            if norm > 0:
                finite_count += 1
            elif norm == 0:
                ideal_count += 1
            else:
                hyperideal_count += 1

    vinberg_match = None
    if vinberg_path is not None:
        vinberg = json.loads(vinberg_path.read_text(encoding="utf-8"))
        vinberg_match = (
            finite_count == int(vinberg["finite_vertex_count"])
            and ideal_count == int(vinberg["ideal_vertex_count"])
        )
        if not vinberg_match:
            raise ArithmeticError("cdd cone endpoints differ from literal Vinberg endpoints")
    return {
        "method": "cddlib GMP rational double description with exact integer replay",
        "dimension": dimension,
        "facet_count": positive.rows,
        "coordinate_basis_facet_positions": [index + 1 for index in basis],
        "inequality_scale": inequality_scale,
        "basis_gram_scale": basis_gram_scale,
        "ray_count": len(rays),
        "finite_ray_count": finite_count,
        "ideal_ray_count": ideal_count,
        "hyperideal_ray_count": hyperideal_count,
        "active_facet_count_histogram": {
            str(key): value for key, value in sorted(active_histogram.items())
        },
        "maximum_primitive_coordinate": coordinate_bound,
        "int64_arithmetic_bound": arithmetic_bound,
        "literal_vinberg_endpoint_counts_match": vinberg_match,
        "finite_volume_certified": hyperideal_count == 0 and vinberg_match is not False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("screen", type=Path)
    parser.add_argument("--vinberg", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.screen, args.vinberg)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
