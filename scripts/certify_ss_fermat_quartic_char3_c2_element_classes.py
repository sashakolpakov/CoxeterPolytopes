"""Reconstruct the full Fermat chamber group and screen high-rank C2 classes."""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction
from functools import reduce
from math import gcd
import json
from pathlib import Path
import re
import subprocess

import numpy as np
import sympy as sp

from scripts.search_ss_fermat_quartic_char3_folds import parse_source


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/external/kondo_shimada_fermat_quartic_char3/compdataAutFQChar3.txt"


def extract(text: str, name: str) -> str:
    marker = text.index(name)
    start = text.index("[", marker)
    depth = 0
    for end in range(start, len(text)):
        if text[end] == "[":
            depth += 1
        elif text[end] == "]":
            depth -= 1
            if depth == 0:
                return text[start : end + 1].replace("\\\n", "")
    raise ValueError(name)


# GF(9) = F_3[i], i^2=2, encoded as a+3*b for a+b*i.
def gf(a: int, b: int = 0) -> int:
    return (a % 3) + 3 * (b % 3)


def add(left: int, right: int) -> int:
    return gf(left % 3 + right % 3, left // 3 + right // 3)


def neg(value: int) -> int:
    return gf(-(value % 3), -(value // 3))


def mul(left: int, right: int) -> int:
    a, b = left % 3, left // 3
    c, d = right % 3, right // 3
    return gf(a * c + 2 * b * d, a * d + b * c)


def power(value: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = mul(result, value)
        value = mul(value, value)
        exponent >>= 1
    return result


def inverse(value: int) -> int:
    if value == 0:
        raise ZeroDivisionError
    return power(value, 7)


def matmul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    out = np.zeros((left.shape[0], right.shape[1]), dtype=np.int8)
    for i in range(left.shape[0]):
        for j in range(right.shape[1]):
            value = 0
            for k in range(left.shape[1]):
                value = add(value, mul(int(left[i, k]), int(right[k, j])))
            out[i, j] = value
    return out


def matrix_inverse(matrix: np.ndarray) -> np.ndarray:
    size = len(matrix)
    augmented = np.concatenate((matrix.copy(), np.eye(size, dtype=np.int8)), axis=1)
    for column in range(size):
        pivot = next(row for row in range(column, size) if augmented[row, column])
        if pivot != column:
            augmented[[column, pivot]] = augmented[[pivot, column]]
        inv = inverse(int(augmented[column, column]))
        for j in range(2 * size):
            augmented[column, j] = mul(int(augmented[column, j]), inv)
        for row in range(size):
            if row == column:
                continue
            factor = int(augmented[row, column])
            if factor:
                for j in range(2 * size):
                    augmented[row, j] = add(
                        int(augmented[row, j]),
                        neg(mul(factor, int(augmented[column, j]))),
                    )
    return augmented[:, size:]


def canonical_subspace(rows: np.ndarray) -> tuple[int, ...]:
    matrix = rows.copy().astype(np.int8)
    rank = 0
    for column in range(matrix.shape[1]):
        pivot = next((row for row in range(rank, len(matrix)) if matrix[row, column]), None)
        if pivot is None:
            continue
        if pivot != rank:
            matrix[[rank, pivot]] = matrix[[pivot, rank]]
        inv = inverse(int(matrix[rank, column]))
        for j in range(matrix.shape[1]):
            matrix[rank, j] = mul(int(matrix[rank, j]), inv)
        for row in range(len(matrix)):
            if row == rank:
                continue
            factor = int(matrix[row, column])
            if factor:
                for j in range(matrix.shape[1]):
                    matrix[row, j] = add(
                        int(matrix[row, j]),
                        neg(mul(factor, int(matrix[rank, j]))),
                    )
        rank += 1
        if rank == len(matrix):
            break
    if rank != 2:
        raise ArithmeticError("line parametrization lost rank")
    return tuple(map(int, matrix.reshape(-1)))


def coefficient_to_gf(value: sp.Expr) -> int:
    radical = sp.sqrt(2)
    value = sp.expand(value)
    b = int(value.coeff(radical))
    a = int(sp.expand(value - b * radical))
    return gf(a, b)


def parse_lines(text: str) -> tuple[np.ndarray, ...]:
    payload = extract(text, "FQparalines")
    entries = re.findall(r"\[([^\[\]]+)\]", payload)
    if len(entries) != 112:
        raise ArithmeticError(f"expected 112 line parametrizations, got {len(entries)}")
    U, V = sp.symbols("U V")
    lines = []
    for entry in entries:
        expressions = [sp.sympify(piece.strip()) for piece in entry.split(",")]
        if len(expressions) != 4:
            raise ArithmeticError("line does not have four coordinates")
        rows = np.zeros((2, 4), dtype=np.int8)
        for column, expression in enumerate(expressions):
            expression = sp.expand(expression)
            rows[0, column] = coefficient_to_gf(expression.coeff(U))
            rows[1, column] = coefficient_to_gf(expression.coeff(V))
            if sp.expand(expression - expression.coeff(U) * U - expression.coeff(V) * V) != 0:
                raise ArithmeticError("nonlinear line parametrization")
        lines.append(rows)
    return tuple(lines)


def generator_matrices() -> tuple[np.ndarray, np.ndarray]:
    # GAP's standard GU(4,3) generators.  Its primitive z satisfies
    # z^2-z-1=0; in the source coordinates we choose z=2+i.
    z = gf(2, 1)
    first_standard = np.diag([z, 1, 1, power(z, 5)]).astype(np.int8)
    second_standard = np.asarray(
        [
            [1, 0, power(z, 2), 0],
            [1, 0, 0, 0],
            [0, power(z, 6), 0, 1],
            [0, power(z, 2), 0, 0],
        ],
        dtype=np.int8,
    )
    # GAP's standard group preserves the anti-diagonal Hermitian form J.
    # Rows of C are a J-orthonormal basis, so C*J*C^(3T)=I.  Conjugating
    # transports the generators to the Fermat form sum x_i^4.
    C = np.asarray(
        [
            [1, 0, 0, 2],
            [z, 0, 0, z],
            [0, 1, 2, 0],
            [0, z, z, 0],
        ],
        dtype=np.int8,
    )
    C_inverse = matrix_inverse(C)
    return (
        matmul(matmul(C, first_standard), C_inverse),
        matmul(matmul(C, second_standard), C_inverse),
    )


def line_permutation(lines: tuple[np.ndarray, ...], generator: np.ndarray) -> tuple[int, ...]:
    lookup = {canonical_subspace(line): index for index, line in enumerate(lines)}
    images = []
    for line in lines:
        image = canonical_subspace(matmul(line, generator))
        if image not in lookup:
            raise ArithmeticError(f"GU generator does not preserve archived lines: {image}")
        images.append(lookup[image])
    if len(set(images)) != 112:
        raise ArithmeticError("line action is not a permutation")
    return tuple(images)


def primitive_positive(row: np.ndarray, ample: np.ndarray) -> tuple[int, ...]:
    values = tuple(map(int, row))
    divisor = reduce(gcd, (abs(value) for value in values if value), 0)
    if not divisor:
        raise ArithmeticError("zero restricted wall")
    primitive = tuple(value // divisor for value in values)
    if sum(left * right for left, right in zip(primitive, ample, strict=True)) < 0:
        primitive = tuple(-value for value in primitive)
    return primitive


def coxeter_label(product: int, left_norm: int, right_norm: int) -> str | None:
    if product < 0:
        return None
    if product == 0:
        return "2"
    ratio = Fraction(product * product, left_norm * right_norm)
    return {
        Fraction(1, 4): "3",
        Fraction(1, 2): "4",
        Fraction(3, 4): "6",
        Fraction(1): "infinity",
    }.get(ratio, "ultraparallel" if ratio > 1 else None)


def screen_involution(
    matrix: np.ndarray,
    walls: np.ndarray,
    wall_types: tuple[str, ...],
    gram: np.ndarray,
    h0: np.ndarray,
) -> dict[str, object]:
    """Run the exact strict-facet Coxeter gate on one involution-fixed cast."""

    identity = np.eye(len(gram), dtype=np.int64)
    if not np.array_equal(matrix @ matrix, identity):
        raise ArithmeticError("high-fixed-rank representative is not an involution")
    if not np.array_equal(matrix @ gram @ matrix.T, gram):
        raise ArithmeticError("representative is not a lattice isometry")
    if not np.array_equal(h0 @ matrix, h0):
        raise ArithmeticError("representative does not fix the chamber interior")

    wall_lookup = {tuple(map(int, wall)): index for index, wall in enumerate(walls)}
    images = walls @ matrix.T
    permutation = tuple(wall_lookup[tuple(map(int, image))] for image in images)
    if any(permutation[permutation[index]] != index for index in range(len(walls))):
        raise ArithmeticError("wall action is not involutory")
    if any(wall_types[index] != wall_types[image] for index, image in enumerate(permutation)):
        raise ArithmeticError("wall action mixes published wall types")

    unseen = set(range(len(walls)))
    orbits: list[tuple[int, ...]] = []
    while unseen:
        index = min(unseen)
        image = permutation[index]
        orbit = (index,) if image == index else (index, image)
        unseen.difference_update(orbit)
        orbits.append(orbit)

    ray_sources: dict[tuple[int, ...], list[int]] = {}
    for orbit_index, orbit in enumerate(orbits):
        ray = primitive_positive(np.sum(walls[list(orbit)], axis=0), h0)
        ray_sources.setdefault(ray, []).append(orbit_index)
    rays = tuple(ray_sources)
    restricted = np.asarray(rays, dtype=np.int64)

    gram_sympy = sp.Matrix(gram.tolist())
    inverse = gram_sympy.inv()
    denominator = int(sp.ilcm(*[sp.Rational(value).q for value in inverse]))
    metric_numerator = np.asarray(inverse * denominator, dtype=np.int64)
    metric_rows = restricted @ metric_numerator
    restricted_gram = metric_rows @ restricted.T
    norms = np.diag(restricted_gram)
    ample_pairings = restricted @ h0

    fixed_rank = len(gram) - sp.Matrix((matrix - identity).tolist()).rank()
    restricted_rank = sp.Matrix(restricted.tolist()).rank()
    if restricted_rank != fixed_rank:
        raise ArithmeticError(
            f"restricted rank {restricted_rank} differs from fixed rank {fixed_rank}"
        )

    strict_facets: list[int] = []
    strict_minima: list[int] = []
    for index, norm_value in enumerate(norms):
        norm = int(norm_value)
        if norm >= 0:
            continue
        margins = int(ample_pairings[index]) * restricted_gram[index] - norm * ample_pairings
        if int(margins[index]) != 0:
            raise ArithmeticError("facet witness has nonzero self margin")
        minimum = int(min(np.min(margins[:index], initial=np.iinfo(np.int64).max), np.min(margins[index + 1 :], initial=np.iinfo(np.int64).max)))
        if minimum > 0:
            strict_facets.append(index)
            strict_minima.append(minimum)

    label_counts: Counter[str] = Counter()
    bad_ratio_counts: Counter[str] = Counter()
    bad_examples: list[dict[str, object]] = []
    for left_position, left in enumerate(strict_facets):
        for right in strict_facets[left_position + 1 :]:
            product = int(restricted_gram[left, right])
            left_norm, right_norm = int(norms[left]), int(norms[right])
            label = coxeter_label(product, left_norm, right_norm)
            if label is not None:
                label_counts[label] += 1
                continue
            ratio = Fraction(product * product, left_norm * right_norm)
            ratio_text = str(ratio)
            bad_ratio_counts[ratio_text] += 1
            if len(bad_examples) < 12:
                bad_examples.append(
                    {
                        "ray_indices": [left + 1, right + 1],
                        "normal_product": str(Fraction(product, denominator)),
                        "normalized_squared_product": ratio_text,
                        "normal_squares": [
                            str(Fraction(left_norm, denominator)),
                            str(Fraction(right_norm, denominator)),
                        ],
                        "source_wall_indices": [
                            [
                                wall + 1
                                for orbit_index in ray_sources[rays[left]]
                                for wall in orbits[orbit_index]
                            ],
                            [
                                wall + 1
                                for orbit_index in ray_sources[rays[right]]
                                for wall in orbits[orbit_index]
                            ],
                        ],
                    }
                )

    orbit_histogram = Counter(
        f"{wall_types[orbit[0]]}:size_{len(orbit)}" for orbit in orbits
    )
    return {
        "fixed_rank": int(fixed_rank),
        "hyperbolic_dimension": int(fixed_rank - 1),
        "wall_orbit_count": len(orbits),
        "wall_orbit_histogram": dict(sorted(orbit_histogram.items())),
        "distinct_restricted_ray_count": len(rays),
        "coincident_restricted_ray_count": len(orbits) - len(rays),
        "restricted_rank": int(restricted_rank),
        "metric_denominator": denominator,
        "strict_facet_count": len(strict_facets),
        "strict_margin_minimum": min(strict_minima) if strict_minima else None,
        "coxeter_label_counts": dict(sorted(label_counts.items())),
        "bad_pair_count": sum(bad_ratio_counts.values()),
        "bad_ratio_counts": dict(sorted(bad_ratio_counts.items())),
        "bad_examples": bad_examples,
        "decision": "not_coxeter" if bad_ratio_counts else "survives_strict_facet_gate",
    }


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    lines = parse_lines(text)
    ns_basis = ast.literal_eval(extract(text, "FQNSbasis"))
    line_classes = np.asarray(ast.literal_eval(extract(text, "FQlineclasses")), dtype=np.int64)
    gram = np.asarray(ast.literal_eval(extract(text, "FQGramS")), dtype=np.int64)
    h0 = np.asarray(ast.literal_eval(extract(text, "FQh0")), dtype=np.int64)
    print("NS basis", ns_basis)
    print("line classes", line_classes.shape)
    permutations = []
    lattice_generators = []
    for number, generator in enumerate(generator_matrices(), start=1):
        permutation = line_permutation(lines, generator)
        matrix = np.asarray(
            [line_classes[permutation[line_index - 1]] for line_index in ns_basis],
            dtype=np.int64,
        )
        print("GEN", number, "line cycle fixed", sum(i == j for i, j in enumerate(permutation)))
        print("isometry", np.array_equal(matrix @ gram @ matrix.T, gram), "h0", np.array_equal(h0 @ matrix, h0), "trace", int(np.trace(matrix)))
        permutations.append(permutation)
        lattice_generators.append(matrix)

    frobenius = np.asarray(ast.literal_eval(extract(text, "FQAF")), dtype=np.int64)
    class_lookup = {tuple(map(int, row)): index for index, row in enumerate(line_classes)}
    frobenius_permutation = tuple(
        class_lookup[tuple(map(int, row @ frobenius))] for row in line_classes
    )
    assert np.array_equal(frobenius @ gram @ frobenius.T, gram)
    assert np.array_equal(h0 @ frobenius, h0)
    permutations.append(frobenius_permutation)
    lattice_generators.append(frobenius)

    gap_permutations = ",".join(
        "PermList([" + ",".join(str(value + 1) for value in permutation) + "])"
        for permutation in permutations
    )
    gap_matrices = ",".join(str(matrix.tolist()).replace(" ", "") for matrix in lattice_generators)
    gap_program = f"""
pgens := [{gap_permutations}];;
G := Group(pgens);;
mgens := [{gap_matrices}];;
MG := Group(mgens);;
hom := GroupHomomorphismByImages(G,MG,pgens,mgens);;
Print("HEADER\\t",Size(G),"\\t",StructureDescription(G),"\\t",Size(Image(hom)),"\\n");
C := ConjugacyClasses(G);;
Print("CLASSES\\t",Length(C),"\\n");
i := 0;;
for c in C do
  i := i+1;;
  g := Representative(c);;
  m := Image(hom,g);;
  fixed := 22-RankMat(m-IdentityMat(22));;
  Print("CLASS\\t",i,"\\t",Size(c),"\\t",Order(g),"\\t",fixed,"\\t",TraceMat(m),"\\n");
  if fixed>=11 and Order(g)>1 then
    Print("HIGHBEGIN\\t",i,"\\n",m,"\\nHIGHEND\\t",i,"\\n");
  fi;;
od;;
QUIT;
"""
    completed = subprocess.run(
        ["gap", "-q"], input=gap_program, text=True, capture_output=True, check=True
    )
    header = re.search(r"^HEADER\t([^\n]+)$", completed.stdout, re.MULTILINE)
    if header is None:
        raise ArithmeticError("missing GAP census header")
    print("HEADER", header.group(1))
    classes = {
        int(match.group(1)): {
            "class_size": int(match.group(2)),
            "order": int(match.group(3)),
            "fixed_rank": int(match.group(4)),
            "trace": int(match.group(5)),
        }
        for match in re.finditer(
            r"^CLASS\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(-?\d+)$",
            completed.stdout,
            re.MULTILINE,
        )
    }
    high = []
    for match in re.finditer(
        r"HIGHBEGIN\t(\d+)\n(.*?)\nHIGHEND\t\1",
        completed.stdout,
        re.DOTALL,
    ):
        class_index = int(match.group(1))
        matrix = np.asarray(ast.literal_eval(match.group(2)), dtype=np.int64)
        high.append((class_index, matrix))
    if len(classes) != 61 or len(high) != 5:
        raise ArithmeticError(
            f"expected 61 element classes and five high casts, got {len(classes)} and {len(high)}"
        )

    published = parse_source(SOURCE)
    walls = np.asarray(published["walls"], dtype=np.int64)
    wall_types = published["wall_types"]
    results = []
    # Highest hyperbolic dimension first; stable class size identifies a class
    # even though GAP's class numbering can change between sessions.
    for class_index, matrix in sorted(
        high,
        key=lambda item: (-classes[item[0]]["fixed_rank"], classes[item[0]]["class_size"]),
    ):
        metadata = classes[class_index]
        print(
            "SCREEN",
            "class_size",
            metadata["class_size"],
            "fixed_rank",
            metadata["fixed_rank"],
            flush=True,
        )
        screen = screen_involution(matrix, walls, wall_types, gram, h0)
        record = {
            "gap_class_index": class_index,
            **metadata,
            "representative_matrix": matrix.tolist(),
            **screen,
        }
        results.append(record)
        print(
            "RESULT",
            {
                key: value
                for key, value in record.items()
                if key not in {"representative_matrix", "bad_examples"}
            },
            flush=True,
        )
    output = ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_element_class_screen.json"
    output.write_text(
        json.dumps(
            {
                "schema": "ss_fermat_quartic_char3_c2_element_class_screen_v1",
                "full_group_order": 26127360,
                "full_group_structure": "PSU(4,3) : D8",
                "element_class_count": len(classes),
                "element_classes": classes,
                "high_fixed_rank_casts": results,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print("WROTE", output)
    if completed.stderr:
        print("GAP STDERR", completed.stderr)


if __name__ == "__main__":
    main()
