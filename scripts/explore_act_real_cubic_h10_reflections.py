#!/usr/bin/env python3
"""Explore real reflection walls in the ACT cubic-threefold period lattice.

This is an exploratory, deliberately non-promotional program.  It fixes one
explicit antiunitary involution of the Allcock--Carlson--Toledo Eisenstein
lattice, diagonalizes its real fixed lattice, derives the exact mod-3 gate for
real reflections to extend to Eisenstein isometries, and runs the initial
levels of Vinberg's algorithm.  ``--cdd`` asks cddlib for the extreme rays of
the *partial* chamber.  Finite volume of a partial chamber is not, by itself,
completeness of the full reflection chamber.

Coordinates use Q = diag(1,...,1,-1).  A primitive Q-root has square 1 or 2.
The ACT extension gate below proves that no square-2 root can extend for this
antiunitary class, so all surviving Coxeter walls are right-angled.
"""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Iterator, Sequence
from fractions import Fraction
import math
import resource
import shutil
import subprocess
import tempfile
import time

import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix


DIMENSION = 10
RANK = DIMENSION + 1


def _act_data() -> dict[str, object]:
    """Return exact ACT, fixed-basis, diagonal, and congruence matrices."""

    field = sp.QQ.algebraic_field(sp.sqrt(-3))
    zero, one = field.zero, field.one
    omega = field.convert((-1 + sp.sqrt(-3)) / 2)
    omega_bar = -1 - omega
    theta = omega - omega_bar

    hermitian = [[zero for _ in range(RANK)] for _ in range(RANK)]
    hermitian[0][0] = field.convert(3)
    for offset in (1, 5):
        for index in range(4):
            hermitian[offset + index][offset + index] = field.convert(3)
        for index in range(3):
            hermitian[offset + index][offset + index + 1] = theta
            hermitian[offset + index + 1][offset + index] = -theta
    hermitian[9][10] = theta
    hermitian[10][9] = -theta

    fixed = [[zero for _ in range(RANK)] for _ in range(RANK)]
    fixed[0][0] = one
    column = 1
    for left, right in ((1, 4), (2, 3), (5, 8), (6, 7), (9, 10)):
        fixed[left][column] = fixed[right][column] = one
        column += 1
        # ACT write their Hermitian form linear in the first variable.  The
        # conjugate fixed generator (omega_bar, omega), rather than
        # (omega, omega_bar), is the convention compatible with the explicit
        # diagonalizer below.  The two choices are unimodularly equivalent
        # bases of the same chi=C*conjugation fixed module.
        fixed[left][column] = omega_bar
        fixed[right][column] = omega
        column += 1
    if column != RANK:
        raise AssertionError("the fixed basis has the wrong rank")

    p4 = [
        [-1, -1, -1, 0],
        [-1, -1, 0, 0],
        [0, 1, 0, 0],
        [-1, -1, -1, -1],
    ]
    diagonalizer = [[0 for _ in range(RANK)] for _ in range(RANK)]
    diagonalizer[0][0] = 1
    for offset in (1, 5):
        for row in range(4):
            for col in range(4):
                diagonalizer[offset + row][offset + col] = p4[row][col]
    for row, values in enumerate(((1, 0), (1, 1))):
        for col, value in enumerate(values):
            diagonalizer[9 + row][9 + col] = value

    frame = [
        [
            sum(
                (
                    fixed[row][middle]
                    * field.convert(diagonalizer[middle][col])
                    for middle in range(RANK)
                ),
                zero,
            )
            for col in range(RANK)
        ]
        for row in range(RANK)
    ]
    inverse = DomainMatrix(frame, (RANK, RANK), field).inv().to_list()

    def conjugate(value: object) -> object:
        expression = field.to_sympy(value)
        return field.convert(sp.conjugate(expression))

    def multiply(left: Sequence[Sequence[object]], right: Sequence[Sequence[object]]):
        return [
            [
                sum(
                    (left[row][middle] * right[middle][col] for middle in range(len(right))),
                    zero,
                )
                for col in range(len(right[0]))
            ]
            for row in range(len(left))
        ]

    # ACT's displayed convention is h(x,y)=x^T H conjugate(y), linear in x.
    fixed_transpose = [
        [fixed[col][row] for col in range(RANK)] for row in range(RANK)
    ]
    fixed_conjugate = [
        [conjugate(fixed[row][col]) for col in range(RANK)]
        for row in range(RANK)
    ]
    raw_fixed_gram = multiply(
        fixed_transpose, multiply(hermitian, fixed_conjugate)
    )
    q_fixed_gram = sp.Matrix(
        [
            [field.to_sympy(field.convert(sp.Rational(2, 3)) * value) for value in row]
            for row in raw_fixed_gram
        ]
    )
    u = sp.Matrix(diagonalizer)
    diagonal_gram = u.T * q_fixed_gram * u
    expected = sp.diag(*([2] * DIMENSION + [-2]))
    if diagonal_gram != expected:
        raise ArithmeticError("the fixed lattice did not diagonalize to 2 I_(10,1)")

    def eisenstein_pair(value: object, multiplier: int = 1) -> tuple[int, int]:
        """Return integral coefficients a,b in a+b*omega."""

        coefficients = (value * field.convert(multiplier)).rep
        if not coefficients:
            return (0, 0)
        if len(coefficients) == 1:
            c1, c0 = field.dom.zero, coefficients[0]
        else:
            c1, c0 = coefficients
        a, b = c0 + c1, 2 * c1
        if a.denominator != 1 or b.denominator != 1:
            raise ArithmeticError("an asserted Eisenstein integer has a denominator")
        return int(a), int(b)

    # In E/(3), use epsilon=theta and epsilon^2=0.  Since
    # omega = 1+2 epsilon mod 3, a+b omega=(a+b)+2b epsilon.
    def epsilon_pair(value: object, multiplier: int = 1) -> tuple[int, int]:
        a, b = eisenstein_pair(value, multiplier)
        return ((a + b) % 3, (2 * b) % 3)

    frame_mod3 = np.array(
        [[epsilon_pair(value) for value in row] for row in frame], dtype=np.int16
    )
    adjoint_mod3 = np.array(
        [[epsilon_pair(value, 3) for value in row] for row in inverse],
        dtype=np.int16,
    )
    return {
        "hermitian": hermitian,
        "fixed_basis": fixed,
        "diagonalizer": diagonalizer,
        "frame": frame,
        "inverse": inverse,
        "q_fixed_gram": q_fixed_gram,
        "diagonal_gram": diagonal_gram,
        "frame_mod3": frame_mod3,
        "adjoint_mod3": adjoint_mod3,
    }


def _extension_mask(vectors: np.ndarray, data: dict[str, object]) -> np.ndarray:
    """Exact iff gate for E-linear extension of the real reflection.

    Put F=BU and J=3F^{-1}.  For a square-one or square-two primitive root r,
    the extension matrix is I-(2/r^2)(Fr)(r^T QF^{-1}).  Since 2/r^2 is a
    unit modulo 3, it is Eisenstein integral iff every coordinate product
    L_i R_j vanishes in E/(3), where L=Fr and R=r^TQJ.  In
    E/(3)=F_3[epsilon]/(epsilon^2), this says L=0, or R=0, or both vectors
    have every coordinate in (epsilon).
    """

    frame = data["frame_mod3"]
    adjoint = data["adjoint_mod3"]
    if not isinstance(frame, np.ndarray) or not isinstance(adjoint, np.ndarray):
        raise TypeError("congruence matrices are not numpy arrays")
    left_constant = (vectors @ frame[:, :, 0].T) % 3
    left_epsilon = (vectors @ frame[:, :, 1].T) % 3
    q_vectors = vectors.copy()
    q_vectors[:, -1] *= -1
    right_constant = (q_vectors @ adjoint[:, :, 0]) % 3
    right_epsilon = (q_vectors @ adjoint[:, :, 1]) % 3

    left_zero = np.all((left_constant == 0) & (left_epsilon == 0), axis=1)
    right_zero = np.all((right_constant == 0) & (right_epsilon == 0), axis=1)
    left_in_epsilon = np.all(left_constant == 0, axis=1)
    right_in_epsilon = np.all(right_constant == 0, axis=1)
    return left_zero | right_zero | (left_in_epsilon & right_in_epsilon)


def _sum_of_squares(dimension: int, square: int, prefix: tuple[int, ...] = ()) -> Iterator[tuple[int, ...]]:
    if dimension == 0:
        if square == 0:
            yield prefix
        return
    bound = math.isqrt(square)
    for coordinate in range(-bound, bound + 1):
        yield from _sum_of_squares(
            dimension - 1, square - coordinate * coordinate, prefix + (coordinate,)
        )


def _lorentz_products(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return left[:, :DIMENSION] @ right[:, :DIMENSION].T - np.outer(
        left[:, -1], right[:, -1]
    )


def _chunks(rows: Iterable[tuple[int, ...]], size: int) -> Iterator[np.ndarray]:
    batch: list[tuple[int, ...]] = []
    for row in rows:
        batch.append(row)
        if len(batch) == size:
            yield np.array(batch, dtype=np.int32)
            batch = []
    if batch:
        yield np.array(batch, dtype=np.int32)


def _vinberg_levels(
    data: dict[str, object], max_height: int, chunk_size: int
) -> tuple[np.ndarray, list[dict[str, int | float]]]:
    walls = np.array([[1] + [0] * DIMENSION], dtype=np.int32)
    ledger: list[dict[str, int | float]] = [
        {"height": 0, "raw": 2, "dominant": 1, "extended": 1, "walls": 1}
    ]
    for height in range(1, max_height + 1):
        started = time.monotonic()
        raw_count = dominant_count = 0
        candidates: list[np.ndarray] = []
        rows = _sum_of_squares(DIMENSION, height * height + 1)
        for positive in _chunks(rows, chunk_size):
            raw_count += len(positive)
            batch = np.column_stack(
                (positive, np.full(len(positive), height, dtype=np.int32))
            )
            dominant = np.all(_lorentz_products(batch, walls) <= 0, axis=1)
            batch = batch[dominant]
            dominant_count += len(batch)
            if len(batch):
                batch = batch[_extension_mask(batch, data)]
                if len(batch):
                    candidates.append(batch)
        level = (
            np.vstack(candidates)
            if candidates
            else np.empty((0, RANK), dtype=np.int32)
        )
        if len(level):
            # Completing the height-zero chamber should make equal-height
            # candidates pairwise nonpositive.  Check in blocks to avoid a
            # large permanent Gram matrix.
            for start in range(0, len(level), chunk_size):
                products = _lorentz_products(level[start : start + chunk_size], level)
                for local in range(min(chunk_size, len(level) - start)):
                    products[local, start + local] = 0
                if np.any(products > 0):
                    raise ArithmeticError(
                        f"positive equal-height pairing at level {height}; "
                        "the Vinberg selection rule needs refinement"
                    )
            walls = np.vstack((walls, level))
        ledger.append(
            {
                "height": height,
                "raw": raw_count,
                "dominant": dominant_count,
                "extended": len(level),
                "walls": len(walls),
                "seconds": round(time.monotonic() - started, 6),
            }
        )
    return walls, ledger


def _primitive_ray(tokens: Sequence[str]) -> tuple[int, ...]:
    values = [Fraction(token) for token in tokens]
    denominator = math.lcm(*(value.denominator for value in values))
    integers = [
        value.numerator * (denominator // value.denominator) for value in values
    ]
    divisor = math.gcd(*(abs(value) for value in integers))
    if divisor == 0:
        raise ArithmeticError("cddlib emitted a zero ray")
    return tuple(value // divisor for value in integers)


def _cdd_audit(
    walls: np.ndarray, executable: str, timeout: int, memory_gib: float
) -> dict[str, int | str]:
    inequalities = np.column_stack((-walls[:, :DIMENSION], walls[:, -1]))
    lines = [
        "H-representation",
        "begin",
        f"{len(inequalities)} {RANK + 1} rational",
    ]
    lines.extend("0 " + " ".join(map(str, row)) for row in inequalities.tolist())
    lines.append("end")
    memory_bytes = int(memory_gib * (1 << 30))

    def limit_child_memory() -> None:
        # Linux enforces this address-space ceiling in the child.  Some macOS
        # kernels reject RLIMIT_AS outright; the parent-side RSS watchdog below
        # is therefore mandatory rather than a cosmetic fallback.
        try:
            resource.setrlimit(resource.RLIMIT_AS, (memory_bytes, memory_bytes))
        except (OSError, ValueError):
            pass

    def resident_tree_bytes(pid: int) -> int:
        try:
            import psutil  # type: ignore[import-not-found]

            process = psutil.Process(pid)
            resident = process.memory_info().rss
            try:
                descendants = process.children(recursive=True)
            except (PermissionError, psutil.AccessDenied):
                # Managed macOS blocks psutil's all-PID scan.  Here ``pid`` is
                # cddexec_gmp itself, so direct RSS is the required ceiling.
                return resident
            for child in descendants:
                try:
                    resident += child.memory_info().rss
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            return resident
        except ImportError:
            # cddexec_gmp does not normally fork.  Keep a stdlib-only fallback
            # for environments that do not provide psutil.
            measured = subprocess.run(
                ["/bin/ps", "-o", "rss=", "-p", str(pid)],
                text=True,
                capture_output=True,
                check=False,
            )
            try:
                return int(measured.stdout.strip()) * 1024
            except ValueError:
                return 0
        except Exception:
            # The process may have exited between poll and measurement.
            return 0

    payload = "\n".join(lines) + "\n"
    with (
        tempfile.TemporaryFile(mode="w+") as input_file,
        tempfile.TemporaryFile(mode="w+") as output_file,
        tempfile.TemporaryFile(mode="w+") as error_file,
    ):
        input_file.write(payload)
        input_file.seek(0)
        process = subprocess.Popen(
            [executable, "--rep"],
            stdin=input_file,
            stdout=output_file,
            stderr=error_file,
            text=True,
            preexec_fn=limit_child_memory,
        )
        started = time.monotonic()
        peak_resident_bytes = 0
        killed_for_memory = False
        killed_for_timeout = False
        while process.poll() is None:
            peak_resident_bytes = max(
                peak_resident_bytes, resident_tree_bytes(process.pid)
            )
            if peak_resident_bytes > memory_bytes:
                killed_for_memory = True
                process.kill()
                break
            if time.monotonic() - started > timeout:
                killed_for_timeout = True
                process.kill()
                break
            time.sleep(0.05)
        process.wait()
        output_file.seek(0)
        stdout = output_file.read()
        error_file.seek(0)
        stderr = error_file.read()
    if killed_for_memory:
        raise MemoryError(
            f"cddlib exceeded its {memory_gib:g} GiB RSS ceiling; "
            f"peak sampled RSS was {peak_resident_bytes / (1 << 30):.3f} GiB"
        )
    if killed_for_timeout:
        raise TimeoutError(f"cddlib exceeded its {timeout}-second timeout")
    if process.returncode:
        raise RuntimeError(stderr[-4000:] or stdout[-4000:])
    output = stdout.splitlines()
    begin = output.index("begin")
    count, width, kind = output[begin + 1].split()
    if int(width) != RANK + 1 or kind not in {"rational", "real"}:
        raise ArithmeticError("cddlib returned an unexpected V-representation")
    rays = []
    for line in output[begin + 2 : begin + 2 + int(count)]:
        tokens = line.split()
        if Fraction(tokens[0]) != 0:
            raise ArithmeticError("the homogeneous chamber acquired a vertex row")
        rays.append(_primitive_ray(tokens[1:]))
    ray_matrix = np.array(rays, dtype=object)
    squares = np.sum(ray_matrix[:, :DIMENSION] ** 2, axis=1) - ray_matrix[:, -1] ** 2
    return {
        "engine": executable,
        "partial_wall_count": len(walls),
        "peak_sampled_rss_bytes": peak_resident_bytes,
        "extreme_ray_count": len(rays),
        "finite_timelike": int(np.sum(squares < 0)),
        "ideal_lightlike": int(np.sum(squares == 0)),
        "hyperideal_spacelike": int(np.sum(squares > 0)),
    }


def _stream_gram_summary(
    walls: np.ndarray, block_size: int
) -> dict[str, int | list[int]]:
    """Summarize the wall Gram matrix without allocating its square."""

    norms = sorted(
        set(
            map(
                int,
                np.sum(walls[:, :DIMENSION] ** 2, axis=1)
                - walls[:, -1] ** 2,
            )
        )
    )
    positive_offdiagonal = 0
    zero_offdiagonal = 0
    minimum_pairing: int | None = None
    for start in range(0, len(walls), block_size):
        stop = min(start + block_size, len(walls))
        products = _lorentz_products(walls[start:stop], walls)
        for local, absolute in enumerate(range(start, stop)):
            products[local, absolute] = 0
        positive_offdiagonal += int(np.sum(products > 0))
        zero_offdiagonal += int(np.sum(products == 0)) - (stop - start)
        offdiagonal = products[
            np.arange(stop - start)[:, None] != np.arange(len(walls))[None, :]
        ]
        if len(offdiagonal):
            block_minimum = int(np.min(offdiagonal))
            minimum_pairing = (
                block_minimum
                if minimum_pairing is None
                else min(minimum_pairing, block_minimum)
            )
    return {
        "norms": norms,
        "positive_offdiagonal": positive_offdiagonal,
        "zero_offdiagonal_pairs": zero_offdiagonal // 2,
        "minimum_pairing": minimum_pairing if minimum_pairing is not None else 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-height", type=int, default=3)
    parser.add_argument("--chunk-size", type=int, default=20_000)
    parser.add_argument("--cdd", action="store_true")
    parser.add_argument("--cdd-timeout", type=int, default=1800)
    parser.add_argument("--max-cdd-walls", type=int, default=1_000)
    parser.add_argument("--cdd-memory-gib", type=float, default=2.0)
    parser.add_argument("--gram-block-size", type=int, default=512)
    arguments = parser.parse_args()

    data = _act_data()
    residues = np.array(
        list(__import__("itertools").product(range(3), repeat=RANK)), dtype=np.int16
    )
    extension = _extension_mask(residues, data)
    q_mod3 = (
        np.sum(residues[:, :DIMENSION] ** 2, axis=1) - residues[:, -1] ** 2
    ) % 3
    print(
        "fixed_gram=2*I_(10,1)",
        f"extendable_residues={int(np.sum(extension))}/{3**RANK}",
        f"extendable_q_mod3={np.unique(q_mod3[extension]).tolist()}",
        "norm2_reflections=0",
    )
    walls, ledger = _vinberg_levels(data, arguments.max_height, arguments.chunk_size)
    for row in ledger:
        print(row)
    gram_summary = _stream_gram_summary(walls, arguments.gram_block_size)
    print(
        f"partial_walls={len(walls)}",
        f"norms={gram_summary['norms']}",
        f"positive_offdiagonal={gram_summary['positive_offdiagonal']}",
        f"zero_offdiagonal_pairs={gram_summary['zero_offdiagonal_pairs']}",
        f"minimum_pairing={gram_summary['minimum_pairing']}",
    )
    if arguments.cdd:
        if len(walls) > arguments.max_cdd_walls:
            raise RuntimeError(
                f"refusing an unbounded cdd run on {len(walls)} walls; "
                f"the safety ceiling is {arguments.max_cdd_walls}.  First "
                "compress by symmetry or raise --max-cdd-walls explicitly; "
                "the child still retains --cdd-memory-gib."
            )
        executable = shutil.which("cddexec_gmp")
        if executable is None:
            raise RuntimeError("cddexec_gmp is required for --cdd")
        print(
            _cdd_audit(
                walls,
                executable,
                arguments.cdd_timeout,
                arguments.cdd_memory_gib,
            )
        )


if __name__ == "__main__":
    main()
