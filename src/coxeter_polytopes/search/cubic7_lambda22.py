"""Exact bounded Vinberg tools for ``Lambda_22`` over cubic seven.

The positive lattice is the 22-dimensional laminated lattice from the
Nebe--Sloane catalogue.  Only its authoritative integral Gram matrix is
vendored; no rounded Euclidean basis is imported.  The catalogue gives
determinant 12, minimum 4, and kissing number 49896.

This module deliberately separates three assertions:

* exact integral matrix invariants;
* independent shortest-vector enumeration;
* any later orbit analysis of the minimal shell.

In particular, the shell count is never replaced by an unproved
transitivity assumption.  The first implemented coefficient is the unit
``a=theta``.  Its complete primitive mirror-norm alphabet is
``1,2,4,3,6,12``: the exponent two at the inert prime above two is essential
because the Smith invariants of ``Lambda_22`` are ``1^20,2,6``.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
from math import gcd, isqrt, prod
from pathlib import Path
from typing import Iterable, Sequence

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import (
    hermite_normal_form,
    smith_normal_decomp,
    smith_normal_form,
)

from coxeter_polytopes.exact.cubic7 import CUBIC7_ZERO, THETA7, Cubic7
from coxeter_polytopes.search.silver_d24_plus import _integer_determinant


SOURCE_URL = "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/LAMBDA22.html"
SOURCE_HTML_SHA256 = "b4e4f08d2c4a4926599b6551ab852b578877f7c3cbebbfed57ea8b8924e3aab4"
GRAM_FILE_SHA256 = "4f3dbdc96e7aead7fae09a834abdbe2d5319620ece24656720f13184dd85b2f1"
GRAM_STREAM_SHA256 = "53e1adff6866521330edfeb17d79e9813398c9f3dc73503713761e4153f46495"
GRAM_PATH = (
    Path(__file__).resolve().parents[3]
    / "data/lattices/lambda22_nebe_sloane_gram.txt"
)

# Tr(theta^(i+j)) in the integral power basis 1,theta,theta^2.
TRACE_GRAM = (
    (3, -1, 5),
    (-1, 5, -4),
    (5, -4, 13),
)

# Exact PARI/GP qfminim certificates in the canonical HNF bases returned by
# lambda22_pairing_radical_basis.  Counts include both signs; PARI stores one
# representative of each sign pair in the returned matrix.
PARI_QFMINIM_CERTIFICATES = {
    1: {
        "index": 1,
        "ordinary_minimum": 4,
        "ordinary_kissing_number": 49_896,
        "ordinary_representatives_mod_sign": 24_948,
        "ordinary_vector_stream_bytes": 2_798_502,
        "ordinary_vector_stream_sha256": (
            "9b752381c953619264cf7eb1aba9b73480e6de94553a68db65eec8bfda30f428"
        ),
        "trace_minimum": 12,
        "trace_kissing_number": 49_896,
        "trace_representatives_mod_sign": 24_948,
        "trace_vector_stream_bytes": 6_091_638,
        "trace_vector_stream_sha256": (
            "68f28fe35482a24842236914e18161c8bcd7fbc325db3f7efec9f6cee64bb0d7"
        ),
    },
    2: {
        "index": 2**20,
        "ordinary_minimum": 12,
        "ordinary_kissing_number": 8_448,
        "ordinary_representatives_mod_sign": 4_224,
        "ordinary_vector_stream_bytes": 460_105,
        "ordinary_vector_stream_sha256": (
            "7e1a12da7f909431322b7d18ec7fb04b0823d94262ddc352a7f357b8a98de35d"
        ),
        "trace_minimum": 36,
        "trace_kissing_number": 8_448,
        "trace_representatives_mod_sign": 4_224,
        "trace_vector_stream_bytes": 1_017_673,
        "trace_vector_stream_sha256": (
            "85033270f30c8cadb7d4721a31845746c55a395e06d59d526003f9310da963ee"
        ),
    },
    3: {
        "index": 3**21,
        "ordinary_minimum": 24,
        "ordinary_kissing_number": 1_782,
        "ordinary_representatives_mod_sign": 891,
        "ordinary_vector_stream_bytes": 100_647,
        "ordinary_vector_stream_sha256": (
            "5444d0e47077d2c0d102cc997e74303a09df6dbc9e94b45b56a9aef9012fd1fd"
        ),
        "trace_minimum": 72,
        "trace_kissing_number": 1_782,
        "trace_representatives_mod_sign": 891,
        "trace_vector_stream_bytes": 218_259,
        "trace_vector_stream_sha256": (
            "896e61ddda939787709edd6927179a929360cb2e3c3935cf9168cc48dfed8884"
        ),
    },
    6: {
        "index": 2**20 * 3**21,
        "ordinary_minimum": 96,
        "ordinary_kissing_number": 1_782,
        "ordinary_representatives_mod_sign": 891,
        "ordinary_vector_stream_bytes": 100_791,
        "ordinary_vector_stream_sha256": (
            "80136af8b02b3aedc70110d897fdc6d47f6e278e214c73316105129161028be5"
        ),
        "trace_minimum": 288,
        "trace_kissing_number": 1_782,
        "trace_representatives_mod_sign": 891,
        "trace_vector_stream_bytes": 218_403,
        "trace_vector_stream_sha256": (
            "0fc0d0fd3366b45f35689e1066e2efdf03f6bea25c8dfbe3023bc083fb602ebd"
        ),
    },
}


@lru_cache(maxsize=1)
def lambda22_gram() -> tuple[tuple[int, ...], ...]:
    rows: list[list[int]] = []
    dimension = None
    for raw in GRAM_PATH.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        values = [int(value) for value in line.split()]
        if dimension is None:
            if values != [22, 0]:
                raise ValueError("the Lambda22 Gram header must be '22 0'")
            dimension = 22
            continue
        rows.append(values)
    if dimension != 22 or len(rows) != 22:
        raise ValueError("the Lambda22 lower triangle has the wrong dimension")
    if any(len(row) != index + 1 for index, row in enumerate(rows)):
        raise ValueError("the Lambda22 lower triangle has a malformed row")
    return tuple(
        tuple(rows[max(i, j)][min(i, j)] for j in range(22))
        for i in range(22)
    )


def lambda22_gram_stream_sha256() -> str:
    digest = sha256()
    for row in lambda22_gram():
        digest.update((" ".join(map(str, row)) + "\n").encode("ascii"))
    return digest.hexdigest()


def lambda22_gram_file_sha256() -> str:
    return sha256(GRAM_PATH.read_bytes()).hexdigest()


def lambda22_basis_audit() -> bool:
    gram = lambda22_gram()
    rational_gram = tuple(tuple(Fraction(value) for value in row) for row in gram)
    return (
        len(gram) == 22
        and all(len(row) == 22 for row in gram)
        and all(gram[i][j] == gram[j][i] for i in range(22) for j in range(22))
        and all(gram[i][i] % 2 == 0 for i in range(22))
        and _integer_determinant(rational_gram) == 12
        and lambda22_gram_file_sha256() == GRAM_FILE_SHA256
        and lambda22_gram_stream_sha256() == GRAM_STREAM_SHA256
    )


@lru_cache(maxsize=1)
def lambda22_smith_invariants() -> tuple[int, ...]:
    smith = smith_normal_form(Matrix(lambda22_gram()), domain=ZZ)
    return tuple(abs(int(smith[index, index])) for index in range(22))


@lru_cache(maxsize=None)
def lambda22_pairing_radical_basis(
    modulus: int,
) -> tuple[tuple[int, ...], ...]:
    """Return a canonical column basis for ``{y : G*y=0 mod modulus}``.

    Smith decomposition first gives an exact basis, and column Hermite normal
    form makes the result deterministic and small.  Rows of the returned
    tuple are the rows of the 22-by-22 basis matrix; its columns are vectors
    in the catalogue basis.
    """

    if modulus not in (1, 2, 3, 6):
        raise ValueError("the implemented Lambda22 radical moduli are 1,2,3,6")
    if modulus == 1:
        return tuple(
            tuple(int(i == j) for j in range(22)) for i in range(22)
        )
    gram = Matrix(lambda22_gram())
    diagonal, _left, right = smith_normal_decomp(gram, domain=ZZ)
    factors = tuple(
        modulus // gcd(modulus, abs(int(diagonal[index, index])))
        for index in range(22)
    )
    basis = hermite_normal_form(right * Matrix.diag(*factors))
    expected_index = prod(factors)
    if abs(int(basis.det())) != expected_index:
        raise AssertionError("the Lambda22 radical basis has the wrong index")
    if any(int(value) % modulus for value in gram * basis):
        raise AssertionError("the Lambda22 radical basis violates its pairing gate")
    return tuple(
        tuple(int(basis[i, j]) for j in range(22)) for i in range(22)
    )


@lru_cache(maxsize=None)
def lambda22_pairing_radical_gram(
    modulus: int,
) -> tuple[tuple[int, ...], ...]:
    gram = Matrix(lambda22_gram())
    basis = Matrix(lambda22_pairing_radical_basis(modulus))
    radical = basis.T * gram * basis
    return tuple(
        tuple(int(radical[i, j]) for j in range(22)) for i in range(22)
    )


def lambda22_pair(
    left: Sequence[Cubic7 | int], right: Sequence[Cubic7 | int]
) -> Cubic7:
    if len(left) != 22 or len(right) != 22:
        raise ValueError("Lambda22 basis vectors have length 22")
    gram = lambda22_gram()
    return sum(
        (
            Cubic7(gram[i][j])
            * Cubic7.coerce(left[i])
            * Cubic7.coerce(right[j])
            for i in range(22)
            for j in range(22)
            if gram[i][j] and left[i] and right[j]
        ),
        start=CUBIC7_ZERO,
    )


@dataclass(frozen=True, slots=True)
class Lambda22ShellAudit:
    determinant: int
    minimum: int
    kissing_number: int
    representative_count_mod_sign: int
    vector_stream_sha256: str


def parse_pari_qfminim_output(output: str) -> Lambda22ShellAudit:
    """Parse a deliberately tiny stable line emitted by the audit script."""

    fields = dict(item.split("=", 1) for item in output.strip().split())
    return Lambda22ShellAudit(
        determinant=int(fields["det"]),
        minimum=int(fields["minimum"]),
        kissing_number=int(fields["kissing"]),
        representative_count_mod_sign=int(fields["representatives"]),
        vector_stream_sha256=fields["sha256"],
    )


@dataclass(frozen=True, slots=True)
class Cubic7Lambda22Seed:
    seed_id: str
    timelike_coefficient: Cubic7
    support: str
    canonical_mirror_norms: tuple[Cubic7, ...]
    local_scope: str

    @property
    def admissible(self) -> bool:
        return self.timelike_coefficient.sign(0) > 0 and all(
            self.timelike_coefficient.sign(place) < 0 for place in (1, 2)
        )


LAMBDA22_UNIT_SEED = Cubic7Lambda22Seed(
    seed_id="cubic7-H22-Lambda22-a-theta",
    timelike_coefficient=THETA7,
    support="unit coefficient; spatial discriminant primes (2),(3)",
    canonical_mirror_norms=(
        Cubic7(1),
        Cubic7(2),
        Cubic7(4),
        Cubic7(3),
        Cubic7(6),
        Cubic7(12),
    ),
    local_scope=(
        "For primitive r=(x,y), put I=(theta*x,(y,Lambda22)). Reflection "
        "integrality gives (q(r))|2I. Away from the inert primes (2),(3), "
        "local unimodularity and primitivity force v(q)=0. The Smith diagonal "
        "1^20,2,6 has exponent one at both discriminant primes. At (3), this "
        "gives v_3(q)<=1. At (2), an odd x makes q odd by evenness; when x is "
        "even and y is primitive, v_2((y,Lambda22))<=1 and the factor 2 in "
        "2I permits v_2(q)<=2. Thus v_2(q) is 0,1,2 and v_3(q) is 0,1. "
        "Class number one and U^+=U^2 leave the six complete canonical norms "
        "1,2,4,3,6,12."
    ),
)


@dataclass(frozen=True, slots=True, order=True)
class Lambda22NormGate:
    norm: Cubic7
    pairing_modulus: int
    coefficient_lattice_minimum: int
    coefficient_lattice_index: int


@lru_cache(maxsize=None)
def lambda22_norm_gate(norm: Cubic7) -> Lambda22NormGate:
    norm = Cubic7.coerce(norm)
    table = {
        Cubic7(1): (1, 4),
        Cubic7(2): (1, 4),
        Cubic7(4): (2, 12),
        Cubic7(3): (3, 24),
        Cubic7(6): (3, 24),
        Cubic7(12): (6, 96),
    }
    try:
        modulus, minimum = table[norm]
    except KeyError as error:
        raise ValueError("the norm is not in the complete Lambda22 alphabet") from error
    certificate = PARI_QFMINIM_CERTIFICATES[modulus]
    return Lambda22NormGate(norm, modulus, minimum, int(certificate["index"]))


@dataclass(frozen=True, slots=True, order=True)
class Lambda22ControlLevel:
    height: Cubic7
    norm: Cubic7
    control: Cubic7
    target: Cubic7


@lru_cache(maxsize=None)
def lambda22_control_levels(height_bound: int) -> tuple[Lambda22ControlLevel, ...]:
    """Exhaust every reflective control through exact Vinberg height ``H``."""

    from coxeter_polytopes.search.cubic7_root_blocks import _bounded_for_caps

    if height_bound < 0:
        raise ValueError("height_bound must be nonnegative")
    seed = LAMBDA22_UNIT_SEED
    coefficient = seed.timelike_coefficient
    levels: set[Lambda22ControlLevel] = set()
    for norm in seed.canonical_mirror_norms:
        caps = (
            norm * height_bound / coefficient,
            norm / (-coefficient),
            norm / (-coefficient),
        )
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0:
                continue
            if not (2 * coefficient * control / norm).is_integral:
                continue
            target = norm + coefficient * control * control
            if any(target.sign(place) < 0 for place in (1, 2)):
                continue
            height = coefficient * control * control / norm
            if (Cubic7(height_bound) - height).sign(0) < 0:
                continue
            levels.add(Lambda22ControlLevel(height, norm, control, target))
    return tuple(sorted(levels))


@dataclass(frozen=True, slots=True, order=True)
class Lambda22GramType:
    """Gram data for ``y=u+theta*v+theta^2*w`` in a spatial sublattice."""

    u_square: int
    v_square: int
    w_square: int
    u_dot_v: int
    u_dot_w: int
    v_dot_w: int
    determinant: int

    @property
    def matrix(self) -> tuple[tuple[int, int, int], ...]:
        return (
            (self.u_square, self.u_dot_v, self.u_dot_w),
            (self.u_dot_v, self.v_square, self.v_dot_w),
            (self.u_dot_w, self.v_dot_w, self.w_square),
        )

    @property
    def rank(self) -> int:
        return int(Matrix(self.matrix).rank())


def _rank3_determinant(values: tuple[int, int, int, int, int, int]) -> int:
    u, v, w, uv, uw, vw = values
    return u * v * w + 2 * uv * uw * vw - u * vw * vw - v * uw * uw - w * uv * uv


def _quotient_gram(
    matrix: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, ...], ...]:
    raw = Matrix(matrix)
    if raw.det():
        return matrix
    diagonal, left, right = smith_normal_decomp(raw, domain=ZZ)
    if left * raw * right != diagonal:
        raise AssertionError("SymPy's Smith transport convention changed")
    if abs(int(left.det())) != 1 or abs(int(right.det())) != 1:
        raise AssertionError("Smith transport matrices are not unimodular")
    nonzero = tuple(
        index for index in range(raw.rows) if diagonal[index, index] != 0
    )
    if not nonzero:
        return ()
    zero = tuple(index for index in range(raw.rows) if index not in nonzero)
    if zero and raw * right[:, zero] != Matrix.zeros(raw.rows, len(zero)):
        raise AssertionError("the Smith zero columns do not span the kernel")
    if len(nonzero) != raw.rank():
        raise AssertionError("the Smith quotient has the wrong rank")
    complement = right[:, nonzero]
    quotient = complement.T * raw * complement
    if quotient.det() <= 0:
        raise AssertionError("the transported Gram quotient is not positive definite")
    return tuple(
        tuple(int(quotient[i, j]) for j in range(len(nonzero)))
        for i in range(len(nonzero))
    )


def lambda22_integer_gram_minimum_at_least(
    matrix: tuple[tuple[int, ...], ...], minimum: int
) -> bool:
    """Decide exactly whether the positive quotient has minimum >= ``minimum``."""

    quotient = _quotient_gram(matrix)
    if not quotient:
        return False
    raw = Matrix(quotient)
    inverse = raw.inv()
    bounds = []
    for index in range(raw.rows):
        diagonal = inverse[index, index]
        cap = Fraction(int(diagonal.p), int(diagonal.q)) * (minimum - 1)
        bounds.append(isqrt(cap.numerator // cap.denominator))
    for vector in product(*(range(-bound, bound + 1) for bound in bounds)):
        if not any(vector):
            continue
        value = sum(
            quotient[i][j] * vector[i] * vector[j]
            for i in range(len(vector))
            for j in range(len(vector))
        )
        if 0 < value < minimum:
            return False
    return True


@lru_cache(maxsize=None)
def lambda22_rank3_gram_types(
    target: Cubic7,
    coefficient_lattice_minimum: int = 4,
    pairing_modulus: int = 1,
) -> tuple[Lambda22GramType, ...]:
    """Enumerate every necessary coefficient Gram type for ``y^2=target``.

    Writing ``y=u+theta*v+theta^2*w`` gives

    ``A=U+2R-W, B=2P+4R-W, C=2Q+V-2R+3W``.

    The inverse trace Gram has diagonal ``1,2/7,2/7``.  Thus, for
    ``L=Tr(target)``, completeness follows from ``U<=L`` and
    ``V,W<=2L/7``.  Positive-semidefinite types of rank one and two are
    retained; their integral kernel is quotiented exactly before the minimum
    test.
    """

    target = Cubic7.coerce(target)
    if not target.is_integral:
        raise ValueError("the Lambda22 Gram-type target must be integral")
    if coefficient_lattice_minimum <= 0 or pairing_modulus <= 0:
        raise ValueError("minimum and pairing modulus must be positive")
    a, b, c = (int(target.a), int(target.b), int(target.c))
    trace = int(target.trace())
    if trace < 0:
        return ()
    vw_bound = 2 * trace // 7
    output = []
    for v_square in range(0, vw_bound + 1, 2):
        if 0 < v_square < coefficient_lattice_minimum:
            continue
        for w_square in range(0, vw_bound + 1, 2):
            if 0 < w_square < coefficient_lattice_minimum:
                continue
            cauchy = isqrt(v_square * w_square)
            lower_from_u = (a + w_square - trace + 1) // 2
            upper_from_u = (a + w_square) // 2
            lower = max(-cauchy, lower_from_u)
            upper = min(cauchy, upper_from_u)
            for v_dot_w in range(lower, upper + 1):
                u_square = a - 2 * v_dot_w + w_square
                if u_square < 0 or u_square % 2 or u_square > trace:
                    continue
                if 0 < u_square < coefficient_lattice_minimum:
                    continue
                twice_u_dot_v = b + w_square - 4 * v_dot_w
                twice_u_dot_w = c - v_square - 3 * w_square + 2 * v_dot_w
                if twice_u_dot_v % 2 or twice_u_dot_w % 2:
                    continue
                u_dot_v = twice_u_dot_v // 2
                u_dot_w = twice_u_dot_w // 2
                if (
                    u_dot_v * u_dot_v > u_square * v_square
                    or u_dot_w * u_dot_w > u_square * w_square
                ):
                    continue
                values = (
                    u_square,
                    v_square,
                    w_square,
                    u_dot_v,
                    u_dot_w,
                    v_dot_w,
                )
                determinant = _rank3_determinant(values)
                if determinant < 0:
                    continue
                if any(value % pairing_modulus for value in values):
                    continue
                gram = (
                    (u_square, u_dot_v, u_dot_w),
                    (u_dot_v, v_square, v_dot_w),
                    (u_dot_w, v_dot_w, w_square),
                )
                if not lambda22_integer_gram_minimum_at_least(
                    gram, coefficient_lattice_minimum
                ):
                    continue
                output.append(Lambda22GramType(*values, determinant))
    return tuple(sorted(output))


@dataclass(frozen=True, slots=True)
class Lambda22LevelAudit:
    level: Lambda22ControlLevel
    ambient_gram_types: tuple[Lambda22GramType, ...]
    congruence_gram_types: tuple[Lambda22GramType, ...]
    gated_gram_types: tuple[Lambda22GramType, ...]


@lru_cache(maxsize=None)
def lambda22_level_audit(level: Lambda22ControlLevel) -> Lambda22LevelAudit:
    gate = lambda22_norm_gate(level.norm)
    ambient = lambda22_rank3_gram_types(level.target, 4, 1)
    if gate.pairing_modulus == 1:
        congruence = ambient
    else:
        congruence = lambda22_rank3_gram_types(
            level.target, 4, gate.pairing_modulus
        )
    if (
        gate.pairing_modulus == 1
        and gate.coefficient_lattice_minimum == 4
    ):
        gated = ambient
    else:
        gated = lambda22_rank3_gram_types(
            level.target,
            gate.coefficient_lattice_minimum,
            gate.pairing_modulus,
        )
    return Lambda22LevelAudit(level, ambient, congruence, gated)


def lambda22_lorentz_pair(
    seed: Cubic7Lambda22Seed,
    left: Sequence[Cubic7 | int],
    right: Sequence[Cubic7 | int],
) -> Cubic7:
    if len(left) != 23 or len(right) != 23:
        raise ValueError("a Lambda22 Lorentz vector has rank 23")
    return (
        -seed.timelike_coefficient
        * Cubic7.coerce(left[0])
        * Cubic7.coerce(right[0])
        + lambda22_pair(left[1:], right[1:])
    )


def lambda22_reflection_integral(
    seed: Cubic7Lambda22Seed,
    norm: Cubic7,
    control: Cubic7,
    spatial: Sequence[Cubic7],
) -> bool:
    if len(spatial) != 22:
        raise ValueError("a Lambda22 spatial vector has length 22")
    if not (2 * seed.timelike_coefficient * control / norm).is_integral:
        return False
    gram = lambda22_gram()
    pairings = tuple(
        sum(
            (
                Cubic7(gram[i][j]) * spatial[j]
                for j in range(22)
                if gram[i][j] and spatial[j]
            ),
            start=CUBIC7_ZERO,
        )
        for i in range(22)
    )
    return all((2 * value / norm).is_integral for value in pairings)


__all__ = [
    "GRAM_PATH",
    "GRAM_FILE_SHA256",
    "GRAM_STREAM_SHA256",
    "LAMBDA22_UNIT_SEED",
    "PARI_QFMINIM_CERTIFICATES",
    "SOURCE_HTML_SHA256",
    "SOURCE_URL",
    "TRACE_GRAM",
    "Cubic7Lambda22Seed",
    "Lambda22ControlLevel",
    "Lambda22GramType",
    "Lambda22LevelAudit",
    "Lambda22NormGate",
    "Lambda22ShellAudit",
    "lambda22_basis_audit",
    "lambda22_control_levels",
    "lambda22_gram",
    "lambda22_gram_file_sha256",
    "lambda22_gram_stream_sha256",
    "lambda22_integer_gram_minimum_at_least",
    "lambda22_level_audit",
    "lambda22_lorentz_pair",
    "lambda22_norm_gate",
    "lambda22_pair",
    "lambda22_pairing_radical_basis",
    "lambda22_pairing_radical_gram",
    "lambda22_rank3_gram_types",
    "lambda22_reflection_integral",
    "lambda22_smith_invariants",
    "parse_pari_qfminim_output",
]
