"""Small pure Fraction polynomial/field/Sturm primitives; no driver or I/O.

The caller must supply a PROVED irreducible monic minimal polynomial.
Sturm validation certifies the real embeddings, not irreducibility or an
integral basis. Integer power coordinates characterize integrality only
when the caller separately supplies O_K=Z[theta]. No numerical library
or floating-point comparison occurs here.

New source, prepared without execution. A later producer and independent
replay must name any sharing of this arithmetic dependency explicitly.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
from math import gcd, lcm
from typing import Callable, Iterable, Sequence

Poly = tuple[Q, ...]
Guard = Callable[[], None] | None


def _guard(guard: Guard) -> None:
    if guard is not None:
        guard()


def trim(values: Iterable[Q | int]) -> Poly:
    out = [Q(v) for v in values]
    if not out:
        return (Q(0),)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def padd(a: Poly, b: Poly) -> Poly:
    return trim((a[j] if j < len(a) else 0) +
                (b[j] if j < len(b) else 0)
                for j in range(max(len(a), len(b))))


def pscale(a: Poly, value: Q | int) -> Poly:
    return trim(Q(value) * v for v in a)


def psub(a: Poly, b: Poly) -> Poly:
    return padd(a, pscale(b, -1))


def pmul(a: Poly, b: Poly) -> Poly:
    out = [Q(0)] * (len(a) + len(b) - 1)
    for j, x in enumerate(a):
        for k, y in enumerate(b):
            out[j + k] += x * y
    return trim(out)


def pdivmod(a: Poly, b: Poly) -> tuple[Poly, Poly]:
    a, b = trim(a), trim(b)
    if b == (Q(0),):
        raise ZeroDivisionError("polynomial division by zero")
    q = [Q(0)] * max(1, len(a) - len(b) + 1)
    remainder = list(a)
    while trim(remainder) != (Q(0),) and len(remainder) >= len(b):
        shift = len(remainder) - len(b)
        scale = remainder[-1] / b[-1]
        q[shift] += scale
        for j, v in enumerate(b):
            remainder[shift + j] -= scale * v
        remainder = list(trim(remainder))
    return trim(q), trim(remainder)


def peval(a: Poly, x: Q) -> Q:
    out = Q(0)
    for coefficient in reversed(a):
        out = out * x + coefficient
    return out


def derivative(a: Poly) -> Poly:
    return trim(j * a[j] for j in range(1, len(a)))


def positive_primitive(a: Poly) -> Poly:
    """Multiply by a positive rational to remove content; preserve sign."""
    a = trim(a)
    denominator = lcm(*(v.denominator for v in a))
    integers = [int(v * denominator) for v in a]
    content = gcd(*(abs(v) for v in integers))
    if not content:
        return (Q(0),)
    return tuple(Q(v // content) for v in integers)


def sturm_chain(f: Poly) -> tuple[Poly, ...]:
    f = positive_primitive(f)
    if len(f) < 2:
        raise ValueError("Sturm input must be nonconstant")
    out = [f, positive_primitive(derivative(f))]
    while True:
        _, remainder = pdivmod(out[-2], out[-1])
        if remainder == (Q(0),):
            return tuple(out)
        out.append(positive_primitive(pscale(remainder, -1)))


def variation(chain: Sequence[Poly], x: Q) -> int:
    signs = []
    for p in chain:
        value = peval(p, x)
        if value:
            signs.append(1 if value > 0 else -1)
    return sum(a != b for a, b in zip(signs, signs[1:]))


def root_count(chain: Sequence[Poly], lo: Q, hi: Q) -> int:
    if not lo < hi or peval(chain[0], lo) == 0 or peval(chain[0], hi) == 0:
        raise ValueError("root-count endpoints must be ordered nonroots")
    return variation(chain, lo) - variation(chain, hi)


@dataclass(frozen=True)
class Interval:
    lo: Q
    hi: Q

    def __post_init__(self) -> None:
        object.__setattr__(self, "lo", Q(self.lo))
        object.__setattr__(self, "hi", Q(self.hi))
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    def __add__(self, other: Interval) -> Interval:
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __neg__(self) -> Interval:
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other: Interval) -> Interval:
        return self + (-other)

    def __mul__(self, other: Interval) -> Interval:
        products = (self.lo * other.lo, self.lo * other.hi,
                    self.hi * other.lo, self.hi * other.hi)
        return Interval(min(products), max(products))

    def scaled(self, value: Q | int) -> Interval:
        return self * Interval(Q(value), Q(value))

    def serial(self) -> list[str]:
        return [str(self.lo), str(self.hi)]


def interval_eval(p: Poly, interval: Interval) -> Interval:
    out = Interval(Q(0), Q(0))
    for coefficient in reversed(p):
        out = out * interval + Interval(coefficient, coefficient)
    return out


def matrix_inverse(matrix: Sequence[Sequence[Q]]) -> list[list[Q]]:
    n = len(matrix)
    if not n or any(len(row) != n for row in matrix):
        raise ValueError("inverse requires a nonempty square matrix")
    work = [[Q(v) for v in row] + [Q(i == j) for j in range(n)]
            for i, row in enumerate(matrix)]
    for column in range(n):
        pivot = next((i for i in range(column, n)
                      if work[i][column]), None)
        if pivot is None:
            raise ValueError("singular rational matrix")
        work[column], work[pivot] = work[pivot], work[column]
        value = work[column][column]
        work[column] = [v / value for v in work[column]]
        for i in range(n):
            if i != column and work[i][column]:
                value = work[i][column]
                work[i] = [x - value * y
                           for x, y in zip(work[i], work[column])]
    return [row[n:] for row in work]


def rational_rank(matrix: Sequence[Sequence[Q]]) -> int:
    if not matrix:
        return 0
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError("ragged rational matrix")
    work = [[Q(v) for v in row] for row in matrix]
    rank = 0
    for column in range(width):
        pivot = next((i for i in range(rank, len(work))
                      if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        value = work[rank][column]
        work[rank] = [v / value for v in work[rank]]
        for i in range(rank + 1, len(work)):
            value = work[i][column]
            if value:
                work[i] = [x - value * y
                           for x, y in zip(work[i], work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


class RealField:
    """Exact quotient by a declared minimal polynomial, all roots real.

    Distinct rational nonroot intervals must isolate ALL d roots in
    increasing order. Their Sturm counts are checked on construction.
    Refinements preserve the certified root, and signs use closed outer
    interval arithmetic until zero is excluded. A guard may interrupt
    any refinement; callers must record that as incomplete, never absent.
    """

    def __init__(self, polynomial: Iterable[Q | int],
                 roots: Sequence[Interval], guard: Guard = None):
        self.f = trim(polynomial)
        self.degree = len(self.f) - 1
        if self.degree < 2 or self.f[-1] != 1:
            raise ValueError("need a monic minimal polynomial of degree>=2")
        self.guard = guard
        self.chain = sturm_chain(self.f)
        self.roots = list(roots)
        if len(self.roots) != self.degree:
            raise ValueError("all field embeddings must be supplied")
        for i, interval in enumerate(self.roots):
            _guard(guard)
            if root_count(self.chain, interval.lo, interval.hi) != 1:
                raise ValueError("interval does not isolate one real root")
            if i and self.roots[i - 1].hi >= interval.lo:
                raise ValueError("root intervals must be disjoint and ordered")
        if root_count(self.chain, self.roots[0].lo,
                      self.roots[-1].hi) != self.degree:
            raise ValueError("root intervals do not cover the full degree")
        self.zero = self.element((0,))
        self.one = self.element((1,))
        self.theta = self.element((0, 1))
        self.trace_powers = self._newton_sums(2 * self.degree - 2)

    def _newton_sums(self, maximum: int) -> list[Q]:
        d = self.degree
        powers = [Q(d)]
        for k in range(1, maximum + 1):
            if k <= d:
                total = sum((self.f[d - j] * powers[k - j]
                             for j in range(1, k)), Q(0))
                total += k * self.f[d - k]
            else:
                total = sum((self.f[d - j] * powers[k - j]
                             for j in range(1, d + 1)), Q(0))
            powers.append(-total)
        return powers

    def element(self, polynomial: Iterable[Q | int]) -> Poly:
        _, remainder = pdivmod(trim(polynomial), self.f)
        return remainder + (Q(0),) * (self.degree - len(remainder))

    def add(self, a: Poly, b: Poly) -> Poly:
        return self.element(padd(a, b))

    def sub(self, a: Poly, b: Poly) -> Poly:
        return self.element(psub(a, b))

    def scale(self, a: Poly, value: Q | int) -> Poly:
        return self.element(pscale(a, value))

    def mul(self, a: Poly, b: Poly) -> Poly:
        return self.element(pmul(a, b))

    def inverse(self, a: Poly) -> Poly:
        a = self.element(a)
        if a == self.zero:
            raise ZeroDivisionError("field inverse of zero")
        old_r, remainder = self.f, trim(a)
        old_t, coefficient = (Q(0),), (Q(1),)
        while remainder != (Q(0),):
            quotient, new_r = pdivmod(old_r, remainder)
            old_r, remainder = remainder, new_r
            old_t, coefficient = coefficient, psub(old_t, pmul(quotient, coefficient))
        if len(old_r) != 1:
            raise ValueError("nonfield presentation or invalid minimal polynomial")
        result = self.element(pscale(old_t, 1 / old_r[0]))
        if self.mul(a, result) != self.one:
            raise ArithmeticError("exact inverse identity failed")
        return result

    def power(self, a: Poly, exponent: int) -> Poly:
        if exponent < 0:
            return self.power(self.inverse(a), -exponent)
        out, base = self.one, a
        while exponent:
            if exponent & 1:
                out = self.mul(out, base)
            base = self.mul(base, base)
            exponent //= 2
        return out

    def trace(self, a: Poly) -> Q:
        a = self.element(a)
        return sum((a[j] * self.trace_powers[j]
                    for j in range(self.degree)), Q(0))

    def trace_dual(self) -> tuple[list[list[Q]], list[Poly]]:
        d = self.degree
        matrix = [[self.trace_powers[j + k] for k in range(d)]
                  for j in range(d)]
        inverse = matrix_inverse(matrix)
        dual = [self.element(row) for row in inverse]
        for i in range(d):
            beta = self.element((0,) * i + (1,))
            for j in range(d):
                if self.trace(self.mul(beta, dual[j])) != Q(i == j):
                    raise ArithmeticError("trace-dual identity failed")
        return matrix, dual

    def enclose(self, a: Poly, embedding: int) -> Interval:
        return interval_eval(trim(self.element(a)), self.roots[embedding])

    def refine(self, embedding: int) -> None:
        _guard(self.guard)
        interval = self.roots[embedding]
        midpoint = (interval.lo + interval.hi) / 2
        if peval(self.f, midpoint) == 0:
            raise ValueError("unexpected rational root of declared minimal polynomial")
        left = root_count(self.chain, interval.lo, midpoint)
        if left == 1:
            self.roots[embedding] = Interval(interval.lo, midpoint)
        elif left == 0 and root_count(self.chain, midpoint, interval.hi) == 1:
            self.roots[embedding] = Interval(midpoint, interval.hi)
        else:
            raise ArithmeticError("Sturm refinement lost its unique root")

    def refine_width(self, width: Q) -> None:
        if width <= 0:
            raise ValueError("root width must be positive")
        for i in range(self.degree):
            while self.roots[i].hi - self.roots[i].lo > width:
                self.refine(i)

    def sign(self, a: Poly, embedding: int) -> int:
        a = self.element(a)
        if a == self.zero:
            return 0
        while True:
            _guard(self.guard)
            interval = self.enclose(a, embedding)
            if interval.lo > 0:
                return 1
            if interval.hi < 0:
                return -1
            self.refine(embedding)

    def is_integral_power_basis(self, a: Poly) -> bool:
        return all(v.denominator == 1 for v in self.element(a))

    def generated_degree(self, a: Poly) -> int:
        powers = [self.power(a, j) for j in range(self.degree)]
        return rational_rank(powers)

    def serial_roots(self) -> list[list[str]]:
        return [interval.serial() for interval in self.roots]


def serial_element(a: Poly) -> list[str]:
    return [str(v) for v in a]
