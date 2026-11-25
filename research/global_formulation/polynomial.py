"""Small exact sparse polynomials; symbols are local to an occupation measure."""
from fractions import Fraction as F
from itertools import combinations_with_replacement


class Poly(dict):
    def __init__(self, terms=()):
        super().__init__({tuple(sorted(m)): F(c) for m, c in dict(terms).items() if c})

    def __add__(self, other):
        other = as_poly(other)
        result = dict(self)
        for m, c in other.items():
            result[m] = result.get(m, F(0))+c
        return Poly(result)

    __radd__ = __add__

    def __neg__(self):
        return Poly({m: -c for m, c in self.items()})

    def __sub__(self, other):
        return self+-as_poly(other)

    def __rsub__(self, other):
        return as_poly(other)+-self

    def __mul__(self, other):
        result = {}
        for a, x in self.items():
            for b, y in as_poly(other).items():
                m = tuple(sorted(a+b))
                result[m] = result.get(m, F(0))+x*y
        return Poly(result)

    __rmul__ = __mul__

    def __truediv__(self, number):
        return self*F(1, number)

    def __pow__(self, power):
        if type(power) is not int or power < 0:
            raise ValueError("Nonnegative integer polynomial power required")
        result = as_poly(1)
        for _ in range(power):
            result = result*self
        return result

    def rename(self, names):
        result = Poly()
        for m, c in self.items():
            result += Poly({tuple(sorted(names.get(v, v) for v in m)): c})
        return result

    @property
    def support(self):
        return set(v for m in self for v in m)


def as_poly(value):
    return value if isinstance(value, Poly) else Poly({(): F(value)})


def var(name):
    return Poly({(name,): 1})


def monomials(names, degree):
    names = sorted(set(names))
    return [Poly({m: 1}) for d in range(degree+1) for m in combinations_with_replacement(names, d)]


def determinant(matrix):
    if len(matrix) == 1:
        return as_poly(matrix[0][0])
    return sum(((-1)**j*as_poly(matrix[0][j])*determinant(
        [[row[k] for k in range(len(matrix)) if k != j] for row in matrix[1:]])
        for j in range(len(matrix))), as_poly(0))
