"""Reusable Pell-boost criterion for arithmetic reflection chambers.

Let ``M=<-a> + L`` over a real quadratic integer ring and let ``d in L``
have positive rational square ``m``.  Integral elements ``A,U`` satisfying

``A^2 - a*m*U^2 = 1`` and ``m | A-1``

define an integral Lorentz boost.  A finite strict-ray certificate through
height ``A^2-1`` then promotes to global nonreflectivity, provided the
basepoint is regular and the mirror norm alphabet is complete.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from coxeter_polytopes.exact.silver import SILVER_ONE, Silver


@dataclass(frozen=True, slots=True)
class PellBoostCriterion:
    coefficient: Silver
    direction_norm: int
    pell_A: Silver
    pell_U: Silver
    witness_scale: Fraction

    @property
    def spatial_from_time(self) -> Silver:
        """Coefficient of ``x*d`` in the spatial component."""

        return self.coefficient * self.pell_U

    @property
    def time_from_direction(self) -> Silver:
        """Coefficient ``U`` of ``(d,y)`` in the time component."""

        return self.pell_U

    @property
    def direction_from_direction(self) -> Silver:
        return (self.pell_A - SILVER_ONE) / self.direction_norm

    @property
    def separator_direction_slope(self) -> Silver:
        """Lower bound for ``(d,y)/x`` for a separating oriented mirror."""

        return self.coefficient * self.pell_A / self.spatial_from_time

    @property
    def control_square_bound(self) -> Silver:
        """The exact Cauchy bound ``x^2/q`` for a separating root."""

        slope = self.separator_direction_slope
        return Silver(self.direction_norm) / (
            slope * slope - self.direction_norm * self.coefficient
        )

    @property
    def separator_height_bound(self) -> Silver:
        """The Vinberg-height cutoff, identically ``A^2-1``."""

        return self.coefficient * self.control_square_bound

    @property
    def witness_separator_coefficient(self) -> Silver:
        """Lower coefficient in ``(z,r)>=coefficient*x``, ``z=(1,t*d)``."""

        return (
            Silver(self.witness_scale) * self.separator_direction_slope
            - self.coefficient
        )

    def audit(self) -> bool:
        A = self.pell_A
        U = self.pell_U
        m = self.direction_norm
        return (
            m > 0
            and self.coefficient.sign(0) > 0
            and (A - SILVER_ONE).sign(0) > 0
            and U.sign(0) > 0
            and A * A - self.coefficient * m * U * U == SILVER_ONE
            and self.direction_from_direction.is_integral
            and self.spatial_from_time.is_integral
            and self.separator_height_bound == A * A - SILVER_ONE
            and self.witness_separator_coefficient.sign(0) > 0
        )


CRITERION_THEOREM = (
    "Let M=<-a>+L, d^2=m, and suppose A,U are integral, A>1, U>0, "
    "A^2-a*m*U^2=1, a*U is integral, and (A-1)/m is integral. Then "
    "g(x,y)=(A*x+U*(d,y), y+(a*U*x+((A-1)/m)*(d,y))*d) is an integral "
    "isometry with integral inverse obtained by reversing the two mixed "
    "signs. If r=(x,y), x>0, separates v0=(1,0) from g*v0, then "
    "(d,y)>=(a*A/(a*U))*x. Cauchy and r^2=q give H=a*x^2/q<=A^2-1. "
    "For z=(1,t*d), the same inequality gives "
    "(z,r)>=(t*a*A/(a*U)-a)*x. Thus, if this coefficient is positive and "
    "a complete strict-ray certificate excludes every mirror through "
    "H=A^2-1, no mirror separates v0 and g*v0. If v0 is regular, g "
    "stabilizes the complete reflection chamber C. Since A>1 makes g an "
    "infinite-order hyperbolic translation, its powers give distinct cosets "
    "of the full reflection subgroup W in the integral isometry group Gamma: "
    "W is normal because Gamma permutes all reflective mirrors, while the "
    "W-stabilizer of C is trivial. Thus [Gamma:W] is infinite. Whenever the "
    "host is admissible, Gamma is an arithmetic lattice, so C has infinite "
    "volume."
)


__all__ = ["CRITERION_THEOREM", "PellBoostCriterion"]
