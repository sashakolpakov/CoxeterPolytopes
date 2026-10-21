# Complete arithmetic distance domains for the new B1 dotted cores

The [exact source](../scripts/fresh_sqrt5_global_b1_dotted_parameters_20261016.py)
and [certificate](fresh_sqrt5_global_b1_dotted_parameters_20261016.json) cover
all 82 minimal baseline-one dotted classes in the inherited complete
[adjacent-discount frontier](fresh_sqrt5_adjacent_discount_frontier_20261016.md).
There are 81 one-dot classes and one two-dot class. Four one-dot classes fail
a conjugate principal restriction. The remaining domains contain 377 necessary
arithmetic parameter seeds: 335 nonsingular and 42 singular. This is a finite
exact distance enumeration, not a sample of distances or an attachment closure.
Generation and a fresh saved replay passed. A separate
[arithmetic review](fresh_sqrt5_global_b1_dotted_parameters_20261016_review.md)
also reconstructs every parameter domain.

Write phi=(1+sqrt(5))/2. For either dotted edge, its two-edge cyclic invariant
is t=4cosh(d)^2. Arithmeticity requires t in Z[phi], t>4, and 0<=t'<=4 at the
other real embedding. The field isomorphism makes the endpoints t'=0,4
impossible for a dotted parameter, but retaining them in the enumeration box
does no harm.

For a one-dot core, deleting both endpoints leaves a nonsingular Lorentzian
block. Eliminating that block gives a residual matrix with diagonal entries
a,b and off-diagonal -x-c, where x is the magnitude of the missing Gram entry.
The unique negative direction is already in the eliminated block, so the
residual must be positive semidefinite: (x+c)^2<=ab. Hence

    t <= 8(ab+c^2)/(d0*d2).

When the finite-edge graph disconnects the dotted endpoints, c=0 and the
bound sharpens to 4ab/(d0*d2). The coefficient data and principal restrictions
are saved for every class, including nonsingular cores without a field-valued
zero-further-contact completion. Those cores must not be discarded merely
because a singular stopping value is absent.

The two-dot core has four black roots, no white roots, and only one finite
nonorthogonal bond, a 3-edge between roots 2 and 3. Its dotted pairs are
(0,2) and (1,3). Let their invariants be t and u. Eliminating roots 0 and 1
from the Gram matrix gives the exact necessary inequalities

    (t-4)(u-4) <= 4,
    (4-t')(4-u') >= 4.

The smallest algebraic integer t>4 with 0<=t'<=4 is 3+phi. Indeed, if
t=A+B*phi then t-t'=B*sqrt(5) forces B>0; B>=3 is already too large to
improve that minimum, and the cases B=1,2 are finite. Thus each parameter
is at most 4+4phi. There are twelve field integers in this box. All 144 ordered
pairs are checked, leaving seventeen arithmetic seeds, four singular.
No requirement that sqrt(t) or sqrt(u) lie in the ground field is imposed:
a tree scaling gives a field-valued Gram matrix without that restriction.

All accepted seeds pass exact physical and conjugate inertia checks and
integrality of every cyclic product. Every singular seed has the alternating
black kernel signs of a quadrilateral. The two-dot inertia decisions are
also checked against the two scalar product inequalities above.

The enumeration uses exact integer intervals: B*sqrt(5)<=upper_bound and
-B*phi'<=A<=4-B*phi'. Both floor and ceiling are obtained with rational
bounds and exact field comparisons. There are no floating-point gates.

Further touching white roots, isomorphisms with dotted parameters, all
singular stopping choices, and common spectator payments still need their
own complete checks. No improved dimension bound or compact ambient
realization is claimed.

```sh
python3 scripts/fresh_sqrt5_global_b1_dotted_parameters_20261016.py --verify
```
