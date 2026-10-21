# Whole A18: an exact isolated-cap domain and a necessary time-norm screen

21 October2026. Algorithm specification and supporting mathematical
argument by root, pending outside-author review and bounded execution.
There is no computed output or whole-A18 exclusion asserted in this
source. This is supporting material for the existing arithmetic paper.

## Mathematical domain

Use the accepted fresh-rank theorem
`results/arithmetic_whole_a18_point52_actual_field_degree_reduction_20261021.md`,
SHA256 `b57aa75adbde5c1a1b6f78354c58bba47d03b2047a923438b133be9aee5f2707`.
At a whole-A18 farthest vertex in compact arithmetic first-type
dimension18, every .52-weak isolated dotted opposite has primitive
integral beta=4gamma² in the ACTUAL field K and

    4<beta<40/3,    0<sigma(beta)<delta=19/35

at every other embedding. The accepted source gives degree at most3.
To exclude degree1, use the standard rational first-type compactness
obstruction stated at lines1050–1052 of
`manuscripts/compact_arithmetic_dimension_bound.tex` (SHA256
`1f09df2976085fd1ad1aa56ed35f1ab182ff95a8332dcfa0b0b80d4849b82203`).
The full old-plus-isolated-cap Lorentzian Gram has entries in Q(beta),
so strict conjugate positivity forces ACTUAL K=Q(beta). If beta is
rational, its rational indefinite ambient form in19 variables is
isotropic by Meyer's theorem (indefinite rational forms in at least
five variables are isotropic), incompatible with compactness of its
first-type arithmetic quotient. This is an explicit external rational
isotropy input, not part of the coefficient enumeration. Hence degree
is2 or3. In degree3 every weak
central cap is isolated dotted; in degree2 other finite/helper profiles
must be handled separately. We do not claim to enumerate those here.

For degree3 the accepted source proves the complete coefficient box

    p=X³-TX²+SX-M, T=5..14, S=1..14, M=1..3.

This contains420 integer triples. For degree2 write p=X²-TX+M.
The same strict embedding bounds give

    4<T<40/3+19/35<14,    0<M<(40/3)(19/35)<8.

Thus T=5..13 and M=1..7 is a complete63-pair box. The two boxes
are necessary domains for actual weak caps, not existence domains.

Keep exactly the irreducible polynomials with d-1 roots in(0,19/35)
and one root in(4,40/3). For every retained polynomial consider
every index i=5..14, put k=i(19-i), and retain its index cell if
the d-1 small roots are all in(0,38/k). This is the complete
conjugate Schur condition for the full old-A18-plus-isolated-cap
Gram. No actual weak cap can be lost by these tests. A retained cell
need not satisfy every actual-wall, height or compatibility condition.

## Exact root-count and arithmetic algorithms

All polynomials and evaluations use rational coefficients. A monic
quadratic is irreducible here exactly when its discriminant is positive
and nonsquare. A monic cubic is irreducible exactly when none of the
positive or negative integer divisors of M is a root. A cubic with three
distinct real roots has positive discriminant; the root-count test
establishes the required locations rather than relying on that sign.

Use the Sturm sequence p,p', followed by negative Euclidean remainders.
Discard zero final remainders; no division by zero occurs. Let V(a) be
the number of sign changes after zero entries are removed from the
sequence evaluated at rational a. Then V(a)-V(b) counts roots in(a,b)
when p(a)p(b) is nonzero. All the endpoints above are rational and
the retained degree>=2 polynomials are irreducible, so p never vanishes
there. We supply the standard Sturm theorem as an external univariate
algebra fact; the worker implements only rational Euclidean arithmetic.

Each remainder has smaller degree. All outer coefficient boxes and
index loops are finite. Thus this root-count algorithm terminates.
The rational-root test, Sturm theorem and unique coefficient tuple
prove that it enumerates the stated algebraic domain completely, once
per minimal polynomial and once per retained index. It does not identify
different generators of an isomorphic field as the same polynomial.

## The raw time and its independent norm check

Normalize the old A18 roots to square2 over K. The isolated root can
be scaled OVER K to norm2beta with sole old product-beta: if its
original product is a, scale by-beta/a, using
a²=(beta/2)*norm(original root). No square root is adjoined.
If Omega_i is the raw old dual, Omega_i²=k/19. Its raw projection
is beta*Omega_i, hence its time square is

    F_(p,i)=beta*(k*beta-38)/19.                       (1)

It is positive physically and negative at every other embedding.
Its field norm is therefore positive in degree3 and negative in
degree2. Exact resultant evaluation gives

    degree2: Norm(F)=M*(M k²-38T k+1444)/19²;
    degree3: Norm(F)=M*(M k³-38S k²+1444T k-54872)/19³. (2)

The worker independently constructs (1) in Q[X]/(p), computes the
matrix of multiplication by F in the basis1,beta,...,beta^(d-1),
and takes its determinant. It requires exact equality with (2).
This check tests the displayed norm formula, not a floating-point
approximation of its physical root.

## Necessary pair groups

If two such actual caps coexist, their fields have the same degree.
The generator discriminants differ from their common FIELD discriminant
by integer-index squares. Thus their polynomial discriminants have the
same positive rational square class. This is necessary only: we do not
infer field isomorphism from an equal discriminant square class.

The orthogonal complement over K of the FULL old rank18 span is one
dimensional. Its time coordinates for two K-valued actual normals have
a ratio in K. Consequently F_1/F_2 is a K-square and
Norm(F_1)/Norm(F_2) is a positive rational square. For a fixed degree
the denominator19^d cancels in (2), so the absolute integer numerators
must have the same squarefree part.

Group index cells by

    (degree, squarefree_part(disc(p)),
             squarefree_part(abs(numerator Norm(F)))).           (3)

Every actual pair must belong to one group. Different groups exclude
coexistence. Membership in a group asserts ONLY these two necessary
norm/discriminant conditions, not a field isomorphism, a K-square,
actual facet compatibility, or a geometric realization.

Squarefree extraction uses trial division with exact integer remainder.
The worker records a witness n=d*s² for each extracted part. Trial
division terminates when the divisor squared exceeds the remaining
factor. This invariant and the final prime factor prove its output.
All integers are bounded by the displayed finite coefficient/index boxes.
There is no external algebra package or finite field list.

## A direct obstruction to a reflected same-polynomial pair

Suppose the two actual caps have the SAME polynomial p and reflected
indices i,j=19-i. Each parameter is the unique root of p above4, so
their physical beta values agree. Their positive time components agree
by (1). Their normalized raw product is therefore

    c_i.c_j=2beta-beta²*i(19-2i)/19.                  (4)

The Cartan cycle through c_i, the old chain from i to j, and c_j
has product plus or minus c_i.c_j: the old chain entries are-1,
the two isolated contacts are-beta, the new norms are2beta,
and the factors cancel exactly. Arithmetic cyclic integrality makes
(4) an algebraic integer. Since2beta is integral,
beta²*i(19-2i)/19 must also be integral. But its norm has absolute value

    M²*abs(i(19-2i))^d/19^d,

which is not an integer: i=5..9 gives the nonmultiples45,42,35,24,9,
and M is at most7 in degree2 or3 in degree3. Compact arithmeticity
therefore excludes the reflected pair. This uses the normalization
and actual time equality, not an assumption that the roots are native.

If every computed degree3 group contains only ONE polynomial and ONE
value k, every two distinct weak cubic indices would be such a reflected
pair. In that EVENT the computation and (4) prove at most one weak cubic
index. The worker reports this predicate rather than assuming it true.
It asserts no energy consequence at an omitted endpoint. Other quadratic
weak types and any mixed groups still require separate arguments.

## Output and reproduction contract

The worker `scripts/check_a18_isolated_cap_domain_and_norm_screen_20261021.py`
will record the complete kept polynomials, every index cell and all groups,
including groups with several polynomial/k types. It also records counts
for all420 and/or63 visited input rows, exact rejection counts, positive
squarefree witnesses and the norm-formula equality checks. There is no
expected emptiness or geometric conclusion hardcoded into the output.

Run only through the existing one-thread59-second/1GiB guarded runner,
with its unchanged40GiB disk/6GiB RAM start and8GiB/3GiB abort gates.
The proof, accepted domain source, rational-isotropy source and worker are hash-bound before and
after execution. The destination and execution receipt must be new.
Only compact JSON and the actual guarded receipt are retained; no
temporary build products or large algebra packages are used.

No computed output is attached to this pre-execution specification.
A separate source-bound execution record supplies actual counts. If
the watchdog stops the worker, that is a failed or deferred execution,
not a partial-domain proof or a reason to reduce the admission gates.
