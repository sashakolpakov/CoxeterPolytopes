# The remaining D18 dotted head has field degree at most six

21 October2026. Root-authored manual supporting candidate, pending full
outside-author review. No mathematical worker, symbolic or numerical
computation, enumeration, Lean, TeX or manuscript change was used. The
argument is a direct algebraic norm bound. It does not exclude the dotted
head, wholeD18 or dimension18.

## Input and result

The sole geometric interface is the accepted source
`results/arithmetic18_whole_d18_long_end_isolation_and_head_field_reduction_20261021.md`,
SHA256 `6ef34bb417647b846c9b26b7c35fb8018e27b999c0395396c33c6d5ba691b3ac`.
For an actual isolated dotted long-end opposite at a whole-D18 farthest
vertex, it supplies the integral own two-cycle U and the ACTUAL field

    K=Q(U), 4<U<28/3,
    0<U_sigma<2 at EVERY nonphysical embedding.          (1)

Only these necessary conditions are used. Elementary algebraic integer
norms, extension of embeddings and differentiation of a cubic are the
other mathematical foundations. No finite field list, discriminant
table, invariant-lattice transfer or polynomial enumeration is assumed.

**Theorem.** Under(1), [K:Q]<=6. If the degree is six, then

    1<=|N(Q(U))|<=5,
    Q(X)=X(2-X)(X-1)^2(2X^2-4X+1),                     (2)

and consequently

    N(U), |N(2-U)|, |N(2U^2-4U+1)| <=5,
    |N(U-1)| in {1,2}.                                 (3)

If |N(U-1)|=2, the three other norm magnitudes in(3) equal1
and |N(Q(U))|=4. These norms are over the actual field K.

## A polynomial small on the entire conjugate interval

Put t=(X-1)^2. The integer polynomial in(2) is exactly

    Q(X)=t(1-t)(2t-1).

For0<X<2 we have0<=t<1. The two nonzero critical points of
g(t)=t(1-t)(2t-1) on[0,1] are

    t=(3-sqrt3)/6, (3+sqrt3)/6.

Indeed g'(t)=-6t^2+6t-1; evaluating at these two points gives
the values -sqrt3/18 and sqrt3/18. At0,1/2,1 its value is zero.
Thus on the FULL interval, not just at a selected set of embeddings,

    |Q(X)| <= c=sqrt3/18 <1.                            (4)

Physically4<U<28/3 implies9<t<625/9. The positive function
t(t-1)(2t-1) increases for t>1, since its derivative is
6t^2-6t+1>0. Therefore

    0<|Q(U)| < A=(625*616*1241)/729
                  =477785000/729 <656000.               (5)

The physical value is nonzero. Accordingly Q(U) is a nonzero algebraic
integer; none of its conjugates can be zero either.

## The integral norm gives the degree bound

Write d=[K:Q]. Equations(1), (4) and(5) give

    1<=|N_K(Q(U))| < A*c^(d-1).                        (6)

If d>=7, c<1 gives an upper bound A*c^6. But

    c^6=1/1259712,
    A*c^6=477785000/918330048 <1.

This contradicts(6). Hence d<=6. The proof does not need a lower
bound on a field discriminant or a bound on individual unit-normal
square roots.

If d=6, then c^5=sqrt3/209952. Using sqrt3<7/4 and A<656000,

    |N_K(Q(U))| <656000*7/839808 <6,

because4592000<5038848. It is a positive integer, so it is at most5.
Multiplicativity gives

    |N(Q(U))|
      =N(U)*|N(2-U)|*N(U-1)^2*|N(2U^2-4U+1)|.        (7)

Every displayed factor is a positive integer. None is zero, since
Q(U) is nonzero. Equation(7) gives all of(3), and if the squared
factor is4, every other factor is1 because a second factor>=2
would make the product>=8. This proves the theorem.

## What is finite, and what is not yet proved

Together with(1), the degree bound gives a finite necessary polynomial
domain: the minimal polynomial is monic of degree at most6, and every
root is positive, one in(4,28/3) and all others in(0,2). Its integer
coefficients are bounded by the elementary symmetric functions of
these intervals. Explicitly, for k=1,...,d, its coefficient magnitude
s_k is a positive integer with

    s_k < (28/3)*binom(d-1,k-1)*2^(k-1)
                         +binom(d-1,k)*2^k.             (8)

This is a finite-domain statement, not an executed polynomial census.
At degree6 the stronger value/norm restrictions(2)--(3) also apply.
No list of fields, possible mirror labels, complete opposite-wall
inventory or new compact polytope is produced by this proof. The
surviving dotted heads still require geometric and arithmetic exclusion.
No accepted guide, checkpoint, PDF, script output or Lean declaration
is changed by this candidate.
