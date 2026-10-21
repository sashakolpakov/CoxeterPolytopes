# Squarefree-ideal exclusion of the long-endpoint plus e3 rows

19 October 2026. Argument proposed by the root agent and independently
checked by batch_guard. Paper/source only; no enumeration, scientific
job, Lean or build was run. This is a local necessary-condition
obstruction, not a refund or a dimension-bound claim.

The exact row interfaces come from
`arithmetic_rank20_isolated_long_dotted_height_bounds_20261019.md`,
SHA256 `9fc0583124cba17afed393bb569f83835e8f2c4be8ce5cfa64839f0369254516`.
In particular d remains isolated, c is orthogonal to the spectator,
and this note concerns ONLY a simple long-endpoint contact together
with the e3 contact. It does not cover the endpoint3-only row.

## Algebraic statement

Let K be any totally real number field, 8<=N<=14, and M=N+2.
Suppose A,delta belong to O_K and satisfy the following necessary
conditions of the actual quadrilateral:

* physically A>=2;
* in the far row C=A(NA+M), and C*delta=M;
* in the near row C=(A+2)(NA+N-2), and C*delta=M;
* at every nonphysical embedding,
  Z_i^2<2M, where Z=2NA+M in the far row and
  Z=2NA+3N-2 in the near row.

If the principal ideal (M) is squarefree in O_K, neither row is
possible. Neither the degree nor the eight-field list is needed
for this algebraic statement.

## Proof

Set U=A in the far row, U=A+2 in the near row. In both cases write

    C=U*B,
    B=NU+M  (far),
    B=NU-M  (near).

All these elements are integral. Since C*delta=M is nonzero, no
factor is zero, and the principal ideal (C) divides (M).

Suppose a nonzero prime ideal p divided (U). Because (U) divides
(C), it would divide (M) as well. But B=NU+/-M then belongs to p,
so p divides both factors U and B. Consequently v_p(C)>=2,
contrary to v_p(C)<=v_p(M)=1. Thus no prime ideal divides (U):
U is a unit.

It follows that B divides M. It also divides NU=B-/+M, so because
U is a unit it divides N. Rational Bezout then gives

    (B) divides (gcd(N,M)) divides (2).

Hence abs(Norm(B))<=2^[K:Q]. We now contradict that bound at
EVERY embedding, without assuming conjugates of U have a sign.

Physically,

    B=NA+M>=3N+2>2                  (far),
    B=N(A+2)-M>=3N-2>2             (near).

At a nonphysical embedding, B=(Z+M)/2 in the far case and
B=(Z-M)/2 in the near case. Therefore

    abs(B_i)>M/2-sqrt(M/2)>2.

The last inequality holds for M>=10: with t=M/2>=5,
(t-2)^2-t=(t-1)(t-4)>0 and t-2>0. Multiplying the strict
absolute-value inequalities over all embeddings yields
abs(Norm(B))>2^[K:Q], a contradiction.

The number-theoretic input here is ordinary prime-ideal factorization
and valuations in the full integer ring of a number field. There is
no assumption that a selected power order is already the full ring.

## Conservative coverage within the eight fields

For a squarefree rational integer M, if every prime dividing M is
unramified in K, then (M) is a squarefree ideal. We use only this
sufficient test, not its converse. The external ramification inputs
are the quadratic discriminants12,24 for Q(sqrt3),Q(sqrt6), and
the fact that a maximal-real subfield of Q(zeta_m) is unramified
away from primes dividing m.

|N|M|fields covered by this sufficient test|
|---|---|---|
|8|10|K9,K11|
|9|11|Q(sqrt3),Q(sqrt6),K9,K15,K16,K20,K24|
|11|13|all eight|
|12|14|K9,K11,K15|
|13|15|K11,K16|

Here Km=Q(2cos(2pi/m)), and K24=Q(sqrt2,sqrt3). The relevant
possible ramification-prime sets are {2,3} for the quadratic fields,
{3} for K9, {11} for K11, {3,5} for K15, {2} for K16,
{2,5} for K20, and {2,3} for K24. These give exactly the
conservative entries above. M=12 and16 at N=10 and14 are not
squarefree rational integers, so this sufficient test says nothing
there. Missing entries are not claimed to be realizable or to have
escaped other obstructions.

This theorem supplies named analytic omissions for a future exact
enumeration manifest. The separately prepared K9,N=9,endpoint3-far
batch concerns a different contact type and is unchanged. No new
simultaneous consumer payment, rank20 allowance or compact arithmetic
dimension ceiling follows without the remaining geometric composition.
