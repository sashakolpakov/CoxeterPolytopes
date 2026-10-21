# Dotted long contacts with an isolated opposite wall

19 October 2026. New root-agent hand arguments, independently checked
and written out by batch_guard. Paper only: no enumeration, scientific
job, calibration, Lean, TeX or build was run. All earlier source-bound
proofs and worker sources are preserved unchanged.

The new arguments below exclude both endpoint3+e3 rows without a
squarefree-ideal hypothesis, and close the last quadratic far
endpoint3-only possibility. Together with the named earlier norm
lemmas and the explicit E8 check, this closes the isolated-d,
spectator-orthogonal D_r long-contact inventory for dotted x.
This is a local theorem, not a global refund or dimension bound.

## 1. Actual-quad scope and inventory

Use the old entire finite vertex D_r+T_s of the ambient20 frontier,
r+s=20,10<=r<=16. Let a be the D_r fork, b,e its short leaves,
and l_1,...,l_L its long arm, with N=r-2, L=N-1. The complete
white support of the actual quadrilateral a,b,c,d is

    W=A1(e)+A_L+T_s.

The hypotheses in this note are:

* d contacts ONLY b among the old vertex roots, with
  beta=-b.d>1; in particular it is orthogonal to a,W;
* c is orthogonal to b,d and the ENTIRE spectator T_s;
* c contacts at least one long-arm root;
* x=-a.c>=1 is an EXPLICIT hypothesis;
* the inherited actual-cycle integrality, corank, and nonphysical
  proper-block positivity interfaces hold in the full actual
  totally real ground field K.

In particular the ambient products of arbitrary opposite walls are
not presumed dotted. The entire spectator remains orthogonal, not
an omitted subdiagram. Compactness gives positive-definite W+c
at the actual finite quad corners.

The complete finite W+c long-contact inventory, using the finite
spherical Coxeter diagram classification, is:

* one simple contact at k=1,2,N-2,N-1, with no e contact;
* for N=8 only, one simple contact at k=3 or5, with no e contact;
* one endpoint4 contact, with no e contact;
* one simple endpoint contact together with e3 or e4.

This is the explicit inventory in
`arithmetic_rank20_isolated_long_finite_product_20261019.md`, SHA256
`9644c76993b26d97803686787c5a87c284051184aa08cd8d2d7d6250a06e1ece`.
Two long contacts create a cycle; a4 bond is terminal in these
finite path ranks; adding e to a next-to-end fork or to a long4
bond is not finite. The two exceptional A7-to-E8 contacts at k=3,5
are retained, not suppressed by a classical asymptotic inventory.

The displayed arithmetic formulas below are those of
`arithmetic_rank20_isolated_long_dotted_height_bounds_20261019.md`,
SHA256 `9fc0583124cba17afed393bb569f83835e8f2c4be8ce5cfa64839f0369254516`.
The new algebraic arguments work for any N>=8 once those interfaces
hold. The current geometric application is the stated rank20 range.
No field list or explicit integral basis is used in this closure.

## 2. A product-norm fact

If P,Q are nonzero elements of O_K and

    abs(Norm(PQ))<2,

then both are units. Indeed the two norms are nonzero integers,
whose product has absolute value1. Each norm is therefore +/-1,
which is the ordinary unit criterion for algebraic integers.
No classification of units, power order, or primitivity is involved.
All norms below use EVERY embedding of the full ground field K,
with multiplicities if a selected element lies in a proper subfield.

## 3. Both endpoint3+e3 rows over every totally real field

Put M=N+2>=10 and A=2x in O_K. The actual simple-long cycle gives
this integrality. The opposite product gives
delta=4beta^2-4 in O_K. The two rows are

    far:  Z=2NA+M,
          C=A(NA+M)=(Z^2-M^2)/(4N);

    near: U=A+2,
          Z=2NU-M,
          C=U(NU-M)=(Z^2-M^2)/(4N).

In both cases C*delta=M and every nonphysical Z_i satisfies
Z_i^2<2M. Hence every nonphysical C_i<-M/4.

For degree d=1, physically C>=6N+4 in the far case and
C>=12N-8 in the near case. Both exceed M, contradicting the
positive integral identity C*delta=M. Suppose henceforth d>=2.
The nonzero integral norm of delta gives the strict bound

    0<C_0<4^(d-1)*M.                               (2)

### Far row: A becomes a unit

Physically A_0>=2. Nonphysically,

    A_i=(-M+Z_i)/(2N),       -1<A_i<0.

The lower bound uses M+sqrt(2M)<2N, which follows from
(N-2)^2-(2N+4)=N(N-6)>0 for N>=8. The upper bound uses M>2.
Thus every abs(A_i(A_i+1))<=1/4. Physically,

    C_0=N*A_0*(A_0+1)+2A_0>N*A_0*(A_0+1).

Combining with(2),

    abs(Norm(A(A+1)))
       <[4^(d-1)*M/N]*(1/4)^(d-1)
       =M/N<=5/4<2.

The factors A,A+1 are nonzero, so both are units by Section2.
Set B=NA+M. Since C*delta=A*B*delta=M, B divides M.
It also divides B-M=NA; because A is a unit, B divides N.
Rational Bezout gives

    B divides gcd(N,M)=gcd(N,2), hence B divides2.   (3)

Yet B_0>=3N+2>2, while every nonphysical conjugate satisfies

    B_i=(M+Z_i)/2>M/2-sqrt(M/2)>2.                 (4)

For the last inequality let t=M/2>=5. Then
(t-2)^2-t=(t-1)(t-4)>0 and t-2>0, so t-sqrt(t)>2.
Thus abs(Norm(B))>2^d, contrary to(3).

### Near row: U becomes a unit

Here U_0=A_0+2>=4 and

    U_i=(M+Z_i)/(2N),       0<U_i<1.

The interval follows from the same comparisons just used. Thus
every abs(U_i(U_i-1))<=1/4. Physically,

    C_0=N*U_0*(U_0-1)-2U_0
       >=(N-2/3)*U_0*(U_0-1),

because the difference is (2/3)*U_0*(U_0-4)>=0. Equality at
U_0=4 causes no issue because(2) is strict. Hence

    abs(Norm(U(U-1)))
       <M/(N-2/3)<=15/11<2.

Again both factors are nonzero integral elements; U and U-1 are
units. Put B=NU-M. From U*B*delta=M we get B divides M;
from B+M=NU and the unit U we get B divides N. Thus B divides2.
Physically B_0>=3N-2>2, whereas

    B_i=(Z_i-M)/2<-M/2+sqrt(M/2)<-2

at every other embedding. The same norm contradiction excludes
the near row.

Both endpoint3+e3 dotted rows are therefore impossible over ANY
totally real K under the stated interfaces. The earlier
squarefree-M hypothesis was a sufficient device for forcing a
unit; the product norm proves that unit property without it.

## 4. Closing the far quadratic endpoint3-only possibility

The earlier theorem
`arithmetic_rank20_endpoint3_basis_free_norm_obstruction_20261019.md`,
SHA256 `03809292c144e54ed85ddbdf91e0796426116641e9b00bb34a44222f89f02f39`,
excludes endpoint3 alone in every degree and role except the
necessary far degree2 possibility K=Q(sqrt2), A_0=1+sqrt2.
It does not claim that this last possibility actually occurs.

Substitute it in the far endpoint3-only identity. One gets

    C=2NA^2+4A-N-3
       =5N+1+4(N+1)*sqrt2,
    Norm(C)=-(7N^2+54N+31)=-D.

Since C*delta=4(N+1) with integral delta, the positive number

    abs(Norm(delta))=16(N+1)^2/D                   (5)

must be an integer. For N>=8,

    16(N+1)^2-D=(N-3)(9N+5)>0,
    D>7(N+1)^2.

Thus(5) lies strictly between1 and16/7<3. Its only possible
integer value is2. Equality to2 would imply

    N^2-38N-23=0,
    (N-19)^2=384.

But 19^2<384<20^2, so no integer N solves this equation. The
last far possibility is excluded. Endpoint3 alone is now
impossible in BOTH roles over every totally real ground field,
without an explicit integral-basis assumption.

## 5. Exceptional E8 contacts retained explicitly

At N=8 the simple k=3,5 contacts in the inventory have
kappa=5/16,3/16 respectively and t=1/16. Put A=2x in O_K
and let k=5 or3 accordingly. Their exact identities are

    Z=8A+k,
    C=(Z^2-5)/4=16A^2+4kA+(k^2-5)/4,
    C*delta=4,
    Z_i^2<1 at every nonphysical embedding.

C is integral and congruent to1 modulo2. Its principal ideal is
therefore coprime to(2); since C divides4, it is a unit. But
A_0>=2 makes C_0>1, and -5/4<C_i<-1 at every other embedding.
Thus abs(Norm(C))>1, contradicting that it is a unit. This
argument also works in degree1 and needs no field list.

## 6. Complete local long-contact composition

The remaining three classical types are excluded over every
totally real ground field by
`arithmetic_rank20_three_dotted_long_rows_norm_obstruction_20261019.md`,
SHA256 `76c2a552bf8ab1d898469bf06bd4ba0f8dc07a0222532205fe6abd5c78472cf0`.
Their actual integral variables reduce C*delta=8 to the far form
NV^2+4V-2 or the near form NV^2-4V-2:

|type|actual integral variable before the near shift|
|---|---|
|next-to-end3|2x; near shift+1|
|endpoint4|2sqrt2*x; near shift+2|
|endpoint3+e4|2(2x)+sqrt2; near shift+2|

That proof keeps the actual cycle distinctions: endpoint4 does
not separately require sqrt2 in K; in the e4 row the two actual
cycles prove sqrt2 in K and hence its integrality.

The inventory in Section1 is now covered without omission:
endpoint3-only by Section4 and its named earlier theorem;
next3, endpoint4, endpoint3+e4 by the common norm lemma;
endpoint3+e3 by Section3; and the two exceptional E8 contacts
by Section5. Therefore the actual quad in Section1 is impossible
whenever x>=1, c contacts the long arm, d is isolated as stated,
and c is orthogonal to the entire spectator. The exclusion is
independent of the eight-field list and of explicit integral bases.

The still-open obligations are not hidden in this composition:
c contacting the spectator, other nonisolated-d attachment types,
zero or finite x outside their separate proofs, all required
triangle consumers, and simultaneous modified-angle payment across
all whole-vertex types still need their own complete arguments and
reviews. This local closure alone does not authorize a new global
refund or lower a dimension ceiling.
