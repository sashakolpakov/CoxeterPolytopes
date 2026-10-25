# Cubic B14+I: even terminal norm and exclusion of foreign angle pi/4

23 October 2026. **Complete scoped written candidate, pending a full
different-author review and root acceptance.** Root derived the new cubic
interval obstruction, the norm-two cutoff and their application below.
Replay supplied the preceding terminal setup and identified odd norm as
forcing V/2 integral. Those inherited contributions are disclosed; they
are not an independent review of this new proof. No worker, enumeration,
symbolic calculation, test, Lean or TeX process supplied the argument.
Every terminating decimal denotes an exact rational number.

## Statement and exact boundary

Consider the original farthest vertex with full forest B14 plus an actual
short isolated component I in a compact arithmetic Coxeter 15-polytope.
Assume that the actual field K has degree three and that the actual
opposite of the short B14 wall is isolated within B14 and has dotted
own contact. Assume also that the head has one of the two residual
types described in H below: pure dotted, or native-three contact to old
I with zero mutual product to its opposite.

Write B=4 gamma^2 for the short cap's own parameter and Q for its
squared contact with old I. Then the following necessary conclusions
hold:

* N_K(B) is even, and in fact is either 2 or 4.
* The short cap's old-I contact is either right-angled or native three;
  native four and every higher finite label are excluded.
* If N_K(B)=2, then B>31 in the physical embedding.

This does not exclude either remaining terminal norm or both residual
head families. In particular, it neither deletes group19 nor proves a
global dimension14 bound. No polytope is asserted to realize the data.

The literal interfaces, all under results/, are:

|ID|Source|SHA256|Used interface|
|---|---|---|---|
|R|arithmetic15_b14_actual_i_strong_five_and_complete_nonisolated_head_exclusion_candidate_20261023.md|9646d72541fb450b4c3346ca9b9dfdc0edf4773369cdbcc799d3a29a83863e40|The original strong5 through12 floors and full-I cyclic cancellation|
|RA|review_arithmetic15_b14_actual_i_strong_five_nonisolated_head_root_full_source_acceptance_20261023.md|666bebcc762e6e295df9ebf2db7c79f25047b5f5bccccfb024db6f7f590cd9b2|Accepted complete scope of R|
|H|arithmetic15_b14_actual_i_all_finite_pure_heads_and_positive_foreign_head_exclusion_candidate_20261023.md|6328c4d018d7a7e1da01863d946071f060e43738d4efce5bb17caff9e1d0858b|The two stated residual head hypotheses, without importing conditional strong13|
|HA|review_arithmetic15_b14_all_finite_pure_heads_positive_foreign_root_acceptance_20261023.md|dbd62d35c1c4fe3e24325752f5d913dc7a9f7161ab511907513af2a931fb5d49|Accepted complete head reduction|
|T|arithmetic15_b14_actual_i_complete_short_white_and_finite_own_isolated_dotted_reduction_candidate_20261023.md|4ad318cf28f7bccd9531de153bf63d064b16aa3b4cda5a85a42604e6104d8ca4|Only Section2's explicit cleared-head point, prices and short-length bound, reproduced below; its whole terminal classification is not presumed accepted at this freeze|
|S|classical_dotted_head_scalar_295_and_cubic_norm_eight_pressure_candidate_20261023.md|1c28409b9e14f51cc150eed1cf9e6a3bd43c9392d92f63f9a8d653140fcde307|The displayed scalar F(u)>.295 on0<u<1|
|SA|review_classical_dotted_head_295_and_cubic_norm_eight_root_acceptance_20261023.md|c00c37c9566ba8cabb935388103e21c012cbb9bde52aade08318265a4063683e|Accepted scalar scope only|

The usual arithmetic reflection-form model, integral cyclic products,
positive definiteness at every nonphysical embedding, finite native
angles and the original farthest-point comparison are explicit external
foundations. No face is assumed to inherit arithmeticity. No field list,
restricted suborder, principal-ideal hypothesis or disconnected equality
of root lengths is used.

## 1. A cubic integer cannot fit the I half-strip

**Lemma.** There is no algebraic integer of degree three whose conjugates
are v,-a,-b with 0<v<1 and 0<a,b<2.

Its trace is an integer -t, with t in {0,1,2,3}. Thus a+b=v+t.
Its norm N=vab is a positive integer, strictly less than4. Its second
elementary symmetric coefficient is the integer

    S=ab-v(a+b)<(v+t)^2/4-v(v+t)<t^2/4.

The monic polynomial is

    p(X)=X^3+t X^2+S X-N.

For t=0 or1, a+b=v+t<2, so ab<1 and N<1, impossible.
For t=2, N<v(v+2)^2/4<9/4, hence N is1 or2, and S<=0.
Since all three roots are greater than-2, p(-2)<0. At t=2 this
is -2S-N<0, which forces S>=0 in both norm cases. Thus S=0.
At N=1 the polynomial factors as (X+1)(X^2+X-1), contradicting
degree three. At N=2 its discriminant is-44, contradicting three
real roots.

For t=3 we have S<=2. Now p(-2)=4-2S-N<0. The complete
integer possibilities and their exact discriminants are:

|N|S|Discriminant of X^3+3X^2+S X-N|
|---:|---:|---:|
|1|2|-23|
|2|2|-104|
|3|1|-76|
|3|2|-239|

They are all negative. The table is exhaustive from N in{1,2,3},
S<=2 and 2S>4-N; no coefficient enumeration is left to a program.
The discriminants follow from

    t^2 S^2-4S^3+4t^3 N-27N^2-18tSN.

This proves the lemma using only integer coefficients and the stated
root intervals.

## 2. The actual terminal norm must be even

Use the original unit old-I normal e and its actual opposite

    d=sqrt(h)p-sqrt(1+h)e,
    V=4h in O_K, 0<V<2, -4<V_sigma<0.

If V/2 were integral, it could not be rational: a rational algebraic
integer strictly between0 and1 is impossible. Since [K:Q]=3, it would
therefore have degree three. Its conjugates would be in the forbidden
intervals of Section1. Consequently

    V/2 is NOT in O_K.                                      (1)

The isolated short cap has unit form

    c=A p-gamma P14-beta e,
    B=4gamma^2 in O_K, Q=beta^2,
    A0^2=A^2-Q=14gamma^2-1=(7B-2)/2.

All nonphysical embeddings satisfy 0<B_sigma<2/7. This follows from
the original positive-definite B14/cap principal and its Schur complement;
the external beta coordinate is not silently removed from the full form.

For completeness, let lambda=-c.d. The original two/three cycles give
the integral elements

    L=4lambda^2, W=4beta^2,
    T0=8beta sqrt(1+h)lambda.

The exact relation lambda=A sqrt(h)-beta sqrt(1+h) gives

    V A0^2=L+T0+W in O_K.                              (2)

This retains every finite old-I contact, including beta>0. No individual
Cartan entry is inferred to be integral from unit normals.

Suppose N_K(B) were odd. Multiplication by B on the free integer
module O_K has determinant N_K(B), so its reduction modulo2 is
invertible. The element 7B-2 has the same reduction as B, and is
therefore a unit of O_K/2O_K. Equation(2) says

    V(7B-2) is in 2O_K.

Multiplying by the inverse residue class gives V in2O_K, contradicting
(1). Thus N_K(B) is even. This uses the full maximal order as a free
integer module; it does not require any splitting classification of2.

## 3. A norm-two cubic terminal has B>31

**Lemma.** If B,a,b are the three conjugates of a cubic algebraic
integer, with B>4, 0<a,b<2/7 and Bab=2, then B>31.

Suppose B<=31. Since ab<4/49, B>49/2. Also ab=2/B>=2/31>1/16,
so a+b>=2sqrt(ab)>1/2. The integer trace T=B+a+b therefore
satisfies 26<=T<=31. Write the monic polynomial as

    X^3-T X^2+S X-2.

Because both small roots are below2/7 and B is above it, evaluating
at2/7 gives

    S<(14T+339)/49.

On the other hand B>T-4/7, so

    S=B(a+b)+2/B>2sqrt(2(T-4/7)).

The complete six-trace contradiction is:

|T|Largest permitted integer S|7[8(T-4/7)-S_max^2]|
|---:|---:|---:|
|26|14|52|
|27|14|108|
|28|14|164|
|29|15|17|
|30|15|73|
|31|15|129|

Each positive entry says the strict lower bound exceeds S_max. This
exhausts all traces under B<=31 and proves the lemma.

## 4. The one-vector geometric estimate and the norm range

Here are the exact local geometric statements used from T Section2.
Keep ONE original farthest vector, written x in B14 and y in the
orthogonal old-I direction, with ||x||^2+y^2<1. Let C_H,C_S be the
actual head and short edge cosh lengths. The strong5 through12 floors
give a positive prefix reference with E12>.86989, first-five mean
b<.361 and separate available P5 and P12 prices greater than.14.

For either residual head write Q_H=0 or1/4, U=4g_H^2, and
u=(1-Q_H)/g_H^2 in(0,1). The positive actual-wall combination

    xi=c_H+beta_H e+((1-Q_H)/g_H)n1

is future, has square -g_H^2 u(1-u), and spatial part
-sqrt2 g_H Z with Z=(1-u/2)e1+(u/2)e2. Every original wall has
nonpositive product with xi. Hence the SAME original farthest comparison
gives x.Z>=H(u), where

    H(u)=[sqrt(2-u)-sqrt(u(1-u))]/sqrt2,
    R(u)=||Z||^2=1-u+u^2/2.

By S/SA, .8H-.16R>.295. The single legal residual .4(Z-P5/5)
spends.08 of P5 and has zero cross with the original reference. Thus
the rebuilt B14 estimate is

    e_B>.86989+.295-.8*.361+.032=.90809.

Append.14(P14-P12). It spends the separate P12 price, has square2
and is orthogonal to that whole rebuilt reference. With

    T_S=sqrt14 sqrt((C_S-1)/(C_S+1)), f12<3.034,

the ONE vector satisfies

    ||x||^2+y^2>
    .90809+.28(T_S-f12)-.0392+q(h),
    q(h)=(sqrt(1+h)-1)/(sqrt(1+h)+1).                (3)

This is one funded reference and one orthogonal I coordinate, not an
addition of independently optimized head/body/terminal energies.
If C_S>=16, then T_S>=sqrt(210/17)>3.514, and(3), even dropping
q(h), exceeds1.00329. Consequently C_S<16.

For an isolated short cap the actual edge identity is

    C_S^2=7B/[2(1-Q)],
    4<B<512(1-Q)/7.

The nonphysical bounds therefore give

    0<N_K(B)<[2048/343](1-Q)<6.

Together with Section2 this leaves precisely N_K(B) in{2,4}.
If the norm is4 then Q<169/512; if the norm is2 then Q<681/1024.
In particular all native labels at least6, whose Q is at least3/4,
are absent. Native5 is also absent: its cycle puts sqrt5 in the
actual degree-three field. Before the final geometric exclusion the
only possible old-I contacts are therefore 0, native3 or native4.

## 5. Native-four old-I contact exceeds the original budget

We first record a lower bound on h shared by both residual head types.
The ungrafted E12 reference with a legal first-coordinate graft at.684
would exceed1 at

    .86989+(5/4)(.684-.361)^2=1.00030125.

Therefore x1<.684. The actual head edge then has C_H<14/5, since
at C_H=14/5 its required x1^2 is9/19>.684^2.

For a pure head, U<392/25 and the actual positive mutual with d
has lambda_H>=1/2. Since lambda_H^2=h(U/2-1),

    h>25/684>1/30.

For the quarter-foreign orthogonal head, C_H^2=2U/3<196/25
and V(U-2)=2. Hence U<294/25 and

    h=1/[2(U-2)]>25/488>1/30.

The rational seal 31*62^2-30*63^2=94>0 shows q(1/30)>.008.

Now suppose Q=1/2, namely native-four contact to old I. The norm
upper bound is less than3, so its only even possibility is2. Section3
gives B>31. Thus C_S^2=7B>217, so C_S>14.7. Consequently

    T_S^2>14*(13.7/15.7)=1918/157>3.49^2.

Use(3) with the SAME rebuilt reference and h>1/30. It gives

    ||x||^2+y^2>
    .90809+.28*(3.49-3.034)-.0392+.008
    =1.00457>1.

This contradicts the original norm budget and excludes the entire cubic
native-four terminal-contact branch, for both residual head types.

## Remaining work

The cubic terminal may still have norm2 or4 and old-I contact0 or
native3. Those cases require additional pair geometry or arithmetic;
neither has been excluded here. Quadratic and quartic fields have
different statements and are not covered by this note. The use of T's
local Section2 is explicit; this source does not turn its pending whole
terminal theorem into an accepted classification. All new statements
still await full independent source review and root acceptance.

No script execution, successful Lean declaration, numerical dimension
improvement, new example, manuscript revision or publication is claimed.
