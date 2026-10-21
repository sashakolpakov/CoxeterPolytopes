# Whole D18: the cubic dotted long head is impossible

21 October 2026. New manual supporting proof, pending full root and
outside-author review. No mathematical/scientific/symbolic worker,
enumeration, import, optimizer, Lean or TeX ran. All small polynomial and
discriminant tables below were derived by hand. No frozen paper, source,
guide or receipt is edited. Root supplied the sharpened head-norm product
bound and independently confirmed the stronger two-cycle identity; the
author derived that identity, the local cubic reduction and the final
seven-row norm/field obstruction. Earlier actual-foot/prefix foundations
have shared authorship.

**Theorem.** At an actual farthest vertex of a compact arithmetic
first-type Coxeter18-polytope with WHOLE diagram D18, an isolated dotted
long head cannot have actual field degree three.

The generic two-cycle lemma in Section2 also applies in other degrees,
but their exclusions are NOT claimed here. This is not a whole-D18
exclusion, arithmetic17 theorem, new example or field-list assumption.

## 1. Exact accepted interfaces and conventions

Use the accepted source
`results/arithmetic18_whole_d18_dotted_head_degree_four_and_isolated_forks_20261021.md`,
SHA `0590f2b9d8747f6edbe0cafd74f2909c54c87f62e5ebe7d52081abbdfe35fa16`.
Its source input
`results/arithmetic18_whole_d18_long_end_isolation_and_head_field_reduction_20261021.md`
has SHA `6ef34bb417647b846c9b26b7c35fb8018e27b999c0395396c33c6d5ba691b3ac`.
Only the following actual cubic conclusions and arithmetic foundations
are used:

    K=Q(U), U integral, degree K=3,
    4<U<28/3, 0<U_sigma<2 at both other embeddings;
    U/2 is NOT integral;
    either actual fork opposite is isolated and dotted;
    K=Q(V), V integral, 4<V<676/9,
    0<V_sigma<4/9, 1<=N(V)<=14;
    UV/2 and (2-U)V/2 are integral.                    (1)

The raw old D18 roots alpha_i have square2. Their rational Cartan
inverse has (1,1)=1, fork diagonal9/2, and (1,fork)=1/2.
The actual head root r0 is normalized OVER K to square2U and old
product-U only at1; the actual fork root rf is normalized to square2V
and old product-V only at that fork. Their raw time squares are

    H0^2=U(U-2), Hf^2=V(9V-4)/2.                       (2)

The orthogonal time residuals belong to one K-valued line, so
Hf^2/H0^2 is a square in K. Its rational norm is consequently a rational
square. This is not a norm-two normalization of the fork cap.

First-type cyclic integrality is used for the actual closed head/fork/
old-chain cycle and for the actual head/fork two-cycle. The head and
fork walls ARE distinct: the former has nonzero old contact only at1,
the latter only at its fork. No distinctness of the two fork opposites
is needed, and no pair of independently optimized fork energies is added.

We use elementary ideal valuations and coprime Hensel factor lifting at2.
They classify only degree-one/two LOCAL components, not global fields.
For a prime over2, normalize v_P(P)=1 and v_P(2)=e_P; its contribution
to v_2 of a norm is f_P*v_P. Squares have even valuations. These
statements apply in the maximal integer ring, with no assumed integral
power basis. A monic cubic's reduction into coprime factors lifts to the
corresponding product of local algebras even if Z[U] has even index.

## 2. A stronger actual two-cycle: UV/4 is integral

Let g=(r0,rf). The accepted actual closed old-chain cycle makes g
integral. Spatial projection gives product UV/2; write the time product

    T=UV/2-g, T^2=U(U-2)V(9V-4)/2.

By (1), T is integral. The actual TWO-cycle is

    L=g^2/(UV),

and cyclic integrality makes L integral as well. Direct expansion yields

    L+T=UV/4+(U-2)(9V-4)/2
       =19UV/4-2U-9V+4.                               (3)

Thus 19UV/4 is integral. Since 3*19-14*4=1 and UV is integral,

    UV/4=3*(19UV/4)-14UV is integral.                  (4)

No unit-mod2 assumption on U enters. This lemma uses actual walls and
both actual cyclic integrality conclusions; a field-square norm alone
would not supply it.

Put

    A=N(U)>0, B=N(U-2)=|N(2-U)|>0, n=N(V)>0.

Taking norms in (4) and (1) gives

    64 divides A*n, 8 divides B*n.                    (5)

## 3. Only one pair of head norms survives

The elementary strips give A<=37 and B<=29. There is a stronger
product bound. Put s=U_phys-1 in (3,25/3), and let r1,r2 be the
absolute values of the other two U-1 conjugates, each in (0,1).
The nonzero integer norm of U-1 gives s*r1*r2>=1. Therefore

    A*B=(s^2-1)(1-r1^2)(1-r2^2)
       <=(s^2-1)(1-r1*r2)^2
       <=(s^2-1)(1-1/s)^2
       <(616/9)*(22/25)^2=298144/5625<54.               (6)

Here the first inequality is r1^2+r2^2>=2r1r2. Both factors in the
penultimate upper bound increase for s>1, so the strict endpoint bound
is valid. No conjugate was set equal to its interval endpoint.

If n is odd, (5) forces A>=64, impossible. If v_2(n)=1, (5) forces
A>=32 and B>=4, contradicting (6). Since n<=14, only n=4,8,12 remain.
At n4/12, A is a multiple of16 and B a multiple of2, leaving only
(A,B)=(16,2). At n8, A is a multiple of8. Moreover A and B have the
same parity: the integer minimal polynomial p_U satisfies
p_U(0)=-A and p_U(2)=-B. Thus B is even. Applying (6) leaves EXACTLY

    (A,B)=(8,2),(8,4),(8,6),(16,2),(24,2).              (7)

Write p_U(X)=X^3-TX^2+S X-A. Then

    5<=T<=13, S=2T+(A-B)/2-4.                         (8)

All pairs in (7) except (8,4) have S odd. We exclude them uniformly.
For even A and odd S, p_U modulo2 is x times a unit-constant quadratic;
these factors are coprime. Over Q_2 there is a degree-one component
with U even, and a complementary degree-two algebra on which U is a
unit. By (4), V/4 is integral on that degree-two component. Its local
norm is therefore divisible by 4^2=16; the norm from the other component
is integral. Consequently 16 divides n, contrary to n<=14. This argument
does NOT assume that p_U modulo2 describes the maximal order at a
repeated factor: invertibility of U on the entire lifted unit component
is all that is used.

It follows that

    A=8, B=4, n=8,
    p_U(X)=X^3-TX^2+(2T-2)X-8.                        (9)

## 4. Exactly two remaining head polynomials, one excluded locally

Set y=X-1 and t=T-3. The shifted minimal polynomial is

    q_t(y)=y^3-t*y^2+y+t-6, 2<=t<=10.                 (10)

Its physical root is in (3,25/3) and the other roots in (-1,1).
For t<=5 and |y|<1, q_t(y)<=q_5(y). For y<=0, q_5(y)<0 termwise;
for 0<=y<1, y^3<=y^2 gives

    q_5(y)<=-4y^2+y-1=-4(y-1/8)^2-15/16<0.

There would be no roots in (-1,1). At t6 the polynomial has root0,
contrary to cubic minimality. At t>=9,

    q_t(25/3)<=q_9(25/3)=-944/27<0,

whereas all three roots would be smaller than25/3, forcing a positive
value. Hence ONLY t7/8, or T10/11, can occur.

First take T11:

    p_U=X^3-11X^2+20X-8.

Modulo2 this is x^2(x+1), with coprime factors. It splits locally into
a degree-one unit-U component and a degree-two even-U component.
The latter has U-norm valuation3 and (U-2)-norm valuation2, since the
unit component contributes zero to both, by A8/B4.

The degree-two algebra cannot be unramified quadratic: its U-norm
valuation would be even. It cannot be Q_2 times Q_2 either: two even
U values with norm valuation3 have valuations1 and2, respectively;
their U-2 valuations are at least2 and exactly1, whose sum is>=3,
not2. Thus the component is ramified quadratic, with v(2)=2,

    v(U)=3, v(U-2)=2, v(H0^2)=5.                      (11)

On the degree-one unit-U component, (4) forces v(V)>=2. Since n=8
has total norm valuation3, the ramified component has v(V)<=1. There
v(4)=4, so v(9V-4)=v(V) and

    v(Hf^2)=2v(V)-2 is even.

Equation (11) makes v(Hf^2/H0^2) odd. It cannot be a square in K.
This excludes T11 without a global field classification.

We are left with T10:

    p_U=X^3-10X^2+18X-8.                               (12)

Put delta=2/(U-2). Substitution in (12) gives

    delta^3+5delta^2+4delta-2=0.                        (13)

Thus delta is integral and K=Q(delta). This minimal cubic has
discriminant316=4*79. The integer index squared divides316; its only
possible nontrivial value is4. Index2 would give field discriminant79,
which is3 modulo4, impossible for a number-field discriminant. Hence

    O_K=Z[delta], Disc(K)=316.                         (14)

No compatible lattice or other cubic field has been discarded by choice.

## 5. The fork time norm leaves exactly seven polynomial rows

Let v be the physical V and a,b its other conjugates. Write its cubic
minimal polynomial as

    p_V=X^3-R X^2+S_V X-8.

Since ab=8/v<(4/9)^2, v>81/2. Thus the integer trace satisfies
41<=R<=75. Define the positive integer

    m=N(9V-4)=5768+144R-324S_V.                        (15)

The actual time-norm ratio, with the RAW fork factor1/2 retained, is

    N(Hf^2)/N(H0^2)=(n*m/8)/(A*B)=m/32.               (16)

It must be a rational square. Since m is a positive integer, its square
class is2, so m=2k^2 for an integer k. Formula (15) makes m divisible
by4; hence k is even and

    m=8l^2, l a positive integer.                     (17)

There is a strong complete upper bound. AM--GM gives

    (4-9a)(4-9b)<=[4-9sqrt(ab)]^2.

Here sqrt(ab)=sqrt(8/v)>3sqrt2/13, while 9v-4<672. Using the strict
rational bound sqrt2>24/17 gives

    0<m<672*(236/221)^2=37427712/48841<800.             (18)

The square test for sqrt2 is 576<578; the final comparison is
37427712<39072800. Thus 1<=l<=9. Reducing (15)/(17) modulo9 gives
l^2=1 modulo9, so ONLY l1 or8 remains. Solving (15) with these values
and 41<=R<=75 yields the following COMPLETE seven possibilities:

| l | R | S_V | Disc(p_V) |
|---:|---:|---:|---:|
|1|41|36|-2704|
|1|50|40|30272|
|1|59|44|198448|
|1|68|48|617792|
|8|49|38|-250812|
|8|58|42|-256784|
|8|67|46|-72956|

For l1 the equation is S_V=(160+4R)/9 and R=5 modulo9; for l8 it
is S_V=(146+4R)/9 and R=4 modulo9. These two arithmetic progressions
in the complete trace interval give exactly the table. Its discriminants
come from the single integer formula

    R^2*S_V^2-4*S_V^3-32R^3-1728+144R*S_V.

Every negative row contradicts total reality. The remaining actual V
generates K, so its discriminant must be316 times an INTEGER SQUARE.
But 30272 and617792 have remainders15 and12 modulo79, respectively.
The last value198448 divided by316 is628, strictly between25^2 and26^2.
It is not a square either. All seven rows are impossible.

## 6. Conclusion and precise remaining scope

Every actual cubic dotted-head case supplies an isolated dotted fork
with (1). The actual two-cycle, norm bounds and complete local argument
reduce it to the unique head field (14); the complete fork time-norm
table then supplies the contradiction. Therefore the cubic dotted long
head is excluded in WHOLE arithmetic D18.

This proof uses one fork at a time, and neither assumes distinct fork
caps nor supplies a geometric completion from formal Gram data. It
invokes no global admissible-field list, numerical search or saved
survivor census. Quadratic and quartic cases are not decided by this
source; any separate quartic consequence of Section2 needs its own
reviewed composition. Whole D18, an arithmetic17 ceiling, new example
and geometric Lean remain outside this theorem. All six papers and
the accepted9353 checkpoint remain unchanged while this proof is reviewed.
