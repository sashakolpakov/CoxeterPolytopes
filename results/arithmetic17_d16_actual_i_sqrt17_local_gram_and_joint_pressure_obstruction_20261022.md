# D16 plus I: a local sqrt17 model and the limit of separate gains

22 October 2026. **Manual supporting candidate, pending root and other-author review.** No worker, enumeration, numerical/symbolic process, scientific import, test, Lean or TeX ran. This is a local arithmetic Gram construction and an explicit relaxation of listed pressure constraints. It is NOT an actual compact Coxeter polytope, a full-arrangement chamber, a farthest vertex, or a counterexample to a possible D16+I exclusion.

The purpose is precise: the accepted head foot, ordinary I barriers and individual loss-.61 floors5..14 cannot finish the norm contradiction by themselves. Their first-five gains overlap. The displayed local model even respects native labels and integral reflection coefficients. A simultaneous finite projection onto the two cap walls supplies an additional constraint which the relaxation violates; thus further genuine geometric pressure remains available.

## 1. Input boundary and contributions

The restricted D16+actual-I reduction is
`results/arithmetic17_d16_actual_i_strong_five_finite_head_and_height_reduction_candidate_20261022.md`, SHA `71b49bf5272c31374702a91c271e6f203f2abc7a773b7f609b53f3a39a5c0412`.
It supplies strong5..13, an isolated-within-D dotted head,0<h<.23 and the actual pure-D head foot. The new centre14 candidate is
`results/arithmetic17_d16_actual_i_unconditional_centre_via_field_strips_candidate_20261022.md`, SHA `bf01ed7f20334c2d5fb1d17d00d4dadb09b390e45e7bf938a395da7b5c194523`.
For the relaxation below we impose its centre floor as an additional constraint, independently of its review status. Its full outside report is843a272a. The actual finite-projection mechanism is
`results/compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`, SHA `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`.

I contributed earlier I/point/paid-reference foundations and reviewed both D16 sources. Root proposed testing the final simultaneous head/I budget. The explicit sqrt17 configuration and relaxation here are my derivation. No earlier source's conclusion is changed.

## 2. A full-rank admissible local Gram over Q(sqrt17)

Put s=sqrt17, K=Q(s), O=Z[(1+s)/2], and

    U=(7+s)/2, J=(5+s)/2,
    (U-2)(J-4)=2, UJ=13+3s.                       (1)

Both U,J are integral, with norms8 and2. Under the other embedding,

    0<U'=(7-s)/2<2, 0<J'=(5-s)/2<4.

Take the rational D16 Cartan space with roots alpha_i of square2,
native bonds-1 and all other products0. Let omega be its first dual,
so omega.alpha_i=delta_(1i), omega²=1. Adjoin orthogonal vectors f,z,
also orthogonal to D16, with

    f²=2UJ, z²=U(2-U).

The resulting K-form has physical signature(17,1): D16 is positive,
f²>0 and z²<0. At the other embedding all three summands are
positive definite, since2U'J'>0 and U'(2-U')>0. In particular this
is an admissible, K-anisotropic local form; its anisotropy follows
already from the positive definite conjugate. No finite-volume
reflection chamber is inferred.

Define two further root vectors

    a=z-U*omega,
    c=z/(U-2)-f/2.                                (2)

They lie in this K-space, and direct expansion using(1) gives

    a²=c²=2U, a.c=-U,
    a.alpha1=-U, a.alpha_i=0 (i!=1), a.f=0,
    c.alpha_i=0, c.f=-UJ.                         (3)

Thus the nineteen vectors alpha1,...,alpha16,f,a,c span rank18;
the first seventeen form exactly the orthogonal finite forest D16+A1.
All distinct displayed products are nonpositive physically. The
mutual a,c label is EXACTLY3. Their other positive contacts are
dotted, because the invariant squared unit contacts are U/4>1 and
J/4>1. All remaining bonds are the old D16 native3 bonds or zero.

Let M be their finite O-span. It is a full lattice in this K-space:
the independent set D16,f,a spans the space, and the additional
vector c has fixed K denominators. Reflection in any displayed root
preserves M. Indeed the nonzero coefficients `2(v.r)/r²` between
generators are precisely-1 on old D bonds; -U or-1 on alpha1/a;
-1 in either direction on a/c; and-J or-1 on c/f. All belong to O.
This is why f was given square2UJ: using the shorter vector f/2 as
an O-generator would leave a coefficient-J/2. No root-length scale
between disconnected components has been assumed.

The cap Schur matrix after the old D16+f block is

    [[ U(2-U),       -U ],
     [    -U, U(4-J)/2 ]].                        (4)

Its determinant is0 by(1). Physically it has one negative and one
zero eigenvalue; at the other embedding it has one positive and one
zero eigenvalue. The zero simply records the one linear relation
among nineteen vectors in dimension18. The independent ambient
form remains nondegenerate, with the signatures proved above.

Choose z future and p=z/sqrt(U(U-2)). In outward UNIT normalization,
with n_i=alpha_i/sqrt2, e=f/sqrt(2UJ), the caps are

    a_unit=A p-gamma*omega1_unit,
    c_unit=sqrt(h)p-sqrt(1+h)e,
    gamma²=U/4, A²=(U-2)/2,
    h=(J-4)/4=(s-3)/8 in(0,.23).                 (5)

Here omega1_unit=sqrt2*omega. Both time coefficients are positive.
Moreover A²h=1/4, so the mutual product is indeed-1/2. The actual
cycle parameters would be U,J,1, all integral; their conjugate
strips are compatible with definiteness. This is a concrete local
arithmetic model, not a proof that these roots are all the actual
simple walls of any full integral chamber.

## 3. An explicit vector satisfying the separate available constraints

Work in physical D16 orthonormal coordinates, and put

    f_i=sqrt(i)-.61/sqrt2, b=f5/5,
    w=(b,b,b,b,b,sqrt6-sqrt5,...,sqrt14-sqrt13,0,0),
    u=gamma^-2=4/U=(7-s)/4,
    Z=(1-u/2)e1+(u/2)e2,
    Y=Z-P5/5, mu=9/25=.36,
    x=w+mu*Y, y=sqrt(q(h)),
    q(h)=(sqrt(1+h)-1)/(sqrt(1+h)+1).              (6)

These are explicit radical quantities. The first five coordinates
of Y sum to0. Hence S_i=x.P_i=f_i for ALL5<=i<=14. The head
foot norm and pressure are

    Q=Y²=4/5-u+u²/2=(89-15s)/80,
    H=[sqrt(2-u)-sqrt(u(1-u))]/sqrt2.

The rational square bounds4.1231<s<4.1232 give

    .7192<u<.719225, .3394<Q<.34,
    .3609<b<.361.

The bound on b follows directly from2.236<sqrt5 and
.61/sqrt2<.4314, while b<.361 is the earlier elementary seal.
Consequently

    x.Z=b+mu*Q>.3609+.36*.3394=.483084.

On the other hand `sqrt(2-u)<1.13173` and
`sqrt(u(1-u))>.44935`: square the positive bounds, using that
u(1-u) decreases for u>.5. Thus H<.4827. The actual normalized
head-foot inequality `x.Z>=H` is therefore satisfied strictly.

The ordinary actual-I endpoint barrier also holds strictly:

    x1=b+.36*(.8-u/2)>.5193,
    t(h)=sqrt(1+4h)-2sqrt(h)<.502.

For the latter use h>.14, monotonicity of t, sqrt1.56<1.25 and
sqrt.56>.748. The other normalized D dual pressures exceed t as
well. At indices2,3,4 the prefix sums are2b+.216,3b+.144,4b+.072;
after normalization they exceed.66,.70,.75. At5..14 they are
f_i/sqrt(i)>.80. Both fork pressures are S14/4>.82 before the
small chamber perturbation below. The external I energy is exactly
its ordinary endpoint floor y²=q(h).

There is nevertheless no norm contradiction. Since h<.141 and
q is increasing, q(h)<.034. This last inequality follows from
`.141 < 4*.034/(1-.034)²`, the inverse relation for q. Also

    E14=w²<.918.

For a fully rational upper proof, f5²/5<.652 and every increment
square `(sqrt(i+1)-sqrt(i))²<1/(4i+1)`. The latter holds because
the squared denominator `(sqrt(i+1)+sqrt(i))²` exceeds4i+1.
For i5..13 these nine rational upper bounds sum to less than
.265326, giving E14<.917326<.918. Since w.Y=0,

    ||x||²+y²=E14+mu²*Q+q(h)
              <.918+.1296*.34+.034=.996064<1.    (7)

Initially x3=x4=x5 and x15=x16=0. Arbitrarily small positive
decreasing perturbations of the first fourteen coordinates, with
x15 positive and x16=0, give the strict D chamber and all strict
floors5..14 while preserving(7), both strict pressure margins and
y>0. Increasing y slightly also makes its endpoint inequality
strict. The norms and inequalities are continuous and every used
strict margin is positive; no limiting equality earns actual
farthest-point status.

Thus the listed individual floors, the head foot, ALL normalized
old-D I direction barriers, the I energy, native a/c label and the
arithmetic Gram restrictions above are insufficient by themselves.
This forbids adding the two independent first-five gains without
their cross term. It does not show that the full actual-point
inequality system is feasible.

## 4. A simultaneous actual projection explicitly rejects this relaxation

In fact this same model exposes a useful missing pressure. If a_unit
and c_unit were actual walls, their positive definite Gram is A2,
with diagonal1 and off-diagonal-1/2. The positive acute projection
from p onto these TWO walls has coefficients

    lambda_a=(4/3)(A+sqrt(h)/2),
    lambda_c=(4/3)(sqrt(h)+A/2),
    T=1+lambda_a*A+lambda_c*sqrt(h)
       =(13+3s)/6.

Its spatial displacement is

    v=lambda_a*sqrt2*gamma*e1
                  +lambda_c*sqrt(1+h)*e.

The projected point is future of square-1 after normalization,
and the farthest comparison would require

    lambda_a*sqrt2*gamma*x1
      +lambda_c*sqrt(1+h)*y >= T-sqrt(T).         (8)

The coefficients are positive, so actual-wall acuteness would make
this point feasible against every unselected wall. This is a single
simultaneous pressure, not a sum of two separately optimized gains.

The vector(6) violates(8), with ample rational margins. From the
same square brackets, A<1.3345, sqrt(h)<.3748,
sqrt2*gamma<1.668, sqrt(1+h)<1.068, x1<.52, y<.185.
Hence lambda_a<2.0292 and lambda_c<1.3894, and the left side of(8)
is less than

    2.0292*1.668*.52+1.3894*1.068*.185<2.035.

Meanwhile T>4.2282 and sqrt(T)<2.057, so its right side exceeds
2.1712. The strict perturbations in Section3 may be chosen small
enough that this violation persists. Thus this vector is explicitly
NOT an actual farthest vector when both caps are actual walls.

## 5. Precise conclusion

The proven obstruction concerns only the final proof strategy based
on the separate old head/I pressures and floors5..14. The local
sqrt17 arithmetic model shows that those data do not already imply
energy>=1, even after imposing their native mutual contact and
integral reflection coefficients. The genuine two-cap finite
projection supplies additional geometry which can overcome that
relaxation; whether it closes the entire remaining family requires
a complete proof over every other mutual parameter and field.

No actual chamber, finite-volume completion, new example, entire
D16+I exclusion or improved global bound is claimed. The six papers,
accepted sources and prepared norm379 worker/spec remain unchanged.
