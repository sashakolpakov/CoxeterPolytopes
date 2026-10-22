# Arithmetic D16 plus actual I: nonisolated forks force a quadratic actual field

22 October2026. Complete manual mathematical candidate, **pending full
outside-author review and root acceptance**. No worker, mathematical
import, enumeration, numerical/symbolic process, test, Lean or TeX ran.
No frozen source, paper, guide, public bound or earlier receipt changed.

**Theorem.** At an actual farthest vertex of a compact arithmetic
Coxeter17-polytope with ENTIRE finite forest D16 plus actual isolated
short I, any fork opposite nonisolated within D16 forces the ACTUAL
defining field to be quadratic. In particular, in actual degree3 both
fork opposites must be isolated within D16 and dotted. Their old-I
contacts are not asserted zero in that isolated case.

This excludes the COMPLETE cubic nonisolated branch, separately for
either fork. It does not exclude quadratic nonisolated forks, any
isolated dotted fork, the whole forest, or arithmetic dimension17.
Possible contact with the I OPPOSITE is retained and supplies the
new integral cycle below. No relative D/I length squareclass is assumed.

The author contributed inherited actual-foot and field-normalization
foundations and their reviews. Root proposed the use of the actual
fork/I-opposite cycle and reviewed related quadratic leads. The present
half-integrality lemma, complete cubic coefficient proof and final
cubic-field strip are newly written here; independent review is required.

## 1. Exact accepted input interfaces

All paths below are under `results/`. The current direct hashes match.

| Input | SHA-256 | Exact use |
|---|---|---|
| `arithmetic17_d16_actual_i_all_forks_dotted_and_nonisolated_degree_four_candidate_20261022.md` | `6cb2a4662222dd7984f8e29f8ed532eee22456d2726034e4f115e49b65252387` | Complete fork list, dotted own, nonisolated helper3 at j1/j15, integral t and ACTUAL K=Q(t), exact conjugate and edge formulas |
| `arithmetic17_d16_actual_i_quartic_nonisolated_forks_exclusion_candidate_20261022.md` | `651bb43b1062034eb7b132143f54e6e292af2ff0b1e7daa9b1f0aeb82210a698` | Every nonisolated fork has actual degree at most3; no isolated-fork conclusion used |
| `arithmetic17_d16_actual_i_nonisolated_forks_orthogonal_to_i_candidate_20261022.md` | `ccfec9e21cb75aa51d6b5e9e21de0dfb67a25ec592d7878898ca5eafdd6be48f` | Such a nonisolated fork has beta0 against OLD I, without removing its I-opposite contact |
| `arithmetic17_d16_actual_i_joint_head_height_and_shorter_forks_candidate_20261022.md` | `3a5c21e33622c62bb77a81d8cdbe394e1179001cc7a19a984cbbf6ebd6030da4` | BOTH actual fork edges have C<37/4 |
| `arithmetic17_d16_actual_i_unconditional_centre_via_field_strips_candidate_20261022.md` | `bf01ed7f20334c2d5fb1d17d00d4dadb09b390e45e7bf938a395da7b5c194523` | Actual-I invariant J=4(1+h) integral, 4<J<123/25 and 0<J_sigma<4 at EVERY other embedding |
| `arithmetic_uniform_classical_weak_helper_obstruction_20261021.md` | `535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905` | Actual arithmetic root lines, invariant-cycle integrality and strict nonphysical definiteness |

These are earlier independently reviewed statements; this conclusion
is not used upstream. Their original geometric, admissible-form and
actual farthest-point hypotheses remain in force. Arithmetic reflection
lattices require no extra first/simplest-type restriction. Elementary
integer traces/norms, minimal-polynomial separability, the discriminant
index formula and polynomial interpolation are used explicitly below.
No list of admissible fields or finite computational output is imported.

Release fork16 and retain fork15, in the original unit D16 coordinates
of the inputs. Every nonisolated opposite is an ACTUAL wall

    c=A p-(t/2)*omega16-omega_j/2, j=1 or15, t>2,
    c_I=sqrt(h)*p-sqrt(1+h)*e, h>0.

The old-I contact of c is zero by ccfe. The two walls c,c_I are distinct:
the former has nonzero contacts with old D walls and the latter has none.
The other original fork is covered by changing the sign of e16; no
symmetry or equality of its opposite wall is assumed.

Put v=t at j1 and v=t+1 at j15. The exact formulas and C<37/4 imply

    j1: 2<v<4.65; all other v_sigma in(-.641,.391);
    j15: 3<v<4.9; all other v_sigma in(-.391,.641).       (1)

For example C=(8t+1)/sqrt17 or(8t+7)/sqrt17; using sqrt17<4.124
gives the stronger v<4.644/4.894 respectively. The strict intervals
are the original quadratic-Schur intervals, enlarged only outward.
In both roles v is integral and K=Q(v) is the ACTUAL field. Its degree
is at most3. Degree1 is impossible since the integral J would be a
rational integer in(4,123/25). It remains to exclude degree3.

## 2. A small actual-I parameter cannot be twice an integer

Write y=J-4. Its physical value is in(0,23/25); every other value
is in(-4,0). In an actual field of degree2 or3,

    y/2 is NOT an algebraic integer.                    (2)

Indeed put x=y/2 and suppose it integral. In degree2 its nonzero
norm has magnitude less than (23/50)*2=23/25<1. In degree3 write
the other values as -b,-c with 0<b,c<2. Then

    0<N(x)=x*b*c<(23/50)*4<2,
    N(x)=1, bc=1/x.

AMGM gives b+c>=2/sqrt(x)>2.94. Thus the integral trace x-b-c
is strictly between-4 and-2 and equals-3. Consequently b+c=x+3.
Since (b+c)^2>=4bc, one has x*(x+3)^2>=4. The left side increases
for x>0 and its value at1/3 is100/27<4, so x>1/3.

The second elementary coefficient is the integer

    E2=bc-x*(b+c)=1/x-x^2-3x.

It is positive because x<23/50 and the decreasing expression at
23/50 is still positive. It is less than2 because x>1/3 gives
1/x<3 and x^2+3x>10/9. Therefore E2=1. Its proposed polynomial is

    X^3+3X^2+X-1=(X+1)*(X^2+2X-1),

which is reducible. The physical value x is not rational: a rational
value would remain positive at the other embeddings. In a cubic field
it therefore has degree3, contradicting this reducibility. This proves(2).

## 3. The actual two-cycle gives a polynomial parity obstruction

The unit time squares of c are

    A^2=q(v)/2,
    q(v)=4v^2+v-1 at j1; q(v)=4v^2-v-1 at j15.         (3)

The ACTUAL product c.c_I=-A*sqrt(h) is nonzero. Its invariant
two-cycle, regardless of whether that contact is finite or dotted, gives

    4(c.c_I)^2=A^2*(J-4)=q(v)*y/2 in O_K.              (4)

This uses the invariant product, not membership of a disconnected unit
length or of sqrt(h) in K. In particular (4) does not set this actual
contact to zero or discard an inconvenient native value.

If N(q(v)) were odd, multiplication by the integral cofactor
N(q(v))/q(v), followed by integer Bezout with2, would turn(4) into
y/2 in O_K, contrary to(2). Thus N(q(v)) is even. Since

    q(v) = v-1 modulo 2O_K

in BOTH roles, the determinant definition of the integral norm gives
N(q(v)) congruent to N(v-1) modulo2. If P is the monic CUBIC
minimal polynomial of v, then N(v-1)=-P(1). Consequently

    P(1) is EVEN.                                      (5)

This necessary parity is valid in the full integer ring; no assumed
power integral basis is needed for it.

## 4. Complete cubic windows in the far-long role j1

Suppose K has degree3. Zero conjugates are impossible. If both small
conjugates were positive, |N(v)|<4.65*.391^2<1. If both are negative,
write them -a,-b, with 0<a,b<.641. Their positive integral norm is
less than4.65*.641^2<2 and hence is1. Thus vab=1.

Here v>2.43 since2.43*.641^2<1, and a+b>=2/sqrt(v)>.92.
The integral trace T=v-a-b is therefore2 or3. Write

    S=ab-v*(a+b)=1/v-v*(v-T).

If T2, then v<3.282. The inequality v*(v-2)^2>=4 implies
v>3.1, since3.1*1.1^2<4. These bounds give -5<S<-3, so S=-4.
The polynomial is X^3-2X^2-4X-1=(X+1)*(X^2-3X-1), reducible.

If T3, AMGM gives v*(v-3)^2>=4, hence v>=4, while v<4.282.
Then -6<S<-3.75 and S is-5 or-4. For S=-5 the polynomial
X^3-3X^2-5X-1=(X+1)*(X^2-4X-1) is reducible. For S=-4,
P(1)=-7 is odd, contradicting(5). This exhausts equal negative signs.

For mixed signs write the small roots a,-b, with a<.391,b<.641.
Now the negative norm has magnitude less than4.65*.391*.641<2,
so vab=1. It follows v>3.98 because3.98*.391*.641<1. Also
b>1/(4.65*.391)>.55>a. Thus the integral trace T=v+a-b is
strictly between3 and5 and equals4. In particular v>4. The integer

    S=-1/v+v*(4-v)

lies between-4 and0: its lower bound exceeds-.25-4.65*.65>-4.
Thus S=-3,-2 or-1. Since P(1)=S-2 must be even, only S=-2 remains.
But P(X)=X^3-4X^2-2X+1 has no root in(-.641,0). Indeed

    P(-b)=1+2b-4b^2-b^3

is strictly concave on[0,.641], and both endpoint values are positive
(the second is .375101279). Concavity keeps it positive throughout.
This contradicts the required negative root. Every cubic j1 sign
pattern is therefore impossible.

## 5. Complete cubic windows in the near-fork role j15

Both small conjugates negative would give |N(v)|<4.9*.391^2<1.
If both are positive a,b<.641, the positive norm is1 or2.

Norm2 would imply v>4.86 since4.86*.641^2<2. Moreover
a+b>=2sqrt(2/v)>1.27 because4.9*1.27^2<8. Its integral trace
would then lie strictly between6.13 and4.9+1.282<7, impossible.

For norm1, vab=1. The trace T=v+a+b is greater than4: AMGM gives
T>=v+2/sqrt(v), which increases for v>1 and at3 is greater than4.
Also T<7, so T5 or6. If T6, a+b=6-v>1.1 and
v>6-1.282>4. As(.641-a)(.641-b)>0,

    ab>.641*(a+b)-.641^2>.294219,

contradicting vab=1. If T5, then v>3.718. AMGM gives
v*(5-v)^2>=4, so v<=4. The integer S=1/v+v*(5-v) satisfies
4<S<6 (use1/v<.27 and v*(5-v)<3.718*1.282<4.767).
Thus S=5 and P=X^3-5X^2+5X-1=(X-1)*(X^2-4X+1), reducible.

Only mixed signs remain. Write a,-b with a<.641,b<.391. Again
vab=1 because4.9*.641*.391<2. The inequalities

    v>3.98, a>1/(4.9*.391)>.52

place the integral trace v+a-b strictly between4 and6, so T=5.
In particular v>4.359. The integer S=v*(5-v)-1/v is positive
(greater than.49-.23) and less than4.359*.641<3. Therefore S=1
or2. Now P(1)=S-3 is even by(5), forcing S=1. The ONLY remaining
cubic possibility is

    P(X)=X^3-5X^2+X+1.                                (6)

No polynomial enumeration or prescribed field list underlies this
exhaustion; each coefficient follows from the stated root windows,
integer norms/traces, AMGM and the actual cycle parity.

## 6. The full integer ring of(6)

Let v be the physical root of(6). The discriminant is404=4*101.
Hence the index [O_K:Z[v]] is1 or2 by the discriminant index formula.
If the index were2, an integral element outside Z[v] could be written

    alpha=(A+Bv+Cv^2)/2, A,B,C in{0,1}, not all zero.

Since Tr(v)=5 and Tr(v^2)=23, integral trace requires A+B+C even.
There are exactly three nonzero possibilities. But

    N((1+v)/2)=6/8=3/4,
    N((1+v^2)/2)=36/8=9/2,
    N((v+v^2)/2)=-6/8=-3/4,

none an integer. Here N(v)=-1 and N(1+v)=-P(-1)=6.
Also N(1+v^2)=1+23+11+1=36, because the sum of pairwise
squares is S^2-2*T*N=1-2*5*(-1)=11. Thus

    O_K=Z[v].                                         (7)

This index argument checks the full integer ring; replacing it by an
unjustified monogenic-order assumption would not suffice below.

## 7. The actual-I strip is empty in this last field

Direct positive rational sign evaluations of(6) place its three roots in

    4.74<v<4.75, .60<a<.61, -.35<b<-.34.               (8)

Indeed its values at4.74/4.75 are-.101576/.109375; at.60/.61
they are.016/-.023519; at-.35/-.34 they are-.005375/.042696.
These disjoint intervals contain three distinct real roots and hence
all roots. The physical root is the one greater than3.

Suppose the required actual integral J exists. By(7) write
J=A+Bv+Cv^2 with integers A,B,C. Its other embeddings J_a,J_b
satisfy0<J_a,J_b<4, whereas4<J_v<4.92. Lagrange interpolation
expresses the quadratic coefficient as

    C=J_v/[(v-a)(v-b)]
      -J_a/[(v-a)(a-b)]+J_b/[(v-b)(a-b)].             (9)

The intervals(8) give v-a>4.13, v-b>5.08, a-b>.94.
The positive terms in(9) are less than1/4 and1 respectively; the
negative term is greater than-2. Hence -2<C<2 and

    C in{-1,0,1}.

The small-embedding difference gives

    |B+C*(a+b)|=|J_a-J_b|/(a-b)<4/.94<4.26.

Since .25<a+b<.26, this implies |B|<4.52 and thus -4<=B<=4.
Each of the three coefficient possibilities now contradicts the strip:

* C=-1 gives J_v-J_b=(v-b)*(B-v-b)<0, since v+b>4.39>B.
  But J_v>4>J_b requires this difference positive.
* C=0 requires B>0 by that same difference, hence B>=1. Then
  J_v-J_b=B*(v-b)>5.08, whereas J_v-J_b<4.92.
* C=1 gives J_v-J_a=(v-a)*(B+v+a)>
  4.13*(-4+4.74+.60)=5.5342, again greater than4.92.

Thus NO integral J has the actual-I embedding strip in the final
cubic field. This excludes(6) and completes the cubic branch.

## Conclusion and exact remaining scope

Every nonisolated fork has actual degree at most3 by651, old-I
contact zero by ccfe, no rational actual field by the J strip, and no
cubic actual field by Sections2-7. Its actual field is therefore
quadratic. This statement includes BOTH original outgoing forks and
ALL dotted own values permitted by the prior actual length bound.

The new proof uses no independent energy at a replacement vertex,
no addition of distinct optimizers and no further bootstrap of centre
or head floors. All geometric work is the exact previously accepted
original-vector interface. The extra actual c/c_I cycle is used as
an integral scalar invariant, never as a claim that sqrt(h), a unit
normal or a disconnected relative length lies in K.

Quadratic nonisolated forks remain OPEN in this source. It also does
not remove foreign old-I contacts from isolated forks, identify the
two fork parameters/walls, or exclude an isolated dotted pair. No
whole-forest exclusion, arithmetic16 ceiling, field-wide nonexistence,
example, runtime or Lean result is claimed. Root and full other-author
review must accept this new theorem before it changes any public status.
