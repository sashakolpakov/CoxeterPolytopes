# Arithmetic D16 plus actual I: a quartic field forces isolated forks

22 October2026. Root-authored complete manual candidate, **pending full
outside-author mathematical review**. No mathematical worker/import,
enumeration, numerical/symbolic process, test, Lean or TeX ran. All
accepted sources, six manuscripts and global ceilings are unchanged.

**Theorem.** In the actual farthest ENTIRE D16+actual-I configuration
of a compact arithmetic Coxeter17-polytope, any nonisolated fork
opposite forces the ACTUAL field degree at most3. Consequently, if
the actual field has degree4, BOTH fork opposites are isolated within
D16 and have dotted own contact. Their foreign old-I contacts are NOT
asserted zero. This does not exclude the quartic isolated-fork case,
the entire forest, degrees2/3 or arithmetic dimension17.

## Exact inputs

All paths in the table are under `results/`.

| Source | SHA-256 | Restricted interface |
|---|---|---|
| `arithmetic17_d16_actual_i_all_forks_dotted_and_nonisolated_degree_four_candidate_20261022.md` | `6cb2a4662222dd7984f8e29f8ed532eee22456d2726034e4f115e49b65252387` | Complete retained-D fork roles; dotted owns; nonisolated t=2gamma integral with ACTUAL K=Q(t), degree<=4 and all conjugate strips; exact actual C formulas |
| `arithmetic17_d16_actual_i_joint_head_height_and_shorter_forks_candidate_20261022.md` | `3a5c21e33622c62bb77a81d8cdbe394e1179001cc7a19a984cbbf6ebd6030da4` | BOTH actual fork lengths C<37/4; not its restricted foreign3 degree corollary |
| `review_arithmetic17_d16_actual_i_dotted_forks_height_and_field_four_root_acceptance_20261022.md` | `50e65569883340da5a5269679d1a53e645a576d34a8bceca433f437ccd830ddd` | Accepted acyclic scope of these two interfaces, including full source-bound other-author reports |

The new proof below uses elementary integral minimal-polynomial
coefficients, norm/trace and AMGM. No field list, finite enumeration,
relative D/I root-length squareclass or new geometric point is assumed.
All arithmetic and actual-wall foundations remain explicit in the two
accepted sources. Arithmetic reflection groups need no extra type
restriction.

## 1. Exact quartic windows and the signs of ALL conjugates

Suppose a fork is nonisolated. Its sole retained-D contact is a simple3
at j1 or at the other fork15, and its old-I contact is beta0 or1/2.
Its own is dotted. Put t=2gamma. The accepted raw normalization gives
ACTUAL K=Q(t), t integral and degree at most4. It suffices to rule out
degree4.

Set v=t for j1 and v=t+1 for j15. The exact conjugate intervals are

    j1:  (-1-sqrt17)/8<v_sigma<(-1+sqrt17)/8;
    j15: (1-sqrt17)/8<v_sigma<(1+sqrt17)/8.

Using sqrt17<4.124, these lie respectively in(-.641,.391) and
(-.391,.641). The exact actual edge formula and beta>=0 give

    j1:  0<v<(C sqrt17-1)/8<4.644;
    j15: 0<v<(C sqrt17+1)/8<4.894<4.9.             (1)

Here C<37/4. These are ACTUAL field generators and all three
nonidentity embeddings are included. If any small conjugate had the
wrong sign (positive for j1, negative for j15), then

    0<|N(v)|<4.9*.391*.641²<1,

impossible for a nonzero integral generator. Zero conjugates are also
impossible. Thus all three small conjugates are NEGATIVE for j1 and
POSITIVE for j15. Write their absolute values x1,x2,x3, each in(0,.641).
Moreover |N(v)|<4.9*.641³<2, so its positive absolute integer value
is1. Consequently

    v*x1*x2*x3=1,
    v>3.79, s:=x1+x2+x3<1.923.                    (2)

The lower bound follows from3.79*.641³<1. Since v<4.9<1.7³,
AMGM also gives

    s>=3/v^(1/3)>30/17.                           (3)

No sign is assigned to a conjugate merely by the physical sign of v;
the integral norm establishes the complete sign pattern first.

## 2. The near-fork role has an impossible third coefficient

For j15 all three small conjugates are positive. The integral trace
is v+s. Equations(1)-(3) place it strictly between5 and7, so

    Tr(v)=6, s=6-v, v>6-1.923=4.077.              (4)

AMGM gives v+3/v^(1/3)<=6. This function increases for v>1. At
v=21/5 it is strictly greater than6, because (21/5)^(1/3)<5/3.
Therefore v<21/5=4.2.

Let p=x1*x2+x1*x3+x2*x3 and q=x1*x2*x3=1/v. The third elementary
symmetric coefficient of v,x1,x2,x3 is the integer E3=v*p+q.
AMGM gives p>=3q^(2/3). Since v>4.077>1.59³,

    E3>=3v^(1/3)+1/v>4.77+5/21>5.

On the other hand p<3*.641²=1.232643, v<4.2 and1/v<1/4, so

    E3<4.2*1.232643+1/4=5.4271006<6.

This integer lies strictly between5 and6, impossible.

## 3. The far-long role has an impossible second coefficient

For j1 the small conjugates are-x1,-x2,-x3. The integral trace
is v-s. Equations(1)-(3) place it strictly between1 and3, so

    Tr(v)=2, s=v-2, v<2+1.923=3.923.              (5)

AMGM requires v-2>=3/v^(1/3). The difference of these two sides
increases for positive v. At v=3.9 it is negative, since

    3.9*1.9³=26.7501<27.

Hence v>3.9. With p as above, the second elementary symmetric
coefficient of v,-x1,-x2,-x3 is the integer E2=p-v*s. It satisfies

    E2<3*.641²-3.9*1.9=-6.177357<-6.

Also v<4 gives q=1/v>1/4. AMGM on x1*x2,x1*x3,x2*x3 gives
p³>=27q²>27/16>1, so p>1. Thus

    E2>1-3.923*1.923=-6.543929>-7.

This integer lies strictly between-7 and-6, another contradiction.

## Conclusion and scope

Both nonisolated roles therefore exclude degree4, with EVERY permitted
old-I contact included. Their prior degree<=4 result gives degree<=3.
The complete fork list and accepted dotted-own theorem imply that a
quartic actual field forces both forks to be isolated within D16 and
dotted. The proof neither equates their parameters nor asserts their
foreign contacts vanish. It does not classify fields or produce a new
global bound, example, runtime or Lean verification. Full outside-author
review and root acceptance remain required for this new refinement.
