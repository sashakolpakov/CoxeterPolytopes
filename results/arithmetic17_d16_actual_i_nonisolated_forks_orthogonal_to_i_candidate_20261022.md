# Arithmetic D16 plus actual I: nonisolated forks are orthogonal to old I

22 October2026. Root-authored complete manual candidate, **pending full
outside-author mathematical review**. No worker, mathematical import,
enumeration, numerical/symbolic process, test, Lean or TeX ran. Earlier
accepted sources, six manuscripts and global ceilings are unchanged.

**Theorem.** In the actual farthest ENTIRE D16+actual-I configuration of
a compact arithmetic Coxeter17-polytope, every fork opposite that is
nonisolated within D16 is orthogonal to the old I wall. This does NOT
remove the fork's possible contact with the I opposite, exclude an
isolated dotted fork, or prove whole-forest nonexistence.

## Exact premises

All paths in the table are under `results/`.

| Source | SHA-256 | Restricted interface |
|---|---|---|
| `arithmetic17_d16_actual_i_all_forks_dotted_and_nonisolated_degree_four_candidate_20261022.md` | `6cb2a4662222dd7984f8e29f8ed532eee22456d2726034e4f115e49b65252387` | Nonisolated roles j1/j15; dotted own; old-I contact beta0 or1/2; integral t=2gamma, ACTUAL K=Q(t), strict conjugate intervals and exact C formulas |
| `arithmetic17_d16_actual_i_joint_head_height_and_shorter_forks_candidate_20261022.md` | `3a5c21e33622c62bb77a81d8cdbe394e1179001cc7a19a984cbbf6ebd6030da4` | Both ACTUAL C<37/4; a nonisolated fork with beta1/2 has ACTUAL field degree<=3 |
| `arithmetic17_d16_actual_i_unconditional_centre_via_field_strips_candidate_20261022.md` | `bf01ed7f20334c2d5fb1d17d00d4dadb09b390e45e7bf938a395da7b5c194523` | Actual-I height0<h<23/100; J=4(1+h) integral, with0<J_sigma<4 at every nonidentity place |
| `arithmetic_uniform_classical_weak_helper_obstruction_20261021.md` | `535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905` | Integral invariant cycles and ACTUAL K Gram normalization; in particular squared physical unit products belong to K |

These sources have separate full outside and root acceptances. This new
proof is independent of the later quartic-nonisolated coefficient lemma:
the beta1/2 degree<=3 interface was already proved by the independent
head-height/fork-length argument. Standard external algebraic facts are
integer minimal-polynomial coefficients, norm/trace, AMGM, and the
explicit quadratic integer rings used below. No script or field census
is assumed, and no relative D/I root-length squareclass is imposed.

## 1. Assume the remaining positive old-I contact

The complete accepted list leaves beta0 or1/2. Suppose beta1/2. Set
t=2gamma, where gamma>1 is the dotted own contact. The exact actual
edge formulas are

    j1: C=(8t+1)/3; j15: C=(8t+7)/3.

Consequently t<107/32 in the far-long role j1 and t<83/32 in the
near-fork role j15. The ACTUAL field is K=Q(t), t is integral and its
degree is at most3. Degree1 is impossible: the integral actual-I
invariant J would be a rational integer strictly between4 and123/25.

For convenience put v=t at j1 and v=t+1 at j15. The physical upper
bounds are107/32 and115/32 respectively. Their nonidentity conjugates
lie in(-.641,.391) at j1 and(-.391,.641) at j15.

## 2. Both cubic possibilities are impossible

Suppose the actual degree is3. If either small conjugate has the wrong
sign, its absolute value is below.391, while the other is below.641.
Then |N(v)|<(115/32)*.391*.641<1, impossible. Thus both are negative
at j1 and positive at j15. Since |N(v)|<(115/32)*.641²<2, its
positive value is1. Write the two small magnitudes x,y in(0,.641).
Then

    v*x*y=1, x+y<1.282, v>2.43.

The last inequality follows from2.43*.641²<1. Also v<4, so AMGM
gives x+y>=2/sqrt(v)>1.

For j15, v=t+1>3. Its integral trace v+x+y is therefore greater
than4 and less than115/32+1.282<5. This is impossible.

For j1, its integral trace v-x-y is greater than2.43-1.282>1 and
less than107/32-1<3. Thus it equals2, giving v=2+x+y<3.282.
AMGM requires v-2>=2/sqrt(v). The difference increases for v>0;
at3.1 it is negative since3.1*1.1²<4. Hence v>3.1.

The second elementary coefficient is the integer

    E2=x*y-v*(x+y)=1/v-v*(v-2).

Its upper bound is1/3.1-3.1*1.1<-3. Its lower bound is
-3.282*1.282=-4.207524>-5. Therefore E2=-4. The norm is1 and
trace2, so the proposed cubic minimal polynomial would be

    P(X)=X³-2X²-4X-1=(X+1)*(X²-3X-1),

which is reducible. This contradicts actual degree3. Only quadratic
actual fields remain under the beta1/2 assumption.

## 3. The complete quadratic generator table

Use the same strict conjugate windows and physical bounds. At j15,
v=t+1 is between3 and115/32. Its integral trace is3 or4. At trace3
the norm is negative with magnitude below2, so it is-1, giving
v=(3+sqrt13)/2. At trace4, its positive small conjugate is greater
than13/32 and below.641; the integral norm lies strictly between1
and3, hence is2, giving v=2+sqrt2.

At j1, v=t lies in(2,107/32) and its trace is2 or3. At trace2 its
negative norm has magnitude below2, giving t=1+sqrt2. At trace3,
a positive norm lies in(0,2) and hence is1, giving t=(3+sqrt5)/2;
a negative norm has magnitude below2 and hence is-1, giving
t=(3+sqrt13)/2. Zero norm is impossible. These interval arguments
exhaust the quadratic possibilities, with no prescribed field list:

| Role | t | Actual K |
|---|---|---|
|j15|`(1+sqrt13)/2`|`Q(sqrt13)`|
|j15|`1+sqrt2`|`Q(sqrt2)`|
|j1|`1+sqrt2`|`Q(sqrt2)`|
|j1|`(3+sqrt5)/2=1+phi`|`Q(phi)`|
|j1|`(3+sqrt13)/2`|`Q(sqrt13)`|

For detail on the strict norm bounds: near trace3 implies v<3.391,
so v(v-3)<2; near trace4 gives v(4-v)>3*(13/32)>1 and less than
(115/32)*.641<3. Far trace2 implies v<2.641, so v(v-2)<2.
Far trace3 with positive conjugate has v>2.609, giving
v(3-v)<2.609*.391<2; with negative conjugate it has
v(v-3)<(107/32)*(11/32)<2.

## 4. The complete actual-I strips in those fields

The integral J satisfies4<J<123/25 and0<J'<4. In each quadratic
integer ring write J=a+bw, with w=phi, sqrt2 or(1+sqrt13)/2.
The embedding difference forces b>0. For sqrt2 or sqrt13 the upper
bound on that difference forces b1; the physical interval then fixes
a3 or a2 respectively. For phi it permits b1/2, but b2 and J'>0
force a>=2 and J>123/25, so again only b1 remains, with a3. Thus

    Q(sqrt2): J=3+sqrt2, h=(sqrt2-1)/4;
    Q(phi): J=3+phi, h=(phi-1)/4;
    Q(sqrt13): J=(5+sqrt13)/2, h=(sqrt13-3)/8.

Each displayed J has the required conjugate. No additional row in any
of these ACTUAL fields survives the strict integer intervals.

## 5. An actual two-cycle excludes all five rows

The fork wall c and the actual I opposite c_I are distinct by their
old-D contacts. Their unit descriptions give

    c.c_I=-F, F=A*sqrt(h)-(1/2)*sqrt(1+h),
    A²=2t²+t/2-1/4 for j1,
    A²=2t²+7t/2+5/4 for j15.

All these squared quantities belong to the ACTUAL K. The invariant
two-cycle gives4F² in O_K, in particular F² in K. But

    F²=h*A²+(1+h)/4-sqrt(h*(1+h)*A²).

It follows that h*(1+h)*A² must be a square in K. On each COMPLETE
row above this quantity is Z/64, where:

| Role and generator | Z | Rational norm N(Z) |
|---|---|---:|
|j15, t=(1+sqrt13)/2|`29t+37`|-81|
|j15, t=1+sqrt2|`77+56sqrt2`|-343|
|j1, t=1+sqrt2|`47+32sqrt2`|161|
|j1, t=1+phi|`77t-33`|-605|
|j1, t=(3+sqrt13)/2|`33t+12`|243|

These expansions use respectively t²=t+3, sqrt2²2, t²=3t-1
and t²=3t+1. The rational norm of a square in a quadratic field is
a nonnegative rational square. Each displayed norm is negative or
a nonsquare positive integer (161=7*23 and243=3^5). Since64 is
itself a rational square, none of Z/64 can be a K square. This
contradicts the ACTUAL two-cycle in every row.

Thus beta1/2 is impossible in every possible actual degree1/2/3,
exhausting its prior degree<=3 interface. All nonisolated fork
opposites must have beta0, which proves the stated orthogonality.

The proof never identifies disconnected physical unit lengths with
K squares. It uses the actual invariant squared product, including
its mixed square-root term. The two fork parameters and opposite
walls need not coincide. All remaining zero-old-I-contact roles,
isolated dotted forks and the entire D16+I forest remain open. Full
outside-author review and root acceptance are required before
promoting this new conclusion; no runtime or Lean credit is claimed.
