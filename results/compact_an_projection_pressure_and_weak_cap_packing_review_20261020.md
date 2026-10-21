# Actual-face pressure and two native cap packets at a whole A_n corner

20 October 2026. Paper only. No worker, computation process, optimization,
dimension claim, manuscript or Lean change. The projection review below
is independent of its author; the supporting cap calculations were
developed in discussion with root and guard.

## 1. Source-bound projection review: PASS

Full source reviewed:
`compact_farthest_spherical_face_projection_pressure_20261020.md`,
SHA256 `abfcb96b6eb74aefdc6498cda362d81b1d4a1f9c32ad17fe3b464cff26fba565`.
Its imported actual-projection input is
`compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`,
SHA256 `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`,
previously reviewed in full. All paths are under `results/`.

I rederived the following. For UNIT normals, `c=a*p-v`, a finite
old set J has `w=sum(B_J^-1 beta)_j*n_j` and
`delta=1-beta^T B_J^-1 beta>0`. Its actual projection has displacement
`a(v-w)/(a^2+delta)` and displacement square `a^2/(a^2+delta)`.
Comparison with the farthest vertex gives

    x.(v-w) >= a/[1+sqrt(delta/(a^2+delta))].          (1)

The positive inverse and actual halfspaces prove feasibility, not an
intrinsic-Coxeter assertion about a metric face. Strict old interiority
gives `x.w>0` for a nonzero contact vector. Together with `x.v<a`,
(1) also gives the useful strict pressure

    x.w < a*sqrt(delta)/(sqrt(a^2+delta)+sqrt(delta)) < sqrt(delta).

For an isolated contact, the printed square of the ridge floor is
correct. Its comparison with the corrected long floor reduces to
`tan(theta)<F(q)`, q=sqrt(D), with the printed F. Clearing the
positive denominator in `(F^2)'` gives

    2q(2q+1)(2q^2+8q+5)
      =2q(4q^3+18q^2+18q+5)>0.

At q=sqrt(6), comparison with8/15 is exactly `331>116sqrt(6)`;
`8/15>tan^2(pi/5)` follows from `sqrt(5)>67/30`.
No square-root floor was confused with its square.

The H4 inverse-end coefficients `(3/2,1,1/2)`,
`delta_5=(7-3sqrt(5))/16`, and both averaging identities check.
In particular `w_+.w_-=(3/4)beta^2` gives
`||v-(w_++w_-)/2||^2=a^2+1-(15/8)beta^2`.
The strict inequality in its equation(3) follows from `||x||<1`.
All other contacts remain in v and a. No false monotonicity under an
arithmetic embedding, or elimination of near4 profiles, is asserted.

## 2. Far5/simple3 cap contact indices are three-separated

This is a geometric packing lemma for the following LITERAL profiles.
Use raw old A_n roots of norm2, N=n+1>=15, duals omega_i and

    C_ij=min(i,j)-ij/N.

Suppose an actual outside wall c_i has just two nonzero old products:
`c_i.alpha_i=-phi` and `c_i.alpha_e=-1`, where e=1 or n and
`5<=i<=N-5`. Its norm can be chosen2 since the far contact is3.
Write `c_i=a_i*p-v_i`,

    v_i=phi*omega_i+omega_e, H_i=a_i^2=v_i^2-2>0.

For two such walls at i<j put `S=v_i.v_j>0` and `R=c_i.c_j`.
They are distinct, since their exact old contacts differ. Then

    R=S-sqrt(H_i H_j)
     =2-||v_i-v_j||^2/2+(sqrt(H_i)-sqrt(H_j))^2/2.   (2)

All actual-wall products are nonpositive. If `R>-1`, the native
angle/disjointness gap forces R=0. If `R>-phi`, the only candidates
are0,-1,-sqrt(2): a finite m-bond gives `R=-2cos(pi/m)`;
m>=5 already has R<=-phi, while a dotted/tangent pair has R<=-2.
The candidate -sqrt(2) cannot occur here. Indeed S,H_i,H_j belong
to Q(phi), and squaring `sqrt(H_i H_j)=S+sqrt(2)` would put
sqrt(2) in Q(phi), since S>0. No arithmeticity assumption is needed.

Here are all distance-one/two cases. Reversing the chain handles
both same-end choices. BOTH opposite-end orientations are treated
explicitly below; exchanging the roots alone would not justify that.

* Same far end, j=i+1:
  `||v_i-v_j||^2=phi^2*(N-1)/N<4`, so R>0, impossible.
* Same far end, j=i+2:
  `||v_i-v_j||^2=2phi^2*(N-2)/N`; hence R>2-phi^2>-1,
  so R=0. Put `P=(i+1)(N-i-1)`, `B=N-2i-1`.
  Expanding `S^2=H_i H_j` and comparing the rational and phi
  coefficients gives

      B=0,             P=2(N+1).

  Thus `i=(N-1)/2` and `N^2-8N-9=0`, giving N9 or-1.
  Neither is in the stated range.
* Opposite far ends: put d=j-i and epsilon=-1 for far ends
  `(e_i,e_j)=(1,n)`, epsilon=+1 for `(n,1)`. Then

      S=phi^2*C_ij+phi*(1+epsilon*d/N)+1/N,
      ||v_i-v_j||^2=phi^2*d(N-d)/N+2(N-2)/N
                    -epsilon*4phi*d/N.

  At d1 this is less than6 in either orientation, so R=0.
  At d2 it is less than `2phi+4` in either orientation, so R>-phi;
  after excluding -sqrt(2), only R=r=0 or-1 remains. Put
  `U=H_i-S,V=H_j-S`. The equation for R=r is

      N^2*[S*(U+V+2r)+U*V-r^2]=0.

  The rational coefficients of its left side are as follows;
  `P=i(N-i-2)+1` is used in the d2 rows.

  |d|epsilon|r|rational coefficient|
  |---|---|---|---|
  |1|-1|0|N(N+6)|
  |1|+1|0|-N^2|
  |2|-1|0|2N(P+N+3)|
  |2|+1|0|2N(P-N-1)|
  |2|-1|-1|N(N+6)|
  |2|+1|-1|-N(3N+2)|

  All are nonzero. For the only row needing the central range,
  `5<=i<=N-7` gives `P-N-1=i(N-i-2)-N>=4N-35>0`.

For transparency, in the same-end distance-two row the same expansion
is a nonzero common factor times `B+[P-2(N+1)]phi`.
The opposite-end d2 rows can be checked from

    epsilon=-1: NS=P+(P+N-3)phi,
       NU=(2i-N-2)+(N+2)phi,
       NV=(N-2i-6)+(N+2)phi;
    epsilon=+1: NS=P+(P+N+1)phi,
       NU=(2i-N-2)+(4i-N-2)phi,
       NV=(N-2i-6)+(3N-4i-10)phi.

For d1 use `P=i(N-i-1)+1`. At epsilon=-1,
`NS=P+(P+N-2)phi`, `NU=(i-N-2)+(N-i+1)phi`,
`NV=-i-3+(i+2)phi`. At epsilon=+1,
`NS=P+(P+N)phi`, `NU=(i-N-2)+(3i-N-1)phi`,
`NV=-i-3+(2N-3i-4)phi`. These give every table entry directly.

Consequently any two contact indices in this literal far5/simple3
family differ by at least3, even when their far endpoints differ.
For arbitrary nonnegative index weights, their total is bounded by
the original-index recurrence

    F_i=max(F_(i-1), weight_i+F_(i-3)),

with unavailable indices omitted rather than compressed. This counts
these cap contact indices only, not every direction a general wall
could close. Other native profiles require their own constraints.

## 3. A near4/simple3 cap forces a long actual mate

Use UNIT old A_n normals and `D_j=2j(N-j)/N` here. Let an actual
outside wall have exactly the contacts

    beta_i=1/sqrt(2), beta_(i+1)=1/2,
    5<=i<=N-5, N>=15.

This profile is not locally excluded by arithmetic positivity. It
does, however, also close the ACTUAL direction j=i+1. Indeed
`I\{j}+c` is the positive definite whole forest
`B_(i+1)+A_(n-i-1)`. Positive actual-wall projection proves that its
finite vertex is feasible, hence no earlier opposite can truncate
this mate. The Schur complement at this mate is

    delta=1/(i+1),
    sigma=1/2+i/[(i+1)sqrt(2)],
    cosh L_j=sigma*sqrt((i+1)D_j).                  (3)

For j in[6,N-4], the minimum of `j^2(N-j)` is at an endpoint.
The right-minus-left endpoint value is `4(N-7)(N-10)>=0`, so

    jD_j>=72(N-6)/N>=216/5.

Also `sigma>13/12` from `1/sqrt(2)>7/10`. Therefore
`cosh^2 L_j>507/10>49` and

    tanh^2(L_j/2)>3/4.                              (4)

The reflected near profile has the reflected mate. This is an actual
same-wall packet, not an assumed shadow bound or a discarded profile.
If used in a global budget, the mate must be counted only once:
each direction has one actual opposite. The exact two-contact
hypothesis is retained, so no extra spectator is silently suppressed.

These lemmas support a complete weak-profile analysis; neither they
nor the reviewed projection source alone exclude a whole A_n corner
or improve a uniform dimension ceiling.
