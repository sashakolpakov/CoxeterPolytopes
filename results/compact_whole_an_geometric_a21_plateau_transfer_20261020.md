# Whole A21--A29: transfer of the accepted farthest-vertex certificate

20 October 2026. Paper reasoning only. No worker, optimizer, Lean, TeX,
or modification of accepted evidence. Awaiting independent review.

**Theorem.** In a compact hyperbolic Coxeter n-polytope, for
21 <= n <= 29, a vertex farthest from a strict interior point cannot
have WHOLE finite diagram A_n. No arithmeticity is assumed.

Together with the already imported universal compact dimension ceiling29,
this excludes every possible whole A_n farthest vertex with n >= 21.
The argument below itself covers ranks21--29; it does not assert an
unbounded-rank local theorem or a component theorem with spectators.
Other possible vertex diagrams prevent any new dimension ceiling here.

## 1. Accepted inputs and notation

All source paths are under results/. The unchanged accepted whole-A21
proof is `compact_whole_a21_farthest_vertex_exclusion_candidate_20261020.md`,
SHA256 `3abe11bb70ce677069279462a91111ed6b8ab8889e3ec9e656ac71c1094200fb`.
Its complete outside review is
`review_compact_whole_a21_farthest_vertex_batch_guard_20261020.md`,
SHA256 `618e255766d77a00ea4457bf3f5b6d983ebb845bcd8352a259d291503dbf7c3f`.
Use its actual-point projection formulas, weak5 positive-dual table,
and exact reference-floor energy and hole-loss seals. The general-rank
far5/simple3 packing and near4 actual-mate input is
`compact_an_projection_pressure_and_weak_cap_packing_review_20261020.md`,
SHA256 `8b0a59493f1fa5dc7081e5e70beff2839ac5918bfcd50e089e05ae9b5546e43b`.
The A21 adjacent-packet statement is NOT imported at larger rank;
Section3 proves its required extension.

Put N=n+1, so 22 <= N <= 30. Unit old A_(N-1) normals have duals

    D_i=omega_i^2=2i(N-i)/N,
    C_ij=2i(N-j)/N (i<j),
    h_i=x.omega_i, h_0=h_N=0,
    ||x||^2=(1/2)sum(h_(i+1)-h_i)^2<1.

The h_i are strictly concave. Call a central index 5 <= i <= N-5
weak when h_i <= sqrt(D_i)-c, with c=61/100. Central dual norms obey

    85/11 <= D_i <= 15,  5/2 < sqrt(D_i) < 31/8.

Write F(q)=[1+(1-c/q)^2]/[1-(1-c/q)^2]. Direct clearing of positive
denominators shows

    F(q)/q < 761/500,  5/2 <= q <= 31/8.          (1)

The numerator after multiplication by 5,000,000 is the concave quadratic
`-715800q^2+3268319q-1860500`; its endpoint values are positive.

## 2. Complete weak inventory and the uniform contact strip

The accepted central finite-white classification and zero-own-contact
proof extend verbatim: both retained paths have rank at least4; a
two-path completion has connected rank N-1 >= 21, so is A_(N-1).
Thus the actual opposite is isolated, has one white3, one end white4,
or has the two end-white3 contacts. The same affine/nonpositive
subdiagrams exclude all the same released labels. Schur residuals are
unchanged functions of the corresponding path length m.

All the accepted isolated strong-floor estimates need only D >= 85/11,
and remain valid. All other profiles except near4/simple3 and the
literal far-m/simple3 packet have sigma/sqrt(delta)>761/500:
the smallest relevant lower bounds are
`(101/125+2/5)/sqrt(3/5)` for near-white3/released>=5 and
`(101/125+1/5)/sqrt(2/5)` for the D positions. Both exceed761/500;
the E, white4 and two-white3 bounds are larger. Equation(1) therefore
makes these directions strong.

For near4/simple3, the ratio is increasing with m. At m8 it equals
`2sqrt(2)/3+7/12 > 761/500`. At m7 and N22 or23 it exceeds3/2,
while D_7 <= 224/23 < (25/8)^2 and D_7 > 9; the accepted comparison
F(q)/q<3/2 on [3,25/8] applies. Consequently a weak near4 has its
actual strong mate at index4 or5 (or the reversal) for N22/23, and
at index4,5 or6 (or the reversal) for N24--30. Each mate has normalized
height greater than sqrt(3)/2 > 43/50, by the unchanged actual-mate lemma.

Every other weak central opposite has exactly

    c_i=A_i p-v_i, v_i=gamma_i omega_i+omega_e/2,
    e=1 or N-1, c_i^2=1.                       (2)

Its allowed released products are gamma=cos(pi/5), cos(pi/m) for
m>=7, or gamma>1 dotted. Put f=101/125 and g=27/25. Then

    f < gamma < g,   A_i>2.                    (3)

Only the dotted upper bound needs proof. For a far-side path parameter
m>=5, delta=(m+1)/(2m), sigma=gamma+1/(2m). At gamma=g,

    sigma^2/delta >= 4g-2=58/25 > (761/500)^2,

since the numerator of the first difference is
`[2m(g-1)-1]^2`. Equations(1) and the actual edge comparison exclude
gamma>=g. The time bound follows from A_i^2>f^2(85/11)-1>4.

## 3. No adjacent weak literal packets, uniformly in N

Let j=i+1 be central and have profiles(2). Write A=gamma_i,
B=gamma_j, and R=c_i.c_j. Weighted AM--GM gives

    A_i A_j <= [(B/A) A_i^2+(A/B) A_j^2]/2.

Substitution of the exact old Gram entries yields these lower bounds:

    same helper:
    R >= -AB(1-1/N) +/- (A-B)/N +(N-1)/(2N)
                                 +(N+1)(A/B+B/A)/(4N);
    outward opposite helpers (e_i=1,e_j=N-1):
    R >= 1/2-AB+(A-1)(B-1)/N;
    reversed opposite helpers:
    R >= 1/2-AB+(A+1)(B+1)/N.                  (4)

The same-helper bound is greater than1-g^2, hence greater than-1/2.
The two opposite-helper bounds exceed
`1/2-g^2-(1-f)(g-1)/22 > -7/10 > -1/sqrt(2)`.
The actual native gap consequently gives R=0 for a common helper,
and R in {0,-1/2} for opposite helpers. Unlike the former rectangle
estimates, (4) has no growing D_i(A-B)^2 error.

For outward opposite helpers, R=0 is excluded just as in the accepted
packet: eta_l=delta_l/sigma_l^2<4/5 and

    (D_i-4/5)(D_j-4/5)-C_ij^2
       = (2/5)C_ij-24/25+8/(5N) > 0,

using C_ij>=80/11. For the reversed assignment put l=i,r=N-j,
l+r=N-1. Its corresponding determinant F obeys

    NF=4lr A^2B^2-(l+1)(r+1)(A^2+B^2)-2AB
                         -(l+1)A-(r+1)B+(N+2)/4.             (5)

This is strictly increasing in A,B on [4/5,g]. If at least one is
>=9/10, concavity in l reduces its minimum to l5 or r5; substitution
at (4/5,9/10), in either order, gives
`NF >= (509N-9974)/500 >= 306/125 > 0`.
For the derivative, `4lr B^2-(l+1)(r+1) >= (39/25)lr-N`, and
l,r>=5,N>=22 give strict positivity immediately. Both remaining
gammas equal cos(pi/5); the imported general-rank far5 packing excludes
that distance-one pair. Thus opposite helpers have R=-1/2.

The accepted mixed replacement is therefore actual and finite, with
entrywise positive inverse Schur matrix, for EVERY remaining pair.
Retain its notation Z=W_i omega_i+W_j omega_j, W=W_i+W_j>0, and
Dcal=A^T H^(-1)A. The accepted bounds

    sqrt(1+Dcal)/W < 7/16 (same helper),
    sqrt(1+Dcal)/W < 3/10 (opposite helpers)

use only l,r>=5, delta_l in (1/2,3/5], gamma>4/5 and A_i>2,
so remain unchanged. Also, with E0=W_i sqrt(D_i)+W_j sqrt(D_j),

    E0-||Z|| < W/[2(sqrt(85/11)+sqrt(80/11))] < 3W/32.

Here sqrt(D_iD_j)-C_ij <= 1-1/N < 1. The actual mixed-point pressure
now makes the weighted average of the two losses less than
17/32 or63/160, respectively. Both are less than61/100. Two weak
literal directions therefore cannot be adjacent.

Weak near4 packets can now be removed exactly as in the accepted proof.
For N22/23 their mates have index<=5, and an opposite mate or strong
anchor among N-6,N-5 has correlation<1/3. For N24--30 the mate index
is<=6, and the correlation is at most6/(N-6)<=1/3. If there is no
opposite near4 packet, the adjacent-packet result guarantees a strong
anchor among N-6,N-5. The normalized heights exceed43/50 and39/50,
whose sum has square (41/25)^2>8/3. Cauchy contradicts ||x||<1.
Hence ALL weak central indices have form(2), and none are adjacent.

## 4. Exact plateau transfer removes weak5 and weak(N-5)

No new table of numerical certificates is required. View an old
reference vector W in the accepted A21 proof through its potential
t_i=W.omega_i, i0,...,22. Insert N-22 copies of t_11 between the
two old halves. This defines an isometric linear lift, since its energy
`(1/2)sum(t_(i+1)-t_i)^2` acquires only zero differences.

Every accepted weak5 certificate has the form

    W=x0+mu*(u-proj(u)).

The residual potential is supported entirely in the unchanged first6
and last6 indices (last5 for the alternate right anchor). The central
potential is that of x0: increasing to its peak at11 and decreasing
afterwards. Its positive coefficient at11 splits into the two
nonnegative slope drops at11 and N-11. All other certificate coefficients,
including the positive old-anchor coefficients after the correction,
are unchanged. Thus the lifted W has EXACTLY the accepted norm and
nonnegative coefficients on the available old constraints and pressure.

If5 is weak,6 is strong and at least one of N-6,N-5 is strong.
Their current floors are at least the accepted values a,d. Concavity
and nonadjacency imply h_11,h_(N-11) >= b, where the accepted middle
floor is b=sqrt(120/11)-.61: if an anchor is weak, both its neighbors
are strong, and their current dual norms dominate the old norms at10,12.

For each fixed gamma and helper assignment the actual time coefficient
is increasing with N:

    helper N-1: A_N^2=10gamma^2-1/2-(10gamma-1)^2/(2N);
    helper 1:   A_N^2=10gamma^2+2gamma-1/2-(10gamma+1)^2/(2N).

The local selected-face deltas in every accepted pressure row are
unchanged; H(A,delta) and A-sqrt(-delta) increase with A. All lifted
constraint lower bounds therefore dominate the accepted ones. The
strict accepted certificates, including their positive-dual signs,
contradict weak5 with identical norms and no reduced margin.
Reversal excludes weak(N-5).

## 5. The same reference energy closes every remaining rank

N22 is the accepted theorem. Suppose N>=23. At each of the two
reference peak locations11,N-11 the ACTUAL height is at least
`sqrt(11)-.61`, even if weak: its neighboring strong floors have
dual square at least `D_10(N23)=260/23>11`. Strong peak locations
satisfy the bound directly.

Lift the accepted full reference-floor curve to a flat plateau between
11 and N-11. Its energy is the SAME `Efull>1007/1000`, and its old
coefficients are nonnegative. The remaining reference constraints are
at5,...,10 and their reversals; all current strong floors dominate
these fixed old floors. Endpoints5,N-5 have already been proved strong.

Only the five sites6,...,10 and their reversals can now be omitted.
Nonadjacency allows at most three in each flank. Each omission has the
same unchanged local loss <1/1600 as in the accepted proof, and their
segment interiors are disjoint. The surviving reference potential
remains concave with nonnegative coefficients. Consequently

    ||x||^2 >= Efull-6/1600 > 4013/4000 > 1,

the required contradiction. This proves the stated rank21--29 extension
without altering the accepted A21 proof or assuming spectators harmless.
