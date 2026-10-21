# Every A component of rank at least eight has actual projection above 4/7

20 October 2026. Complete paper candidate by resource_cleanup, held for
root and outside-author mathematical review. No scientific computation,
search, optimizer, Lean, TeX, manuscript edit or status promotion.

**Proposed theorem.** At a vertex p farthest from a strict interior point
of a compact hyperbolic Coxeter polytope, every whole A_r component,
r>=8, has ||x_A||^2>4/7, with arbitrary finite spectators. All pressure
inequalities use the SAME actual farthest vector z=E(p-x), ||x||<1.
This is an actual-projection theorem, not an edge-only energy bound.

## 1. Exact accepted inputs and what changes

All paths below are under results/. The r8 case, including its complete
L/H/F feet and five mixed-foot comparisons, is the accepted source
`compact_universal_a8_three_foot_actual_projection_20261020.md`, SHA256
`0f1b3bbe8bfcaba6f53d717664286a7b19dbaaba0911bba8626a50fe91a8800f`.
The complete ALL-RANK inventories and general-N foot formulas are in
`compact_large_an_actual_half_uniform_transfer_20261020.md`, SHA256
`7351b471ba3b0edaf0a524ccf082010e458b7fedd8d063491a75632223bd70e7`.
Their geometric inputs include the corrected A7 source, SHA256
`5ae61f571761f517222bef38f35468462ec8fb9ad9d879850c419a6c316828a5`.

The actual-point, positive acute projection, reciprocal-shadow and
orthogonal-component opposite interfaces are retained exactly as stated
there. In particular, for an actual opposite c=A p-v with foreign cost
Q, each contacted foreign component has one contact with cost q>=1/4,
there is a gap(1/4,1/3), and q<A^2. A positive combination giving a
future negative vector is feasible against EVERY unknown wall, not just
the old walls. These are proved geometric inputs, not new hypotheses.

For r>=9 put N=r+1>=10 and use unit A normals and inverse Gram

    n_i=(e_i-e_(i+1))/sqrt2,
    C_ij=2min(i,j)(N-max(i,j))/N,
    u_i=omega_i/sqrt(C_ii), h_i=x_A.u_i.

Assume ||x_A||^2<=4/7. Actual edges give h_i>=tanh(L_i/2)>0. For
DISTINCT actual opposites, the exact pair law is

    F(h_i)F(h_j)>=u_i.u_j, F(t)=2t/(1+t^2).

We strengthen the uniform-half argument in two places: the paired-end
forcing Gram is monotone, and the L foot's far correlation improves
already at N10. The finite inventory itself needs no enlargement.

## 2. The A8 four-pair certificate is monotone in N

Take pairs1/2,3/4,(r-3)/(r-2),(r-1)/r. Their respective correlations
increase with N. The exact N9 calibration and concavity of log F in the
accepted A8 source therefore give, whenever the pair is distinct,

    h1+h2>41/40, h3+h4>1219/1000,

and the reflected bounds at the other end. Put

    M_i=u_i+u_(N-i) (1<=i<=4),
    W=2M1+2M2+M3+M4.

For 1<=i<=j<=4 direct inverse-Gram substitution gives the exact identity

    M_i.M_j=2N*sqrt(i/[j(N-i)(N-j)]).                    (1)

It includes the diagonal case2N/(N-i). Its logarithmic derivative is

    1/N-1/[2(N-i)]-1/[2(N-j)]<0.

Hence EVERY entry of this paired-end Gram decreases with N>=9, and
all coefficients defining W are positive. The exact N9 seal from A8
therefore applies uniformly:

    W^2<=W_9^2<747/10.

If all four actual pairs are distinct then x_A.W>3269/500. But

    (3269/500)^2-(4/7)*(747/10)=104527/1750000>0,

contradicting Cauchy. Thus, after merely reversing the numbering, a
common actual opposite occurs at1/2 or3/4. No polytope symmetry is
assumed by this relabeling.

## 3. Complete uniform common-pair inventory at the stronger target

The all-rank finite classification in the pinned uniform-half source
is independent of the assumed projection threshold. Here are all its
branches and the two places where the stronger edge cutoff is needed.
Every retained white diagram is an ACTUAL positive replacement.

For common3/4, there is at most one contact on each of the old paths
{1,2,3} and {4,...,r}. For double contacts not exactly at3,4, one
deletion retains both new bonds internally in a connected rank-r finite
diagram. A large finite path has only a terminal4, and D/E have only3.
Thus both contacts must be3, giving an affine cycle, or the white set
is already nonpositive. This includes the rank8 E8 exceptions.

At3,4 the possible labels are3/4/5 at3 and3/4 at4. The3/3 triangle
is affine. The two minimal nontrivial endpoint bounds are

    cosh L3>=sqrt(3/(2N))*[(N-3)+sqrt2*(N-4)],
    cosh L4>=(3sqrt2+4)*sqrt((N-4)/(2N)).

Both increase with N and exceed11/3 already at N9 by the exact A8
seals. Larger bonds and foreign cost only increase the appropriate
endpoint bound. Therefore its half-tangent square exceeds4/7.

Single3 contacts at positions3,...,r-2 contain the principal affine
E8 arms1,2,5 or E7 arms1,3,3, with longer arms truncated. The four
remaining positions1,2,r-1,r have the Q bounds

    Q<5/8, Q<1/2, Q<2/(N-3), Q<1/2+1/[2(N-3)],

respectively. At the ends the height is<1/4; at2,r-1 there is at most
one foreign component and its q=Q>A^2. All are excluded. Single4 is
affine internally; at1 orr its white bound is Q<1/4 orQ<1/(N-3), so
Q=0 and its height is negative. Single5 can occur only at1 or3, and
its H4 white cost>3/4 clears all spectators. At3 it also closes4,
with cosh L4=3phi^3*sqrt((N-4)/N)>11/3, again excluded. At1 it is F.
Zero own contact has zero own-edge derivative and cannot close the pair.

For common1/2, there is at most one contact on2,...,r, and possibly
contact1. The full finite-white list for double positions1,j>=3
has both bonds3, except4@1+3@3 and4@1+3@r. The all3 cases are
affine cycles;4@1+3@r is the exact affine F4 path2-1-c-r-(r-1).
For4@1+3@3 the foreign bound Q<1/[2(r-1)] clears spectators; this
is type H. At adjacent1,2 the bond at2 is3/4, that at1 any finite
native label;3/3 is affine and the other rows are type L.

Single3 internally has the same affine principals. At1,r its height
is<1/4. At r-1 the bound Q<2/r leaves no spectator and negative
height. At2, Q<1/2+1/(2r) permits at most two foreign components:
one contradicts q>A^2, while two have first bonds3 and give the exact
affine two-fork D5 principal c,old2,old1,old3 and the contacted foreign
nodes. Single4 is affine internally; at1 at most one foreign component
has q>A^2, and at r no foreign component is possible and height is
negative. A single label>=5 is possible only at1. Label6 is affine3,6;
for native>=7, Q<1-beta1^2<1/4 forces Q0. Isolated5,Q0 is F, and
isolated5 with Q>0 is L. Thus the surviving types are exactly:

* L: adjacent1,2, isolated finite>=7 at1, or isolated5 at1 with Q>0;
* H:4@1+3@3, Q0;
* F: isolated5 at1, Q0.

All exclusions in this inventory are finite-white, actual affine,
spectator-cost or explicit endpoint statements. The old half-bound
exclusion of H/F is NOT imported at the stronger threshold: we retain
both feet below. Arbitrary foreign components have not been omitted.

## 4. Actual feet, monotone pressures and the improved L correlation

For each type the pinned sources give a positive combination
xi=A p-U, of norm -k<0 and positive time coefficient. Its products
with every actual wall are nonpositive. Thus

    x_A.(U/||U||)>=g, g=(A-sqrt(k))/sqrt(A^2-k),
    g^2=(A-sqrt(k))/(A+sqrt(k)).                          (2)

The raw supports of L,H,F occupy the first3,4,5 coordinates,
respectively, with raw sum sqrt2*B>0. Normalization subtracts their
mean sqrt2*B/N. Their exact general-N formulas imply that A^2 increases
with N, k and B stay fixed, and hence BOTH g^2 and U^2/B^2 increase
with N. For L this is read at fixed native/foreign parameters, whose
domain below is rank-independent. The accepted N9 lower bounds are

| Type | g^2 | U^2/B^2 |
| --- | --- | --- |
| L | >9/20 | >4/9 |
| H | >53/100 | >1/3 |
| F | >11/20 | >2/9 |

For clarity the L identities, after positively clearing ALL foreign
old roots, are B=b+2d, z0=h/B^2, h=4(b^2+bd+d^2)/3+Q-1, and

    8/81<z0<1/3, a_N=2/3-2/N,
    A^2=B^2(z0+a_N), k=B^2(z0-3z0^2),
    U^2/B^2=a_N+3z0^2,
    U_raw=(B/sqrt2)*(2/3+z0,2/3+z0,2/3-2z0,0,...,0).

The explicit positive old-A2 coefficients, future sign, complete wall
products and z0 domain are those proved in the uniform-half source;
none depends on the projection threshold. H has
k=sqrt2-5/4 and A^2=3/2+sqrt2-(11/2+3sqrt2)/N. F has
k=34phi-55 and A^2=phi^2*(N-1)/(2N)-1. Thus their monotonicity is
also immediate, rather than an assumption about different polytopes.

Put W_f=u_(r-1)+u_r. For all three types the raw-sum identity gives

    rho_T=(U.W_f)/||U||
      =[2/sqrt(N(N-2))+sqrt(2/[N(N-1)])]/sqrt(U^2/B^2).    (3)

Its numerator decreases and its denominator increases with N.
Consequently the accepted N9 bounds rho_H<3/4 and rho_F<9/10
remain valid. For L we use N>=10. At N10 the numerator square is
5/36, while the denominator square exceeds

    7/15+64/2187>49/100.

The last inequality has positive cross-product margin3891. Thus

    rho_L^2<125/441<64/225, rho_L<8/15,                 (4)

where the final cross-product margin is99. Finally
W_f^2=2+2sqrt((N-2)/[2(N-1)])<2+sqrt2<7/2.

## 5. Distinct far opposites are impossible

If the two far actual opposites are distinct, the same outer-pair law
gives x_A.W_f>41/40. Write v=x_A.(U/||U||), J=W_f^2, E0=4/7.
Perpendicular Cauchy gives

    x_A.W_f<=rho*v+sqrt((E0-v^2)(J-rho^2)).               (5)

For L, v>=sqrt(9/20). The right side decreases with v on this domain:
the derivative condition is rho^2 E0<v^2 J, and it holds already from
rho^2<64/225, J>=2 and v^2>=9/20. Substitute the lower endpoint for v,
then use(4) and J<7/2 to obtain

    x_A.W_f<(8/15)*sqrt(9/20)+sqrt((17/140)*(7/2))
            <9/25+33/50=51/50<41/40.                    (6)

The rational square seals are sqrt(9/20)<27/40 and
sqrt(17/40)<33/50. This is the only new far-foot estimate needed.

For H and F no derivative estimate is needed: v<=sqrt(E0) and the
strict pressure lower bounds give, using sqrt(E0)<19/25,

    H: x_A.W_f<(3/4)*(19/25)+sqrt(29/200)<.57+.4=.97;
    F: x_A.W_f<(9/10)*(19/25)+sqrt(3/40)<.7+.28=.98.

Both are below41/40. Thus EVERY far-distinct case is excluded.

## 6. Common far opposites and conclusion

If the two far actual opposites are common, Section3 supplies a right
foot L,H or F. Reversing normals uses e_j mapping to -e_(N+1-j),
so its raw sum is -sqrt2*B'. At N>=10 the left/right raw supports
are disjoint for EVERY pair, including F/F. Therefore

    U.U'=2BB'/N,
    corr(U,U')=(2/N)/sqrt((U^2/B^2)(U'^2/B'^2)).           (7)

For each pair its correlation upper bound improves with N, while its
two normalized pressure lower bounds increase. Except for F/F we may
therefore use directly the complete accepted A8 sum-direction table:

| Types | Correlation upper bound | Sum pressure lower bound | Projection lower bound |
| --- | --- | --- | --- |
| L/L | 1/2 | 2sqrt(9/20) | >3/5 |
| L/H | 3/5 | 139/100 | >19321/32000 |
| L/F | 5/7 | 141/100 | >139167/240000 |
| H/H | 2/3 | 36/25 | >3888/6250 |
| H/F | 5/6 | 73/50 | >15987/27500 |

Every last entry exceeds4/7. These are unchanged exact rational
Cauchy seals, not a new search over feet. The disjoint-support formula
already holds at N9 for these five pairs, so its comparison base is
valid even at H/F, whose raw supports together have exactly nine entries.

Finally F/F means two distinct actual isolated end5 walls, with no
foreign contact. Their exact product is

    c.c'=1-phi^2/2+phi^2/N.

For every N>=9 this lies strictly between -1/2 and0: phi^2<3 gives
the lower bound, and its value at N9 is1-7phi^2/18<0, since
phi>8/5. This violates the actual native gap. No accidental sign from
normal reversal or an unsupported small-overlap formula is used.

All possibilities contradict ||x_A||^2<=4/7. The accepted A8 base
supplies r8 and the present transfer supplies every r>=9, proving the
proposed theorem with arbitrary spectators.

## 7. Limited dimension23 consequence and contributions

If independently reviewed, this theorem eliminates the A16+rank7
necessary configuration in the separately held dimension23 frontier.
Even after one short-A1 discard its spectator cost is>6/14, so
4/7+6/14=1 with strict actual inequalities. With one isolated rank-two7
defect the total is>4/7+9/100+5/14>1. Ordinary spectators are stronger.
This consequence does not itself exclude the A19/A20 configurations;
those require separate pressure inputs. No new dimension ceiling is
promoted by this source before its complete mathematical reviews.

Finite57_research proposed checking monotonic paired-end Gram and foot
parameters; root assigned the A16 target. The exact decreasing paired
Gram identity(1), the N10 L-correlation seal and the uniform4/7 assembly
were derived here from their accepted A8/uniform-half inputs. Earlier
actual geometry and foot classifications retain their shared authorship.
All accepted and frozen sources remain unchanged.
