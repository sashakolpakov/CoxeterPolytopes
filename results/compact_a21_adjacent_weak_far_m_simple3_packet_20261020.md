# Adjacent weak far-m/simple3 opposites at a whole A21 maximum

20 October 2026. Own paper proof, no computation or execution. This
completes the finite-angle coverage of the literal two-contact packet
below. It does NOT close the whole A21 inventory or imply a dimension
bound. Isolated caps, other multiwhite contacts and spectators remain
outside the statement.

## 1. Exact statement and inputs

At a farthest vertex p of a compact Coxeter polytope let the WHOLE
finite diagram be A21, with UNIT outward normals n_1,...,n_21. Set
N=22, let omega_i be their duals, and write

    D_i=omega_i^2=2i(22-i)/22,
    C_ij=omega_i.omega_j=2i(22-j)/22, i<j.

The spatial farthest vector x has ||x||<1. Put h_i=x.omega_i and
call i weak if h_i<=sqrt(D_i)-61/100. Consider adjacent central
i,j=i+1 with5<=i<j<=17. Assume their ACTUAL opposite walls are

    c_l=A_l*p-v_l, c_l^2=1,
    v_l=gamma_l*omega_l+(1/2)*omega_(e_l),
    e_l in{1,21}, A_l>0.                         (1)

These are the COMPLETE old contacts: a nonzero native finite m-angle
or a dotted contact at l, one simple3 at a far end, and no others.

**Theorem.** The two indices cannot both be weak.

The actual positive acute projection/farthest comparison input is
`compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`,
SHA256 `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`.
The already proved far5 distance-one packing input is Section2 of
`compact_an_projection_pressure_and_weak_cap_packing_review_20261020.md`,
SHA256 `8b0a59493f1fa5dc7081e5e70beff2839ac5918bfcd50e089e05ae9b5546e43b`.
It is geometric, with no arithmetic hypothesis, and handles BOTH
far-end assignments. Its norm2 convention rescales (1) by sqrt2.

## 2. The complete contact range under weakness

Finite m3 would make the all3 cycle consisting of c_l and the old path
from e_l to l, an affine A diagram. Finite m4 makes the affine star
{c_l,n_l,n_(l-1),n_(l+1)} with incident squared cosines1/2,1/4,1/4.
Finite m6 makes the affine3-node path6,3. Parallel contact is affine
rank2. All are impossible among actual walls of a compact Coxeter
polytope. Thus either gamma=cos(pi/5)=phi/2, gamma=cos(pi/m) with
m>=7, or gamma>1 is dotted. In particular

    f=101/125 < gamma.

For m>=7, gamma>9/10: u=2cos(pi/7) is the root above sqrt3 of
u^3-u^2-2u+1, a strictly increasing polynomial there, and its value
at9/5 is -1/125. The rational lower f uses phi>202/125.

Weakness also gives gamma<g=27/25. For finite gamma this is automatic.
For dotted gamma the exact individual edge formula is, with
m=l for e_l=1 and m=22-l otherwise,

    delta_m=(m+1)/(2m), sigma=gamma+1/(2m),
    A_l^2=D_l*sigma^2-delta_m,
    cosh L_l=sigma*sqrt(D_l)/sqrt(delta_m).        (2)

If gamma>=g, then sigma/sqrt(delta_m)>38/25: clearing denominators
gives28m^2-188m+625>0. For q=sqrt(D_l) in[5/2,10/3], the comparison
with weakness is the positive concave quadratic

    -36400q^2+163602q-93025>0,

with endpoint values88480 and430835/9. This is exactly
[1+(1-(61/100)/q)^2]/[q*(1-(1-(61/100)/q)^2)]<38/25.
Using h_l>=sqrt(D_l)*tanh(L_l/2) proves the claimed upper bound.

Hence all subsequent bounds use f<gamma_i,gamma_j<g. Also

    85/11<=D_i,D_j<=11, |D_i-D_j|<=1,
    C_ij>=80/11, (omega_i-omega_j)^2=21/11,
    A_i,A_j>2.                                   (3)

For the last assertion v_l^2>gamma_l^2 D_l and
f^2*(85/11)-1>4; its rational seal is10201*17>15625*11.

## 3. All native products between the two caps

Let R=c_i.c_j. The exact identity

    R=1-||v_i-v_j||^2/2+(A_i-A_j)^2/2             (4)

applies to these distinct ACTUAL facets.

If their far ends agree, the helper cancels. The quadratic in both
gammas attains its maximum on the four rectangle corners. Equal
corners have square at most2g^2. The mixed corners have square at most

    2fg+11(g-f)^2+f(g-f)<3.

Thus R>-1/2, so the complete native gap forces R=0.

If the far ends differ, put A=gamma_i,B=gamma_j, Delta=D_j-D_i.
For the order e_i=1,e_j=21 the square in (4) is

    Psi=1+2AB-2(A-1)(B-1)/22
                   +(A-B)^2 D_i+(A-B)(1-B)Delta. (5)

Equal rectangle corners give Psi<=1+2g^2. Both mixed corners are at
most

    1+2fg+(1-f)(g-1)/11+11(g-f)^2+(g-f)(1-f)
      =620937/171875 <452/125=2+2f.               (6)

For the REVERSED assignment e_i=21,e_j=1, the exact square is

    Psi_rev=Psi-(4/22)[2A+(A-B)(22-2j)].          (7)

The equal corners are no larger than in (5). At A=g,B=f, the remaining
quadratic is increasing for5<=i<=16 and attains its maximum at16;
at A=f,B=g it is decreasing and attains its maximum at5. Indeed its
derivative has sign(A-B)*[(A-B)(11-i)+1+B], whose second factor is
positive throughout this interval (its endpoint lower bound is56/125).
Both maximum values are608585/171875, also below2+2f. Convexity in
each gamma completes ALL rectangle branches.

Thus opposite ends always have R>-f>-cos(pi/5). The only native
products are0,-1/2,-1/sqrt2; tangent and dotted products are excluded.
The last product is impossible: letting c_L have its simple3 atold1
and c_R atold21, the ACTUAL five-wall path

    n_2 --3-- n_1 --3-- c_L --4-- c_R --3-- n_21

is affine F4 (labels3,3,4,3, not the finite3,4,3). All other products
in this subset are zero by (1). Compactness excludes this affine
diagram.

In fact R=0 is also impossible for opposite ends. For the order in
(5), project off the retained old paths. The cap spatial parts become
sigma_i*omega_i,sigma_j*omega_j, with
eta_l=delta_l/sigma_l^2<1/(2f^2)<4/5. Therefore

    (D_i-eta_i)(D_j-eta_j)-C_ij^2
      >(D_i-4/5)(D_j-4/5)-C_ij^2
      >=556/275>0.                               (8)

Here D_i D_j-C_ij^2=2C_ij and D_i+D_j=2C_ij+21/11.
The spatial cap product is positive, so (8) gives R<0.

For the REVERSED order put l=i,r=22-j, so l+r=21 and l,r>=5.
The retained projections are orthogonal and their residual cap norms
are delta_r=1/2+1/(2r),delta_l=1/2+1/(2l). The remaining spatial
parts are

    z_i=gamma_i*omega_i+omega_j/(2r),
    z_j=gamma_j*omega_j+omega_i/(2l).

For F=A_i^2*A_j^2-(z_i.z_j)^2, direct expansion gives

    22F=4lr A^2 B^2-(l+1)(r+1)(A^2+B^2)
                     -2AB-(l+1)A-(r+1)B+6.      (9)

This polynomial is strictly increasing in either A or B on[4/5,g]:
4lr B^2-(l+1)(r+1)>100 and the corresponding derivative is
2A times this bracket minus[2B+(l+1)], hence positive. If at least one
gamma>=9/10, evaluate (9) at(4/5,9/10) or the reverse. As a function
of l it is concave, so its minimum is at5 or16 and equals306/125>0.
Consequently R<0 here too. The only leftover case has BOTH gammas
equal phi/2, already excluded by the far5 distance-one packing input.

We have proved the exhaustive cap products

    same far end: R=0;
    opposite far ends: R=-1/2.                    (10)

## 4. Actual mixed vertex and full displacement pressure

Remove old i,j. All remaining old facets are positive definite. If
the far ends agree, their common retained helper projection has
square t=(ell-1)/(2ell), with ell=i for a left helper or ell=22-j
for a right helper. Its two-cap Schur matrix is

    H=I-t*ones(2,2), det(H)=1/ell>0.              (11)

If the ends differ, H has diagonal entries delta_l,delta_r in either
order and off-diagonal-1/2. Here l=i,r=22-j and each delta lies in
(1/2,3/5]. Thus det(H)>0 and H^-1 is ENTRYWISE POSITIVE. The same
holds in (11). All chosen walls are actual, so positive acute
projection gives an ACTUAL finite vertex, not a formal Gram model.

Put lambda=H^-1(A_i,A_j)^T>0, Dcal=A.lambda. The actual point and
farthest inequality are

    q=[(1+Dcal)p-Z]/sqrt(1+Dcal),
    Z^2=Dcal(1+Dcal),
    x.Z>=1+Dcal-sqrt(1+Dcal).                     (12)

The helper terms cancel from Z, leaving Z=W_i*omega_i+W_j*omega_j
with positive weights. For a common left end,

    W_i=lambda_i*gamma_i+(lambda_i+lambda_j)/(2i),
    W_j=lambda_j*gamma_j;

the common-right formula is the reversal. For opposite ends in the
order e_i=1,e_j=21,

    W_i=lambda_i*(gamma_i+1/(2i)),
    W_j=lambda_j*(gamma_j+1/[2(22-j)]).

For the REVERSED assignment the correct, mixed weights are

    W_i=lambda_i*gamma_i+lambda_j/(2i),
    W_j=lambda_j*gamma_j+lambda_i/[2(22-j)].       (13)

In every case W=W_i+W_j>(4/5)(lambda_i+lambda_j). Write
E0=W_i*sqrt(D_i)+W_j*sqrt(D_j). Then

    E0-||Z|| <(7/80)W.                           (14)

Indeed sqrt(D_i D_j)-C_ij<=21/22, E0>=W sqrt(85/11), and
||Z||>=W sqrt(80/11). Rationalizing gives the upper coefficient
(21/44)/[sqrt(85/11)+sqrt(80/11)]<7/80, using
sqrt(935)>61/2 and sqrt(880)>59/2.

For a common far end, put Aplus=A_i+A_j. Equation(11) gives

    lambda_i+lambda_j=ell*Aplus,
    Dcal=[ell*Aplus^2+(A_i-A_j)^2]/2,
    W>(4/5)*ell*Aplus, Aplus>4.

Hence sqrt(1+Dcal)/W<7/16: its squared bound is
[(ell+1)/2+1/16]/[(16/25)ell^2], decreasing for ell>=5 and
equal49/256 at5.

For opposite ends, the explicit inverse of H gives
5/6<lambda_i/lambda_j<6/5 for BOTH diagonal assignments. Let
L=lambda_i+lambda_j. Using each delta<=3/5 gives

    Dcal/L^2<=3/55,
    Aplus=(delta_r-1/2)lambda_i
                     +(delta_l-1/2)lambda_j<=L/10

in one assignment, with the diagonal names swapped in the other.
Thus L>40, W>32, and

    sqrt(1+Dcal)/W
       <sqrt(15/176+1/1024)<3/10.                (15)

Now (12), (14), and ||Z||-(1+Dcal-sqrt(1+Dcal))<sqrt(1+Dcal)
give

    W_i*(sqrt(D_i)-h_i)+W_j*(sqrt(D_j)-h_j)
       <(21/40)W      for the same far end,
       <(31/80)W      for opposite far ends.

Both contradict the two losses being at least61/100. This proves
the theorem, including finite far-m, dotted, and mixed finite/dotted
pairs, all native cap products, and both far-end assignments.

## 5. Limit of the result

Within this LITERAL complete two-contact family, weak central indices
have no adjacent pair at a whole A21 maximum. This is a simultaneous
actual-wall constraint, not an individual .61 floor. It does not
exclude nonadjacent weak directions, isolated profiles, other far or
near multiwhite contacts, or arbitrary spectator components. A whole
A_n obstruction still requires complete coverage and an energy
assembly; no uniform dimension ceiling is asserted.
