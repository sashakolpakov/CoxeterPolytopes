# Adjacent weak dotted/far3 opposites at a whole A21 maximum

20 October 2026. Own paper proof, no computation or execution.
This is an actual two-wall contradiction, not a dimension theorem.
It concerns the stated complete contact vectors only. Finite far-m
contacts and other nonisolated profiles are not discarded.

## 1. Statement and normalization

Let p be farthest from a fixed strict interior point of a compact
hyperbolic Coxeter polytope, and suppose its WHOLE finite diagram is
A21. Use unit outward old normals n_1,...,n_21 and N=22. Their inverse
Gram entries are

    D_i=2i(22-i)/22,
    C_ij=2i(22-j)/22, i<j.

Write the spatial farthest vector as x, with ||x||<1, and put
`h_i=x.omega_i`, where omega_i is the old dual. Its actual edge has
`h_i>=sqrt(D_i)*tanh(L_i/2)`.

Call i weak here if `h_i<=sqrt(D_i)-61/100`. Consider adjacent
central i,j=i+1 with5<=i<j<=17. Suppose their actual opposites are

    c_l=A_l*p-v_l,  c_l^2=1,
    v_l=gamma_l*omega_l+(1/2)*omega_(e_l),
    gamma_l>1, e_l in{1,21}, l=i,j.                 (1)

Thus each COMPLETE old contact vector has one dotted contact at its
own closing direction and one label3 at a far chain end. No hidden
additional contact or spectator is permitted in (1).

**Claim.** The two indices cannot both be weak, regardless of whether
their far ends agree. The proof constructs a genuine finite mixed
replacement using BOTH actual opposites, and compares it with p.

The actual projection and farthest comparison input is
`compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`,
SHA256 `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`.

## 2. Weakness bounds each dotted contact

Let m=l if e_l=1 and m=22-l otherwise, so m>=5. Projecting the far
contact off the retained white path gives exactly

    delta_m=(m+1)/(2m), sigma_l=gamma_l+1/(2m),
    A_l^2=D_l*sigma_l^2-delta_m,
    cosh L_l=sigma_l*sqrt(D_l)/sqrt(delta_m).       (2)

All central D_l lie in[85/11,11]. Put g=27/25 and c=61/100.
If gamma_l>=g, then sigma_l/sqrt(delta_m)>38/25. Indeed clearing
denominators reduces this to

    28m^2-188m+625>0,

whose discriminant is negative. On the other hand, for
q=sqrt(D_l) in[5/2,10/3],

    F_c(q)/q<38/25,
    F_c(q)=[1+(1-c/q)^2]/[1-(1-c/q)^2].

After positive denominators are cleared, the last comparison is
`-36400q^2+163602q-93025>0`. This concave quadratic is positive
at both enclosing endpoints:88480 and430835/9. Therefore (2) would
force `sqrt(D_l)*tanh(L_l/2)>sqrt(D_l)-c`, contradicting weakness.
Consequently simultaneous weakness implies

    1<gamma_i,gamma_j<g=27/25.                   (3)

## 3. Every native pair branch yields a finite actual replacement

Set R=c_i.c_j. The exact identity

    R=1-||v_i-v_j||^2/2+(A_i-A_j)^2/2             (4)

reduces pair control to a convex quadratic on the rectangle (3).
The old adjacent difference has square21/11<2 and
`|D_j-D_i|<=1`.

If the far ends agree, convexity in each gamma shows
`||v_i-v_j||^2<3`. At the two equal corners this is at most
`2g^2<3`; at a mixed corner it is at most
`2g+(g-1)[(g-1)*11+1]<3`. Hence R>-1/2. Distinct actual outward
Coxeter normals are nonpositive and have no product in(-1/2,0),
so R=0.

If the far ends differ, use e_i=1,e_j=21 first. With A=gamma_i,
B=gamma_j, Delta=D_j-D_i, direct expansion gives

    ||v_i-v_j||^2=1+2AB-2(A-1)(B-1)/22
                    +(A-B)^2 D_i+(A-B)(1-B)Delta.

Its four rectangle-corner bounds are3, `1+2g^2`, and
`1+2g+(g-1)^2*D_i` or the same with D_j. All are bounded by
`1+2g^2=2083/625`. Reversing the two far-end assignments decreases
this norm: the sign-changing cross term is positive because
`A*(4/22)+(A-B)*[2(22-2j)/22]>0`; its lower bound is
`4/22-(2/25)*(24/22)>0`. Thus in either orientation

    R> -833/1250 > -1/sqrt2.

The COMPLETE native alternatives are therefore R=0 or R=-1/2.
Every other finite Coxeter angle, a parallel pair and an ultraparallel
pair are excluded by this strict bound. No arithmetic field restriction
has been used.

Remove old i,j. If the far ends agree, the common retained projection
has square `t=(ell-1)/(2ell)`, where ell=i for the left end or
ell=22-j for the right end. The two-cap Schur matrix is

    H=I-t*ones(2,2), det(H)=1/ell>0.              (5)

If the far ends differ, their retained projections are orthogonal and

    H=[[delta_i,R],[R,delta_j]], delta_i,delta_j>1/2.
                                                               (6)

Both R=0 and R=-1/2 give positive definite H. The retained old Gram
is also positive definite. Thus the two caps and retained old walls
have a positive definite rank21 Gram. They are ALREADY actual facets,
so positive acute projection gives an ACTUAL finite replacement q.
This is not a conjectural diagram realization.

## 4. The far-end helpers cancel from the displacement

Put `lambda=H^-1(A_i,A_j)^T>0` and
`mathcalD=(A_i,A_j).lambda`. The actual projection has

    q=[(1+mathcalD)p-Z]/sqrt(1+mathcalD),
    Z^2=mathcalD*(1+mathcalD),
    x.Z>=1+mathcalD-sqrt(1+mathcalD).              (7)

In every case `Z=W_i*omega_i+W_j*omega_j` with positive weights.
For opposite far ends in the order e_i=1,e_j=21,
`W_l=lambda_l*sigma_l`. For a common left end the weights are

    W_i=lambda_i*gamma_i+(lambda_i+lambda_j)/(2i),
    W_j=lambda_j*gamma_j.                         (8)

The common-right case is the reversal. To check cancellation in (8),
the retained far-end projection is
`(omega_1-omega_i/i)/2`, so all omega_1 terms cancel. The opposite-end
case in this order cancels the two far-end terms separately. For the
reversed order e_i=21,e_j=1, put l=i, r=22-j, k_l=1/(2l), k_r=1/(2r).
Then the exact weights instead are

    W_i=lambda_i*gamma_i+k_l*lambda_j,
    W_j=lambda_j*gamma_j+k_r*lambda_i.             (8r)

Thus all weights stay positive, but the helper cancellation mixes the
two directions. This cancellation is
why retaining the full mixed displacement is stronger than an
individual height or edge floor.

Let W=W_i+W_j and `E0=W_i*sqrt(D_i)+W_j*sqrt(D_j)`.
For all adjacent central i,j,

    E0-||Z|| < (7/80)*W.                          (9)

Indeed `sqrt(D_i D_j)-C_ij<=21/22`, while
`D_i,D_j>=85/11` and `C_ij>=80/11`. Hence

    E0-||Z|| <= W*(21/44)/[sqrt(85/11)+sqrt(80/11)]
              <(7/80)*W.

The last strict inequality uses
`sqrt(935)>61/2`, `sqrt(880)>59/2`.

For a common far end, (5) gives

    lambda_i+lambda_j=ell*(A_i+A_j),
    mathcalD=[ell*(A_i+A_j)^2+(A_i-A_j)^2]/2.

By (8), W>=ell*(A_i+A_j), and
`(A_i+A_j)^2>=4*(85/11-1)=296/11`. Therefore

    sqrt(1+mathcalD)/W <7/20.                    (10)

For an explicit seal, the square of its upper bound is at most
`[(ell+1)/2+11/296]/ell^2`, decreasing for ell>=5, and at5 equals
899/7400<49/400.

For opposite far ends in the order e_i=1,e_j=21, rescale (6) by the
sigmas. Write
`eta_l=delta_l/sigma_l^2` and `a_l=A_l/sigma_l=sqrt(D_l-eta_l)`.
The rescaled H' has nonpositive off-diagonal and

    eta0=625/1458 < eta_l <1/2, H' Wvec=avec.      (11)

For the lower bound in (11), use gamma<g and
`g(g-1)m>1/4` for m>=5. For the upper use gamma>1 in (2).
The explicit two-by-two inverse shows

    max(W_i/W_j,W_j/W_i)
      <=g^2*sqrt[(85/11+1-eta0)/(85/11-1/2)] <5/4.

This last rational-square seal is
`16*531441*133093 <25*390625*115911`.
Consequently

    mathcalD/W^2<=41/162,
    W>=2(a_i+a_j)>=4sqrt(159/22),
    sqrt(1+mathcalD)/W
      <=sqrt(41/162+11/1272)<41/80.               (12)

The first inequality uses the nonpositive off-diagonal in H' and
`(W_i^2+W_j^2)/W^2<=41/81`; the second sums its two equations.
For the last seal, `41/162+11/1272=8989/34344<1681/6400`.

For the REVERSED order, the diagonal entries of H are
delta_r=1/2+k_r and delta_l=1/2+k_l, respectively. First let R=0.
Put s_i=gamma_i+k_r, s_j=gamma_j+k_l, eta_i=delta_r/s_i^2 and
eta_j=delta_l/s_j^2, and a_i=A_i/s_i,a_j=A_j/s_j. The same scalar
calculation as (11) gives eta0<eta_i,eta_j<1/2. If q_i=k_r/s_i,

    a_i^2=D_i-(2i/11)q_i+(21/11)q_i^2-eta_i,

and the j formula has i replaced by r and q_j=k_l/s_j.
Since q_i<=1/(2r+1), q_j<=1/(2l+1), l,r>=5 and l+r=21,

    D_i-1/2-32/121 < a_i^2 < D_i-eta0,
    D_j-1/2-32/121 < a_j^2 < D_j-eta0.

It follows that max(a_i/a_j,a_j/a_i)^2<6/5: substitute
D_min=85/11 and |D_i-D_j|<=1; the cleared positive difference is
`85/11-8-192/121+5*eta0>0`. Also eta_max/eta_min<6/5.
Set t_i=s_i*lambda_i,t_j=s_j*lambda_j. As R=0,
t_i=a_i/eta_i and t_j=a_j/eta_j, so their ratio lies between
25/33 and33/25, hence between3/4 and4/3. By (8r), W=t_i+t_j.
Therefore

    mathcalD/W^2 < (1/2)*(t_i^2+t_j^2)/W^2 <=25/98,
    W>4*sqrt(1685/242),
    sqrt(1+mathcalD)/W
        <sqrt(25/98+121/13480)<13/25.             (12r)

Every seal is rational after squaring; the last comparison is
`25/98+121/13480<169/625`.

For reversed R=-1/2, the inverse of H gives
5/6<lambda_i/lambda_j<6/5. Writing L=lambda_i+lambda_j,

    mathcalD/L^2 <=3/55, W>L,
    A_i+A_j=(delta_r-1/2)lambda_i
                     +(delta_l-1/2)lambda_j <=L/10.

Since A_i,A_j>sqrt(74/11)>5/2, W>50. Consequently
`sqrt(1+mathcalD)/W<sqrt(3/55+1/2500)<1/4`, also below13/25.

## 5. Contradiction and exact remaining scope

Writing `loss_l=sqrt(D_l)-h_l`, equation(7) gives

    W_i*loss_i+W_j*loss_j
      <=E0-[1+mathcalD-sqrt(1+mathcalD)].

Moreover `||Z||-[1+mathcalD-sqrt(1+mathcalD)]<sqrt(1+mathcalD)`.
For a common far end, (9),(10) bound the right side strictly by
`(7/16)W`. For ordered opposite far ends, (9),(12) bound it strictly
by `(3/5)W`; for reversed opposite ends, (9),(12r) bound it by
`(7/80+13/25)W=(243/400)W`. All contradict losses at least61/100,
since the weights
are positive. This proves the claim.

Thus the weak indices in this LITERAL dotted/far3 family have no
adjacent pair at a whole farthest A21 vertex. This does not exclude
all isolated weak directions, finite far5/simple3, finite far-m/simple3
with m>=7, near4 packets, other contact vectors or arbitrary spectators.
It proves neither a whole-A21 exclusion nor a new compact/general
dimension bound. A complete simultaneous coverage argument remains
necessary; no finite profile is silently promoted to a closed inventory.
