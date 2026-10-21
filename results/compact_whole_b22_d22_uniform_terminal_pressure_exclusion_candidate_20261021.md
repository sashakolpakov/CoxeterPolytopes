# Whole B22 and D22: one smaller-loss window and terminal pressure proof

21 October 2026. Written mathematical candidate, assembled by the
composition-audit lane from joint root, sharp-foundations and composition
derivations. Awaiting complete root and outside-author review. No scientific
job, TeX, Lean, manuscript edit or global status promotion is asserted.

**Theorem.** A farthest vertex of a compact hyperbolic Coxeter polytope
cannot have WHOLE finite diagram B22 or D22.

This is a whole-diagram theorem. In particular the new .58 weak inventory
below is NOT an arbitrary-spectator transfer of the old .61 inventory.
The argument completes a single long-prefix window through20, rather than
requiring a classification at the final B long node21 or short node22.

## 1. Exact geometric inputs and normalization

All paths in this table are under results/. Only their geometric
classifications, actual-point constructions and indicated scalar identities
are imported; no whole-B22 or whole-D22 exclusion is an input.

| Input | Source | SHA256 |
| --- | --- | --- |
| Whole classical finite white lists, affine exclusions, mixed replacement, native pressure table | `compact_whole_b_bulk_weak_packet_and_endpoint_pressure_20261020.md` | `e5a3c51858bb7b1e392b2828c874387be8d2044d5ae1be583a546b67eaddb54c` |
| B3/B2 finite continuation and actual B2 mate | `compact_whole_b_short_tail_extension_20261020.md` | `b6adc03a5950a60a1e263aeb79869785268f1bb7905e74ec47a784446c989c5c` |
| D bulk/D4/D3 finite lists and positive finite/negative-core feet | `compact_dn_large_components_short_tail_transfer_20261020.md` | `651a9d047d904f1b6496ba2664df5f313878466f67d011dcda144edef6a5620a` |
| Prefix dual interpolation and old rational increment list | `compact_classical_generic_prefix_energy_20261020.md` | `f637d7294d17e37581330ca508c914c23a08ce259845614b055e45815d2790ef` |
| Full last-D-center finite list, fork mates and H3 foot | `compact_last_d23_center_and_isolated_i_frontier_candidate_20261020.md` | `a77f88f5bd2c2c02588dffae7be4bc9ad96e83421bea2d525a3911c4659d996c` |
| Actual finite projections and positive feasible future points | `compact_farthest_vertex_mixed_replacement_inequalities_20261020.md` | `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84` |

The ordinary finite Coxeter classification and the previously established
compact/farthest-vertex geometric framework are explicit inputs. All new
threshold classifications and all applications at ranks below the former
tail statements are checked below.

Write z=E(p-x), ||x||<1, at a farthest vertex p. In B22 use

    n_i=(e_i-e_(i+1))/sqrt2 (i<22), n22=e22,
    omega_i=sqrt2*P_i (i<22), omega22=P22.

In D22 use the same n_i for i<22, n22=(e21+e22)/sqrt2,
omega_i=sqrt2*P_i for i<=20, and

    omega21=(P21-e22)/sqrt2, omega22=(P21+e22)/sqrt2.

Here P_i=e1+...+e_i. Put S_i=x.P_i and

    c=29/50, f_i=sqrt(i)-c/sqrt2, q_i=sqrt(2i).

Call i strong when S_i>f_i and weak otherwise. The coordinates are
positive and decreasing through21 in B, and satisfy
x1>...>x21>|x22| in D; hence the prefix sequence is strictly concave
through21 in both cases. For long nodes C_ij=2min(i,j).

An actual opposite wall is cwall=A*p-v, A>0, ||v||^2=A^2+1. A
positive finite old-normal correction w, with delta=1-||w||^2>0,
gives the ACTUAL finite-point pressure

    x.(v-w)>=H(A,delta),
    H(A,delta)=A/[1+sqrt(delta/(A^2+delta))]>A-sqrt(delta).       (1)

A positive combination xi=cwall+w with xi^2=-k<0, time coefficient
A>0 and all selected products nonpositive is feasible against EVERY
actual wall and gives x.(v-w)>=A-sqrt(k). Unselected wall products
are nonpositive term by term by acuteness. These are the pinned
actual-point results, not formal Gram completions or Coxeter metric faces.

## 2. Rechecking the smaller-loss bulk window

For an edge with C=cosh L=qR, the condition R>T(q) suffices to make
its direction strong, where

    T(q)=50/29-(q-29/50)/[q(2q-29/50)].                       (2)

T is increasing for q>=sqrt10. The following rational seals suffice:

    T(q)<83/50 for q<7,
    T(q)<801/500 for q<15/4,
    T(q)<159/100 for q<16/5,
    T(q)<33/20 for q<32/5.                                  (3)

These follow by substitution at the displayed right endpoints and
clearing positive denominators. We use q_i<=sqrt40<32/5 in the final
window. The first three slightly looser comparisons simplify the early
finite endpoints.

For B22, ordinary bulk is5<=i<=18; for D22 it is5<=i<=17. Deleting i
retains A_(i-1) and respectively B_(22-i) or D_(22-i). Their finite
one-node white extensions are exactly the pinned lists: on the left,
endpoint3/4, next-to-end3, or the admissible E positions3; on the right,
long-leaf3 and the finite D fork3 exceptional extensions; a two-piece
bridge is the full B22 or D22 endpoint bridge. The old zero-own affine
principals and impossible positive-index extensions are unchanged.
In particular no threshold change creates a new finite white diagram.
Because the diagram is WHOLE, there is no foreign inverse-Gram cost Q.

Here are new checks of all estimates that used the former loss .61.
For left near3 and own gamma>=phi/2, with phi=(1+sqrt5)/2, the least
edge ratio is

    R_i=[i(phi+1)-1]/sqrt(2i(i+1)),

which increases with i. At i6, R6=(6phi+5)/sqrt84>801/500; this
handles i6,7 by q<15/4. For i>=8, R_i>=R8=(8phi+7)/12>83/50.
These seals use phi>809/500. At i5, own>=7 gives R>5/3>83/50;
the only newly admitted near3 role is own5, treated in Section3.

For D positions,

    R_i>=(i*phi/2+1)/sqrt(2i),

which increases for i>=5. At i5 this is>159/100, using
gamma>809/1000 and sqrt10<317/100. At i>=6 it is at least
2927/1750>83/50, using sqrt12<7/2. E positions have residual<=1/4,
j>=3 and i<=8, so R>=2gamma+3/8>83/50. Far4/own3 remains the
exact affine F4 principal; near4/own3 has
C=(i-1)+i/sqrt2 and is stronger than (3). Its own>=5 roles are
stronger still. Right long3 and the extra B3 short3 give
R>=1+1/sqrt2>83/50. For D_s fork3 with s>=5,
R>(.8+.25)/sqrt(3/8)>83/50. Two-piece endpoint bridges have
smaller residuals and larger numerators. Thus all those rows are strong.

The isolated finite rows also need new loss seals. Finite5 has loss
<25/76+5/24<.58 by the pinned H4 foot. Finite7 has loss
<1000/4779+9/25<.57 by its old left-A2 negative core. Finite>=8
has worst value at i5,m8,

    sqrt10+sqrt2-4<.578<.58.

The exact singleton projection formula makes this monotone in the
remaining finite native interval. Isolated3,4,6 and the parallel null
case retain their exact affine/null exclusions.

For an isolated dotted cap put u=gamma^-2 in(0,1). Its actual foot
gives

    h_i=x.omega_i>sqrt(2i-u)-sqrt(u(1-u)).                    (4)

For i>=6 its loss is at most

    a_i*u+sqrt(u(1-u)),
    a_i=1/[sqrt(2i)+sqrt(2i-1)]<37/250.

Indeed sqrt12+sqrt11>6.76>250/37. Since
sqrt(1+(37/250)^2)<1011/1000, the maximum is<.5795<.58.
Only the isolated dotted profile at i5 is additionally admitted.

Consequently, except for that one i5 role and the near5 role excluded
next, every weak ordinary bulk opposite is precisely the PURE helper

    v=gamma*omega_i+omega1/2,
    A^2=2i*gamma^2+2gamma-1/2,
    own native5, finite>=7, or dotted.                       (5)

This is a new .58 statement. No foreign-far .61 estimate is transported
to .58, and no arbitrary-spectator .58 assertion is made.

## 3. The newly weak near5 is impossible independently

Suppose v=gamma*omega5+omega4/2, gamma=phi/2. Select the RIGHT old
H4 chain5,6,7 and put

    w=gamma*(3n5/2+n6+n7/2), delta=(5-3phi)/8>0,
    Y=v-w=(1/2+3gamma/4)*omega4+(gamma/4)*omega8.

The positive combination xi=cwall+w is orthogonal to those selected
old walls.
In particular its extra product with old4 is negative, not a hidden
chord in the selected H4 Gram. The finite projection is actual by (1).
Direct calculation gives

    A^2=(7+13phi)/2, Y^2=(33+49phi)/8,
    sum Y_j=sqrt2*(2+5phi/2).

Using 809/500<phi<81/50 and sqrt2<71/50,

    x.Y>18/5, Y^2<281/20, sum Y_j<43/5.                    (6)

For the demand, A>187/50 and sqrt(delta)<17/125 suffice.
The OLD .61 bulk helper classification/nonadjacency already supplies
a strong K11 or12, independently of index5 and all new .58 work.
This is the geometric part of the pinned generic M11 lemma, whose
rank hypotheses B>=16,D>=17 hold here. Its old floors satisfy
f11(old)>72/25 and f12(old)>303/100. Set Z=Y+(3/4)P_K.
If K11, x.Z>144/25 and Z^2<33.1375, with squared-demand excess
.0401. If K12, x.Z>2349/400 and Z^2<33.7, with excess .78625625.
Cauchy and ||x||<1 contradict both. Thus this new near5 role is
impossible outright; no new strong5 premise is used.

## 4. New helper nonadjacency and actual helper pressure

Every weak helper still has phi/2<=gamma<6/5. The upper cutoff is
reproved from (3): its far-path ratio squared at gamma>=6/5 is
at least4gamma-2>=14/5>(83/50)^2. In particular A>2.
For adjacent helpers i,j=i+1, their exact time-coefficient estimate is

    c_i.c_j>=1/2-gamma_i*gamma_j
                   +(gamma_i/gamma_j+gamma_j/gamma_i)/4
             >=1-gamma_i*gamma_j> -11/25> -1/2.

The actual native gap forces orthogonality. The pinned actual two-cap
replacement, retaining all old walls except i,j, has positive Schur
inverse and positive displacement. Its two weighted losses are

    <7/16+1/(4sqrt10)<.58.                                 (7)

Only i>=5, gamma>.8 and A>2 enter that computation. Hence two NEW
weak helpers cannot be adjacent. This does not assume that they were
old-.61-weak. In particular one of6/7 and one of15/16 is strong,
independently of the isolated dotted5 possibility.

The actual helper table gives

    x1>sqrt2*[T_native-gamma*(q_i-.58)],                    (8)

where T_native is H(A,2/3-gamma^2) for5,
A-sqrt(gamma^2-3/4) for7/8, H(A,1-gamma^2) for finite>=9,
and A-sqrt(1-gamma^-2) for dotted. These are the pinned positive
finite/negative-core feet, with every unknown wall retained.
Their vector identities use actual simple normals; at the D center
the identity omega21*=omega21+omega22=sqrt2*P21 is merely a spatial
prefix identity, with concavity through21, not a new old normal.

Three explicit consequences of (8) will be used:

    i5 weak helper:                 x1>.5455;
    any weak helper i<=20:          x1>.378;
    any weak helper i<=9:           x1>.45.                (9)

For the first, the old point-table expression at .61 is>.5965 on
the SAME gamma interval, so subtracting .03sqrt2*gamma leaves
>.5965-.03*(99/70)*(6/5)>.5455. This uses its scalar expression,
not inclusion in the former weak class.

For the second, finite5/7/8/>=9 respectively use the brackets
.58*.809+.079-.111, .58*.900+.079-.253,
.58*.923+.079-.323, .58*.939+.079-.344; each times sqrt2 is>.378.
For dotted the bracket decreases in gamma and q, so its infimum is

    sqrt(119/2)-(6/5)*sqrt40+.696-sqrt11/6.

The rational seals sqrt(119/2)>7.71362, sqrt40<6.32456,
sqrt11/6<.55278 and sqrt2>1.4142 give a product>.3781118.
For the third, at q<=sqrt18 the finite rows have A-gamma*q>.15:
this difference increases with gamma and decreases with q, and at
gamma=.8 its numerator is1.1 and denominator<7. The same finite
brackets with .15 in place of .079 exceed .45 after multiplication
by sqrt2. The dotted endpoint uses sqrt(1391/50)>5.274,
sqrt18<4.243, sqrt11/6<.553 and sqrt2>1.4, giving>.45556.

## 5. Early references without assuming strong5

Let d_i=sqrt(i+1)-sqrt(i) and

    E_K(c)=f5^2/5+sum_(i=5..K-1)d_i^2.

The old exact rational list gives E15(.61)>.92564,
E16(.61)>.94176, E19(.61)>.9847 and E20(.61)>.99752.
Changing the initial block to c=.58 adds

    .03sqrt(2/5)-.00357>.0154,                            (10)

using sqrt10>1581/500. A nonadjacent omitted floor j costs less
than1/[32(j-1)^3]. Thus all holes j>=6 cost altogether

    <(1/32)*integral_(4..infinity)t^-3 dt=1/1024.           (11)

We require the following reference contract for a strong K15 or16:
a vector W supported in1..K, with sum W_j=f_K; a certified pressure
x.W>A_W made from positive prefix and actual-point duals; endpoint-K
prefix coefficient>1/9; and effective energy

    L_W=2A_W-||W||^2>q_K,
    q15=.94006, q16=.95618.                              (12)

Effective energy is deliberately distinguished from ||W||^2. There
are three exhaustive cases for index5.

If5 is strong, take the concave prefix reference through all strong
floors5..K. Its squared norm is>E_K(.58)-1/1024, giving (12).
Its endpoint-K coefficient is at least d_(K-1)>1/8>1/9.

If5 is a weak helper,6 is strong. Interpolate through1,6,K with
values .5455,f6,fK. The slopes are positive and decreasing, and the
first-six energy is

    .5455^2+(f6-.5455)^2/5>.74385.

The tail energies for K15,K16,K19 are respectively>.225,>.2404,
>.2804. Hence the first two reference energies are>.96885,.98425,
well above (12). Endpoint-K coefficient is1/(sqrtK+sqrt6)>1/9.

Finally suppose5 is isolated dotted and weak. Set u=gamma^-2 in(0,1).
The positive feasible point cwall+gamma^-1*n5 gives the ACTUAL pressure

    x.U'>=T', U'=P5+(u/2)*(e6-e5),
    T'=[sqrt(10-u)-sqrt(u(1-u))]/sqrt2.                   (13)

Writing a=sqrt10-3<13/80, its loss is at most
a*u+sqrt(u(1-u))<(a+sqrt(1+a^2))/2<.588, since
sqrt(1+a^2)<2027/2000. Weakness forces u>.4: for u<=.4,
the coarser upper bound a*u+.5<.565 is already<.58.
Also T'>1.8185, using sqrt10>3.162 and sqrt2<99/70.

Choose strong J6 or7 and strong K15 or16. Let W0 have constant first-J
slope fJ/J and constant tail slope
t_JK=(sqrtK-sqrtJ)/(K-J). Put

    V=U'-(5/J)*P_J, W=W0+mu*V,
    mu=.21 for J6, mu=.2 for J7.

V has sum0 and support<=J, hence W0.V=0 and sum W_j=fK. Its squared
norm is at most77/150 for J6 and194/175 for J7, because
||V||^2=5-u+u^2/2-25/J and u>.4. The residual demands
Delta=T'-(5/J)fJ exceed .119 and .221 respectively. Subtracting
mu*(5/J) from the endpoint-J positive dual is legal: for J6 this
payment is .175, while f6/6-t_6K>.1808 even at K15; for J7 its
payment is1/7, below f7/7-t_7K>.165. Consequently

    L_W>||W0||^2+2mu*Delta-mu^2*||V||^2.

The respective gains exceed .027342 and .044057. The following
table records lower bounds for the first-J base plus gain and for
the tail; every entry is a rational square comparison.

| J | First-J base plus gain | Tail15 | Tail16 | Tail19 |
| --- | --- | --- | --- | --- |
| 6 | .69312+.027342=.720462 | .225 | .2404 | .2804 |
| 7 | .7139+.044057=.757957 | .188 | .2037 | .2445 |

Both early columns imply (12). The endpoint-K coefficient remains
t_JK>1/9; the residual did not change it. The same construction is
legal for K19 or20, since its endpoint-J payment only decreases.
The tail grows with K, so its last column will later exclude this
entire weak5 case. No strong5 or late endpoint was presumed to obtain
the early references.

## 6. One flat-prefix pressure comparison

Suppose an actual pressure has a normalized vector Z whose first16
coordinates are1, x.Z>G and ||Z||^2<=R. For either strong K15 or16,
put W*=W+(Z-P_K)/9, using Section5's applicable reference.
The added residual vanishes on support of W, so its cross product
with W is exactly0 and its norm squared is at most R-K. In the pressure
dual reduce the endpoint-K coefficient by1/9, then add x.Z/9.
This payment is legitimate by (12), and gives

    ||x||^2>L_W+2(G-fK)/9-(R-K)/81.                     (14)

One must NOT infer x.(Z-PK)>=G-fK merely from SK>=fK; (14) uses
the positive-dual payment just described. Nor is L_W silently
identified with ||W||^2 in the isolated5 case.

Use f15<3.463 and f16<3.590. Direct rational substitution shows that
EACH following row contradicts ||x||<1 in BOTH K cases:

| G | R | Minimum excess above1 in (14) |
| --- | --- | --- |
| 3.956 | 19 | >.00023 |
| 4 | 19.7 | >.0013 |
| 4.05 | 20 | >.0087 |
| 4.1 | 20 | >.0198 |
| 4.15 | 21 | >.0186 |
| 4.4 | 22 | >.0618 |

All applications of (14) below have first16 coordinates constant
alpha>0; Z is their displacement divided by alpha. The D4 own5
row instead uses the already pinned OLD first12-coordinate certificate,
as stated explicitly in Section8. Some rows have strictly smaller R
or larger G, and therefore satisfy the displayed table as well.

## 7. B22 continuation through20

At i19 the right B3 white extensions are long3 (old B bulk) or short3
making F4. The latter has residual1/4 and numerator gamma+1/(2sqrt2).
Own0 is the exact affine F4 principal; otherwise its ratio is at
least1+1/sqrt2>83/50. A left bridge is full B22 and strong. Thus the
weak inventory at19 is (5).

At i20 the right B2 extensions are simple3 at either endpoint.
The long-end branch is already covered. For the new short-end branch,
residual=1/2, numerator=gamma+1/(2sqrt2). Own0/4/6 have the pinned
affine F4, (4,3,3)-star and6,3 principals. Own>=7 is .58-strong.
The two remaining own values are3 and5. A left bridge is full B22
with residual1/(2i) and numerator at least
gamma+1/(2i)+1/(2sqrt2), directly strong even at own0.

For own5, put gamma=phi/2 and select the LEFT H4 old path20,19,18.
Its positive correction has delta5=(5-3phi)/8 and, with
alpha=sqrt2*gamma+1/2,

    A^2=20alpha^2-1/2, Y^2=A^2+delta5<20alpha^2.

A/alpha>4.45 and sqrt(delta5)/alpha<.083 give G>4.3,R20,
contradicting (14).

For own3 the cap is v=(omega20+omega22)/2. Replacing old21 gives
a positive definite full A22 actual endpoint. Its exact mate cosine is

    C=(21+20sqrt2)/sqrt23>245/24.

It forces x.P21>=sqrt21*sqrt((C-1)/(C+1))>83/20=4.15.
For the last comparison 1511*(245/24)>15289 suffices. Thus
(G,R)=(4.15,21), again excluded by (14). We did not call this mate
.58-strong; its actual pressure is what was used.

Every weak B22 index5..20 is therefore a pure helper, apart from
the still-permitted isolated dotted5 case.

## 8. D22 continuation through18 and19

At i18 the right D4 extensions are simple3 at any leaf. Long-leaf3
and full endpoint bridges are directly strong. Fork3 has residual1/2
and numerator gamma+1/4; own0/3/4/6 are the old exact affine
principals. Own>=8 is .58-strong, since
sqrt2*(cos(pi/8)+1/4)>33/20>T(q18). Retain BOTH own5 and newly
weak own7 rather than importing the old native>=7 strong claim.

For own5, the LEFT H4 correction of Section7 has now
alpha=sqrt2*(gamma+1/4), alpha^2=phi+5/8 and
A^2=18alpha^2-1/2. Its residual ratio is<.01, so
A/alpha>4.2 and sqrt(delta5)/alpha<.1 give G>4.1,R<18.
Here the left correction reaches coordinate16, so it would be wrong
to call the vector flat through16 and apply (14) at K16. Instead use
the pinned OLD early-prefix certificate (last-D source Section2):
one old strong K11/12 has respectively q>.8472225/.8689525 and
fK<2.886/3.034, and its endpoint coefficient exceeds1/7. Its old
strong5 bootstrap uses only the complete old bulk through17, with
no new tail premise. The first12 coordinates here are flat. Therefore
q+2(4.1-fK)/7-(20-K)/49>1.0102 in both cases. This excludes own5
independently of all new .58 references.
For own7, .9<gamma<.902, use the positive future point

    xi=cwall+gamma^-1*n18+(2gamma)^-1*n17,
    k=1-3/(4gamma^2)>0.

The cwall and n17 products are0; the n18 product is
-gamma+3/(4gamma)<0. With the same alpha,
A^2=18alpha^2-1/2, Y^2=A^2-k<18alpha^2.
A/alpha>4.22 and sqrt(k)/alpha<.18 give G>4.04,R18.
The own7 pressure dominates rows of (14). Hence18 has only (5).

At i19 the right D3=A3 extensions are middle3, end3, end4, end5.
Middle3 and end4 retain their affine/strong alternatives with ratio
at least1+1/sqrt2. End5 is sealed H4 and has a very large ratio even
at own0; it is not discarded as affine. Full left bridges have
residual1/(2i), or1/8+1/(2i) in the end3 A22 case, and are directly
strong after the old own0 affine cycle is removed. No foreign cost
exists. Thus only end3 with no left contact remains, with

    residual=5/8, numerator=gamma+1/4.

Own0/3/4/6 retain the exact affine principals. Own5 is the left H4
case, A^2=19alpha^2-5/8, giving G>4.1,R19. For all remaining finite
or dotted contacts, gamma>.9. Values gamma>=53/50 are .58-strong:
their ratio squared is at least (131/100)^2*(8/5)>(33/20)^2.

For .9<gamma<53/50, select each of the two old A2 sides:

    w_+-=gamma*(4n19/3+2n_side/3),
    side=18 or20, d=4gamma^2/3-1>0.

Each cwall+w_+- has old selected products0, cwall product=-d,
norm=-d and future time A. All coefficients are positive; the
fork helper's additional product is negative. Average their actual
pressures to obtain Y=v-(w_++w_-)/2 and

    x.Y>=A-sqrt(d), alpha=sqrt2*(gamma+1/4),
    A^2=19alpha^2-5/8,
    Y^2=A^2+1-14gamma^2/9<19alpha^2.

The ratio d/alpha^2 increases with gamma. If gamma<=26/25,
A/alpha>4.33 and sqrt(d)/alpha<.365. If gamma>=26/25,
A/alpha>4.337 and sqrt(d)/alpha<.381, using gamma<53/50.
Thus in the entire interval G>3.956,R19. These are direct rational
squares at .9,26/25,53/50 and give the first row of (14).
Hence19 too has only (5).

## 9. Complete finite continuation at the D22 center20

Deleting center20 retains A19 and two orthogonal fork nodes21,22.
Let gamma be own contact, beta the one possible left contact at j,
and delta,epsilon the fork contacts. The exact data are

    Delta=1-q_left-delta^2-epsilon^2,
    sigma=gamma+beta*j/20+(delta+epsilon)/2,
    C=sqrt40*sigma/sqrt(Delta).                           (15)

Finite white classification is the pinned last-D-center list with
left rank19 in place of20: endpoint3/4 or next-to-end3 on the left;
with both forks only the full endpoint D22 bridge; with one fork
only a left endpoint3 and fork3/4. Internal attachments have two
forks or long internal4 and are not finite. The bridge residuals
are1/40 or11/40. Their nonaffine ratios are>83/50, including own0
where the endpoint bridge is directly strong. The remaining own0
left/fork3 case is the exact all3 affine cycle; own3 is an affine
triangle. The finite tree list has not changed at this left rank.

With no fork contact the complete Section2 left inventory applies.
For isolated5 at this endpoint, the actual left H4 correction uses
the spatial identity omega21*=omega21+omega22=sqrt2*P21;
concavity through21 proves its old foot estimate. No virtual normal
is introduced. With no left contact, finite white roles are exactly
one fork with any finite native label, or two forks (3,3),(3,4),(3,5),
up to interchange. We now exhaust these roles.

For one fork and own0, fork3/4/6 are the exact affine E8, star and
6,3 principals. For fork5 the actual H4 on f,20,19 has

    alpha=delta/sqrt2, Y^2/alpha^2=19,
    A^2/alpha^2=6+8phi, sqrt(delta5)/alpha<15/56.

Thus G>4.08,R19. The last residual seal uses phi>79/49,
not merely phi>8/5: sqrt5>109/49 follows from
5*49^2-109^2=124. This retains the corrected premise in the input.

For own0 and fork>=7 put b=delta^2 in(.81,1). The positive point

    xi=cwall+delta^-1*n_f+2(delta^-1-delta)*n20

has products0 with cwall,n_f, product(3-4b)/(2delta)<0 with n20,
norm -k with k=(1-b)(4b-3)/b, and positive future time A.
Its displacement has first16 coordinates alpha=delta/sqrt2 and

    R(b)=30-16/b+6/b^2,
    G(b)=sqrt(22-2/b)-sqrt(-8+14/b-6/b^2).                (16)

For b<=.9, R<19.7, first root>4.41 and second root<.41,
giving (G,R)=(4,19.7). For b>=.9, R<=20, first root>4.44 and
second root<.39, giving (4.05,20). The R derivative has numerator
16b-12; the second radicand has derivative numerator12-14b,
with maximum1/6 at b6/7. The sharper b>=.9 bound is4/27.
Both rows contradict (14).

For own3, fork3 is an affine triangle. Fork4/5 makes the SAME cap
an actual contacted-fork edge mate with full positive definite D22
replacement and

    C_f=(20+22delta)/2>17.

Its pressure has Z=sqrt2*omega_f, R22 and
x.Z>sqrt22*sqrt(16/18)>4.4. Fork>=6 is directly strong.
Own4 has the exact affine center/left/otherfork/c principal.

For own>=5 and fork>=4, direct ratios remove every case except
own5/fork4. That exceptional new role has the SAME left-H4 data as
Section7's B2 own5/short3: alpha=sqrt2*gamma+1/2,
A^2=20alpha^2-1/2, G>4.3,R20. It is therefore excluded, rather
than incorrectly importing the former ratio phi/sqrt2+1/2 as
.58-strong. Own6 has the old6,3 affine principal.

For fork3 and own5, the left H4 has
alpha=sqrt2*(gamma+1/4), A^2=20alpha^2-3/4, giving G>4.1,R20.
For all remaining finite>=7 or dotted own values, gamma>=6/5 is
directly strong. For .9<gamma<6/5 use

    xi=cwall+gamma^-1*n20+(2gamma)^-1*n19,
    k=1-3/(4gamma^2), A^2=20alpha^2-3/4,
    Y^2=A^2-k<20alpha^2.

Its selected products are as in Section8. The bounds A/alpha>4.44
and sqrt(k)/alpha<.34 follow from gamma>.9 and
k/alpha^2<=575/5046<.34^2. Hence G>4.1,R20, again excluded.

Finally consider two forks. Own3 gives an affine triangle, and
own>=4 has ratio>=1+1/sqrt2>83/50. At own0, pair3/3 is an affine
square. Pair3/4 has an ACTUAL full B22 mate at the3-contacted fork,
with C_f=(40+22sqrt2)/4>17, hence the same (4.4,22) pressure.
Pair3/5 has the actual finite H3 foot

    alpha=phi^2/(2sqrt2), Y=alpha*P20,
    A^2=20alpha^2-kappa, kappa=(2-phi)/4,
    q=kappa/alpha^2=2phi^-6<1/8.

It gives x.Y/alpha>=sqrt(20-q)/(1+sqrt(q/20))>
4.45/1.08>4.1, with R20. This is the last row family.

Every remaining non-helper role has contradicted (14). Thus the
new weak inventory (5) is complete through20 in D22 as well,
apart from the separately retained isolated dotted5.

## 10. Closing the window without a strong5 assumption

Indices19 and20 are now pure-helper-or-strong, so nonadjacency gives
a strong K19 or20. First eliminate the two possible weak5 profiles.
If5 is a weak helper, Section5's interpolation has energy
>.74385+.2804>1.024. If5 is isolated dotted, Section5's SAME
residual reference with K19 or20 has effective energy greater than

    .720462+.2804=1.000862       (J6),
    .757957+.2445=1.002457       (J7).

Its endpoint-J payment is legal because the tail coefficient decreases
with K, and the tail contribution increases with K. Both contradict
||x||<1. Thus5 is strong, proved only AFTER the terminal inventories.

If20 is strong, (10)--(11) give

    ||x||^2>E20(.58)-1/1024>.99752+.0154-1/1024>1.

Otherwise20 is a weak helper and19 is strong. Here

    E19(.58)>.9847+.0154>1.00010,
    f5/5<.366.                                          (17)

If some weak helper k lies in6..9, (9) gives x1>.45. Replacing the
first-five constant slope by .45 followed by four equal slopes
adds (5/4)(.45-f5/5)^2>.00882. All slopes remain positive and
decreasing at the junction with the later prefix reference. Even the
full loss1/1024 in (11) leaves total energy>1.

If no weak index lies in6..9, all omitted floors are j>=10 and
their total loss is

    <(1/32)*integral_(8..infinity)t^-3 dt=1/4096.

The weak20 helper gives x1>.378, adding
(5/4)(.378-f5/5)^2>.00018 with the same legitimate first-five dual
adjustment. Thus

    ||x||^2>1.00010+.00018-1/4096>1.0000358>1.             (18)

This exhausts every branch and proves the stated WHOLE B22/D22
obstruction. Every correction was charged once on the same vector.
The near5 H4, early isolated-dot residual reference, later flat-prefix
pressure and final first-coordinate gains occur in the disjoint proof
branches explicitly indicated above; no overlapping gains were added.

The application to a new global compact ceiling additionally requires
the separately accepted dimension22 forest exhaustion and both rank21
isolated-I couplings. They are not premises of this whole-diagram
theorem, and no new global ceiling is promoted by this candidate note.
