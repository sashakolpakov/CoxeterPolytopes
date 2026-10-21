# Whole B_n: bulk weak packets and a high-rank exclusion

20 October 2026. Paper reasoning only; no scientific execution.
Supporting geometric results, not a dimension bound. Awaiting review.

Let p be a farthest vertex from a strict interior point of a compact
Coxeter polytope, with WHOLE finite diagram B_n, n>=21. Use unit old
normals (e_i-e_(i+1))/sqrt(2), i<n, and e_n. The farthest vector has

    x_1>...>x_n>0, ||x||<1,
    omega_i=sqrt(2)(e_1+...+e_i), D_i=2i (i<n),
    omega_n=e_1+...+e_n,
    C_ij=2min(i,j) (i,j<n).

Put h_i=x.omega_i and c0=61/100. An index is weak if
h_i<=sqrt(2i)-c0. Here BULK means exactly 5<=i<=n-4.

**Conclusions.** Every weak bulk opposite has precisely the old contacts
own m at i and simple3 at1, with m5, m>=7, or dotted. No two weak bulk
indices are adjacent. In addition, index5 is strong. Consequently a
whole B_n farthest vertex is impossible for every n>=26.
The proof is geometric and permits every actual native angle.

## 1. Complete actual bulk inventory

Deleting i leaves A_(i-1)+B_(n-i), whose ranks are at least4.
Finite classification permits: no white contacts; one contact to the
left A path (end3/4, next-to-end3, or an admissible E position3);
one right contact, necessarily simple3 at i+1; or two simple3 contacts
joining both pieces into B_n. No finite exceptional diagram contains
the right B component and the additional root, since its rank is>=5.
A tree has at most one contact to each old component; a two-component
completion of rank n>=21 is necessarily B_n.

All zero-own-contact branches are excluded. A single left end3 gives
finite B_(n+1); an internal left3 contains the affine fork with two
simple leaves and terminal4 on its long arm. Left end4 gives the
double4 affine path to n. A sole right3 similarly gives an affine fork.
A two-white completion gives an all3 affine cycle. A finite full
extension has impossible positive index n+1.

If c has no right contact, its own4 would give the EXACT affine path
c--4--i--3--...--4--n, irrespective of its other left contacts.
Own6 gives the affine6,3 path using i+1. For a left white3, own3
gives an all3 cycle. Thus only own5, own>=7, or dotted remain in
that case. A left FAR white4 at1 with own3 contains the exact affine
F4 path old3,old2,old1,c,oldi. This is not an admissible exceptional
packet. A left NEAR white4 at i-1 with own3 is allowed, but is strong
by its own edge estimate below; no An mate estimate is imported.

Let gamma=-c.n_i. For a left contact at j<i, the white Schur data are

    white3: delta=1-j(i-j)/(2i), sigma=gamma+j/(2i);
    white4: delta=1/i, sigma=gamma+j/(sqrt(2)i);
    cosh L_i=sqrt(2i)*sigma/sqrt(delta).          (1)

For a right simple3, its white cost is1/2 and its contribution to sigma
is1/2. Its own3 makes an affine triangle, hence gamma>=1/sqrt(2),
and sigma/sqrt(delta)>=1+1/sqrt(2)>100/61. Extra left contacts only
increase this lower bound. Such a direction is strong, since
q(1-tanh(L/2))<q/cosh L.

The same bound handles every left white4 except the already forbidden
far4/own3. The near4/own3 case has the exact value

    cosh L_i=(i-1)+i/sqrt(2),

so is strong for i>=5. For left near3/own>=5, the ratio in(1) is
increasing in i; from i8 onward it exceeds100/61. At i5 it exceeds
38/25 and q=sqrt(10)<10/3. At i6,7 it exceeds8/5 and q<=sqrt(14)<4.
The elementary edge comparisons F(q)/q<38/25 on [5/2,10/3] and
F(q)/q<8/5 on [3,4], where
F(q)=[1+(1-c0/q)^2]/[1-(1-c0/q)^2], verify these finite endpoints.
For D positions use delta=2/i,sigma>=gamma+1/i: i5 uses the first
comparison, i6 the second, and i>=7 already exceeds100/61.
For E positions i<=8, delta<=1/4 and gamma>4/5 give the second
comparison. These inequalities use only cos(pi/5)>101/125.

Isolated bulk opposites are strong by the actual spherical/Lanner-foot
estimates in the accepted A21 proof, which need only D>=85/11 and two
old simple arms of length>=3. Those arms are available here. The exact
same local affine exclusions remove isolated3,4,6 and parallel contact.
Consequently the ONLY weak bulk profile is

    c_i=A_i p-v_i,
    v_i=gamma_i omega_i+omega_1/2,
    A_i^2=2i gamma_i^2+2gamma_i-1/2.             (2)

No right-tail contact has been omitted.

## 2. Native gap and actual two-cap replacement

In(2), put f=101/125, g=6/5. Weakness implies

    f<gamma_i<g, A_i>2.                         (3)

Indeed F(q)/q<1/c0 for q>c0. At gamma=g, the far-path ratio obeys
sigma^2/delta>=4g-2=14/5>(100/61)^2, since its difference numerator
is [2i(g-1)-1]^2. Thus gamma>=g is strong at EVERY rank.

For adjacent weak i,j=i+1, put a=gamma_i,b=gamma_j. Weighted AM--GM
of their time coefficients gives the exact lower bound

    c_i.c_j >= 1/2-ab+(a/b+b/a)/4 >= 1-ab > -11/25 > -1/2.

The actual native gap forces c_i.c_j=0. After retaining all old walls
except i,j, their common helper projection has square t=(i-1)/(2i).
The Schur matrix is I-t*ones(2,2), positive definite with entrywise
positive inverse. Its projection is an ACTUAL finite vertex.

The positive mixed displacement has form Z=W_i omega_i+W_j omega_j.
Writing W=W_i+W_j, the same exact calculation as the accepted adjacent
packet gives sqrt(1+Dcal)/W<7/16; it uses only i>=5, A_i>2,
gamma_i>4/5. Also

    W_i sqrt(D_i)+W_j sqrt(D_j)-||Z|| < W/(4sqrt(10)).

Here C_ij=D_i>=10 and sqrt(D_iD_j)-C_ij<1. The two-loss weighted
average is therefore less than 7/16+1/(4sqrt(10))<61/100,
contradicting both indices weak. This proves nonadjacency.

## 3. Actual core pressure excludes weak5

Suppose5 is weak. Then6 is strong, and at least one of16,17 is strong.
All those indices are bulk, even at n21. Put g5=sqrt(10)-.61.
The actual pressures for the complete native range in(2) are:

|own contact|u with x.u>=T|T|
|---|---|---|
|5|gamma(omega4+omega6)/2+omega3/6|H(A,2/3-gamma^2)|
|7,8|gamma(omega4+omega6)/2+omega2/4|A-sqrt(gamma^2-3/4)|
|finite>=9|gamma(omega4+omega6)/2+omega1/2|H(A,1-gamma^2)|
|dotted|(gamma-gamma^-1)omega5+(omega4+omega6)/(2gamma)+omega1/2|A-sqrt(1-gamma^-2)|

Here A^2=10gamma^2+2gamma-1/2 and
H(A,delta)=A/[1+sqrt(delta/(A^2+delta))]. These are exactly the actual
finite/negative-core and dotted-foot constructions of the accepted
A21 pressure table, with the far helper on the LEFT. Concavity gives
h4+h6<=2h5, h2<=2h1, h3<=3h1. Therefore EVERY row implies

    x1=h1/sqrt(2) > sqrt(2)*(T-gamma*g5).        (4)

The right side exceeds1193/2000=.5965 throughout the complete range.
For finite5,7,8 it in fact exceeds .8,.7,.65 respectively, using
the elementary intervals .809<cos(pi/5)<.81,
.9<cos(pi/7)<.902 and .923<cos(pi/8)<.9239.
For finite>=9 the difference H(A,1-gamma^2)-gamma*g5 is increasing
from gamma=.939: A'>sqrt(10), H_A>9/10, H_delta<0, so its derivative
is positive. At .939, A>3.1929, sqrt(A^2+delta)>3.211 and
sqrt(delta)<.344 give H>2.88 and hence the right side>.65.

For dotted gamma in (1,6/5), the difference in(4) is decreasing:
A'<13/4, the derivative of sqrt(1-gamma^-2) exceeds1, and g5>5/2.
Its infimum is at6/5. There it is exactly

    sqrt(2)*(sqrt(163/10)-sqrt(11)/6-(6/5)(sqrt(10)-.61)) > .5965.

For a rational seal use sqrt(163/10)>4.03732,
sqrt(11)/6<.552771, sqrt(10)<3.162278 and sqrt(2)>1.41421.
Each is verified by squaring positive endpoints.

Let S_k=sum_(j<=k)x_j. Strong6 gives
S6>sqrt(6)-.61/sqrt(2)>2.01815. If16 is strong then
S16>4-.61/sqrt(2), whose difference from the preceding floor is
4-sqrt(6)>1.55051. The positive minimum-energy prefix interpolation
through1,6,16 has slopes .5965, (S6_floor-.5965)/5,
(4-sqrt(6))/10, then0, in strictly decreasing order. Thus

    ||x||^2 > .5965^2+(2.01815-.5965)^2/5+1.55051^2/10
            =1.00043812051>1.

If17 instead is strong, replace the last contribution by
(sqrt(17)-sqrt(6))^2/11, which exceeds
(4-sqrt(6))^2/10. The contradiction persists. This proves5 strong.

## 4. Whole B_n is excluded for n>=26

Now both21,22 are bulk. Choose M21 if21 is strong and M22 otherwise;
nonadjacency makes22 strong in the latter case. Let
f_i=sqrt(i)-.61/sqrt(2), the corresponding lower floors on S_i.
The minimum-energy reference curve through0 and f_5,...,f_M is
linear up to5, concave thereafter, then constant after M. Its norm is

    E_M=f_5^2/5+sum_(i=5)^(M-1)(sqrt(i+1)-sqrt(i))^2
       >=E_21>1009/1000.                         (5)

A rational seal uses f_5^2/5>6514/10000 and
(sqrt(i+1)-sqrt(i))^2 >=1/(4i+2). For i5,...,20 these reciprocals
have respective lower bounds in units1/100000

    4545,3846,3333,2941,2631,2380,2173,2000,
    1851,1724,1612,1515,1428,1351,1282,1219.

Their sum is35831/100000. Thus E_21>1.00971, more than (5) requires.

Both reference endpoints5,M are strong. Omit the weak internal floor
constraints; at most eight are omitted and none are adjacent. Since
(sqrt(t))''=-1/(4t^(3/2)), its absolute curvature on[5,M] is less
than1/40. Removing one isolated floor reduces this prefix-coordinate
energy by less than (1/40)^2/2=1/3200. The affected segment interiors
are disjoint, so losses add. The surviving curve is concave, with
nonnegative dual coefficients; it remains the actual lower-constraint
minimum. Therefore

    ||x||^2 >= E_M-8/3200 > 1009/1000-1/400 >1.

This proves the whole-B_n exclusion for n>=26 without assumptions on
the last four opposite profiles or any arithmetic hypothesis.

## 5. Exact remaining obligation

This does NOT exclude whole B_n for21<=n<=25 (or below21). Even all
scalar floors from5 through n-4 have, at n21, minimum prefix energy
below1; the nonadjacent missing
floors weaken that relaxation further. The still-needed mechanism is
simultaneous pressure from the distinct actual bulk caps and from the
last four directions. The terminal triples are not treated as automatically
finite or as inheriting the A plateau transfer. No complete native cap
configuration or compact polytope is asserted by that scalar relaxation.

Sources: the accepted A21 proof and its actual-projection imports remain
unchanged. Its SHA256 is
`3abe11bb70ce677069279462a91111ed6b8ab8889e3ec9e656ac71c1094200fb`;
the general actual-point input is
`compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`,
SHA256 `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`.
