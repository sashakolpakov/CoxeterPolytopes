# Whole B_n: the two short-tail continuation layers

20 October 2026. Paper only, developed with root. Supporting geometric
theorem, not an ambient dimension bound. No scientific job, manuscript,
Lean change, or accepted-source edit.

**Theorem.** A vertex farthest from a strict interior point of a compact
Coxeter polytope cannot have WHOLE diagram B_n for n>=23.

Use the notation and unit normalization of the accepted source
`compact_whole_b_bulk_weak_packet_and_endpoint_pressure_20261020.md`,
SHA256 `e5a3c51858bb7b1e392b2828c874387be8d2044d5ae1be583a546b67eaddb54c`:

    x1>...>xn>0, ||x||<1, S_i=sum_(j<=i)x_j,
    omega_i=sqrt(2)*prefix_i, D_i=2i (i<n), omega_n=prefix_n,
    h_i=x.omega_i, c0=61/100.

An index i<n is weak if h_i<=sqrt(2i)-c0. The accepted source proves
the complete weak-packet inventory for 5<=i<=n-4, their nonadjacency,
index5 strong, and exclusion n>=26. The actual-projection and finite
classification inputs are the same; no Coxeter property of a metric
facet is assumed. The new ingredients below cover the next two layers.

## 1. Complete continuation at n-3 and n-2

The accepted left-white inventory, its affine exclusions and its edge
estimates remain valid at both layers. Its isolated-cap estimates also
remain valid: only the LEFT H4 segment is needed for label5. It gives

    v-w_minus=beta*(3omega_(i+1)/4+omega_(i-3)/4),

whose x-product is <=beta*h_i by concavity of S. This formula is valid
through i=n-2, since i+1<=n-1 is still a long node. The other isolated
estimates use the old node i, one simple neighbor and/or its positive
interior product; those remain available. Isolated3 is excluded by the
affine fork with two simple leaves and the right arm terminating in4;
isolated4 by the double4 path, and isolated6 by the6,3 path.

At i=n-3 the retained right diagram is B3. Its complete one-node finite
extensions are B4 by a3 at its long endpoint i+1, and F4 by a3 at its
SHORT endpoint n. No other finite extension occurs. If a left contact
also occurs, the connected rank-n white completion must be B_n, so the
F4 branch is sealed. The long-end branch is already in the accepted
inventory. For the extra short-end branch the white Schur data are

    delta=1/4, sigma=gamma+1/(2sqrt(2)), gamma=-c.n_i.

Gamma0 gives the exact affine F4 path i,i+1,n-1,n,c. Therefore
gamma>=1/2 and sigma/sqrt(delta)>=1+1/sqrt(2)>100/61, which is strong.
Consequently the accepted weak inventory and nonadjacency extend to
EVERY 5<=i<=n-3, with no change:

    v=gamma*omega_i+omega1/2,
    A^2=2i*gamma^2+2gamma-1/2,
    own label5, >=7, or dotted; 101/125<gamma<6/5.       (1)

At i=n-2 the retained right diagram is B2. Its one-node finite
extensions are precisely a3 at either endpoint, both B3. The long-end
branch is again accepted. In the new short-end branch with NO left
contact, delta=1/2 and sigma=gamma+1/(2sqrt(2)). Gamma0 gives the affine
F4 path i-1,i,n-1,n,c. Gamma4 gives the affine(4,3,3) star on
c,i,i-1,n-1; gamma6 gives the affine6,3 path on c,i,n-1. Every native
gamma>=cos(pi/5) is strong, because

    sigma/sqrt(delta)=sqrt(2)*gamma+1/2>100/61.

Thus the sole additional possible weak profile is

    gamma=1/2, v=(omega_(n-2)+omega_n)/2.              (2)

If the new short-end branch also contacts the left A_(i-1), finite
classification forces a simple3 at one of that path's endpoints and
no further contact: the full connected white diagram is B_n. Here

    delta=1/(2i),
    sigma>=gamma+1/(2i)+1/(2sqrt(2)),

so sigma/sqrt(delta)>sqrt(i)/2>100/61 for i>=21. This covers the
ZERO-own-contact branch as well; it has not been discarded as an
affine cycle. Therefore, for n23..25, a weak n-2 is either (1) or (2).

## 2. Profile (2) has an actual strong mate

Put i=n-2. Deleting n-1 in the old diagram and adding c from (2) gives
A_n; deleting n gives D_n. Both complete Grams are positive definite.
The acute positive-projection lemma makes these actual edge endpoints,
so the SAME c closes directions n-1 and n. This is not a formal
completion asserted without feasibility.

For its n-1 edge the exact Schur calculation gives

    delta=(i+3)/(4(i+1)),
    sigma=i/(2(i+1))+1/(2sqrt(2)),
    C=cosh L_(n-1)=[(1+sqrt(2))*i+1]/sqrt(i+3).

For 21<=i<=23 set q=sqrt(2i+2). Then C/q>79/50, C>9 and q<7.
For the first bound use 1/sqrt(2)>7/10 and square: the required
positive difference is 984i^2-19014i-17498, positive at21 and
increasing thereafter. Since sqrt(C^2-1)>C-1/C,

    q*(1-tanh(L/2))
      <q/[C+1/2-1/(2C)]
      <1/(79/50+4/63)<61/100.                       (3)

Thus n-1 is strong whenever (2) is needed for n23..25.

## 3. A reusable far-helper demand

Every weak packet (1) with 5<=i<=21 forces

    x1>53/125=.424.                                 (4)

Indeed the actual pressure table in accepted Section3 works at ANY
such i. Its finite5 H4 uses the nodes1,2,i; its7/8 negative core uses
1,i; its other rows use i and the dotted feet. All these nodes and the
two simple neighbors of i remain available through i=n-2. Concavity
therefore gives

    x1>sqrt(2)*[T-gamma*(q-.61)], q=sqrt(2i),
    A^2=q^2*gamma^2+2gamma-1/2,

with T=H(A,2/3-gamma^2) for5, A-sqrt(gamma^2-3/4) for7/8,
H(A,1-gamma^2) for finite>=9, and A-sqrt(1-gamma^-2) for dotted.
For finite labels, A-gamma*q>79/1000: its numerator2gamma-1/2>1.116
and its denominator A+gamma*q<14. Using H(A,delta)>A-sqrt(delta),
the lower bounds for T-gamma*(q-.61) are respectively

    .61*.809+.079-.111,
    .61*.900+.079-.253,
    .61*.923+.079-.323,
    .61*.939+.079-.344.

Each times sqrt(2)>7/5 exceeds.424. The elementary native intervals
and positive square-root bounds here follow by squaring the displayed
endpoints (cos(pi/9)>.939 for the last row).

For dotted gamma in(1,6/5), the bracket is decreasing in gamma and q.
For gamma-monotonicity, A'-q<(q^2+2)/(4q^3)<1/10 for q>=sqrt(10),
whereas (sqrt(1-gamma^-2))'>1 throughout this interval. For
q-monotonicity use A>gamma*q. Its infimum is consequently

    sqrt(2)*[sqrt(3119/50)-(6/5)*sqrt(42)+.732-sqrt(11)/6]>.424.

Positive rational seals are sqrt(3119/50)>7.8981,
sqrt(42)<6.480741, sqrt(11)/6<.552771 and sqrt(2)>1.4142.

## 4. Prefix-energy closure

Keep the accepted reference floors f_j=sqrt(j)-.61/sqrt(2), and its
positive-dual concave prefix reference. The accepted exact bounds give

    E21>1.009,
    E20=f5^2/5+sum_(j=5)^19(sqrt(j+1)-sqrt(j))^2>.99752,
    f5/5<.361.                                      (5)

The second bound uses the same rational list with its last entry
removed: .6514+.34612=.99752. Isolated omitted floors lose less than
1/3200 each. Two consecutive omitted floors lose less than1/800:
their three segment slopes have successive gaps<1/40, so their
squared variance is<2*(1/40)^2. Removed clusters have disjoint segment
interiors; their losses add and the surviving reference is concave.

For n25, both21 and22 lie in the enlarged range (1), so one is strong.
The accepted E21 argument applies unchanged. For n24, either21 is
strong, or22 is strong, or22 has profile (2) and23 is strong by (3).
The remaining alternative, adjacent far-weak21 and22, is impossible.
The same argument works for n23 when21 is strong or has profile (2).
In the latter cases the reference ends at22 or23. Besides isolated
holes, there can be only ONE two-hole terminal cluster: the new
profile (2) and the preceding far-weak index. A deliberately excessive
loss bound9/3200+1/800 is still less than.009, so (5) contradicts
||x||<1.

Finally, at n23 suppose21 is far-weak. Then20 is strong, there are at
most seven nonadjacent holes in6..19, and (4) applies. Fixing x1>.424
adds to the prefix energy through20 precisely
(5/4)*(.424-f5/5)^2, while preserving decreasing positive slopes at
1 and5. Thus the legitimate positive-dual lower bound is

    ||x||^2 > .99752+(5/4)*(.424-.361)^2-7/3200
             =1.00029375>1.

This exhausts n23..25. The unchanged accepted n>=26 theorem supplies
the remainder. No statement about n<=22, proper B components with
spectators, arithmetic lattices, or an ambient dimension ceiling is
made here. The still-open step is a complete usable tail pressure
at those smaller ranks, rather than assuming that the terminal4 is
harmless or importing an A_n plateau transfer.
