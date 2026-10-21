# A positive fraction of large three-face carrier excess in dimension twelve

Date: 2026-10-18. Written Lorentz/Schur argument with a tiny exact rational
interval certificate. No optimizer, sampling proof, Lean verification or
new dimension ceiling.

## Theorem

For any finite-volume right-angled twelve-polytope, let f be its finite
vertices and

    E=sum_non-P3 (nu3-10)
     =B-10a3+10a2+28160c,

where nu3 is facet-plus-cusp sum and B counts canonical P3 faces. Then

    E >=f/5390.

The same method proves that at least f/385 finite vertices are incident
to a three-face with nu3>=11. This is a genuine global excess fraction,
extending the earlier existence and maximum-energy local-defect theorems.
It is too small to bridge the current triangle-density lower and upper
bounds. No dimension reduction is claimed. If f=0 the excess inequality
is trivial. The present certificate is for dimension twelve only.

## Actual moves and small-carrier rows

Fix one future unit timelike u for the whole polytope. At finite vertex q
put H(q)=-u.q>0 and g(q)=1/H(q). Consider the finite undirected weighted
graph consisting of:

* finite/finite edges whose cosh length is sqrt(2), with weight one;
* actual P3 oppositions, with weight 1/20.

Every move has uniformly bounded hyperbolic distance. The two edge types
have distinct distances, and the reverse of every actual P3 move has the
same weight. Write d(q) for weighted degree and

    R(q)=sum_neighbors w(q,q')g(q')/g(q).

Summing reciprocal pairs gives sum_q[R(q)-d(q)]>=0.

Call a finite vertex small if all its incident three-faces have nu3<=10.
The established P3/Q8/R7 classification then applies locally. As in the
conditional dimension-eight theorem, every finite edge at such a vertex
has cosh length sqrt(2), and every triple of its ideal directions is an
actual P3. If k directions end at ideal vertices and m=12-k at finite
vertices, its weighted degree is d=m+binomial(k,3)/20.

In the orthonormal frame at q put p_i=-u.n_i/H(q), so sum p_i^2<1. The
finite and P3 energy ratios are sqrt(2)-p_j and 2-sum_S p_i. Replacing
all p_i by their absolute values can only increase R. Increasing a
nonnegative coordinate until the norm reaches one also increases R,
and the denominators remain positive (at least sqrt(2)-1 and 2-sqrt(3)).
It therefore suffices to prove the uniform row inequality on nonnegative
unit vectors. The certificate below proves

    R(q) <=d(q)-1/8                 at every small vertex.

## Moment bounds and analytic maximization

For any positive list r with mean mu and minimum at least a>0, the exact
identity gives

    average(1/r)=1/mu+average[(r-mu)^2/(mu^2 r)]
                 <=1/mu+variance(r)/(mu^2 a).

For the nonnegative unit vector put x=sum_finite p_j^2 and 1-x for the
ideal squared norm. Let sigma_f and sigma_i be the group sums. Then

    sqrt(x)<=sigma_f<=sqrt(mx),
    sqrt(1-x)<=sigma_i<=sqrt(k(1-x)).

For the finite group the reciprocal bound is

    m/mu_f+[x-sigma_f^2/m]/[a_f mu_f^2],
    mu_f=sqrt(2)-sigma_f/m,       a_f=sqrt(2)-sqrt(x).

For the uniform triples of the ideal group (k>=3), put
C0=3(k-3)/(k(k-1)). Their mean is mu_i=2-3sigma_i/k and variance
C0[(1-x)-sigma_i^2/k]. Thus their contribution, including the move weight,
is at most

    binomial(k,3)/20 * {1/mu_i
       +C0[(1-x)-sigma_i^2/k]/[a_i mu_i^2]},
    a_i=2-sqrt(3(1-x)).

These are rigorous moment bounds; no assumption that coordinates are
equal is made.

For fixed norm data each group's maximum over its allowed sum is found
analytically. For the finite group write y=sigma_f/m, b=sqrt(2),
a=a_f and q=x/m. Its normalized bound is

    F(y)=1/(b-y)+(q-y^2)/[a(b-y)^2].

The derivative has sign ab+2q-y(a+2b), so its maximum on an interval is
at the critical point (ab+2q)/(a+2b), clipped to that interval.
For the ideal group set y=sigma_i/k, z=1-x, a=a_i. Its bound has
derivative sign 6a+6C0z-y(9a+4kC0), so its critical point is
(6a+6C0z)/(9a+4kC0), likewise clipped. This includes k=3, where C0=0.

## Exact interval certificate, rather than point sampling

The norm-split interval [0,1] is divided into sixteen closed rational
intervals [l,h]. At k=0 only x=1 is possible, and at k=12 only x=0 is
possible, giving 178 total boxes over all k.

Use b0=1414213562373095/10^15<sqrt(2), and enclose every other square root
between adjacent multiples of 10^-7 by integer square roots. In each box,
for the finite group replace x by h, set a=b0-upper_sqrt(h), and maximize
the displayed rational F with

    lower_sqrt(l)/m <= y <= upper_sqrt(h/m).

This is an upper bound for every actual point of the box: it inflates the
variance numerator, decreases the positive denominators, and enlarges
the allowed sum interval. For the ideal group use z in [1-h,1-l], inflate
it to 1-l in the variance and a, and maximize over

    lower_sqrt(1-h)/k <= y <=upper_sqrt((1-l)/k).

All parameters and clipped critical points are rational. The derivative
sign formulas prove the maxima over the whole intervals. It is harmless
if an enlarged interval includes values with a negative inflated variance;
the analytic maximum is over the larger rational function domain, which
contains every actual allowed value where the variance is nonnegative.

Every resulting row gap exceeds 1/8; the smallest is greater than
0.125535. The checker verifies the exact rational inequalities and records
the smallest gap and worst box for each k. This is a cover of all norm
and group-sum possibilities, not a numerical sample or optimizer claim.

## Exceptional rows and conversion to an ordinary count

At any vertex, even one with large carriers, a selected finite-edge ratio
contributes at most sqrt(2) above its unit weight. A selected P3 ratio
contributes at most (1+sqrt(3))/20 above its weight. Since there are at most
twelve finite edges and 220 corner triangles,

    R(q)-d(q) <=12sqrt(2)+11(1+sqrt(3)) <48.

The last comparison is checked with rational upper bounds sqrt(2)<10/7
and sqrt(3)<7/4. Let b count nonsmall finite vertices. Combining this bound
with the small-vertex deficit and global reciprocal-pair inequality gives

    0 <=sum[R-d] <=-(f-b)/8+48b,

so b>=f/385.

Finally a large three-face with nu>=11 has at most 14(nu-10) finite
vertices. At nu=11, facets+2cusps>=12 gives facets<=10 and hence
p=4facets-2nu-4<=14. At nu>=12, p<=2nu-4<=14(nu-10).
Every nonsmall vertex is incident to at least one such face, so

    b <=sum_large p <=14E.

Therefore E>=f/(14*385)=f/5390, as asserted. No normalization by energy
remains in the final count.

## Limits and reproduction

This proves a positive excess fraction in dimension twelve, but its
coefficient is far from enough to close the known dimension-twelve
aggregate gap. The same groupwise second-moment bound does not currently
certify the analogous row estimate in dimension eleven near k=6,m=5 and
finite squared norm one; no eleven-dimensional fraction is claimed.

Run `python3 scripts/check_right_angled_mixed_carrier_excess_fraction_research_20261018.py`
for optimizer-free exact replay. Generation with `--write` acquires the
shared resource lock and applies the resource gates; the computation is
only 178 tiny rational boxes. Frozen earlier packets are unchanged.
