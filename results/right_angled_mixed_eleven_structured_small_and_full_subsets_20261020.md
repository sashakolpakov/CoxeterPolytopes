# Structured opposition estimate: seven, eight and eleven directions

20 October 2026. **Written analytic subproblem, with a reproducible
rational interval certificate.** This is a partial estimate towards the
mixed dimension-ten question, not a dimension theorem. The nine- and
ten-direction cases remain separate obligations.

Let p_1,...,p_11 be real, with sum p_i=3 and sum p_i^2=1. For a
three-element subset T write r_T=2-sum_(i in T) p_i. For a set I of
k directions the estimate considered here is

    sum_(T subset I, |T|=3) (r_T^-2-37/100)_+
        <=756/25+(27/100) C(k,3).                 (1)

All sums use distinct triples. No assumption about which triples are
actual graph triangles enters this lemma. The proof below establishes
(1) for k=7,8,11 if its finite rational certificate replays successfully.
It makes no claim that checking these cases proves (1) at other k.

## Common ratio and moment identities

Center the full frame at3/11. Its squared deviation is2/11, so Cauchy
gives the full three-subset range

    (13-4sqrt3)/11 <= r_T <= (13+4sqrt3)/11.

In particular every r_T>11/20 and r_T<20/11; the lower comparison
follows by squaring 139/80>sqrt3. This bound holds before restricting I.

For k<11 put S=sum_I p_i and V=sum_I(p_i-S/k)^2. Cauchy on the
complement gives

    0<=V<=v_k(S)=1-S^2/k-(3-S)^2/(11-k).

Its center is S=3k/11, its maximum is2/11, and its allowed interval
has endpoints [3k +/- sqrt(2k(11-k))]/11. For k=7 it is contained
in [27/22,57/22]; for k=8 it is contained in [17/11,31/11].

Let y_i=p_i-S/k, mu=2-3S/k, d_T=-sum_T y_i, T3=sum y_i^3 and
T4=sum y_i^4. Averaging uniformly over the C(k,3) triples gives

    E d=0,
    E d^2=alpha2 V,
    E d^3=-alpha3 T3,
    E d^4=alpha4 T4+b4 V^2,

where

    alpha2=3(k-3)/[k(k-1)],
    alpha3=3(k-3)(k-6)/[k(k-1)(k-2)],
    alpha4=3(k^2-17k+54)/[k(k-1)(k-2)],
    b4=18(k-4)/[k(k-1)(k-2)].

These identities follow by expanding indicators of distinct indices;
their one-, two- and three-index inclusion probabilities are respectively
3/k, 6/[k(k-1)] and 6/[k(k-1)(k-2)]. At the sizes used here,
alpha3>0 and alpha4<0. Cauchy gives T4>=V^2/k, while
|T3|<=V^(3/2) follows from |y_i|<=sqrt(V).

## A rigorous reciprocal and positive-part envelope

The exact identity

    (mu+d)^-2=mu^-2-2d/mu^3+3d^2/mu^4-4d^3/mu^5
                  +d^4*(4(mu+d)+mu)/[mu^5*(mu+d)^2]

has a remainder factor decreasing in the positive ratio mu+d.
Thus, with a=11/20 and K(mu)=(4a+mu)/(mu^5*a^2), the average
reciprocal square is at most

    J=mu^-2+3alpha2 V/mu^4+4alpha3 V*sqrt(V)/mu^5
                    +K(mu)*(b4+alpha4/k)*V^2.     (2)

Every term on this upper bound is nonnegative, decreases with mu and
increases with V. In particular replacing mu by a positive lower bound
and V by an upper bound preserves the inequality. The coefficient
b4+alpha4/k is positive at k=7,8.

Put t=37/100 and b=t^-1/2-mu. The positive-part average equals
E(r^-2)-t+E(t-r^-2)_+. Its last term is at most27/400, by r<20/11.
If an interval has upper mean mu_hi<18/11, we also have

    E(t-r^-2)_+
       <=(113/500)*alpha2 V/[2*(18/11-mu_hi)].      (3)

Indeed t^(3/2)<113/500 and t^-1/2>18/11. The function
(t-r^-2)_+ is at most 2t^(3/2)*(d-b)_+, and
(d-b)_+<=d^2/(4b), since (d-2b)^2>=0 when d>b.
When an interval additionally satisfies mu_hi+sqrt(3(k-3)V/k)<=18/11,
the correction is exactly zero. Taking the minimum of these justified
bounds is allowed.

## Complete rational interval proof for k=7,8

The checker splits the enclosing S intervals above into consecutive
closed intervals of length1/1100. There are1500 and1400 intervals.
On each it takes the exact maximum of the concave quadratic v_k(S),
at the nearest point to3k/11. An interval whose maximum is negative
contains no admissible frame and is skipped. It uses the exact lower
and upper means at the interval endpoints.

Square roots are rounded UP by exact integer-square comparisons at
denominator10^6. Substitution in (2), and the best correction from
(3) or the preceding range bounds, gives a rational upper average.
The checker proves this is strictly less than

    37/100+[756/25+(27/100)C(k,3)]/C(k,3)

on every nonempty interval. The closed intervals cover all admissible
S, so these finite rational inequalities imply (1). This is an interval
proof, not sampled evidence. The certificate stores the count and least
gap; replay recomputes every interval without trusting a stored verdict.

## The full eleven-direction case

Here S=3, mu=13/11 and V=2/11 exactly. Write T3=sum y_i^3.
The sharper Cauchy identity

    T4>=V^2/11+T3^2/V

follows by projecting the vector of squares onto the orthogonal vectors
1 and y. The cubic bound is |T3|<=18/(121sqrt5)<1/15. For completeness,
with M=sqrt(10V/11), z=-M/10, sum(M-y_i)(y_i-z)^2>=0 gives its
upper half; replace y by -y for the lower half. Expanding uses only
sum y=0 and sum y^2=V.

For a=11/20 and mu=13/11 the rational K(mu)<5. Hence the average
reciprocal square is at most

    mu^-2+3*(24/605)/mu^4+5*(15/121)*V^2
                 +16*T3/(33mu^5)-T3^2.

The function C*T3-T3^2, C=16/(33mu^5), is increasing on
[-1/15,1/15], as C>2/15. Replace its last two terms by C/15-1/225.
The positive-part correction is at most

    (113/500)*(24/605)/(2*(5/11)),

because t^-1/2-mu>5/11. Every remaining number is rational. The
checker verifies their sum minus t is below
[756/25+(27/100)*165]/165, establishing (1) for k=11.

The checker also reconstructs the balanced ordinary incidence row
3722B>159195f+2794E from its exact identities. That is a separate
reviewed necessary-count theorem. The proposed universal opposition
upper bound still requires (1) at every remaining subset size and
its written transport argument; it is not assumed in these checks.
