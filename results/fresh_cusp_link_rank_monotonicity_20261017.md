# Product-link curvature increases with each sufficiently large factor

For a product of simplex cusp links with factor ranks r1,...,rs, put d=sum(r_i), V=product(r_i+1), and let kappa be the exact link curvature in the [shape formula](fresh_cusp_link_shape_curvature_20261017.md). If every r_i>=17, increasing any one rank by1 increases kappa/V by at least7/76. Consequently, when n=d+1>=30,

    kappa/V >= 3s,

with equality only for the single factor of rank29. This is an all-rank analytic inequality. It does not use the existing ambient bound163, a finite partition cutoff, or an optimizer.

If the proposed fine face payments of cost3n at finite vertices and3sV at ideal vertices were complete under this cusp-rank hypothesis, this inequality and the positive additive shape slack would imply n<=29. The all-finite-quadrilateral transfer is still open under that weaker hypothesis, so no such ambient theorem is asserted here.

## A pointwise comparison of centered discrete uniforms

Let X_r be uniform on {-r/2,1-r/2,...,r/2} and define f_r(y)=E|y+X_r| for real y. Both f_r and f_(r+1) are even and piecewise linear. Their combined nonnegative breakpoints are integers and half-integers, and their difference is zero outside the larger support interval. Thus it is enough to compare those breakpoints exactly.

For r=2k, put D=(2k+1)(2k+2). At integer j and half-integer j+1/2 in the respective support ranges,

    f_(2k+1)(j)-f_(2k)(j) = ((k+1)^2-j^2)/D,
    f_(2k+1)(j+1/2)-f_(2k)(j+1/2) = (k(k+1)-j(j+1))/D.

The ranges are 0<=j<=k+1 for the first expression and0<=j<=k for the second. The maximum is attained at y=0 and equals (r+2)/(4(r+1)).

For r=2k+1, put D=(2k+2)(2k+3). The analogous expressions are

    f_(2k+2)(j)-f_(2k+1)(j) = ((k+1)^2-j^2)/D,
    f_(2k+2)(j+1/2)-f_(2k+1)(j+1/2) = ((k+1)(k+2)-j(j+1))/D.

Both ranges are 0<=j<=k+1. The maximum is attained at y=1/2 and equals (r+3)/(4(r+2)). The displayed numerators also verify nonnegativity on their ranges. Linear interpolation and the zero tails complete the pointwise proof for every real y.

For every r>=17, either maximum is at most5/19. The estimate is uniform in y, so it remains valid after averaging over an arbitrary independent random variable Y. There is no assumption about Y's lattice parity or number of summands.

## Curvature monotonicity

The exact normalized shape identity is

    kappa/V = 3d/4 - (5/2) E|X_(r1)+...+X_(rs)|.

Increasing one rank changes d by1. The preceding pointwise estimate, averaged over the sum of the other centered uniforms, therefore gives

    Delta(kappa/V) >= 3/4 - (5/2)*(5/19) = 7/76 > 0.

In particular, for a fixed number s of factors with all ranks at least17, the normalized curvature is minimized when every rank is17.

For s>=2, the variance formula and Cauchy--Schwarz at that minimal profile give

    kappa/V >= (51/4)s - (5/2)*sqrt(323s/12) > 3s.

Both sides in the final comparison are positive before squaring. Its squared form is 4563s>8075, and at s=2 its integer margin is1051; the margin increases with s. This proves the result for every number of factors, without enumerating partitions.

For s=1, the simplex-link identity gives kappa/V=beta_n/n, where beta_n=n(n-6)/8 for even n and (n-1)(n-5)/8 for odd n. For n>=30 this is at least3, with equality only at n=30. The exact numerator identities are

    8(beta_(30+2t)-3(30+2t)) = 4t^2+60t,
    8(beta_(31+2t)-3(31+2t)) = 4t^2+64t+36.

## Conditional geometric consequence and the remaining obligation

The shape theorem gives

    D5 >= beta_n*v_f + sum_v kappa(L_v) + 5(n-1).

For n>=30 and all factor ranks at least17, the right side is strictly larger than 3n*v_f+3*sum_v s_v V_v. Thus a complete family of face payments with those vertex costs would be impossible. This is a conditional implication, not a replacement for checking every face type. The new large-factor proof treats ideal-corner quadrilaterals; the preliminary private-diagram steps for four-finite-corner quadrilaterals still require work.

Equivalently, after triangle and ideal-corner payments are established, let Q_bad count the four-finite-corner quadrilaterals receiving no fine payment. The same accounting would force

    Q_bad >= (beta_n-3n)*v_f + sum_v(kappa(L_v)-3s_v V_v) + 5(n-1).

Under the present rank and dimension assumptions this is at least5(n-1). This identifies the remaining geometric incidence obligation instead of silently discarding it.

## Verification

The [checker](../scripts/fresh_cusp_link_rank_monotonicity_20261017.py) and [data](fresh_cusp_link_rank_monotonicity_20261017.json) verify all four rational parity identities, the exact squared comparison and the single-factor parity polynomials. Sixty-one literal discrete-uniform calculations corroborate the formulas; they do not supply their unbounded scope. Generation and a fresh saved replay pass. Syntax, all four direct source hashes and local links pass. A narrow root mathematical/source review passes the four pointwise parity formulas, uniform increment bound, variance comparison, one-factor thresholds and conditional incidence accounting. That review does not claim a duplicate arithmetic replay. The source, data and report are final and frozen.
