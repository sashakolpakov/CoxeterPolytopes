# Root review of the last ten-direction analytic case

20 October 2026. **PASS: the ten-direction structured subset contract
is proved by the written reductions and exact finite certificates.**
The full geometric/counting composition is a separate review obligation.
This record does not claim geometric Lean verification or a polytope census.

I read both complete written sources and the entire198-line worker:

* `right_angled_mixed_k10_subset_moment_reduction_20261020.md`,
  SHA256 `8c2051b3be0f692e5097eedd81cf494f15313037d36b0f6625ca78cd01aabd34`;
* `right_angled_mixed_k10_rational_envelopes_20261020.md`,
  SHA256 `7fdcad68b9b4ae113c482e9e17266a6f5e9aea0d21409eb7db332ec915627e45`;
* `scripts/check_mixed_k10_rational_envelopes_20261020.py`,
  SHA256 `d0520c374cb3e8689dedec78284dc50ea8d7d9c182b47d71ab2e81430bcf99fb`.

## Independent derivation of the reductions

For the omitted coordinate v, the mean and centered norm are
mu=(11+3v)/10 and V=(1+6v-11v^2)/10. Their feasible domain is
exactly [(3-2sqrt5)/11,(3+2sqrt5)/11]. The centered triple moments
are 7V/30, -7T/60, and 3V^2/20-U/15. I rederived these from
the one-, two- and three-index inclusion probabilities.

The elementary zero-sum cubic bound and U>=V^2/10+T^2/V apply
jointly. The exact reciprocal remainder therefore produces

    mu^-2+7V/(10mu^4)+7T/(15mu^5)
                   +K(43V^2/300-T^2/(15V)).

For every feasible T the bracket is nonnegative:
43/300-(64/90)/15=259/2700>0. This justifies replacing the WHOLE
coefficient K by7. It is not legitimate to replace its positive and
negative occurrences independently; the source does not do that.

The K<7 proof is valid on each half of the true domain. For v<=0,
V<=1/10 and V'>=3/5 imply a'<0, so a>16/25, while mu>21/20.
For v>=0 use a>11/20 and mu>=11/10. Both rational upper constants
in the source equal the corresponding decreasing remainder factors
and are strictly below7. The full-frame lower a bound applies to
the extremal centered triple as well, since the centered deviation
vectors with the specified v are feasible full frames.

For v<=0 the concave maximum after setting K=7 is at its upper
cubic endpoint. Indeed its vertex is V/(2mu^5), and the stated
strict endpoint comparison is valid. Substituting gamma<17/20 and
sqrt(V)<=1/10+(5/2)V gives exactly Jneg, including1813/2700.
Here mu+w<39/25<10/sqrt37, so the clipping correction is zero.

For v>=0, completing the square gives7V/(60mu^10), which yields
Jpos. On [0,1/20], the upper bounds mu<=223/200 and
w<13/25 keep all ratios below the clipping threshold. On the
remaining true domain the tangent/second-moment correction is at most
791V/[30000(41/25-mu)], with a positive denominator. I independently
checked the direction of both rational radical replacements.

Consequently clearing the positive denominators gives exactly the
three source polynomials Pneg, Pzero and Pcor. Checking their larger
rational intervals is valid even where V would be negative: those
extra points are only part of the polynomial certificate, never
asserted to be geometric or moment data.

## Exact coefficient proof and execution

The worker constructs the three polynomials by Fraction arithmetic;
it does not copy saved polynomial coefficients into its calculation.
It applies the affine substitution on each closed grid interval and
the power-to-Bernstein coefficient identity. The Bernstein basis is
nonnegative and sums to one, so strict positivity of every coefficient
proves full-interval positivity, rather than positivity at grid points.

Generation and a separate read-only replay both passed. The complete
840 intervals have12550 strictly positive coefficients. The certificate
also stores the exact global minima, all coefficient digests and one
positive lower floor per interval. The bounded ledgers are
`mixed_eleven_k10_envelope_generation_20261020.json` and
`mixed_eleven_k10_envelope_replay_20261020.json`. Each execution took
less than half a second and used less than56MB process-tree RSS.
The shared lock,59-second/1GiB limits and normal disk/RAM gates were
in force. No build, Lean process or optimizer ran.

Thus the corrected reciprocal average is strictly below223/250;
subtracting37/100 gives the required261/500 positive-part average.
Multiplication by120 proves the full k10 contract. No symmetry
assumption or substitution of a uniform local maximum is used.

I also checked the complete elementary k<=6 proof in Section7 of the
moment source: the full-frame lower ratio and complementary-triple
convexity estimates suffice. Together with the other separately
reviewed subset sizes, these statements supply the analytic input
to the proposed dimension-ten assembly. That final passage still
requires the actual opposition graph, ordinary counts and inherited
three-face geometry; coefficient positivity alone proves none of them.
