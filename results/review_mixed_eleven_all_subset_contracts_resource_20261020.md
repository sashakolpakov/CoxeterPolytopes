# Written review of every eleven-frame subset contract

20 October 2026. **PASS: no flaw found in the analytic reductions or
exact-certificate methods.** I independently read the complete sources
and workers listed below and rederived their moments, clipping estimates,
domain reductions and polynomial envelopes. I did not run a scientific
program, import a proof source, replay a certificate, optimize, enumerate
frames, build TeX or run Lean. Source hashing and reading the small existing
execution ledgers were the only execution in this review.

One independence distinction matters: I authored the new global k9
proof. Its independent written review belongs to the coordinator. I
separately audited the earlier k9 interval proof and worker, which I did
not author; the other subset proofs and all k10 source derivations below
were also independently reviewed here. The execution statements below
are attributed to the existing generation/replay ledgers, not to a new
execution by this reviewer.

## Exact statement and scope

For eleven arbitrary real coordinates satisfying sum p_i=3 and
sum p_i^2=1, an index set I of size k, and
r_T=2-sum_(i in T)p_i, the reviewed contract is

    sum_(T subset I, |T|=3) (r_T^-2-37/100)_+
       <=756/25+(27/100) C(k,3),       0<=k<=11.              (1)

The triples are formal distinct subsets of I, not necessarily actual
two-face or graph objects. No proof assumes nonnegative individual p_i,
a uniform selected vector, a maximal-height actual subset, or global
pair-balancing monotonicity. This review concerns(1); it does not itself
establish the actual-frame construction, the unweighted actual-subset
summation, the balanced geometric lower bound, or a dimension theorem.

## Source pins

|Source|SHA256|
|---|---|
|`results/right_angled_mixed_k10_subset_moment_reduction_20261020.md`|`8c2051b3be0f692e5097eedd81cf494f15313037d36b0f6625ca78cd01aabd34`|
|`results/right_angled_mixed_k10_rational_envelopes_20261020.md`|`7fdcad68b9b4ae113c482e9e17266a6f5e9aea0d21409eb7db332ec915627e45`|
|`scripts/check_mixed_k10_rational_envelopes_20261020.py`|`d0520c374cb3e8689dedec78284dc50ea8d7d9c182b47d71ab2e81430bcf99fb`|
|`results/right_angled_mixed_eleven_structured_small_and_full_subsets_20261020.md`|`04d3b8fa433365210600234709894d777803e7db8591cbb32876a0796e87a1ce`|
|`scripts/check_mixed_eleven_balanced_and_subset_rows_20261020.py`|`e403456a9f880af672407173e51b18e4e3b03b3cb474cdb1030ce54ba381f830`|
|`results/right_angled_mixed_eleven_nine_direction_subset_20261020.md`|`651b2b90fe0e1b4df52cd03e99a8c53ce0f92af38f89c424f1e539df5c2ebe59`|
|`scripts/check_mixed_eleven_nine_direction_subset_20261020.py`|`f73ca6c70a354244e7c9768e46ad9295fc10f35fc36925aee36a3eea3f9ef68a`|
|`results/right_angled_mixed_eleven_nine_direction_global_proof_20261020.md`|`bf288985ade5acf93e13faafb25544b0642b38f13e422ad38a3f1223e9706d15`|
|`scripts/check_mixed_eleven_nine_direction_global_20261020.py`|`357070546d69aaa634c01292fb15c174f6d85552e87cb42be6adbf4d2cdb0b62`|

## Common identities checked independently

Centering all eleven p_i at3/11 gives squared norm2/11. The centered
indicator of a triple has squared norm24/11. Thus

    (13-4sqrt3)/11<=r_T<=(13+4sqrt3)/11,
    11/20<r_T<20/11.

The two rational comparisons follow from139^2>3*80^2 and48<49.
This proves positivity before any reciprocal or clipping manipulation.

Within k selected coordinates, center at S/k and put V=sum y_i^2,
T3=sum y_i^3,T4=sum y_i^4,mu=2-3S/k,d=-sum_T y_i. Expanding
with one-, two- and three-index inclusion probabilities
3/k,6/[k(k-1)],6/[k(k-1)(k-2)] reproduces all four moment coefficients
printed in the sources. In particular the quartic coefficients are

    alpha4=3(k^2-17k+54)/[k(k-1)(k-2)],
    b4=18(k-4)/[k(k-1)(k-2)].

The sign alpha4<0 for7<=k<=11 is essential and is used correctly.
The sharp cubic estimate follows by summing
(M-y_i)(y_i+M/(k-1))^2 with M=sqrt((k-1)V/k). The joint bound
T4>=V^2/k+T3^2/V follows by Cauchy after subtracting V/k from
the squared entries. Neither is a critical-point classification.

I checked the exact reciprocal Taylor identity and its remainder
coefficient(4r+mu)/(mu^5 r^2). For positive r it decreases with r;
for a fixed positive lower r bound it also decreases with mu. This
justifies the conservative endpoint replacements used by the workers.

For c=37/100 and rho=c^-1/2, the correction to the full reciprocal
mean is E(c-r^-2)_+. Its global bound27/400 follows from r<20/11.
The variance bound uses the concave tangent at rho and
(d-b)_+<=d^2/(4b), b=rho-mu>0. I checked the factors of2:
the resulting bound is c^(3/2) E d^2/[2b]. The source consistently
uses an upper mean when lowering b, and sqrt bounds when declaring
the correction zero. Thus it never substitutes an uncorrected full
reciprocal estimate for a positive-part estimate.

## Sizes zero through six

The k10 reduction source includes these complete cases, not just its
ten-direction work. The uniform lower r>6/11 gives each summand less
than673/225. For k<=5 at most ten triples suffice, since
6730/225<756/25.

At k6, pair each triple with its complementary triple. Their ratios
sum to C=4-sum_I p>18/11, while each ratio exceeds6/11. The positive-
part function is convex and decreasing. Its maximum on a fixed-sum
pair lies at the interval endpoints, so the pair sum is bounded by
f(6/11)+f(12/11)=12461/3600. Ten pairs give12461/360<891/25,
the required k6 right side. The six-coordinate Cauchy bound and strict
rational comparisons cover the entire real domain.

## Seven and eight directions

Cauchy on the complement gives
V<=1-S^2/k-(3-S)^2/(11-k), with exact center3k/11 and maximum2/11.
The rational enclosing S domains used by the worker contain both true
endpoints. The 1/1100 grids therefore cover the full admissible domain,
including boundary points. Clamping the quadratic vertex to each closed
interval computes its exact maximum. An interval is skipped only if
even that maximum is negative.

The mean decreases with S, and the worker takes its lower value at
the right endpoint. The quartic coefficient b4+alpha4/k is positive,
so the positive moment upper bound increases with V and decreases with
mu. Rounded-up square roots are obtained by integer-square comparisons,
not floating-point evaluation. Every correction uses a justified global,
variance or zero bound. Taking their minimum is valid.

The worker rechecks each rational interval and compares the correct
positive-part target, including subtraction of37/100. Its certificate
reports1498 checked plus2 empty k7 intervals, and1386 checked plus14
empty k8 intervals. These are full interval covers, not sample points.
The stored gaps are positive, and separate generation and replay ledgers
record PASS for these exact source bytes.

## All eleven directions

Here mu=13/11,V=2/11. I independently obtain

    E d^2=24/605,
    E d^4<=15V^2/121-T3^2/5,
    |T3|<=18/(121sqrt5)<1/15.

The global lower ratio gives K<5, so the quartic term becomes
at most5*(15/121)V^2-T3^2. The remaining expression
C*T3-T3^2, C=16/(33mu^5), increases on[-1/15,1/15] because
C>2/15. Thus the source's substitution C/15-1/225 is valid.
Its clipping gap uses rho-mu>5/11. The final rational target is
2493/5500. The worker's strict positive gap checks the full analytical
bound, not an assumption that T3 and T4 are independently extremal.

## Nine directions

For the earlier interval proof I independently checked the exact
delta domain[0,12/11], V bound2delta/3-11delta^2/18, all three
moments and the sharp joint quartic inequality. The small-delta
neighborhood is analytic, covers0<delta<=1/100 and has zero clipping
correction; its mean reciprocal bound1-5delta/36 is valid. At delta0
the uniform nine-entry equality frame is forced by Cauchy.

On the remainder, a0=max(11/20,mu0-sqrt_up(2v0)) is a valid lower
ratio bound. Increasing V to v0 raises the joint fourth-moment bound
while retaining its SAME remainder constant. Negative T3 can be
discarded before replacing the positive cubic coefficient. The worker
then maximizes the concave quadratic over a containing interval by its
endpoint or vertex. Both the root rounding and clipping correction are
in the conservative direction. Its10810 closed intervals cover through
the exactly truncated endpoint12/11.

The separate authored global proof avoids that interval computation.
Its new content is the global bound(c-r^-2)_+<=3(r-mu)^4 and a
single cleared polynomial whose residual is at least635/968. It gives
the stronger mean bound63/100-delta/18. I do not label my own derivation
an independent review; the coordinator supplied the independent check.
The earlier interval proof remains separate evidence for the exact same
k9 contract and equality frame.

## Ten directions: complete envelope reconstruction

I independently derived, for the omitted coordinate v,

    mu=(11+3v)/10, V=(1+6v-11v^2)/10,
    (3-2sqrt5)/11<=v<=(3+2sqrt5)/11.

The centered ten-entry moments are7V/30,-7T3/60 and
3V^2/20-T4/15. With T4>=V^2/10+T3^2/V, the full quartic
bracket is43V^2/300-T3^2/(15V). This remains nonnegative throughout
the sharp cubic range: its lower coefficient is259/2700>0. Increasing
K to seven is consequently justified for the entire bracket; the
negative term has not been changed separately.

For v>=0, mu>=11/10 and the full-frame ratio bound imply
K<=12000000/1771561<7. For v<=0, V<=1/10 and V'>=3/5 imply
a=mu-sqrt(21V/10) decreases up to v0. Hence a>=a(0)>16/25;
mu>21/20 gives K<7220000000/1045529856<7. The continuity case
V0 and all denominator signs are valid.

In the negative branch the quadratic vertex is beyond the sharp
cubic endpoint, because2gamma*mu^5*sqrt(V)<1. Using gamma<17/20
only in the positive cubic term, and the tangent
sqrt(V)<=1/10+5V/2, gives exactly the printed Jneg coefficients
119/3000,119/120 and1813/2700. The negative gamma-square contribution
has correctly been retained at its actual value. The clipping correction
is zero since mu+sqrt(21V/10)<39/25<41/25<rho.

In the nonnegative branch completing the square over all real T3 gives
Jpos with coefficients301/300 and7/60. For0<=v<=1/20 the exact
mean and variance bounds make the correction zero. Elsewhere the bound
791V/[30000*(41/25-mu)] follows from the common variance correction.
That denominator stays positive on the entire enlarged rational interval.

I reconstructed all three cleared polynomials from those formulas:
Pneg, Pzero and Pcor have degrees9,14,15. The source's rational domains
[-7/50,0],[0,1/20],[1/20,7/10] cover the actual domain. Values with
V<0 on an enlarged endpoint are used only as polynomial test points;
no moment interpretation is asserted there. The target is correctly
223/250 before subtracting37/100, giving261/500 afterward.

The source also correctly warns that uniform selected coordinates are
not a global k10 maximum; its displayed two-level frame satisfies the
full-frame sum and norm. The proof does not rely on that symmetry
heuristic.

## Exact finite certificate methodology and recorded executions

The k10 worker first forms exact Fraction polynomials. After substituting
v=left+(right-left)u it converts power coefficients to Bernstein
coefficients by

    b_i=sum_(j<=i) a_j*C(i,j)/C(n,j).

This is the correct exact conversion. Strictly positive coefficients
give strict positivity of the polynomial on each entire closed interval.
There are140 degree9,50 degree14 and650 degree15 intervals, totaling
840 intervals and12550 coefficients. No solver, root finder or sampled
positivity inference is involved. Every coefficient is reconstructed on
replay; the certificate stores polynomial coefficients, a coefficient
digest, exact minimum and a positive rational floor for each interval.

All workers stream their source hashes, use exact square-root rounding
or rational polynomial arithmetic, and exclusively create their initial
certificate. Read-only replay compares the whole reconstructed result.
The inspected ledgers record separate generation/replay PASS. This
review did not execute any of them.

|Certificate|SHA256|
|---|---|
|`results/mixed_eleven_balanced_and_subset_rows_20261020.json`|`eb0b79f20b856e13ec4ec78ba57f147328eb07ad08b846d8c4d87baf2134e257`|
|`results/mixed_eleven_nine_direction_subset_20261020.json`|`18073931501990399f28b58e75dca098af833a18ae5f10e7d140db8e8c8f0c77`|
|`results/right_angled_mixed_k10_rational_envelopes_20261020.json`|`a04a84b1d997b1479d5478408e277684228afd1b1d65e67114bf747c7f4c6163`|

The relevant separate replay ledgers are
`mixed_eleven_balanced_subset_replay_20261020.json` (SHA
`ce89fa94118d7e572a0e0aade9a931fcc869b0ef7e3ae61d15e90dbd25c2eb25`),
`mixed_eleven_nine_direction_replay_20261020.json` (SHA
`829c2a8092ac09a9746d49db883c1f5b6dcddd4cd32c1ae935f8015a64ef06f0`),
and `mixed_eleven_k10_envelope_replay_20261020.json` (SHA
`f2a2f71c2b328bc274780668634cb379e45be60cb1e47406fb66144dfa6f7066`).

The balanced-row arithmetic sharing the first worker does not itself
prove its geometric premises. Those and the final actual opposition
transport require their separate written review. Subject to that scope,
the analytic subset obligations(1) are complete for every k0,...,11.
No dimension promotion or new Lean claim is made in this review.
