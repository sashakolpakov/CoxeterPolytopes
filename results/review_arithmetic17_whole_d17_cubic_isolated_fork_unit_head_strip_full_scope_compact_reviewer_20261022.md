# Full outside-author review: whole-D17 cubic isolated fork

22 October2026. **FULL mathematical SOURCE PASS; no correction requested.**
This is a manual review of the complete new proof, at the explicit accepted
geometric/arithmetic inputs below. No worker was imported, no finite domain
was executed, and no enumeration, symbolic process, test, Lean, TeX, compiler
or Git operation was used. The result is not an execution or external
refereeing claim.

## Frozen source and reviewed inputs

The reviewed source is
`results/arithmetic17_whole_d17_cubic_isolated_fork_unit_and_head_strip_exclusion_candidate_20261022.md`,
SHA256 `851c6fa4f98a2a901e857b1fdb2c230d89cd2d58009d35c5961c488b720ffc97`
(11,620 bytes;268 lines). I read it in full, with its actual pinned input
interfaces, rather than treating earlier PASS labels as a proof of the new
unit, ramification or integral-basis arguments.

All four literal input hashes match current bytes:

| Input under results/ | SHA256 | Scope used |
|---|---|---|
| arithmetic17_whole_d17_dotted_head_degree_three_and_nonisolated_quadratic_candidate_20261022.md | 4baa131a7609f1af7254e154fae99834f0ab38f9e7a31bc2a49a186b0916ceba | Complete actual degree2/3 reduction; in degree3 both actual forks isolated dotted; actual integral head/fork parameters and strips; shared-time cycle identity |
| review_arithmetic17_whole_d17_degree_three_nonisolated_quadratic_root_acceptance_20261022.md | c028c30b686ef85ef68dc6ddf834fbce77c38d073b479da93d0a24b39deb98e5 | Written acceptance and precise open-interface boundary of the preceding source |
| arithmetic17_whole_d17_fresh_long_head_and_all_finite_head_exclusion_candidate_20261022.md | 9e2cf7d75f12c73b69e0bc9c39044acb56c37be36f94fe707efb808ee5ae9eb6 | Actual dotted head/field; Section8's full cyclic-cubic integer basis, rational root bounds and interpolation identities, used as field algebra |
| arithmetic_uniform_classical_weak_helper_obstruction_20261021.md | 535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905 | Actual admissible K-root lines and invariant cyclic integrality; no individual Cartan-entry or disconnected-length integrality assumption |

I previously independently reviewed the two D17 head/degree sources and
have contributed to shared actual-root/cycle foundations in this research.
I did not author this new cubic proof. Those inherited foundations are
disclosed inputs, not independently regenerated literature or geometry.
The new norm-two, dyadic split, ramified-polynomial completeness and both
full integer-basis head-strip exclusions are independently checked here.

## Actual configuration and new norm-one conclusion

Source lines39–54 correctly retain ONE original whole-D17 farthest vertex
and ONE actual field K. In degree3 the accepted interfaces give

    U integral, K=Q(U), 4<U<11, other0<U_sigma<2;
    W integral, K=Q(W), 2<W<676/17,
    other0<W_sigma<4/17, N(W) in{1,2}.

The head and chosen isolated fork have different old-root supports and are
actual root lines in that same space. Their parameters are not selected
from unrelated feasible Grams, and the two fork parameters are not equated.

For P=X³-A X²+B X-n, positivity gives A,B positive integers and A<=40.
The sign P(4/17)<0 is correct: that point lies above both small roots and
below the physical root. Clearing the denominator gives exactly

    64-272A+1156B-4913n<0.

For n2, wab=2 implies w>289/8 and
B=w(a+b)+ab>=2sqrt(2w)+2/w>17. The upper inequality gives
B<20642/1156<18. No integer survives. Consequently W is a unit in
EVERY remaining cubic field, and 289/16<w<676/17 with A19,...,40.
This argument precedes the ramification split and discards no ramified2
case (source lines56–77).

## Both actual cycles and the unramified branch

I rederived the raw normalization and source lines79–103. The old D roots
have square2; the actual head and fork have squares2U and4W, sole old
products-U and-2W, and spatial mutual product UW. The full old-chain
closed Cartan cycle yields integral g, while the two-cycle yields
L=g²/(2UW) integral. Expanding their SAME actual time-line identity gives

    L-g=8UW-17W-2U+4,
    g²=2UW(g+d0), d0=(8U-17)W-2U+4.

At every prime p above2, W being a unit makes d0 a unit. Equation(5)
forces g into p and g+d0 a unit, with g nonzero; hence
2v_p(g)=v_p(2)+v_p(U). If2 is unramified, v_p(2)=1 at EVERY such p,
so v_p(U) is positive odd and U/2 is integral. Other primes present no
denominator. No unramifiedness is presumed in the alternative branch.

For Y=U/2 the exact product norm bound is99/64<2. Its two nonzero
integer norm factors therefore give N(Y)=1 and N(1-Y)=-1. Trace3,...,7
gives precisely X³-TX²+(T-1)X-1. I checked the discriminants
-23,-31,-23,49,257; at T7 the value at6 is-1 and the two other roots are
below1, so its physical root exceeds6. Only T6 remains. Its reciprocal
field has f=X³-5X²+6X-1; square discriminant49 makes it cyclic, allowing
an auxiliary generator with its physical value the largest root.
The 7-Eisenstein shift and the local valuation argument establish the
FULL ring O_K=Z[u], rather than a guessed suborder.

The inherited rational brackets and derivative estimates used at
source lines123–145 are valid field-algebra facts, independent of a
finite-head assumption. Applying the CURRENT positive W strip gives
3<C0<8 and |B0+C0(beta+chi)|<.174. The complete C0=4,5,6,7 check
leaves only C0=4,B0=-7; the other integer distances exceed.23,.476,.264.
At chi,7chi-4chi² lies in(1.229,1.235), so the remaining W strip forces
1.229<A0<1.471, impossible for an integer. Thus the entire unramified
branch, not merely a native7 head, is excluded.

## Complete ramified table and both maximal orders

If2 ramifies, it divides the field discriminant and therefore the
polynomial discriminant. Reduction of X³-A X²+B X-1 modulo2 has a
repeated factor exactly when BOTH A and B are odd: its derivative is
X²+B; neither the B0 case nor A0,B1 has a common root. This is only a
necessary condition, used in the correct direction (source lines150–158).

The two small positive roots give
B²>4w>4(A-8/17), while the cleared endpoint sign gives
B<(272A+4849)/1156. I independently checked all eleven odd traces
19,21,...,39 against those strict inequalities. The respective lower/
upper odd-B rows are

    (9,7),(11,9),(11,9),(11,9),(11,9),(11,11),
    (13,11),(13,11),(13,11),(13,11),(13,13).

Only (A,B)=(29,11),(39,13) remain. Expanding the cubic discriminants
gives4596=4*3*383 and20084=4*5021. I checked the complete stated prime
tests for383 and5021, including all primes through67 in the latter
case; the displayed remainders are correct. Thus only2 could divide
either power-order index. The shifts are exactly

    X³-26X²-44X-18, X³-36X²-62X-26,

both2-Eisenstein. Their root pi is a local uniformizer with v(2)=3.
For a numerator sum a_i*pi^i not all a_i even, the least term valuation
is <=2 and unique modulo3; division by2 cannot be integral. Any index
factor2 would supply such an integral numerator divided by2, so it is
excluded. This proves O_K=Z[W] in BOTH actual cubic candidates, with no
assumption of maximality or reduced reduction at2.

## The two final head strips and composition

I checked all eight exact rational polynomial evaluations in source
lines213–233. They yield the disjoint root brackets

    P29: .14<u<.15, .23<v<.235,28.615<w<28.63;
    P39: .11<u<.12, .21<v<.22,38.66<w<38.68.

The derivative bounds follow from root differences: in P29 the small
derivative magnitudes exceed28.465*.08 and28.38*.08, both>2.27,
and the physical derivative exceeds800. In P39 the small bound is
38.44*.09>3.45 and the physical derivative exceeds1480.

Writing the ACTUAL head U=a0+b0W+c0W² in the now proved FULL ring,
Lagrange interpolation with the physical value in(4,11) and both other
values in(0,2) gives exactly the two c0 intervals at lines243–244.
Both are strictly inside(-1,1), so c0=0. Then
b0=(U(w)-U(u))/(w-u) is strictly between0 and1, contradicting integral
b0. This argument needs neither a time-square comparison nor a native
label list or equality of the two fork parameters.

The norm-one split is exhaustive: unramified2 has the empty K7 fork
strip; ramified2 has exactly the two full-basis head-strip contradictions.
Together with the separately accepted degree2/3 interface, this excludes
the cubic case and proves that the ACTUAL field at a whole-D17 farthest
vertex is quadratic. It does NOT exclude whole D17, does not remove an
arithmetic forest family, and does not lower17 to16. All quadratic fork
alternatives remain open. No numerical checker, new example, literal-host
nonreflectivity or Lean execution is credited. Six papers and all frozen
inputs remain unchanged.
