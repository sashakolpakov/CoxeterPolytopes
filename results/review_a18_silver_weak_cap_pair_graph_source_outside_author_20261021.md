# Whole A18 silver pair graph: full outside-author source/domain review

21 October2026. **PASS BEFORE EXECUTION.** No source correction is required.
I read the complete112-line mathematical specification and266-line worker,
after independently deriving the pair criterion. No worker was imported
or run; no scientific, Lean or TeX job, chamber construction or new bound
is asserted by this review.

The frozen sources are:

* `results/arithmetic_whole_a18_silver_weak_cap_pair_graph_spec_20261021.md`,
  6182bytes, SHA256
  `38f7a3fec648ef3bf6ba274ad700db47694bba079acd5f54059e5118711d58eb`;
* `scripts/check_a18_silver_weak_cap_pair_graph_20261021.py`,9147bytes,
  SHA256 `b2a8750e4ca32a38dfa8ca07d95785fb82a21108bcf53b67b7babb34e80bfc1e`.

Both direct mathematical pins match their current bytes: complete inventory
`4710025219a3d8a2d5a485de27d3a9f6962fa61289ef9e6120ee880ee4eec618`
and actual-field composition
`fb185ba3e154863520a62729348339ebde51443f5f94c77ff0abaa528c96485a`.
The independent derivation is
`results/arithmetic_a18_silver_pair_test_independent_derivation_20261021.md`,
SHA256 `5486b0cb8dee9dd619259e4201521af5149a82b1af7b75ad9033f978b2f69ddf`.

**Domain and normalization.** Spec lines9–29 and worker138–158 include
exactly20 far3 models,4 correctly oriented/reflected near4 models, and
20 isolated-dot models at all central indices5..14. The two isolated
parameters3+2sqrt2 and6+4sqrt2 both satisfy their strict per-index
Schur strips. Raw old roots have square2, inverse
min(i,j)(19-max(i,j))/19. Helper normals have square2 and projection
t*Omega_i+Omega_j; isolated normals have square2beta and projection
beta*Omega_i. Therefore `spatial(projection,projection)-norm` is precisely
the raw time square in EACH case. The physically positive/conjugately
negative time and totally positive raw norm checks have the right signs.
No sqrt(beta) normalization or accidental unit-root inverse is used.
All946 unordered model pairs are generated, including the76 same-index
pairs, which are correctly impossible as two distinct outgoing opposites.

**Field operations and square test.** Worker30–126 implements exact
Q(sqrt2) arithmetic using reduced Fractions. Opposite-sign field values
have sign(a)*sign(a^2-2b^2); all zero and same-sign cases are correct.
The inversion denominator is nonzero for a nonzero field element.
`positive_field_sqrt` tests BOTH signs of the rational square root of
the norm. For nonzero rational u it reconstructs v=b/(2u) and verifies
the full field square. Its u=0 branch includes pure-sqrt2 roots of a
rational value; the other branch includes pure rational roots. The final
physical sign is selected exactly. These equations prove completeness,
not merely necessity of a norm-square test. No floating point or field
extension is introduced.

**Common time and cyclic integrality.** At worker165–181, the ratio
F_b/F_a must be a field square because actual residuals lie on ONE
K-line. The physical positive witness fixes A_a*A_b=F_a*z; its exact
product and squaring are checked before forming R=P_a.P_b-F_a*z.
There is no freedom to use the other time orientation.

Spec39–52 correctly treats MIXED raw norms. Each helper has an outgoing
unit helper contact, each isolated cap an outgoing marked Cartan entry
-1. The old path and the incoming contact at the second cap cancel
its possible beta denominator in the returning new edge. Thus a
longer cyclic product is plus or minus R, even when the chosen old
contacts coincide. The NEW two-cycle is independently
4R^2/(N_a*N_b); its integrality is required as well, including for
dotted pairs. Worker176–181 imposes BOTH tests. This is essential:
R integral alone does not generally cover the nonunit raw beta case.
All other new-edge simple cycles are R times already integral outgoing
contact ratios; cycles avoiding the new edge are individually calibrated.
The criterion is therefore not missing another cyclic denominator.

**Physical angle filters.** Positive R is rejected by actual all-wall
acuteness. Zero R is orthogonal. For negative R, the two-cycle B is
four times the squared normalized cosine. B<4 gives a positive
two-wall Gram; positive-inverse projection proves an actual ridge,
so the native angle filter is legitimate even without first assuming
these two facets adjacent. Q(sqrt2) permits only m3,4,6,8 at nonzero
finite angle, with B=1,2,3,2+sqrt2; m2 is zero. The other conjugate
2-sqrt2 is not such an angle. B=4 is correctly rejected by compactness
using the actual feasible positive null sum. B>4 is dotted with no
artificial finite-label cutoff, and still requires both integral cycles.
No positivity of conjugate mutual PRODUCTS is imposed; only the correct
conjugate Gram interpretation from the common time relation is needed.

**Known survivor.** I independently recalculated far(5,18),far(7,1):
F_b/F_a=32/25=(4sqrt2/5)^2, spatial product(202+141sqrt2)/19,
positive time product(240+160sqrt2)/19, hence R=-2-sqrt2 and
B=6+4sqrt2>4. The IDs[1,4] are exactly these two models. The regression
assertions at237–240 therefore preserve a genuine necessary dotted
pair rather than hardcoding an unjustified empty compatibility graph.

**Finite algorithm/output scope.** The standard-library-only worker has
no imported census or hidden top-level execution. The model list, exact
pair loops and finite graph traversal terminate. Each pair receives
one explicit first-failure reason or the full passing witnesses; reason
counts depend on this declared order, not on additional mathematics.
Components include isolated vertices, and no clique calculation or
inference from connectedness occurs. Pair data and model objects are
not mutated through aliases. The before/after input/source pins, parent
ownership requirement, new results-only destination and exclusive output
creation match the stated contract. Actual resource enforcement remains
the existing guarded runner's duty; its future run is not supplied here.

Passing pairs remain NECESSARY local arithmetic possibilities. The worker
does not prove simultaneous weak heights, a whole chamber, an invariant
lattice or any dimension exclusion. Those distinctions and the pending
execution are explicit in its output and specification. Any later norm
argument must retain all weak indices in one paid reference; the known
pair cannot be discarded. The separate independent derivation records
possible actual-foot inputs, but no energy closure is assumed here.

I did not author this worker or its specification. I supplied the missing
two-cycle observation before the frozen source was finalized and have
earlier shared actual-foot/reference work; those contributions are
disclosed. The full domain, exact square algorithm, native filters and
control/output composition were independently checked for this review.
No review of my separate pending B18 inventory is claimed.
