# B18 short-anchor/prefix signature census: full pre-execution source review

21 October 2026. **PASS at the stated necessary-screen boundary.** No
blocking mathematical, completeness or implementation defect was found.
This is a manual source review, not a census execution, output verdict,
whole-B18 exclusion or dimension-seventeen claim. No worker was imported
or run; no scientific, Lean or TeX job was launched.

## Reviewed sources and contribution boundary

| Source | SHA256 |
|---|---|
| scripts/check_b18_short_anchor_prefix_signature_census_20261021.py | 6ae6944b6c3300cfec600bbbd956904838f59dc971fe41e6e42c64ca11bb9b1a |
| results/arithmetic_whole_b18_short_anchor_prefix_signature_census_spec_20261021.md | a51427783ba85b47f1088b420db3325db877259095626f856c27831f70a80aad |
| scripts/check_a18_isolated_cap_domain_and_norm_screen_20261021.py | 29da7434634c80fb001784d9d5f0d69c3d6247a7b6de1e55e2a0eb6b7ff22932 |
| results/arithmetic_whole_b18_short_cap_inventory_and_degree_three_reduction_20261021.md | 0d9cffcc58700e38985e31865380457bbbf885d0c329cb6c53ea69c089f675e1 |
| results/arithmetic_whole_b18_point48_weak_prefix_inventory_and_cap_census_domain_20261021.md | 04dfdc05a18b6aa2d7d7ba71a0e6b3824601dce14a797e113a038f209ff003cb |

I read the complete 425-line worker and 254-line specification, all
invoked local helper functions and their import/main boundaries, and the
relevant geometric normalization/domain proofs. I authored the accepted
04df inventory/domain note and supplied independent coefficient-box and
norm derivations used here; I also contributed to inherited arithmetic
normalization/actual-foot inputs. I did **not** author this worker or
specification. This review independently checks their implementation and
logical interface; it does not represent a second outside-author review
of my own earlier geometric classification.

## Domain and exact algebraic decisions

Worker lines31–100 and189–275 implement the stated domains without an
unproved truncation. The anchor loops have 299 quadratic and 924 cubic
coefficient entries, plus the finite silver anchor. Their signed N/M
ranges include every negative constant from the two helper families.
The isolated prefix loops have 420 quadratic and 516 cubic entries; the
strict bound N<32/i is implemented by floor(31/i), including i8 and16.
All twelve typed rows are retained, including the changed B2 native12
row at index16. There is no cubic typed row hidden by a field-list test.
These are raw domain counts, not predicted retained counts.

I checked the coefficient bounds from trace, symmetric sums and norm.
For example the isolated anchor has T<38, N<676/81<9 in degree two,
and M<1352/729<2 and S<1356/81<17 in degree three. For far t and near z,
the physical and small-conjugate bounds give exactly the safe boxes
displayed in specification lines42–46; allowing negative N/M is essential.
The polynomial convention is X²−TX+N and X³−TX²+SX−M throughout.

The cubic rational-root test (worker42–49) checks both signs of every
divisor of abs(M), so a negative M cannot produce an empty divisor loop.
For degrees two and three this proves irreducibility after the real-root
filters. Irreducible characteristic-zero polynomials are squarefree and
have no rational roots. The imported Sturm construction uses strictly
decreasing Euclidean degrees and exact rational signs; interval endpoints
and every bisection point are therefore safe from polynomial-root equality.

The two root counts in worker88–100 cover all d roots in disjoint physical
and nonphysical intervals. I checked the irrational anchor tests directly:

* far: t>2 and t²>8 choose t>2sqrt(2); the inequality
  81t²+18t−3210<0 is exactly (9t+1)²<169·19 on this branch;
* near: z>4 and z²−4z−4>0 choose z>2+2sqrt(2); the inequality
  81z²−18z−3210<0 is exactly (9z−1)²<169·19;
* the nonphysical tests (9/2)t²+t−1<0 and (9/2)z²−z−1<0
  have respectively the exact roots (−b,a) and (−a,b) from the specification.

Thus the rational outer intervals introduce no retained extraneous branch
in these exact strip tests and omit no required root. Isolated roots use
their exact rational intervals directly. Polynomial reduction detects
equality before refinement: a zero remainder is rejected by strict tests;
a nonzero remainder of degree<d cannot vanish at the selected root.
Rational interval Horner bounds converge under bisection. Root separation
and continuity prove termination of both isolation and sign refinement;
no timeout or tolerance is interpreted as a negative decision.

## Raw norms, field signatures and quadratic squares

Worker121–186 recomputes every time as v²−N0 with raw B dual products
Omega_i·Omega_j=min(i,j), then compares the resulting multiplication
determinant with an independent permutation determinant. I rederived the
closed formulas in specification151–160 and worker148–162, including the
signed cubic terms for negative M. In particular:

* isolated prefix: U=alpha(i alpha−2), raw norm2alpha;
* isolated short anchor: U=B(9B−2)/2, raw normB;
* far short anchor: U=(9/2)t²+t−1, raw norm2;
* near short anchor: v=(z/2−1)Omega18+Omega17 gives
  U=(9/2)z²−z−1, not the unshifted far polynomial;
* finite short silver: v=Omega17+(sqrt(2)/2)Omega18 gives
  U=24+17sqrt(2), raw norm2.

All twelve typed prefix projections agree with their displayed times.
In particular native12 at index16 has raw norm1 and
v=((1+sqrt(3))/2)Omega16+Omega18/2, hence
U=(55+32sqrt(3))/2 and norm−47/4. The code does not replace its raw
normalization by an externally adjoined square root.

Physical time is positive and every other conjugate negative, so its
nonzero field norm is negative in degree two and positive in degree
three. For a reduced rational norm n/d the code uses
sign(n)·squarefree(abs(n)·d). This is the correct rational squareclass,
including the denominator8 in cubic anchor norms. The recorded integer
square factor k means n/d=sign(n)·part·(k/d)².

Equal ACTUAL fields imply equal degrees and polynomial-discriminant
squareclasses, since an integral generator's power-order discriminant
differs from the field discriminant by an index square. A K-square time
ratio implies equal signed norm squareclasses. Worker362–402 uses only
these necessary tests in degree three and retains **every** matching
prefix with an explicit unresolved field-map/K-square status. Equal
cubic discriminant squareclasses are never promoted to field isomorphism.

In degree two, the common positive squarefree discriminant D identifies
the field. The image (T+k sqrt(D))/2 in worker352–358 is the unique larger
physical root in every domain, including all three typed generators.
The exact K arithmetic/sign implementation and complete square decision
in worker279–359 are correct: both signs of sqrt(Norm(x)) are tried;
u²=(a+r)/2, v=b/(2u) covers u!=0; the b=0, u=0 branch tests
v²=a/D. This includes pure-radical roots. Exact resquaring and a positive
physical sign check certify every accepted root, and the displayed
necessary equations prove that no K-square is rejected.

The fixture assertion at worker403–408 is valid and does not assume an
empty census. By hand, for s=sqrt(2),

    (16+12s)/(24+17s)=12−8s=[−2+2s]²,

with −2+2s>0. Both models are present and have matching signatures.
The mandatory output root ["-2","2"] is therefore correct. This checks
an arithmetic pair only, without asserting a realizable chamber.

## Imported closure, output and acceptance limits

The worker invokes only the helper's digest, trim, value, remainder,
sturm, roots, squarefree, reduce_mod and determinant, including their
internal calls. I checked all of them. The helper has only standard-library
imports, no dynamic/local further imports, and a guarded main. Its old A18
domain loops, CLI, proof pins and scientific census are not invoked by
this worker. Its exact elimination copies its input matrix; no model or
domain aliases are mutated by the norm calculation.

Worker366–425 hashes both premises, the imported helper, specification
and executing source before and after the bounded work, checks all fixed
pins, preserves every retained model and every anchor's target/rejection/
retention rows, and creates a previously absent result under results/.
It writes no frozen input. Degree selection is explicit in the output;
the default both is required for the complete two-degree screen. The
specification hash is recorded rather than hardcoded, so post-run
acceptance must compare it with the reviewed a5142778 pin above.

The source is suitable for the proposed single guarded execution, subject
to the unchanged parent resource gates. **No execution or retained-count
conclusion is supplied by this review.** Later acceptance must bind the
actual terminal source/spec/helper/output hashes and read every survivor.
Even then quadratic matches certify only same-field time squares;
cubic matches remain necessary signatures. Mutual Gram/native-angle
conditions, actual geometric coexistence and a same-original-vector
argument for all surviving weak sets remain separate obligations.
