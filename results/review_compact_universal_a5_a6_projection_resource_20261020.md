# Outside review: universal actual A5 and A6 projections

20 October 2026. Full paper-only review by resource_cleanup. No
scientific execution, imports, enumeration, optimization, Lean or TeX.
No accepted source or status was changed. These are next-batch review
items, separate from the already reviewed large-component transfers.

**Verdict: PASS for both stated actual-projection bounds.**

* `results/compact_universal_a5_actual_projection_three_eighths_20261020.md`,
  SHA256 `fd6f024d975e525184de05961d4d6943f1b0bc8eef4acc279dde98cbc1bc4456`:
  every whole A5 component at a farthest vertex has squared actual
  projection greater than3/8, with arbitrary finite spectators.
* `results/compact_universal_a6_actual_projection_three_sevenths_20261020.md`,
  SHA256 `022c0cc77512576dfc7bed74be372a94956618498c9ed74508aff53c40a50ae2`:
  the corresponding A6 projection is greater than3/7, again with
  arbitrary finite spectators.

No correction or unresolved mathematical gap was found. These are
actual-vector estimates using farthest geometry and actual walls;
neither is promoted to a theorem about arbitrary edge-cover vectors
or a new unconditional dimension ceiling.

## Shared geometry and spectator bookkeeping

I independently checked the displayed inverse A5 matrix and the A6
inverse formula, the normalized correlations, and the conventions
for full versus half edge tangents. The actual distinct-opposite
law and strict concavity of log F have the correct directions.
All contact bonds are finite and native because each old node is
retained in at least one of the two given actual replacements.

Every retained spectator component has at most one contact, hence
cost at least1/4, with no cost strictly between1/4 and1/3. A
nontrivial connected component has inverse unit diagonal at least
4/3; the next A1 native square after1/4 is1/2. If a contacted
component cost q is at least the full squared wall height, reciprocal
shadow gives respectively an actual orthogonal-component endpoint
or a feasible null projection. Both contradict the stated farthest
or compactness hypotheses. Thus the strict condition q<A^2 is
justified for every contacted component; it is not an assumption
that foreign contacts disappear harmlessly.

The shadow, actual-point and orthogonal-component input pins agree
with the accepted sources previously audited here: `c913ec3b...e8badc`,
`3732b826...f6ed84` and `4216d06c...2f23d`.
The reviewer contributed to checking some of these older interfaces,
but did not author either new A5/A6 proof. The new classifications,
dual constants and actual foot were independently rederived in this
audit.

## A5: complete common-pair classification

The central correlations are exactly1/sqrt2,1/sqrt2 and1/2. The
rational F-seal gives sum>109/100 in each distinct central pair.
For the combined vector, its exact norm is7+4sqrt2; the positive
Cauchy margin209/35000 is correct. Under projection at most3/8,
there is therefore an actual common2/3 or, after numbering reversal,
the same case with all assumptions intact.

I checked all six possible two-contact positions. Actual deletion2
or3 gives the stated rank-five paths or D5, forcing all3 except
4@1+3@3 and the adjacent2,3 profile. Their all3 alternatives contain
genuine affine cycles. The remaining 4@1+3@3 has white cost7/8,
so Q=0. For adjacent2,3, the native pairs3/4/5 and the exact length
formula are complete; its minimum outside3,3 is3,4, with the
reversed4,3 larger. In both cases the actual cosh bound is
`(3+sqrt2)/2>11/5`, giving squared half tangent>3/8.

The single-contact list is also exhaustive. The retained old3
neighbor rules out labels>=6; the middle of the retained A3 allows
only3. For single3, every own cost and Q bound in the five-row
table checks. With two spectators in row2, q<A^2 forces both
costs>1/3, impossible under Q<2/3. In row3 both costs must be
>1/4, so the universal cost gap gives each at least1/3, again
impossible. The other rows have too small a height or only one
possible foreign component. No spectator rank/type is dropped.

Single4 at an internal node contains the affine3,3,4 star; an
endpoint has insufficient height for any contacted spectator.
For single5, all listed retained H3/H4 costs force Q=0. I rederived
all four length formulas from the inverse Gram. In particular the
new direction4 for contact1 is recognized using the complete PD
H4+A1 replacement of actual walls, not by assuming a formal new
wall exists. Each listed length exceeds11/5; the phi^3 estimate
uses the exact identity phi^3=2phi+1. This exhausts the common pair
and proves the A5 bound.

The A5 forest consequence correctly adds its ACTUAL projection to
the previous orthogonal-component accounting. The exceptional
short-A1/rank-two7 cases are retained, so the stated dominant-D
alternative follows without an ambient-ceiling assertion.

## A6: pair reduction and every finite profile

I rederived the three relevant correlations and the three pair-sum
seals. The weighted vector has exact squared norm
`84/5+8sqrt(8/15)+8sqrt(3/10)` and the source's rational upper
bound33781/1250. The Cauchy margin2181583/700000000 is positive,
forcing a common2/3 or3/4 pair after reversal.

The double-contact tables cover all eight position pairs for2/3
and all nine for3/4. The E6 exceptions at1,4 and1,5 in the former
case are correctly retained and force all3. The D6 rows and the
rank-six paths likewise force all3, leaving only the adjacent
profiles and4@1+3@3. Every discarded all3 profile has the asserted
affine cycle in the actual full diagram. The adjacent length
formulas, monotonicity in both contacts and Q, and threshold5/2
check. The last double profile has exact white cost9/10, so Q=0.

All single3 own costs and Q bounds check in both tables. The
two-spectator repair is essential and valid: in the remaining
internal rows Q<5/8 forces BOTH first bonds to be3. The six actual
nodes are exactly the affine two-fork D5 diagram, irrespective of
the rest of either spectator component. This principal diagram
gives the required compactness contradiction. The other rows
vanish directly by the strict cost/height condition or its cost-gap
consequence. Single4 is excluded by the affine star or the endpoint
height bound. The complete single5 positions are exactly those
listed; their H3/H4 white costs force Q=0. Each selected endpoint
is actual by positive acute projection, and each displayed cosh
formula exceeds5/2.

## A6: the remaining ACTUAL positive foot

For4@1+3@3, I multiplied the three old-root corrections explicitly.
All coefficients in
`xi=c+(sqrt2-1/2)n1+(sqrt2-1)n2+(sqrt2/2)n3`
are positive. The products with c,n1,n3 are zero, that with n2
is(sqrt2-3)/4<0, and xi^2=5/4-sqrt2<0. Every unknown wall pairs
nonpositively with every term; spectators pair zero. Its unchanged
positive time coefficient makes xi future. It is therefore an
ACTUAL feasible point in the full polytope.

The source's explicit U in the sum-zero R7 model and
`A^2=(5+4sqrt2)/7`, `U^2=A^2-k`, `k=sqrt2-5/4` all check.
Farthest comparison and Cauchy consequently give the stated ratio
`(A-sqrtk)/(A+sqrtk)>1/2`. The positive seal
`A^2-9k=(335-236sqrt2)/28>0` is correct, including its integer
square comparison. This genuinely uses an actual point, not only
edge-floor information, and closes the sole remaining A6 branch.

The final profile-specific extension to r>=6 has the correct
increasing height formula. Its scope remains that EXACT profile,
including zero foreign contact; it is not a general larger-A
projection theorem. A7/A8 rates and the unclassified large-A
mechanism remain outside these two reviewed results.
