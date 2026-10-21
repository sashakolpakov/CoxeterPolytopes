# D18 dotted degree-four theorem and finite-head precision delta

21 October 2026. **FULL OUTSIDE-AUTHOR MANUAL MATHEMATICAL REVIEW: PASS**
for the current dotted-head degree-four source, and **PRECISION-DELTA PASS**
for the current finite-head successor. No scientific/symbolic computation,
enumeration, worker import, Lean or TeX ran; no reviewed file was edited.
The only programmatic verification here was metadata-only hashing,
gzip decoding and exact text comparison.

Full reviewed source:
`results/arithmetic18_whole_d18_dotted_head_degree_four_and_isolated_forks_20261021.md`,
SHA `0590f2b9d8747f6edbe0cafd74f2909c54c87f62e5ebe7d52081abbdfe35fa16`
(12,648 bytes, 292 lines). I did not author this source. I contributed to
older actual-foot/prefix work and reviewed its direct predecessors; that
inherited involvement is disclosed. I independently rederived the new
degree/divisibility and cubic/quartic arguments without relying on their
previous PASS labels.

## Full mathematical review of the current dotted-head source

* **Inputs and generic fork interface (lines 22–155).** All three direct
  source pins match current bytes: 6ef34bb4, 026e39cd and 48089176.
  The .61 floors and actual U normalization in 6ef apply to dotted heads,
  not only to finite heads. The C<13 paid certificate uses those floors
  alone. I rechecked the complete A17 finite-white extension list, all six
  raw H^2/C formulas, their actual-field equality and the four generically
  excluded roles. Unit normals in the affine/null vectors are now
  explicitly distinguished from the raw K-valued cap normalizations.
  The two possible end helpers have K=Q(t) or Q(z), with nonphysical
  magnitudes <3/5 and physical upper bounds25/4 or13/2. Their cycles
  give Ut/2 or Uz/2 integral, without requiring U to be a unit at2.
  An isolated fork has raw norm2V, time V(9V-4)/2, K=Q(V), and
  0<V_sigma<4/9. Its actual head/fork cycle gives UV/2 integral;
  subtraction from V also gives (2-U)V/2 integral. No finite-head
  field/unit assumption, rational-pair-integrality shortcut or larger
  whole-diagram exclusion is imported.

* **Degree5/6 obstruction (lines 157–182).** The imported Q bound at d5
  is |N(Q(U))|<656000/11664<57, hence at most56; at d6 it is at most5.
  In these degrees the two possible helper norms are <.81 and<.8424,
  respectively, so both actual forks must be isolated. The V norm bounds
  are exactly173056/59049<3 and692224/531441<2. Thus N(V)=1 or2
  at d5, and1 at d6. Taking norms of BOTH integral cycle quantities
  gives the stated divisibilities by2^d. At d5 the nonzero first two
  Q-norm factors are each multiples of16, forcing product>=256; at d6
  each is a multiple of64, forcing product>=4096. These contradict
  the respective Q bounds. All norms are over the actual field, and
  all factor signs are handled by absolute value. One fork suffices;
  neither its energy nor a distinct second opposite is assumed.

* **U/2 cannot be integral (lines 184–205).** If W=U/2 is integral,
  then W(1-W) has physical magnitude<154/9 and other magnitudes<=1/4.
  This excludes every d>=4. At d3 its norm is exactly -1, since only
  the physical factor is negative. Hence N(W)=1, N(1-W)=-1 and its
  minimal polynomial is X^3-TX^2+(T-1)X-1 with T=3,4,5,6.
  I recomputed the discriminants -23,-31,-23,49. At T6 the value at5
  is -1; the monic polynomial has a root>5, incompatible with the
  sole physical root<14/3 and two roots in(0,1). This closes all four
  cases, rather than assuming a unit classification.

* **Quartic helpers (lines 207–218).** Their nonzero norm magnitudes
  are <27/20 and<351/250, hence1. The integral parameter is a unit.
  Multiplication of Ut/2 or Uz/2 by its integral inverse would make
  U/2 integral, already impossible. Thus every quartic fork is isolated.

* **Complete cubic nonunit helpers (lines 220–267).** Their norm
  magnitudes can only be1 or2; the unit cases are eliminated as above.
  For the long helper, one positive nonphysical conjugate would give
  magnitude<3/2, so both are negative and N(t)=2. The trace is5 and
  t is in(50/9,25/4). The decreasing second coefficient
  2/t+5t-t^2 lies strictly between -2997/400 and -5521/2025.
  Thus the ONLY integer possibilities are -3,-4,-5,-6,-7. I manually
  recomputed their complete discriminant list
  -1315,-1172,-883,-424,229. In the last case the values -1 and23/125
  at -1 and -3/5 force a forbidden nonphysical root in(-1,-.6).
  For the other helper, a negative nonphysical conjugate would give
  magnitude<39/25, so both are positive with N(z)=2. Their product
  >4/13 gives sum>1 and the integral trace is7. The norm constraint
  excludes z>=6; hence29/5<z<6. Its decreasing second coefficient
  z(7-z)+2/z is between19/3 and5296/725, leaving only7. The resulting
  cubic has discriminant -59. Every coefficient case is covered by
  these small hand derivations; no polynomial census is assumed.

* **Final fork conditions and scope (lines 269–292).** In d3/4 the
  isolated finite labels5/8/12 would force a quadratic actual field,
  while3/4/6 and parallel have the actual affine/null exclusions.
  Therefore both forks are dotted. The displayed V norm bounds give
  integer N(V)<=14 or<=6 and retain BOTH cycle divisibilities. This
  does not classify the quadratic case, exclude the surviving cubic/
  quartic heads, prove whole D18 or lower the dimension ceiling.

No mathematical correction is required. The full current theorem is
degree<=4, U/2 nonintegral at d>=3, and isolated dotted fork opposites
in the cubic/quartic cases, with exactly the written residual conditions.

## Preserved predecessors and finite-head precision successor

Preservation manifest
`results/arithmetic18_d18_unit_null_vector_clarification_20261021.json`,
SHA `645949e03378ee5e389636ebaf78c1cdac0745a078c5197713ee31ed497084ca`,
was independently checked against compressed bytes, decoded bytes and
both current sources:

| Source | Preserved decoded SHA / bytes | Gzip SHA / bytes | Current SHA / bytes |
|---|---|---|---|
| finite-head source | fc32360368d7136de045ae542de75323e523ec1036f936eb25a9e651e14106e1 /16635 | 1e2448f19d46bba69d94199df8b683c9b7db29f6b5ec2ab4472f370e11119be3 /7575 | 48089176640b49f414e6a0fb56a12cfc1658873c4cf70816e9f8d43506b8986a /16828 |
| dotted degree-four source | f84c0823cc4bb7f70c00e7b1d62f94c83a7178bf5418d773b511a8b90efa8591 /12455 | 71d667e73ea6ae24ee3a2df4253663f3d2c78456f23c4bed2c7a3ea993cce658 /5652 | 0590f2b9d8747f6edbe0cafd74f2909c54c87f62e5ebe7d52081abbdfe35fa16 /12648 |

The exact diffs contain ONLY the same 193-byte explicit convention
sentence in each source, plus the degree-four source's dependent finite
SHA update. No scalar, polynomial, norm, theorem hypothesis or conclusion
changed. The new convention is correct: the displayed positive affine/
null combinations require UNIT cap c; the time/cycle tables use their
separately stated raw normalizations. This matches the interpretation
used in my complete finite-head review
`results/review_arithmetic18_d18_finite_head_fork_exclusion_full_scope_domain_reviewer_20261021.md`,
SHA `58b9888a8906f96486b2893936813896d42dc121b3964e086d27e97e1620c760`,
which remains immutable and bound to fc323603. Its full mathematical
verdict therefore carries to current finite-head source48089176 by this
checked precision delta, not by silently replacing its historical pin.
No execution or new bound is credited by either verdict.
