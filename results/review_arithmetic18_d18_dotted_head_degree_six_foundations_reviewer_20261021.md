# Full written review: the D18 dotted-head degree-six norm bound

21 October2026. Reviewer: sharp_ra_foundations_audit. **PASS at the exact
scope below; no correction requested.** Manual source and mathematical
review only: no scientific, symbolic, enumeration, import, Lean or TeX
job. I authored the sole geometric input6ef34bb4 and therefore do not
claim an independent audit of that foundation here. I did not author
the reviewed root polynomial/norm argument and rederived all its steps.

Reviewed source:
`results/arithmetic18_whole_d18_dotted_head_degree_six_norm_obstruction_20261021.md`,
SHA256 `026e39cd20e780c8fbc5f554615ed58bacf7d6f5dc0829aecd40005fb33d220a`,
4,258 bytes/113 lines. Its sole direct source pin6ef34bb4 matches
`arithmetic18_whole_d18_long_end_isolation_and_head_field_reduction_20261021.md`
at full SHA256
`6ef34bb417647b846c9b26b7c35fb8018e27b999c0395396c33c6d5ba691b3ac`.

The imported interface is exactly an ACTUAL integral two-cycle U,
K=Q(U), 4<U<28/3 and every other U_sigma in(0,2). No quadratic-field
claim, normalized root integrality or larger-rank endpoint is imported.
At these inputs the new polynomial proof is complete:

* With t=(X-1)^2, Q(X)=X(2-X)(X-1)^2(2X^2-4X+1)
  equals t(1-t)(2t-1). For X in(0,2), t is in[0,1).
  Differentiating gives -6t^2+6t-1; the two critical points
  (3 plus or minus sqrt3)/6 give values plus or minus sqrt3/18.
  Checking the endpoints proves the claimed bound on the ENTIRE
  interval. The point X=1, where Q=0, causes no norm exception:
  Q(U) is nonzero physically and embeddings are injective.
* Physically 9<t<625/9, and t(t-1)(2t-1) increases for t>1.
  The upper numerator is exactly625*616*1241=477785000,
  with denominator729. This is a strict upper bound on |Q(U)|.
* c^6=1/1259712 and729*1259712=918330048. Thus degree>=7
  forces the nonzero integral norm below477785000/918330048<1.
  No restriction on unlisted embeddings or field discriminants is used.
* At degree6, c^5=sqrt3/209952. The displayed conservative bounds
  A<656000 and sqrt3<7/4 give norm<4592000/839808<6.
  The norm magnitude is consequently an integer1..5.
* Multiplicativity has precisely the squared U-1 factor. U is totally
  positive and all factors are nonzero algebraic integers, so every
  magnitude is a positive integer. Therefore the other three factors
  are at most5 and |N(U-1)| is1 or2. In the latter case its squared
  factor4 forces all other factors1 and total norm4. No sign was
  dropped where positivity is needed.
* The final coefficient box follows by separating elementary symmetric
  terms that contain the one physical root from terms involving only
  the d-1 roots in(0,2). Both strict interval bounds are valid, including
  k=d with the second binomial coefficient zero. It is a finite necessary
  domain, not a claim that any polynomial or field list was generated.

The result proves degree<=6 and the additional degree6 norm conditions
for the remaining actual dotted long head. It neither excludes that
head nor whole D18, and proves no arithmetic17 ceiling, new example,
wall inventory, computation or geometric Lean theorem. In particular
the separately discussed stronger degree<=4 lead is NOT part of this
review. Frozen papers, sources, guides and receipts were not changed.
