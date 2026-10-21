# Root review: complete isolated-side dotted long-contact exclusion

19 October 2026. **PASS, full written-source review.** This records
written mathematics, not a numerical run or a Lean verification.

Sources read in full:

* `arithmetic_rank20_endpoint3_basis_free_norm_obstruction_20261019.md`,
  SHA256 `03809292c144e54ed85ddbdf91e0796426116641e9b00bb34a44222f89f02f39`;
* `arithmetic_rank20_isolated_d_dotted_long_unit_closure_20261019.md`,
  SHA256 `5c583ce01db154908d5d0e139c61f4986cf1491cd4f3f991f0c7e26c8a5a92c6`.

The finite inventory, height interfaces and common three-type norm
obstruction have separate root reviews. The new arguments were
developed by root and independently checked by batch_guard; this
record is an integration check, not a claim of two independent
discoveries.

## Endpoint3 alone

The integer norm of delta and strict nonphysical bounds give the
physical bounds `A<3*2^(d-2)` in the far role and
`U<3*2^(d-2)+1/4` in the near role. I checked both comparisons at
N8 and their monotonicity in N. Sign exclusion and integral norms
then force A, respectively U, to be a unit. The resulting conjugate
intervals make `2A+1` or `2U-1` have nonzero norm of absolute value
below one in degrees at least four. The special degree-four near
estimate is correctly separated from the simpler degree-at-least-five
bound. The cubic trace, reciprocal-trace and endpoint-sign constraints
give exactly the three stated polynomials of negative discriminant.
These are characteristic polynomials over the full field, so proper
subfields and repeated embeddings are covered.

The far quadratic remainder is A=1+sqrt2. Substitution gives
`Norm(C)=-(7N^2+54N+31)`. Its forced positive integer norm quotient
lies strictly between one and sixteen sevenths, hence must equal
two. That would require `(N-19)^2=384`, impossible for integer N.
This closes the quadratic remainder without a field restriction.

## Endpoint3 plus e3

Put M=N+2. In the far role the product-norm bound makes A and A+1
units; in the near role it makes U and U-1 units. The near estimate
uses U>=4 and `C>0`, with the difference
`(2/3)U(U-4)>=0`; the strict bound on C supplies the strict final
norm bound even at U4. The norm products are below two, so their
nonzero integer norms are exactly one in absolute value.

In each role the remaining factor B divides M and N because the
other factor is a unit. Thus B divides two. Every one of its full
field conjugates has absolute value greater than two, including
the physical one. This contradicts the integer norm of a divisor
of two. The comparison `t-sqrt(t)>2` is valid for t>=5.

The exceptional E8 k3/k5 cases retain their actual integral cycles:
C is odd, divides four, and is therefore a unit, while every
conjugate has absolute value greater than one. No E8 contact is
discarded from the inventory.

Together with the earlier common-form lemma, this excludes the
complete long-contact inventory for an isolated d, c orthogonal
to the entire spectator, and explicitly dotted x>=1, over every
totally real ground field. Zero/finite x and spectator contacts
require their separate proofs. No other first-d type, triangle
consumer, simultaneous refund, or new dimension ceiling is
established by this composition alone.
