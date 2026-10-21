# Independent review: the isolated-wall arithmetic D_n quad obstruction

18 October 2026. **PASS for the stated local obstruction, relative to
the explicit published arithmetic inputs.** This is a read-only
mathematical review of
[the complete focused report](arithmetic_dn_isolated_quad_obstruction_20261018.md).
No enumeration, numerical calculation, proof script, Lean compilation
or TeX build was executed. The author's report and all frozen evidence
were left unchanged.

The statement concerns a full D_n vertex in ambient dimension n≥15,
with all three listed contact hypotheses, including the explicit
ambient product a.c≤−1. It does not classify every quadrilateral at a
fork/whisker angle or establish any new charge refund or dimension bound.

## 1. The complete white support and the Schur identity

The white support really has rank n−2 and type A1+A_(n−3). This uses
the full D_n vertex hypothesis. It would be invalid to use the same
matrix unchanged for a D_n factor with additional old components to
which c or d could attach.

For unit normals, the endpoint inverse entry of an A_k chain is
2k/(k+1). Therefore eliminating the two white branches at the fork
gives

    s=1−1/4−(n−3)/(2(n−2))=n/(4(n−2)).

The actual finite corner W+b+c gives t>0, while orthogonal projection
onto the positive-definite white span gives t≤1. With a_W,c_W
entrywise nonpositive and H^(-1) entrywise nonnegative, the residual
cross product is −X with

    X=−a.c+a_W H^(-1)c_W≥−a.c≥1.

The assumed zero contacts leave exactly the displayed four-node
residual matrix. The n+2 actual roots lie in a vector space of dimension
n+1; hence its determinant vanishes after elimination of the invertible
white block. Eliminating the c and d diagonal entries gives exactly

    (X²/t−s)(y²−1)=1/4.

There is no unnoticed full-star rank assumption here: the complete
white support and all four actual side normals have been retained.

Since s<1 and X²/t≥1, this identity itself implies y²>1. Acuteness
gives y≥0, so y>1. Thus this part can be obtained without asserting
that arbitrary opposite ambient facets are disjoint. The report's
separate geometric argument for the unchanged b,d projections is also
valid: the positive two-normal projection uses an inverse M-matrix,
preserves feasibility against every other projected wall, and would
produce an actual b/d corner if those opposite lines intersected.
A positive null sum excludes the parallel case by compactness.

The upper bound follows with the correct direction:

    y²≤1+1/(4(1−s))=(4n−10)/(3n−8).

Thus 4<z≤U=92/17 is a valid loose uniform bound. Using a non-strict
upper endpoint avoids any need to discuss separately the impossible
compact parallel pair a.c=−1. Most importantly, the proof does not
infer a.c≤−1 from opposition within the two-face; that remains an
explicit hypothesis.

## 2. Arithmetic conjugation, positivity, and the ground field

The arithmetic input is the cyclic-product criterion in the physical
ground field K, together with the positive conjugate defining forms.
With unit normals the integral two-cycle is 4(b.d)²=z. No assertion
that the individual normalized product y belongs to K is needed.

The D_n inverse entry at either short whisker is n/2. It follows that
the old D_n block together with d has Schur determinant ratio

    1−(n/2)y²=1−nz/8.

This is physically negative and nonzero, so these n+1 roots are an
independent full set. In a simplest-type arithmetic realization one
can work with unnormalized roots over K. At each nonidentity embedding,
the defining form is positive definite and this independent Gram
matrix remains positive definite. Positive diagonal normalization then
restores the D_n unit-root ratios and leaves the stated determinant
ratio unchanged.

For clarity, positivity of sigma(z) has a little more content than
injectivity alone. In unnormalized arithmetic roots it is four times
the squared embedded b/d inner product divided by the two positive
embedded root norms. Hence it is nonnegative; since z is nonzero,
injectivity excludes zero. The determinant ratio gives the strict upper
bound. Therefore every nonidentity embedding satisfies

    0<sigma(z)<8/n≤4/7.

Even if only positive semidefiniteness were used initially, the weak
upper endpoint could not be attained: sigma(z)=8/n, a rational number,
would imply z=8/n by injectivity, contradicting z>4.

The primitivity argument is valid. If Q(z) were a proper subfield of
the totally real K, separability supplies a nonidentity real embedding
of K extending the identity on Q(z). That embedding fixes z>4 and
contradicts its small-conjugate bound. Thus K=Q(z), including the
rational case. This is the required bridge between the algebraic degree
of z and the field restrictions applied at the end.

## 3. Independent check of every interval-lemma degree

I checked the exact rational endpoints and all signs by hand.

For degree at least five, total positivity makes the norm a positive
integer, but

    Uv^4=23552/40817<1,  v=4/7,

excludes it. Larger degrees only insert further factors less than one.

For degree four, Uv^3=5888/5831<2 forces norm one. Each small conjugate
exceeds 1/(Uv²)=833/1472, while the large one exceeds 1/v³=343/64.
Hence the integer trace is greater than 2597/368>7 and less than
848/119<8. Both strict inequalities are justified by the strict small-
conjugate interval even though z≤U is allowed.

For degree three, Uv²=1472/833<2 again forces norm one. If the minimal
polynomial is T³−A T²+B T−1, the trace bound gives A≤6. Its other
two roots r,s have rs=1/z, so

    B=z(r+s)+rs≥2sqrt(z)+1/z>4.

Thus B≥5. Evaluation at v gives p(v)≥29/343>0, whereas the factorization
(v−z)(v−r)(v−s) is strictly negative. This excludes the whole cubic
case without any irreducibility or trace subcases.

For degree two, the trace is strictly between four and six and so is
five. The positive integer norm is at most three. Norm three gives
p(v)=23/49>0, contradicting the negative two-root factorization. Norms
one and two give exactly the two possible large roots
(5+sqrt21)/2 and (5+sqrt17)/2. Degree one permits only z=5.

The interval lemma therefore proves precisely the field alternatives
claimed; it has no unexamined degree or sign branch.

## 4. Published composition and exact limits

The source's n≥15 range safely avoids ambiguity at the endpoint of
Vinberg's published high-dimensional field theorem. In that range its
quadratic alternatives are Q(sqrt2), Q(sqrt3), Q(sqrt5), Q(sqrt6), so
neither quadratic field surviving the interval lemma is allowed. This
uses the published field classification; it is not derived from the
new local calculation. The earlier project arithmetic proof map lists
the same explicit cyclic, simplest-type, and Meyer inputs in
[its input register](../docs/compact_arithmetic21_proof_map.md).

If K is rational, simplest-type arithmeticity and Meyer isotropy in
at least five variables exclude compactness. Here n≥15 is well inside
that range. Thus the rational z=5 interval possibility does not evade
the geometric conclusion.

The primary source is
[Vinberg, the published field theorem and arithmetic discussion](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).
The review retains those published results as inputs and independently
checks the new geometric and algebraic composition above. It does not
reprove the field theorem, the arithmetic criterion, or Meyer isotropy.

No flaw was found in the stated local obstruction. The finite ambient
opposite-product branch a.c>−1, additional contacts of d, vertices
other than full D_n, and the global consumer/budget argument remain
outside its conclusion. In particular it does not change the current
compact arithmetic ceiling20 or the separate unconditional compact29
frontier.
