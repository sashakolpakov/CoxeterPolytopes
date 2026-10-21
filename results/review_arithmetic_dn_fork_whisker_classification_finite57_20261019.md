# Independent written review of the D19/D20 isolated-d classification

**PASS within the stated hypotheses and inherited dotted-branch input.**
I found no omitted finite continuation or invalid field/trace step in
[the classification report](arithmetic_dn_fork_whisker_quad_attachment_classification_20261019.md),
SHA256 `46fcb083a3ad783544335ac9920bea00f7a7cfc3e1838a1ce1f171a71f243386`.
The report excludes the specified actual quadrilateral with isolated d and
c perpendicular to b,d. It does not establish those contacts for arbitrary
unpaid quadrilaterals, install a global refund, or lower a dimension bound.

This is a full source/paper mathematical read. No matrix calculation,
symbolic program, enumeration, checker, Lean process or build ran. The
identities discussed below were rederived on paper. Hashing only binds
the reviewed bytes and is not mathematical verification.

## Finite continuation and compactness

The actual finite corner W+b+c makes W+c positive definite. Its long
A_L chain has L=16 or17, so any finite irreducible containing it is A, B
or D. The exceptional E,F,H diagrams are too small, and finite diagrams
have no cycles. Thus c has at most one long-arm contact. For a path that
contact is at an endpoint, with label3 or4; adding e beyond c forces the
long bond to3 and the e bond to3 or4. For an internal contact, the D fork
needs two length-one arms, giving exactly l2 or l_(L-1), label3, with no e
contact. The isolated and dihedral-e-only possibilities are also retained.
This proves the complete table without a finite-label cutoff.

The b,d projections are unchanged because both are perpendicular to W.
The projected polygon is acute: eliminating a positive white M-matrix
only decreases off-diagonal products. If the b,d block were positive
definite, its nonnegative inverse would project a strict timelike polygon
point to a feasible common corner of these opposite sides. The positive
null combination excludes the parallel case by compactness. Hence y>1.
This reasoning does not impose ambient disjointness on a,c.

The zero a.c cases all close. The simple internal attachment gives the
stated affine-D positive-null block, while extending the far endpoint
gives an impossible positive-definite D_(n+1). The endpoint4 attachments
give the stated affine-B block. Both-simple contacts give a chord-free
affine-A cycle. The far long/simple plus e/4 case has the five-node
affine-B block. All positive-null certificates satisfy the stated neighbor
equations and are genuine induced blocks; unused long-arm vertices create
no extra edges inside them.

For the remaining near long/simple plus e/4 case, the old cycle puts
sqrt2 in K. At its negative embedding the projection cost is
n/4+(n-3)/2-(n-3)/sqrt2. Its excess over1 is positive exactly when
(3n-10)^2>8(n-3)^2, with both quantities before squaring positive.
The difference n^2-12n+28 is positive for n>=9. The physical projection
cost exceeds1, so the independent-block admissibility test is legitimate.

## Schur formulas and embeddings

I independently eliminated W and the old D_n block. With N=n-2 the
remaining a,b Gram has diagonal s=n/(4N),1 and product-1/2. The c residual
has square t and product-X with a, where X=x+kappa. It follows that

    q_c=1-t+2N*X^2,
    (X^2/t-s)*(y^2-1)=1/4,
    z=4+t/(X^2-s*t).

The denominator cannot vanish. Since x>=1/2 in every finite nonzero
a.c branch, q_c>1 in all applications. Thus the D_n+c roots are independent,
and conjugate positive definiteness yields the strict inequalities used.

The report properly distinguishes normalized contacts from ground-field
elements. Simple old-path cycles put x in K; a long4 cycle instead puts
2sqrt2*x in its integer ring. In the joint e/4-long/3 case, the two cycles
put both x and sqrt2 in K. In the sole-fork and fork/e cases it uses squared
contacts and their actual cycle product, rather than conjugating a missing
square root. These are determinant-ratio arguments with positive diagonal
scalings at the conjugate embeddings.

## Simple long contacts and their last quadratic traces

The cosine trace formulas in Section5 are correct, including the different
odd/even m correction in Tr(x^2). For0<=kappa<=3/4 their combined correction
is at least-1/2. Comparing at degree4 gives exactly
B=1/2+3kappa^2-2kappa-3r^2. All six table values in Section6 follow from
this expression and are positive. The degree<=3 cosine list is exactly
3,4,5,6,7,9. The stated elementary conjugates exclude every case except
the near m4 and far m5 cases with both simple contacts.

For near m4 I obtain
X^2-s*t=(N-1)/N+kappa*sqrt2, whose quadratic norm is
-(n^2-8)/(8N^2). Its trace yields

    Tr(z)=4+4(3n-8)/(n^2-8).

For far m5 one has kappa=t=s, x^2=x/2+1/4 and
X^2-s*t=((n-1)/N)*x+1/4. Its trace gives

    Tr(z)=4+2(7n-8)/(n^2+2n-4).

Both lie strictly between4 and5 at n19 and20. A trace in the intrinsic
quadratic subfield must be integral if z is an algebraic integer, even
when the arithmetic ground field is larger. No ground-field list is used.

## Long4 contacts and the CRT case

For the single long4 contact, w=2sqrt2*x is indeed an integral cycle product.
At the far endpoint its nonphysical conjugates have absolute value<4/N,
and its physical value is in[sqrt2,2sqrt2). The norm bound forces degree1,
w=2. At the near endpoint A=w+2 has physical value in[2+sqrt2,2+2sqrt2)
and other conjugates in(0,4/N). Degree>=3 is impossible. In degree2 its
norm is1 and its trace4 or5; the reported resulting cosine-squares3/8 and
(11+sqrt21)/16 are correct. Their exclusion uses the elementary degree1/2
cosine-square lists, not an arithmetic field classification. Degree1 has
A=4. Both remaining z formulas are rational nonintegers as stated.

For a long/simple and e/4 contact, the CRT exponent1 modulo2m and5 modulo8
exists precisely in the asserted branch4 not dividing m: the common
2-part is at most4, and5=1 modulo that common part. It is a unit modulo
the common conductor, fixes x and negates sqrt2. Its restriction extends
to K and is nonphysical. If4 divides m, the exponent m-1 is coprime to2m
and sends x to-x; it negates sqrt2 when m=4 modulo8 and fixes it when
8 divides m. The numerical rational brackets in the latter case have
strict margins: cos(pi/8)>9/10, sqrt2/4<5/14 and3/70>1/34. They exclude
every label in (12), without a hidden finite scan.

## No long-arm contact: squared fields and the final cases

In the sole-fork case, the integer v=4x^2 has physical value in[1,4) and
every nonphysical conjugate in(0,2/N). The integer norm excludes degree>=2.
Substitution for v1,2,3 gives the three stated negative/noninteger z values.
It is not necessary to assume x belongs to K.

With contacts on a,e, completing the beta square gives the invariant lower
bound q_c>=(4N/n)*sigma(u). Since u,v are totally positive squared finite
cosines and p^2=u*v, the signed real square-root interpretation at each
embedding is valid. For m>=4, any nonidentity embedding fixing u would
give q_c>=2N/n>1. Total reality therefore forces K=Q(u); no Galois
assumption on K is required.

The cosine-square trace bound gives degree<=4, hence exactly the displayed
phi(m)<=8 list. Each listed table exponent is coprime to m, nonidentity
modulo sign, and at most3m/10; its conjugate exceeds2/7. The surviving
rational and quadratic cases are therefore exhaustive.

The remaining norm/square-class exclusions are correct. In particular
for m8 the three unwanted norms are1/128,1/32,9/128, none a rational square.
For m12, x lies in sqrt2*Q(sqrt3), leaving only ell4 or12. Equal labels
5,8,12 have p=u and q_c'=(9n-16)u'/2>1, as claimed. I independently obtain
the last traces

    (m,ell)=(5,3):
      8-36(n-2)(n-3)/(11n^2-26n-1),
    (m,ell)=(12,4):
      4+4(7n-15)/(3n^2-6n-1).

The first is4504/869 and2224/431 at n19 and20, strictly between5 and6;
the second is strictly between4 and5. Their integrality contradictions hold.

For m3, the actual cycle puts beta in K. The factorization of q_c-1 is
correct and forces every nonphysical beta below-(n-4)/n. Thus A=2(beta+1)
has degree at most2 by its positive integer norm. At degree2 its trace
is strictly greater than3 and less than5, hence4, and its norm is1; this
gives ell6. Degree1 gives ell3, the affine triangle. The ell6 trace is
the same noninteger quadratic expression. This closes the final branch.

## Dependencies and remaining boundary

The inherited reports read for this review were

* [double-whisker obstruction](arithmetic_dn_double_whisker_obstruction_20261018.md),
  SHA256 `6432c844bb7b79a5e2df71f89e4e895ec540d763b0b77c6cd04dcfd4fb88bd69`;
* [isolated dotted/parallel obstruction](arithmetic_dn_isolated_quad_obstruction_20261018.md),
  SHA256 `aac7e6251b5d2b85ad0687cca39c2cb84bca581df5d863d8f41023aa65be9a95`.

The new finite-product proof uses no published field list. Its combination
with the inherited a.c<=-1 branch still depends on that report's cited
high-dimensional arithmetic field theorem and rational isotropy input.
Those published inputs were not re-established in this review.

The completed statement is local to a whole D_n finite vertex, n19 or20,
an actual quadrilateral, isolated d, and c perpendicular to b,d. Nonisolated
d and different c cross contacts remain open. No universal facet-payment
refund, simultaneous current consumer audit, full forest bound, executable
certificate or kernel proof follows from this written PASS alone. The
compact arithmetic ceiling remains20.
