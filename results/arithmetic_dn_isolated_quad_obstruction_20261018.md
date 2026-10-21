# An arithmetic obstruction to one isolated-wall D_n quadrilateral branch

18 October2026. **New paper argument, awaiting independent review.** No
enumeration, numerical solve, script, Lean compilation or TeX build was
run for this argument. It does not lower the current compact arithmetic
ceiling20. Its opposite-wall hypothesis is explicit and essential.

## Statement and boundary of the claim

Let P be a compact arithmetic hyperbolic Coxeter polytope of dimension
n≥15, with the usual simplest-type arithmetic criterion and published
compact ground-field restriction. Suppose a finite vertex has whole
diagram D_n. Let a be its fork root, b one of the two short whiskers,
and W its other n−2 roots. Normalize wall normals to square1.

Consider an actual quadrilateral face at the marked angle(a,b), with
opposite side pairs(a,c) and(b,d). Suppose

1. d is orthogonal to every old D_n root except b;
2. c is orthogonal to b and d;
3. **a.c≤−1** in the ambient Gram matrix.

Then this configuration is impossible.

The third assumption is NOT inferred just from opposite sides of the
two-face: ambient facets can meet away from that two-face. This result
only removes the indicated ambient dotted/parallel branch. It does not
remove finite ambient products a.c=−cos(pi/m), prove every fork angle
refundable, or establish a new dimension bound.

At an unpaid quadrilateral's isolated-A1 adjacent corner, assumptions1–2
are the relevant zero-contact pattern: an adjacent marked pair cannot
have zero Vinberg weight, so zero weights at(b,c) and(c,d) force the
stated orthogonalities. Assumption3 must still be tested separately.

## 1. Exact geometric Schur relation

The white Gram matrix H is A1+A_(n−3), with diagonal1 and simple
edges−1/2. It is positive definite and H^(-1) is entrywise nonnegative.
Project a,c orthogonally off W. Their squared norms and product are

    s = n/[4(n−2)],     0<t≤1,     a_bar.c_bar=−X, X≥1.

Indeed a has product−1/2 with the A1 node and the near endpoint of
the A_(n−3) chain. For a length-k A-chain with unit-root Gram, the
inverse endpoint entry is2k/(k+1). Thus

    s=1−1/4−(n−3)/[2(n−2)]=n/[4(n−2)].

The actual finite corner W+b+c gives t>0; orthogonal projection gives
t≤1. Acuteness and H^(-1)≥0 give

    X=−a.c + a_W H^(-1)c_W ≥ −a.c ≥1.

Write b.d=−y. Since b and d are both perpendicular to W, their
projections are unchanged. The two-face is acute: projection off a
positive-definite white M-matrix only decreases off-diagonal products.
Its opposite b/d lines cannot intersect. Otherwise the positive-definite
two-wall projection, whose inverse has nonnegative entries, produces an
actual common vertex satisfying all other strict inequalities. Parallel
opposite lines similarly produce a feasible ideal point by their positive
null sum. Both contradict a compact quadrilateral. Hence y>1.

Eliminating H leaves, in order a,b,c,d, the matrix

    [[s,   −1/2, −X,  0],
     [−1/2, 1,    0, −y],
     [−X,   0,    t,  0],
     [0,   −y,    0,  1]].

There are n+2 actual normals in dimension n+1, so its determinant
vanishes. Eliminating the positive entries t and1 gives

    (X²/t−s)(y²−1)=1/4.

Since X²/t≥1,

    1<y²≤(4n−10)/(3n−8).

Consequently z:=4y² satisfies

    4<z≤4(4n−10)/(3n−8)≤92/17.

The last, slightly loose endpoint is the n=14 value; the expression
decreases for n≥14.

## 2. Arithmetic conjugates of the isolated wall product

The length-two cyclic product of the normalized Gram matrix gives
z∈O_K. This is the standard arithmetic cyclic-product input, rather
than an assumption that the unsquared normalized product y belongs to K.

The inverse D_n Gram entry at a whisker is n/2. One direct derivation
uses ordinary D_n roots e_i−e_(i+1), e_(n−1)+e_n: the spinor fundamental
weight has squared norm n/4. Unit normalization multiplies the dual
squared norm by2, giving n/2.

The old D_n roots together with d have determinant ratio

    1−(n/2)y²=1−nz/8.

It is physically negative and nonzero. At every nonidentity real
embedding of K the arithmetic conjugate form is positive definite.
Conjugating this determinant ratio therefore gives

    0<sigma(z)<8/n≤4/7.

Strict positivity follows also from z≠0 and injectivity of the field
embedding. This argument is invariant under the positive diagonal
normalizations of the arithmetic roots at each embedding.

In particular z generates K: if it lay in a proper subfield, a
nonidentity real embedding of the totally real K fixing that subfield
would fix z>4, contradicting the displayed bound.

## 3. A totally positive algebraic-integer interval lemma

Put U=92/17 and v=4/7. Let z be a totally positive algebraic integer
with one conjugate in(4,U] and all other conjugates in(0,v). The only
possibilities are degree1 with z=5, or degree2 with

    z=(5+sqrt(17))/2 or (5+sqrt(21))/2.

Here is a proof with exact bounds.

**Degree at least5.** Its positive integer norm is less than

    U v^4=23552/40817<1,

which is impossible. Further factors only decrease the bound.

**Degree4.** The norm is an integer below

    U v^3=5888/5831<2,

so it is1. Writing the three small conjugates r,s,t gives

    z>1/v^3=343/64,
    r,s,t>1/(U v²)=833/1472.

Its integer trace is therefore strictly greater than

    343/64+3*(833/1472)=2597/368>7,

and strictly less than U+3v=848/119<8, a contradiction.

**Degree3.** The norm is an integer below U v²=1472/833<2,
so it is1. Write its integer polynomial as

    p(T)=T³−A T²+B T−1.

The trace satisfies A<U+2v=780/119<7, hence A≤6. If the two small
conjugates are r,s, then

    B=z(r+s)+rs≥2sqrt(z)+1/z>4,

so B≥5. But

    p(v)≥v³−6v²+5v−1=29/343>0,

whereas (v−z)(v−r)(v−s)<0. This is impossible.

**Degree2.** The integer trace lies in(4,U+v), with U+v=712/119<6,
so it is5. The positive integer norm is below Uv=368/119<4, hence
is1,2 or3. Norm3 is impossible because its polynomial at v is

    v²−5v+3=23/49>0,

instead of the negative product(v−z)(v−r). Norm1 and2 give exactly
the two displayed quadratic values, in fields Q(sqrt21) and Q(sqrt17).

**Degree1.** The integer in(4,U] can only be5.

## 4. Field composition and use in the current search

For the compact arithmetic dimensions in question, the published field
restriction permits quadratic ground fields only Q(sqrt2), Q(sqrt3),
Q(sqrt5), Q(sqrt6). It excludes Q(sqrt17), Q(sqrt21). The interval lemma
therefore excludes every nonrational K. The rational case is separately
impossible for a compact simplest-type arithmetic chamber in dimension
at least4, by rational isotropy of the defining indefinite form.

This proves the stated obstruction, in particular at n=19 and20.
The field input is Vinberg's published high-dimensional compact field
restriction (Theorem3 of the primary paper), already used at dimension21
in the current compact20 package. The positive-conjugate Gram and
cyclic integrality conditions are the arithmetic criterion used in those
same proof packets.

The result removes an infinite ambient dotted-parameter branch without
using a numerical angle cutoff. The finite-angle branch a.c>−1 remains
open. Other non-isolated d continuations remain open too. The complete
current compact20 certificate is unchanged, and no new weight or total
budget has been installed.

Primary reference: E.B.Vinberg, *Absence of crystallographic groups of
reflections in Lobachevsky spaces of large dimension*, published Theorem3
and the arithmetic criterion; primary PDF:
https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt.
These new inequalities are our deduction from those inputs, not a quoted
published obstruction or worldwide-priority claim.
