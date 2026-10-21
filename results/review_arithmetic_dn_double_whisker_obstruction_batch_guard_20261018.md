# Independent review: arithmetic double-whisker wall obstruction

19 October 2026. **PASS for the exact local theorem stated in the
[author's report](arithmetic_dn_double_whisker_obstruction_20261018.md).**
I read that report in full and independently checked its Gram, trace,
and integrality calculations by hand. No proof program, enumeration,
numerical solve, Lean compilation or TeX build was executed. The
author's source and all earlier evidence were left unchanged.

The conclusion concerns a full D_n finite vertex in ambient dimension
n≥16, with two additional actual walls having exactly the stated
isolated-whisker contacts and c.d=0. It gives no global refund and no
new dimension ceiling. No published ground-field classification is
used.

## 1. Unit normalization and the rank identity

For the ordinary square-two D_n root system, the two spinor fundamental
weights have squared norms n/4 and mutual product (n−2)/4. Replacing
the simple roots by unit roots multiplies these dual weights by sqrt2.
Thus the two needed inverse Gram diagonals are n/2, and their cross
entry is (n−2)/2, with the signs used in the report.

The c-to-old contact vector has its sole entry −sqrt(u) at e; the
d-to-old vector has its sole entry −y at b. The old inverse therefore
gives the residual matrix

    [[1−nu/2,          −sqrt(u)y(n−2)/2],
     [−sqrt(u)y(n−2)/2, 1−ny²/2]].

Expanding its determinant yields

    1−n(u+y²)/2+(n−1)u*y².

It must vanish because the complete n+2 actual normals lie in Lorentz
dimension n+1, and the old D_n block is invertible. This is a legitimate
rank equation for the full stated wall set. It does not apply a rank
condition to a smaller incomplete face star.

Solving gives exactly

    z=4y²=4(nu−2)/((2n−2)u−n).

If the denominator vanished, the original equation would also require
nu=2, forcing (n−2)²=0. Hence there is no missed denominator-zero
branch in the stated range. The physical assumption y>1 implies z>4.
It is retained as a hypothesis, not inferred from two-face opposition.

## 2. What is conjugated in the arithmetic argument

The arithmetic inputs are the usual simplest-type admissibility and
Vinberg cyclic-integrality conditions. The two-cycles put 4u and z in
the arithmetic ground field, and make z an algebraic integer. This
does not require either unit-normal square root to belong to that
field. Since rational scaling does not change a generated field,
Q(u) is a subfield of the ground field.

The D_n+c determinant ratio is 1−nu/2. It is negative and nonzero:
u≥1/4 and n≥16 give nu/2≥2. Thus the n+1 roots in this principal
block are independent. Using arithmetic root representatives over
the defining field, their independence survives every field embedding.
The admissible nonphysical defining form is positive definite, so
their embedded Gram is positive definite.

Positive normalization at that embedding gives the same inverse
whisker diagonal n/2 and the squared contact sigma(u). Its Schur
complement is therefore 1−n*sigma(u)/2>0. The squared embedded
inner product divided by positive embedded root norms is nonnegative;
injectivity, since u≠0, excludes zero. Consequently

    0<sigma(u)<2/n≤1/8.

This supplies the missing sign information which would not follow
from a bare determinant inequality alone. Sign changes of individual
embedded root representatives do not affect the squared contact or
the inverse diagonal.

Every nonidentity real embedding of Q(u) extends to an embedding of
the totally real ground field. Its extension is nonphysical because
its restriction to Q(u) is nonidentity. Thus the bound applies to
every nonphysical conjugate used in the subsequent trace argument,
whether or not the ambient ground field is larger than Q(u).

## 3. The cyclotomic trace excludes all larger degrees

For every m≥3, u=(2+zeta_m+zeta_m^(-1))/4 and
Q(u)=Q(zeta_m+zeta_m^(-1)) has degree d=phi(m)/2. Pairing the primitive
mth roots into conjugate pairs gives trace mu(m) for the real generator,
so

    Tr(u)=d/2+mu(m)/4≥d/2−1/4.

The argument uses zeta_m, not zeta_(2m), and remains correct for even m.
The physical value u<1 and the strict conjugate bounds give

    Tr(u)<1+(d−1)/8=(d+7)/8.

For d≥3 the lower endpoint (4d−2)/8 is at least this upper endpoint.
The strict inequality gives a contradiction, including at d=3 where
the rational endpoints coincide. Thus d≤2.

The list m=3,4,5,6,8,10,12 is complete. A prime divisor must be 2,3
or5; the constraints on their powers and products under phi(m)≤4
give exactly these values after m=1,2 are removed. There is no finite
label cutoff imposed by a search.

## 4. Independent check of the seven remaining labels

The rational substitutions give:

* m=3: z=2(8−n)/(n+1)<0;
* m=4: z=−2(n−4)<0;
* m=6: z=6+2/(n−3), rational and strictly between 6 and 7.

The first two contradict z>4. The last contradicts integrality of a
rational algebraic integer.

For m=8, the other conjugate (2−sqrt2)/4 exceeds 1/8 because sqrt2<3/2.
For m=10, (5−sqrt5)/8 exceeds 1/8 because sqrt5<4. These violate the
conjugate bound.

For m=5 I independently used the identities u+u'=3/4 and uu'=1/16
in the fractional-linear expression for z. With D=n²−4n−1, they give

    Tr(z)=4(2n²−3n−12)/D=8+20(n−2)/D,
    N(z)=−4(n²−24n+64)/D=−4+20(4n−13)/D.

Thus 4Tr(z)−N(z)=36+100/D exactly. D≥191 and is increasing for n≥16,
so the last fraction is strictly between zero and one. Both the trace
and norm in Q(sqrt5) would be integers if z were integral. Integrality
descends to this subfield even if the chamber's ground field is larger.
This excludes m=5 without any field list.

For m=12, the identities u+u'=1 and uu'=1/16 give

    Tr(z)=4(3n²+n−16)/(3n²−6n−1)
         =4+4(7n−15)/(3n²−6n−1).

The added fraction is positive and less than one. The relevant
difference is 3n²−34n+59, positive at n=16 and increasing thereafter.
Hence the trace is strictly between four and five and cannot be an
integer. This excludes the last label.

## Scope of the PASS

All Gram and algebraic steps in the stated theorem check. The external
inputs retained are the simplest-type arithmetic realization, positive
definiteness at the other embeddings, and cyclic-product integrality.
The proof does not use a ground-field catalogue, a chamber census, or
an assertion that cos(pi/m) itself belongs to the ground field.

The full D_n vertex and the complete zero-contact pattern are essential
to this calculation. Long-arm contacts, extra old components, different
placements of c, nonzero c.d, or removal of the y>1 hypothesis require
new arguments. No assertion about all fork quadrilaterals, a universal
face payment, or an improved dimension bound follows from this local
obstruction alone.
