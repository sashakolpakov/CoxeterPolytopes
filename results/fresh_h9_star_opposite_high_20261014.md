# Opposite single-contact high-label family in star case0

**This family has no compact H9 Coxeter realization for any finite high internal label.** The first three blocks retain both leaves with own contacts `(3,3)`. Leaf11 is right to12 and has π/3 contact to13; leaf17 has the opposite profile. Their missing triples force the last internal label `m>=7`. The literal swap of12 and13 covers the other orientation. All unlisted finite labels and dotted distances remain free.

The [exact checker](../scripts/fresh_h9_star_opposite_high_20261014.py) and [certificate](fresh_h9_star_opposite_high_20261014.json) exhaust47,040 residual cases. Of these,43,412 have positive Bernstein coefficients on a continuous interval containing the entire high-label tail;3,612 have exact nonpositive actual-vertex witnesses;16 use a sharper exact cosine bound. No unresolved case or finite label cutoff remains. The result excludes the [saved105-visit necessary assignment](fresh_h9_star_shared_high_20261014_search_batch_005.json), not the whole incidence.

The reviewed two-dense gates reduce the first three internal labels to4 or5. Actual four-root supports force every cross-block contact right by the relaxed nonpositive path4–3–4. The eight block roots have a positive direct-sum Gram.

The retained-own-root budgets are rebuilt after changing the last contacts: anchor11 retains13, anchor17 retains12, and14 retains its first-block own root. Each anchor has at most one foreign π/3 contact, giving seven patterns. The actual triangle `{12,13,14}` forces14 right to the high block, leaving five patterns for14.

The remaining finite contact domains differ from the shared-root family. The actual triangle `{11,13,17}` excludes a mutual label at least6 by the path3–6, so the mutual label lies in `{2,3,4,5}`. Fresh exact actual-vertex witnesses exclude label4 for11–14 and label5 for14–17; their domains are respectively `{2,3}` and `{2,3,4}`. The witnesses relax the last internal label to4, and all omitted contacts to right. Their failure persists under increasing any cosine. No stronger contact bound is imported from the shared-root theorem.

Let `c=cos(pi/m)` and `D=(1−c²)^{-1}`. Projection off the first six block roots gives the three-by-three residual Gram
\[
R=\begin{pmatrix}
a-D/4&b-cD/4&e\\
b-cD/4&d-D/4&f\\
e&f&q
\end{pmatrix}.
\]
Using `D(1−c²)=1`, its determinant satisfies
\[
(1-c^2)\det R=B(1-c^2)+A+cC,
\]
where
\[
B=\det\begin{pmatrix}a&b&e\\b&d&f\\e&f&q\end{pmatrix},\quad
A=-\frac{(a+d)q}{4}+\frac{e^2+f^2}{4}+\frac q{16},\quad
C=\frac{bq-ef}{2}.
\]
The numerator is a quadratic in c with coefficients in `Q(sqrt2,sqrt5)`. The checker verifies the symbolic identity before the coefficient census. No assumption is made that c belongs to the coefficient field.

Since `m>=7`, we have `c>sqrt3/2>6/7`. There are exactly `8*7*7*5*4*2*3=47,040` combinations of the internal choices, projection patterns and contact labels. In43,412 cases all three exact Bernstein coefficients on `[6/7,1]` are positive. The corresponding residual determinant is nonzero throughout the required interval, giving rank11 in a ten-dimensional ambient space.

In3,612 cases the interval test is replaced by a saved exact nonpositive actual-vertex LDL pivot at the relaxed last internal label4. The finite-diagram classifier only locates the witness; every rejection is directly checked in the exact field.

The remaining16 occurrences have the same numerator
\[
N(c)=\frac{1+\sqrt2}{2}c^2-\frac{2+\sqrt2}{16}c-\frac{11+9\sqrt2}{32}.
\]
They are positive on `[9/10,1]`. Its Bernstein coefficients there are
\[
\frac{54\sqrt2-41}{800},\qquad
\frac{35\sqrt2-4}{320},\qquad
\frac{1+5\sqrt2}{32},
\]
all strictly positive.

The refined bound `cos(pi/7)>9/10` is exact. For `p(x)=8x³−4x²−4x+1`, the identity `T7(x)+1=(x+1)p(x)²` makes `cos(pi/7)` a root. The cubic is strictly increasing on `[3/4,1]`, since its derivative and second derivative are positive there, and `p(9/10)=−1/125`. Together with `cos(pi/7)>sqrt3/2`, this proves the bound. Every larger finite m has an even larger cosine.

Thus all cases are covered continuously. Generation and saved replay are available through:

```sh
python scripts/fresh_h9_star_opposite_high_20261014.py --verify
```

The primary imports the frozen exact field helper; independent implementation scope is stated in the separate review.

The [independent review](fresh_h9_star_opposite_high_20261014_review.md) passes generation and saved replay. Its standalone Fraction implementation imports no producer or project arithmetic. It reconstructs all47,040 residual cases through separate determinant interpolation, including4,363 distinct five-point polynomial checks, verifies all3,612 vertex failures, and checks the fresh literal contact domains and exact cosine bounds.
