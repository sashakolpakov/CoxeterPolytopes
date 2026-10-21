# Independent review of the D4 vertex-frequency obstruction

The [D4 incidence argument](fresh_sqrt5_d4_short_face_incidence_20261017.md)
passes this internal mathematical review. Its scope is a compact arithmetic
Coxeter polytope of simplest type over Q(sqrt5), of dimension d=4k. A special
vertex has full diagram kD4. The conclusion is

    (4k−9) V_D < (12k+2) V_other,  k≥3.

This is a vertex-frequency restriction. It is neither a dimension exclusion
nor a permission to reduce payments at arbitrary D4 corners. The unclassified
mixed D4/D5 vertices remain a substantive obstacle.

The [separate checker](../scripts/fresh_sqrt5_d4_short_face_incidence_review_20261017.py)
imports no producer arithmetic. It reconstructs literal Gram matrices with
SymPy and obtains inertia by exact congruence elimination, including zero
diagonal pivots. Its [saved certificate](fresh_sqrt5_d4_short_face_incidence_review_20261017.json)
contains all four reattachments, 36 triangle cases, 22 Schur cases, the
symbolic residual congruence, and the counting identities for arbitrary k.

## Arithmetic and signatures

An edge joining two special vertices must restore the deleted root's exact
contacts. An intact D4 cannot be joined to another root and remain a component
of a diagram consisting only of D4s. Within the affected component, deletion
of a leaf leaves A3 and deletion of the center leaves three A1s; their unique
D4 reattachments are checked independently.

The old and new roots share a neighbor with both normalized products −1.
Their triangular cyclic product is therefore their mutual product g.
Arithmeticity gives g in the ring of integers of Q(sqrt5). If these two roots
meet in a finite vertex, the finite Coxeter classification and the degree of
2cos(pi/m) leave labels 2, 3, and 5 in this field. The complete exact signature
checks leave only a leaf replacement with label5. At a third finite corner,
the center must be deleted: otherwise a 3–3–5 triangle remains and is
indefinite. Its replacement diagram is 2A1+I2(5), so a triangle cannot have
three special vertices. A triangle with two special vertices uses a
center–leaf pair at each.

For a quadrilateral with four special vertices the two replacements have
identical contacts to the common spherical support and to each other.
Eliminating that support gives the four-root residual matrix in the producer.
Taking differences and sums of opposite roots gives, up to a positive factor,

    diag(lambda_i, lambda_j) ⊕ [[2a−lambda_i, 2c],
                              [2c, 2b−lambda_j]].

Here lambda_i=2−g_i and lambda_j=2−g_j. Physically both lambdas are positive.
The four polygon normals span the Lorentzian three-space orthogonal to the
white support. Thus the residual matrix has rank3 and inertia (2,1,1),
and the displayed final block has rank1 and is negative semidefinite.
Conjugate positivity makes its conjugate positive semidefinite. The two
conjugate lambdas are also strictly positive: zero would give g_i'=2 and,
by injectivity of conjugation, g_i=2, contrary to the physical sign.

For different D4 components, c=0. Each actual hyperbolic edge separately gives
lambda_i>2a and lambda_j>2b, making the residual determinant nonzero. This
contradicts rank3. For two leaves of one D4, a=b=4/3 and c=−2/3. The determinant
and signatures give

    A=−3g_i−2, B=−3g_j−2, AB=16,
    A,B>0 and A',B'<0.

In particular A and B are integral and congruent to1 modulo3. The integer
norm N(A) is negative, divides256, and is1 modulo3. In Z[phi], the norm is
x²+xy−y²; its reduction modulo2 is zero only at (0,0). An even norm therefore
forces both coordinates even, and its 2-adic valuation is even. Its only
negative possibilities dividing256 are −1,−4,−16,−64,−256, all2 modulo3.
This contradiction excludes the leaf–leaf quadrilateral. Only center–leaf
corners can occur in an all-special quadrilateral.

## Global counting

Every special vertex has exactly3k center–leaf corners. Give each triangle
with two special vertices one unit at its other corner and one half at each
special corner. A triangle with at most one special vertex can receive its
two units at two other corners. A quadrilateral with an other vertex can
receive its unit there, while an all-special quadrilateral receives one
quarter at each special corner. Every charged special corner is center–leaf
and receives at most one half. Each other corner receives at most one.
Therefore

    2T+Q ≤ (3k/2)V_D + binom(4k,2)V_other.

Distinct two-faces at a simple vertex occupy distinct corners, so these
charges do not reuse a corner. The symmetric nonnegative h-vector inequality,
including its strict endpoint contribution, gives

    (2k²−3k)(V_D+V_other) < 2T+Q.

The checker verifies the resulting coefficient identities symbolically for
all k. Together these inequalities prove the stated frequency restriction.
No finite rank sample, ambient compact completion of a local example, or
unstated minimum-distance assumption enters this deduction.

Generation and a fresh saved replay are the reproducible arithmetic checks:

```sh
python3 scripts/fresh_sqrt5_d4_short_face_incidence_review_20261017.py
python3 scripts/fresh_sqrt5_d4_short_face_incidence_review_20261017.py --verify
```

The published arithmeticity, Coxeter classification and h-vector facts are
mathematical inputs. This is an internal review, not external peer review or
formal verification.
