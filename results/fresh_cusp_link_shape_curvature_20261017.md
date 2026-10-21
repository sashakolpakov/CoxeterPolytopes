# Exact curvature from the complete cusp-link shape

This is a stronger global counting inequality for finite-volume hyperbolic Coxeter polytopes. It retains all ranks in each cusp's product of simplices, rather than replacing them by a single worst-case average. It supplies no ambient dimension bound on its own. An exact local example below shows why the coarsest distance payment still needs improvement.

Let P have dimension n>=3, let d=n−1, and let D5 be the sum of5 minus the number of sides over its two-faces. Write v_f for its number of finite vertices. At an ideal vertex v let the Euclidean link be

` L_v = product_i Delta^(r_vi),       sum_i r_vi=d. `

Define

` H_v(t)=product_i (1+t+...+t^(r_vi))=sum_(j=0)^d h_vj t^j, `

` V_v=H_v(1),       G_v=sum_j h_vj min(j,d−j), `

` kappa(L_v)=(5/2)G_v−dV_v/2. `

With beta_n=n(n−6)/8 for even n and beta_n=(n−1)(n−5)/8 for odd n, the proved inequality is

**` D5 >= beta_n v_f + sum_(ideal v) kappa(L_v) + 5(n−1). `**

No hypothesis on cusp simplicity, factor count, minimum factor rank, arithmeticity, or induced Coxeter faces is imposed. Finite volume supplies simple vertex figures (products of simplices at ideal vertices) and a bounded projective closure. The geometric facts are the same as in the [factor-payment framework](fresh_cusp_factor_payment_framework_20261017.md).

## Proof retaining the local cut count

Choose a generic linear function on the bounded projective polytope. At every vertex w, form its d-dimensional simple vertex figure L_w and the hyperplane section induced by the level hyperplane through w. Let E_w count all edges of L_w and E_w^0 count edges not met by that section. The exact reduction identities in [Khovanskii, Theorem9, printed p.49](https://www.math.utoronto.ca/askold/1986-Faa-1-english.pdf) are

` I = sum_w E_w,       2 f2(P) = sum_w E_w^0, `

where I is the total number of vertices of all two-faces. Therefore

` D5=sum_w [(5/2)E_w^0−E_w]. `

Here is a direct lower bound for E_w^0 in terms of the complete H-polynomial, rather than a ratio bound. Orient the edges of the simple vertex figure by a generic linear function representing its cutting hyperplane. A vertex of index j has j downward edges and d−j upward edges. Assign each uncut edge below the section to its upper endpoint and each uncut edge above the section to its lower endpoint. Its contribution is j or d−j, respectively, and is at least min(j,d−j). The number of vertices of index j is the H-coefficient h_j, by the simple-polytope index formula (Khovanskii, Theorems2 and6). Consequently

` E_w^0 >= G(L_w)=sum_j h_j min(j,d−j). `

At a global minimum or maximum vertex of P, the section misses the entire vertex figure, so E_w^0=E_w. For any simple d-polytope, symmetry of H and the index formula give

` E_w−G(L_w)=sum_(j>d/2) h_j(2j−d) >= d, `

because h_d=1 and all H-coefficients are nonnegative. The distinct global minimum and maximum thus supply an additional total at least2d in the uncut-edge sum, or5d in D5. This proves the quantitative slack in the displayed theorem. In particular it also proves the strict version without tracking the additive constant.

For a finite vertex, L_w is a d-simplex: H=1+t+...+t^d, G=floor(d^2/4), and E=d(d+1)/2. Substitution gives exactly beta_n. At an ideal vertex, the product formula for face polynomials and the substitution H(t)=F(t−1) give the displayed product H_v; E=dV_v/2. This proves the entire inequality in every dimension. No finite partition census is used.

## A probabilistic expression and a usable relaxation

Let independent X_i be uniformly distributed on the integers0,...,r_i. Since their sum has probability generating polynomial H(t)/V,

` kappa(L)/V = 3d/4 − (5/2) E|sum_i X_i−d/2|. `

Thus the complete link shape enters through a mean absolute deviation, with exact rational value. Its variance is

` Var(sum_i X_i)=sum_i r_i(r_i+2)/12. `

Cauchy–Schwarz gives the convenient rigorous relaxation

` kappa(L)/V >= 3d/4 − (5/2)sqrt(sum_i r_i(r_i+2)/12). `

The exact H-polynomial value is generally stronger. A product of many small factors can therefore have substantially more curvature available than the universal beta_n V/n bound. This is a combinatorial improvement; the corner-payment theorem must still be established for the proposed geometric weights.

## An exact obstruction to one coarse local comparison

Consider the legitimate Euclidean cusp-link type

` Delta^30 × Delta^30 × Delta^30 × Delta^30 × Delta^30 × Delta^13. `

Its link dimension is163, so it is a local cusp model for H164. Choose each irreducible affine diagram to be an A cycle. On the five31-node cycles and the14-node cycle, every node pair has distance at most15. Thus the distance15 rule charges, per link vertex,

` sum_i r_i/2 = 163/2. `

The exact curvature margin is

` kappa(L)/V − 163/2 = −69328961/200404057 < 0. `

Therefore a proposed proof that pays every cusp separately by this distance15 budget cannot be justified by the exact link curvature alone. The global additive slack is positive, but is not a substitute for the failed per-cusp coefficient inequality when cusp counts are unrestricted. This is a scoped obstruction to that payment comparison, not a counterexample to a dimension theorem: no finite-volume H164 completion of this local cusp model is asserted. It is compatible with the independently proved unrestricted bound163.

The smaller distance7 budget is available only under the separate rank restriction of the factor framework. Fine-weight extensions could lower the large-factor cost further, and are a geometric obligation separate from this exact counting lemma.

## Exact checks and status

The [checker](../scripts/fresh_cusp_link_shape_curvature_20261017.py) and [data](fresh_cusp_link_shape_curvature_20261017.json) verify179 simplex-to-beta identities,17 product-link polynomials, symmetry and variance identities, the extremal slack, and the exact negative margin and cycle pair counts. The all-rank proof is the local index argument above; the finite lists corroborate its formulas. Generation and fresh saved replay both terminate with PASS. Syntax, all five direct source hashes and the local report links pass. The root independently read the all-profile index proof, including the quantitative global-extremum slack, and found no gap; see the [root review](fresh_cusp_link_shape_curvature_20261017_root_review.md). That review does not duplicate the finite arithmetic checks. This source, data and report are final/frozen. No optimizer or new dimension claim is part of this package.
