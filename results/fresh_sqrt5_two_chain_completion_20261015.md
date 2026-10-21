# Two disjoint white-chain completions over Q(sqrt5)

All positive integer chain lengths are resolved for the stated two-path construction. The 130 previously open nonsingular, finite, unpaid cores have 458 accepted first-node contact profiles. Pairing these profiles with replacement gives 1,455 unordered trials. Exactly 104 profile pairs give 112 surviving principal-star records, on 59 source cores. These are records with literal labels, not a count of isomorphism classes or new polytopes.

The [source](../scripts/fresh_sqrt5_two_chain_completion_20261015.py) and [exact certificate](fresh_sqrt5_two_chain_completion_20261015.json) are separate from the frozen [single-chain package](fresh_sqrt5_single_chain_completion_20261015.md). Each of the two external components is a simply-laced path; only its first node touches the core, and the two paths have no contacts with each other. Additional components orthogonal to the entire configuration do not change its rank condition. Branching, a component with several nodes touching the core, contacts between the two components, dotted cores and arbitrary further completion are outside this result.

| Exact outcome | Profile pairs |
|---|---:|
| Nonfinite corner already at the two first nodes | 880 |
| Two negative directions already at the first nodes | 283 |
| Conjugate matrix not positive semidefinite | 53 |
| No positive integer length pair | 135 |
| Finite unpaid principal-star completions | 104 |
| Total | 1,455 |

All 112 completed matrices have four positive-definite corner diagrams, physical inertia `(N−2,1,1)` and conjugate inertia `(N−1,0,1)`. They remain unpaid by the distance-two and entire-H4 tests. The inherited crystallographic affine-witness search finds no witness; failure of that search is not asserted to prove the absence of every affine principal subdiagram. Actual future-oriented compact quadrilateral cones and larger polytopes are not certified here.

## Exact all-length equation

Let `H` be the nonsingular core Gram matrix, with first-node contact columns `v,w` and scaled squared norms `d,e`. Put

\[
c_{11}=v^TH^{-1}v,\quad c_{22}=w^TH^{-1}w,\quad c_{12}=v^TH^{-1}w,
\]
\[
a=2c_{11}/d,\qquad b=2c_{22}/e,\qquad
f=4(c_{11}c_{22}-c_{12}^2)/(de).
\]

The chain blocks are `(d/2) Cartan(A_l)` and `(e/2) Cartan(A_m)`. Their first inverse entries are `2l/[d(l+1)]` and `2m/[e(m+1)]`. The determinant lemma therefore gives

\[
\frac{\det G}{\det H\det A_l\det A_m}
=\frac{P(l,m)}{(l+1)(m+1)},\qquad
P(l,m)=(1-a-b+f)lm+(1-a)l+(1-b)m+1.
\]

A full quadrilateral star must be singular, so `P(l,m)=0` is necessary. The equality is in Q(phi), hence enforces both real embeddings simultaneously. The certificate recomputes every core inverse and contact cost, and checks the cleared two-by-two Schur determinant identity at four parameter values for each of the 239 pairs reaching this stage.

Write `P=A lm+B l+C m+1`. When a coefficient is nonrational, its two rational coordinate equations eliminate `m` to the exact quadratic

\[
(A_0B_1-A_1B_0)l^2+
(C_0B_1-A_1-C_1B_0)l-C_1=0.
\]

A nonzero polynomial has only the positive integer roots obtained by the exact rational discriminant test. Each is substituted into both coordinate equations. If the elimination polynomial vanishes identically but the coefficients are not all rational, the original polynomial factors as `(B l+1)(C_0 m+1)`; the nonrational first factor cannot vanish at an integer `l`. This identifies any free-length family explicitly.

When all coefficients are rational, clearing denominators gives integer coefficients `A,B,C,D`. For `A≠0`,

\[
(Al+C)(Am+B)=BC-AD.
\]

Every divisor of either sign is tested. Zero right-hand side gives fixed-length families rather than a finite list. For `A=0`, the linear Diophantine equation is solved by an exact gcd parameterization and positivity bounds. Thus neither the finite solution lists nor the potential infinite families depend on a length cutoff. The implemented branch checks compare all solver branches with 10,800 direct finite evaluations; those checks supplement the algebraic completeness argument.

The 239 equations split into 182 nonrational elimination cases, six nonrational factorization cases and 51 rational bilinear cases. They yield 107 isolated finite solutions and five infinite rank families.

## The five rank families are bounded by actual corners

Each infinite family fixes one chain length to one. At free length two, four families have an affine actual corner and the fifth has an indefinite actual corner. The certificate saves every corner matrix, its exact inertia, and each affine kernel. Every longer chain contains the same failing corner as a principal subdiagram, with the other chain's nodes relabeled as necessary. Therefore all free lengths at least two are excluded, without an upper search cutoff. The remaining length-one member of each family passes the completed-star tests.

| Source class | Source profile indices | Free chain | Failure at length two |
|---|---|---|---|
| 564 | 2, 7 | first | affine corner |
| 564 | 2, 12 | first | affine corner |
| 564 | 3, 6 | second | affine corner |
| 564 | 6, 12 | first | affine corner |
| 674 | 2, 30 | first | indefinite corner |

There are no unresolved integer-length families in this construction. The final 112 records have chain lengths at most five as a consequence of the exact equations and corner bounds. Of these, 96 have lengths `(1,1)`; the other 16 records and their full diagrams are saved individually. The stars have 8–14 nodes. No claim is made that further external attachments are impossible.

Generation and saved replay pass. Frozen source manifests bind the old contact census, its independent review, and both single-chain source and independent review. The next task is containment in already sealed configurations, followed by complete new-white contact tests for genuinely unmatched stars.

```sh
python3 scripts/fresh_sqrt5_two_chain_completion_20261015.py --verify
```
