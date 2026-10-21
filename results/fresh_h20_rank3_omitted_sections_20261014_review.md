# Independent review of the rank-three local section census

The local classification passes. No mathematical issue was found in the conic reduction, removal of zero-coordinate solutions, algebraic root selection, or shift interval. This is a census of **local affine-triangle directions**, not a global Coxeter-section construction or exclusion.

The [separate checker](../scripts/fresh_h20_rank3_omitted_sections_20261014_review.py) imports no producer code. It reconstructs the rational source Gram matrices and checks the [saved census](fresh_h20_rank3_omitted_sections_20261014.json) through a different elimination: it adds an inverse variable `z*x*y=1` to the original conics. The [review certificate](fresh_h20_rank3_omitted_sections_20261014_review.json) binds all inputs and both scripts. Exact generation and saved replay pass.

| Check | Result |
|---|---:|
| Literal rank-three affine factors | 3 |
| Omitted-facet/ordered-target systems | 120 |
| Empty saturated systems | 17 |
| Positive directions | 91 |
| A2 / C2 / G2 directions | 16 / 19 / 56 |
| Rational / quadratic / quartic branches | 33 / 14 / 44 |

The independent checker uses SymPy's rational Groebner bases and real-root isolation, as does the producer. It is a separate reconstruction and formulation, not an independent computer-algebra implementation or formal proof-assistant proof.

## Completeness of the geometric equations

A bounded Euclidean Coxeter triangle has angle labels `(3,3,3)`, `(2,4,4)`, or `(2,3,6)`. These also follow directly from `1/a+1/b+1/c=1`; the standard affine classification is recorded in [Felikson's lecture notes, section 1.2](https://www.maths.dur.ac.uk/users/anna.felikson/talks/cl_cox.pdf). The three target vector triples encode exactly these shapes, with nonpositive pairings and a strictly positive linear dependence. Permuting the roots and identifying equal ordered squared-cosine signatures gives ten target placements. Because all nonzero target pairings are negative, squared cosines do not lose an orientation branch here.

For retained positive-definite Gram `B` and `J=B^-1`, every oriented target Gram can be written `T=XX^T+dYY^T`, where each row has an independent positive length. Positive common scaling permits setting the first length to one, leaving `x,y>0`. A rank-one orthogonal projection of `B` onto this triangle requires the two nonzero eigenvalues of `B^-1/2 T B^-1/2` to agree. This is exactly

```text
X^T J Y = 0,
X^T J X = d Y^T J Y = q > 0.
```

Conversely these equations imply that `qB-T` is positive semidefinite of rank one. The saved normal construction and the identity `norm*(qB-T)=q*a*a^T` are checked for every accepted branch. The target triangle's positive kernel forces `B^-1*a` to have one common strict sign; the code chooses the positive orientation and verifies every coefficient. Thus no extra positive branch is filtered by an unproved geometric heuristic.

## Exact elimination and saturation

Both producer resultants are nonzero. Every strictly positive conic solution therefore lies among their finitely many coordinate roots. Taking squarefree parts preserves their zero sets. Dividing out a zero-coordinate factor discards only `x=0` or `y=0`, outside the intended domain. Adding those resulting polynomials to the conics preserves every positive solution and cannot create a spurious conic solution.

The successful lexicographic Groebner basis has the form `x-X(t), y-Y(t), h(t)` with `t=x+c*y`. This is an exact ideal description; the bound on attempted primitive-element coefficients is not a coverage cutoff because failure would raise an error. All factor roots of the squarefree `h` are isolated over the reals, and exact interval signs select every `x,y>0` branch.

The review independently starts from `[f,g,z*x*y-1,t-x-c*y]`, without using either producer resultant. Some original intersections are tangent and hence nonreduced: for example, system 2 has eliminant `(t-3)^2*(t-11/6)`. In shape position, squarefreeing this univariate polynomial gives the radical ideal after reducing the coordinate polynomials. For every nonempty system, this direct saturated radical agrees with the saved representation. All 91 positive roots are matched bijectively to the archived root intervals; no branch is missing or duplicated within a system.

## Normal and parameter checks

For every direction, the review reconstructs the full 22-coordinate normal, its positive square, the affine null shift, orthogonality to that shift, the projected Gram, and the strict positive kernel. It also verifies the minimum defining `L=min(k_i/m_i)` and its designated witness wall.

With the horosphere normalized by `(D,p)=-1`, put `y_i=-(r_i,p)`. The retained simplex satisfies `y_i>=0` and `sum(k_i*y_i)=-t`; the omitted inequality is `sum(m_i*y_i)<=1`. Maximizing that last expression at the simplex vertices gives precisely `-L<=t<0`. At `t=-L`, an omitted inequality can coincide with a lower-dimensional face and requires separate global endpoint treatment; `t=0` does not give a full-dimensional local section. The census correctly leaves those global questions open.

```sh
python scripts/fresh_h20_rank3_omitted_sections_20261014_review.py --verify
```

The three literal source factors and all four choices of omission are checked against the complete saved cusp census. The 91 count is of accepted omission/target branches; it does not assert 91 globally inequivalent planes, nor that any branch yields a finite-volume Coxeter H20 section.
