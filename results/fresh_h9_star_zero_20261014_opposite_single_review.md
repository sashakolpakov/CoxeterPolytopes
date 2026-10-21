# Independent review of the opposite-single branch

**PASS.** The [independent checker](../scripts/fresh_h9_star_zero_20261014_opposite_single_review.py) reconstructs all 896 determinant polynomials and all 7,488 positive Bernstein coefficients in a separate Fraction implementation of Q(sqrt(2),sqrt(5)). It verifies the 64 remaining actual-vertex obstructions and all new literal domain witnesses. The [review certificate](fresh_h9_star_zero_20261014_opposite_single_review.json) binds the frozen [primary source](../scripts/fresh_h9_star_zero_20261014_opposite_single.py) and [data](fresh_h9_star_zero_20261014_opposite_single.json). Independent generation and saved replay pass.

The scope is the stated family on incidence case 0: the first three two-node blocks have both leaves' own contacts 3/3; the last two leaves have opposite single own contacts, each with finite label at least 4; **all four internal block labels are at least 4**. The last internal label 3 is outside this certificate. This is not a whole-incidence exclusion.

## Domain and literal incidence checks

The checker reads the frozen 378-vertex incidence directly. It verifies all eight relevant missing triples and the complete vertex/missing-face preservation under the last two-root swap. Thus the opposite orientation is justified by a literal automorphism, rather than a diagram resemblance.

The four reductions of an internal label at least 6 use the previously reviewed strong-contact or two-dense theorems. Their exact lower-bound hypotheses are checked against the present fixed contacts, with recursive manifest verification. The proofs of those inherited theorems are not independently repeated here. The surviving internal choices are 4 or 5 in each block.

All 24 prospective cross-block contacts have a literal containing vertex and a relaxed affine path (4,3,4). Its positive null vector `(1,sqrt(2),sqrt(2),1)` is verified directly. Perron monotonicity therefore forces every cross-block contact right for all larger finite labels as well.

Three literal eight-root supports give the foreign-contact budgets. Each last anchor retains its nonright own root with squared cosine at least 1/2; any foreign contact would add at least 1/2, contradicting positive definiteness. The third root retains own cost 1/4, so it may have at most one foreign contact, necessarily pi/3. The literal three-root path bounds the two anchors' mutual label by 3. Four separately reconstructed actual-vertex LDL failures bound each third-anchor label by 3. Unlisted meeting contacts are relaxed to right only in these obstruction tests; increasing their negative Gram magnitudes cannot restore positive definiteness.

For the last missing triple with one right own contact, internal label 4 forces the other own label at least 5, while internal label 5 forces it at least 4. These follow from the exact three-node determinants at the excluded smaller labels. The two own cosines therefore lie in closed rectangles `[cos(pi/5),1]^2` or `[cos(pi/4),1]^2`; the artificial endpoint 1 safely enlarges the finite-label domain.

## Independent arithmetic

The checker uses the basis `1,sqrt(2),sqrt(5),sqrt(10)`, Galois-conjugate inversion and exact signs by squaring. It imports no producer or project arithmetic. For each of the 16 block assignments, seven possible third-root foreign patterns and eight anchor-contact assignments, it builds the full eight-root inverse and the three-root Schur matrix.

Each variable occurs in one added root, so the determinant has degree at most two in that variable. The checker evaluates the Schur determinant by the six signed permutations on the 3-by-3 grid `{0,1,2}^2` and reconstructs its polynomial by tensor interpolation. This exactly matches all 896 saved coefficient tables.

For the 832 positive cases it separately evaluates the determinant at each rectangle's endpoint/midpoint grid. In one variable the degree-two Bernstein middle coefficient is `2*f(midpoint) - (f(left)+f(right))/2`; applying this in each variable recovers all nine Bernstein coefficients. All 7,488 agree with the primary and are strictly positive. This calculation does not use the primary's binomial power-to-Bernstein expansion.

For the other 64 cases, the checker reconstructs the complete saved actual-vertex Gram at the minimum own labels and repeats the exact LDL failure, including its pivot order and value. Perron monotonicity extends these failures throughout the corresponding rectangle. Coverage of the entire 16×7×8 Cartesian product is checked explicitly.

A compact H9 realization would have positive definite S8 and total Gram rank at most 10. Its three-root Schur complement consequently has rank at most 2 and determinant zero. Every case instead has a strictly positive determinant or an already impossible actual vertex. No numerical tolerance, finite high-label cutoff, or metric optimizer is used.
