# Exact rank completions with one dotted edge

The [checker](../scripts/fresh_sqrt5_dotted_completion_20261015.py) and
[certificate](fresh_sqrt5_dotted_completion_20261015.json) examine all 24
unpaid dotted minimal cores and their 369 previously enumerated first-white
contacts. They find **27 necessary arithmetic completions** when every further
white root is orthogonal to the displayed core. These are principal Gram
matrices, not asserted compact polytopes or a dimension theorem.

The unknown entry is the negative inner product of the opposite marked roots
0 and 2. Deleting either endpoint removes it. Exact physical and conjugate
inertia tests on these principal matrices exclude five contact cases;
57 more fail the finite-edge cyclic ground-field test. These exclusions
remain valid with arbitrary additional roots.

For every survivor, the common submatrix B obtained by deleting both
endpoints is nonsingular and has one negative direction. Write its two
Schur diagonal entries as a and b, and the cross projection as c. With dotted
entry -x, the residual determinant is ab-(x+c)^2. If no further white root
touches the core, the quadrilateral-star rank condition therefore forces
x=-c±sqrt(ab). Orthogonal positive-definite spectator blocks do not change
this condition.

If the finite-edge graph connects the dotted endpoints, its cyclic products
force x to lie in Q(sqrt(5)). The exact quadratic-field square-root test
excludes 87 more zero-contact cases. Otherwise c=0 and x^2=ab lies in the
ground field. Scaling the entire finite component of root 0 by sqrt(ab)
then gives a ground-field Gram matrix, with totally positive scaling square
whenever the conjugate completion is admissible.

The remaining 244 cases have all rank roots explicitly enumerated. Forty-one
roots satisfy the strict dotted inequality and the physical rank condition.
All 41 have positive-semidefinite conjugate matrices; fourteen fail an
explicit integral cyclic-product condition, leaving 27. The certificate
saves each full Gram matrix and every simple-cycle product, including
two-edge cycles. General closed walks factor into such cycles.

Finite-corner and short-path domains are inherited from the reviewed input
census. The current package checks rank, field, inertia and integrality;
it does not establish that every passing matrix realizes the claimed
quadrilateral face. It also leaves arbitrary additional touching white
roots open. In particular, a core with no zero-contact completion may still
have a completion after further attachments.

Generation and saved replay pass. The primary computation reuses reviewed
field arithmetic; independent review of this missing-entry reduction is
still pending.

```sh
python3 scripts/fresh_sqrt5_dotted_completion_20261015.py --verify
```
