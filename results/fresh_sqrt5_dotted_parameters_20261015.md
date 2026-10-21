# Complete arithmetic distance domains for the remaining dotted cores

The [exact checker](../scripts/fresh_sqrt5_dotted_parameters_20261015.py)
reduces all 24 unpaid minimal dotted families to **32 discrete parameter
cases**: eighteen singular cores and fourteen nonsingular cores. Its
[certificate](fresh_sqrt5_dotted_parameters_20261015.json) exhausts 532 field
integers in explicit bounded embedding boxes. This parameter reduction
allows arbitrary subsequent attachments; it is not a distance cutoff.

Let the scaled Gram diagonals of the dotted endpoints be d0,d2, and their
inner product -x. In every seed, the common principal matrix B obtained by
deleting these endpoints is nonsingular with one negative direction. The
two-by-two Schur residual has diagonals a,b and off-diagonal -x-c. Any
physical extension with only one negative direction requires this residual
to be positive semidefinite. Consequently (x+c)^2≤ab, and

    t = 4x²/(d0d2) ≤ 8(ab+c²)/(d0d2).

The displayed inequality follows from x²≤2(x+c)²+2c². Exact rounding of
the right side gives upper bounds between 9 and 76. Since the walls are
disjoint, t>4. Arithmeticity requires the two-edge cyclic product t to be
an algebraic integer in Q(sqrt(5)). Conjugate positive semidefiniteness
requires 0≤t'≤4.

Write t=A+Bφ, with φ=(1+sqrt(5))/2 and A,B integers. The embedding bounds
imply B>0, A≥0, and both coefficients are at most the integer upper bound.
The checker enumerates that finite rectangle and tests both embedding
inequalities exactly. The endpoints t'=0,4 cannot actually occur with
t>4, because conjugation is injective and fixes rational numbers.

Of the 532 candidates, 365 fail the dotted cyclic-field condition, 107
force two physical negative directions, and 28 fail conjugate positive
semidefiniteness. All remaining 32 pass every simple-cycle integrality
test. They are necessary principal cores, not asserted realized faces.
The two minimal classes 6 and 7 have no surviving arithmetic distance.

Generation and saved replay pass. Independent review of this new reduction
is pending. Attachment closure and compatible face weights are still needed
before it can contribute to an improved dimension theorem.

```sh
python3 scripts/fresh_sqrt5_dotted_parameters_20261015.py --verify
```
