# Independent review of all2,600 chain triples

**PASS:** any coherent Coxeter cap cutting the fixed H20 weak escape E0 has
at most two nondotted contacts among the26 retained chain walls, and thus
at least24 dotted contacts. This condition applies to an E0-cutting wall;
it is not a statement that every added wall has24 dotted contacts.

The [independent source](../scripts/fresh_h20_arbitrary_triple_escape_gate_20261016_review.py)
and [certificate](fresh_h20_arbitrary_triple_escape_gate_20261016_review.json)
import no producer arithmetic. They reconstruct all300 positive gap pairs
and all2,600 actual triples using SymPy Chebyshev polynomials, polynomial
minors, rational matrix kernels and separate derivative-based enclosures.
They agree with the [primary theorem](fresh_h20_arbitrary_triple_escape_gate_20261016.md).
The physical retained chamber, degree50 after abelian base change, weak E0
feasibility and support-budget theorem are explicitly inherited from the
[completed independent asymmetric review](fresh_h20_asymmetric_arc_escape_gate_20261016_review.md).

## All coefficient identities are reconstructed

For every separation1 through25 the reviewer divides `2*(T_k(C)-C)` by
`1-C` with an independent polynomial library and verifies exact division.
For a triple with gaps p,q it reconstructs all adjugate entries by2x2
minors and the determinant by expansion. The seven coefficient columns of

`q_y*det(G)-f^T*adj(G)*f`

are reconstructed in full, including every zero coefficient. Their maximum
degree is48, below the inherited lower bound50 for any nonzero polynomial
over Q^ab vanishing at C. All300 independently computed kernels and pivots
agree with the primary. Their translation multiplicities sum to2,600.

The checker also verifies the reduction from formal kernel variables to
actual contact squares and products. For242 cases the coefficient matrix
has full rank. In each of the57 remaining nonconsecutive cases, it verifies
that a specified sum of two squared contacts lies in the matrix row space
by pairing that row with every vector in the independently computed kernel.
The two nonnegative squares therefore vanish. Their actual cross products
vanish as well, and the remaining two-column matrix has rank2, forcing
the final contact and q_y to vanish. This separately checks the nonlinear
step instead of relying only on the names of kernel classes.

For the sole consecutive gap pair(1,1), the outer contacts vanish and the
remaining two-column kernel gives `q_y=f_middle^2/2`. Its projection is
`y=-f_middle*n_middle/2`, which has positive pairing with an available
retained wall at separation2 if the middle contact is nonzero. The previously
reviewed24-position check covers every such location. When all three contacts
are zero, their independent normals force y=0; the cap is then a negative
nonnegative combination of support normals and cannot share a strict
interior with the retained chamber. Thus every coherent cap with three
nondotted chain contacts has a dotted support contact.

## Separate sign bounds and the outer-pair estimate

The reviewer uses a rational midpoint and the mean-value theorem to bound
each polynomial: on `0<C<1`, the derivative is bounded by the sum of
`k*abs(coefficient_k)`. This differs from the producer's interval Horner
calculation. It independently proves every triple determinant negative and

`g25 < ... < g3 < g2 < -2`.

Consequently any two nonconsecutive chain normals span a Lorentzian plane.
The reviewer symbolically inverts its2x2 Gram matrix and obtains

`u=-y^2 <= (2f^2+2h^2-2g*f*h)/(g^2-4)`.

The orthogonal complement of this plane in the Lorentz3 chain space is
positive, justifying the inequality. The two partial derivatives of the
right side are nonnegative for `g<-2` and `0<=f,h<=2`. Its maximum is
`8/(-g-2)` at f=h=2. Gap monotonicity then gives `u<=2/C`.

The dotted support contact and minimum inverse diagonal1/2 force u>0.
The outer-root sum is future timelike, and its pairing with y is
`-(f+h)<=0`. A timelike y cannot be orthogonal to a timelike vector, so this
pairing is strictly negative and y is future. This also handles the apparent
zero outer-contact endpoint without an unproved orientation assumption.

## Applying the independently reviewed common escape

The bound `u<=2/C` gives the same total support budget below281/70 as in the
asymmetric consecutive theorem. Every one of the new primary's eight finite
profile costs, escape pairings and margins is compared with the independent
asymmetric review. Both continuous dotted-endpoint inequalities also agree.
That previously complete support reduction proves `r.E0<0` for any cap
having at least three nondotted chain contacts, regardless of which triple
it uses. Its finite-collection perturbation proof supplies a strict spacelike
escape and infinite volume.

Generation and a fresh saved replay pass. The new arbitrary-gap computation
does not repeat the old physical chamber or modular factorization; those
remain pinned independent dependencies. No finite label cutoff, arithmeticity
assumption on the proposed cap, or numerical proposal bank enters this proof.
The result leaves open caps with at most two nondotted chain contacts and
constructions changing retained walls. It gives no new dimension bound and
no finite-volume example.

```sh
python3 scripts/fresh_h20_arbitrary_triple_escape_gate_20261016_review.py --verify
```
