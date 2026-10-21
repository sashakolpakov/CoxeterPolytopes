# All-subset integer packing still leaves an exact obstruction

Testing every subset of the 13 ports, for every nonzero divisor 2, 3 and 4,
strengthens the [optimized partial-model witness](fresh_sqrt5_labelled_packing_full_model_20261017.md)
but does not defeat it. The exact repaired lower bound is

\[
 \frac{2790968068837079040}{217727747504296963}
 =12.818614535026125\ldots>12.
\]

This is a fixed-witness lower bound for an enlarged partial weighting model.
It is not its optimum, a dimension theorem, or a realizable polytope incidence
distribution. No optimizer is used. The new family gives zero payment to
quadrilateral rows.

## Valid charges and exact symmetry averaging

At a finite vertex of a compact Coxeter polytope, distinct maximal simplex
faces have disjoint port blocks C. For any port subset U and integer t≥2,

\[
 \sum_C\left\lfloor |C\cap U|/t\right\rfloor
 \le\left\lfloor |U|/t\right\rfloor.
\]

This is the previously proved [integer packing inequality](fresh_sqrt5_integer_simplex_port_packing_20261017.md),
now tested on every subset, rather than only triples and the full port set.
Compact simplex corners have at most four ports, so t>4 pays nothing.

Let G be the automorphism group of the complete labelled finite vertex
forest. Average a charge over G by replacing its payment on C with

\[
 \frac1{|G|}\sum_{g\in G}\left\lfloor |gC\cap U|/t\right\rfloor.
\]

For each g, all the blocks remain disjoint, so the same vertex capacity
floor(|U|/t) remains valid. Averaging each block uniformly over its orbit
is equivalent by orbit-stabilizer and linearity. This requires no symmetry
of the actual ambient chamber. Consequently every literal U gives a valid
invariant charge, with generally rational face coefficients and integer
vertex capacity.

The saved face masses weight these averaged blocks. For k=|C∩U|≤4,

\[
 \begin{aligned}
 \lfloor k/2\rfloor&=\binom{k}{2}-2\binom{k}{3}+4\binom{k}{4},\\
 \lfloor k/3\rfloor&=\binom{k}{3}-3\binom{k}{4},\\
 \lfloor k/4\rfloor&=\binom{k}{4}.
 \end{aligned}
\]

The source accumulates weighted literal pairs, triples and quadruples, then
uses subset zeta transforms to sum the contained features for all 8,192 U.
All arithmetic uses rational numbers or integers with a common denominator.
The subtraction in these formulas does not make any actual payment negative.

## Complete test on the saved positive simplex support

The [source](../scripts/fresh_sqrt5_all_subset_port_packing_20261017.py)
first rebuilds and verifies the complete saved model and its dual. It then
reconstructs every positive simplex star, using its complete finite
one-deletion vertices and full corner port blocks. Each reconstructed star
must reproduce every inherited ordinary-pair and one-port coefficient.

There are 3,574 positive simplex rows: 451 four-simplexes, 556 tetrahedra,
792 triangles and 1,775 maximal edges. Their 10,405 corners use 6,011 block
orbits at 740 vertex forests. The remaining 18 positive vertex forests have
no positive simplex mass and receive zero load from this family.

For those 740 forests the exact test checks 17,828,080 literal inequalities
with positive capacity. It finds 2,424 violations at 102 forests. Maximizing
subsets include sizes 5, 7 and 9 for t=2; sizes 4, 5, 7, 8 and 11 for t=3;
and sizes 7, 10 and 11 for t=4. Some t=2 violations occur at size 3 because
the new test also reconstructs the formerly unexported positive simplex rows.

For each forest its vertex mass is increased to the largest required ratio
over all U and all three divisors, or retained if already larger. Every face
mass and every actual-provenance mass stays fixed. Thus all inherited signed
loads are unchanged, and every inherited positive-cost inequality remains
valid under the monotone increase in vertex masses. The exact numerator stays
unchanged; the enlarged normalizer gives the bound above.

## Checks, exported data and limits

Every calculated payment is checked nonnegative. Independently of the zeta
placement computation, for each subset size s the sum of all loads is checked
against the elementary identity

\[
 \sum_{|U|=s}\lfloor |C\cap U|/t\rfloor
 =\sum_j\binom{|C|}{j}\binom{13-|C|}{s-j}\lfloor j/t\rfloor.
\]

Every maximizing load is also recomputed directly from literal block
intersections. These are additional arithmetic consistency checks, not an
independently implemented replay of the whole program.

The [saved certificate](fresh_sqrt5_all_subset_port_packing_20261017.json)
contains the full positive simplex block features, literal block orbits,
all per-forest/divisor maxima, the repaired masses and 148 maximizing
violated charges with exact rational coefficients on the exported simplex
classes. Those 148 charges suffice to reproduce this witness repair; they
are not the entire family for a subsequent optimization. A later restricted
primal model may give them zero payment on unexported face classes, but must
state that restriction explicitly.

All subset integer capacities are still only a family of necessary conditions
for joint disjoint-block packing. They do not describe the entire packing
polytope for blocks of size three or four. That stronger local problem is
a separate next step. Full global face and field coverage also remain open.

Generation 97130 and fresh saved replay 28367 both terminated with exit 0
and `PASS_EXACT_ALL_SUBSET_INTEGER_PACKING_SEPARATION`. Syntax, direct source
hashes and local links passed. This is one arithmetic implementation with
the additional checks described above; no independent complete arithmetic
replay or external review is claimed. **FINAL/FROZEN.**

```sh
python3 scripts/fresh_sqrt5_all_subset_port_packing_20261017.py
python3 scripts/fresh_sqrt5_all_subset_port_packing_20261017.py --verify
```
