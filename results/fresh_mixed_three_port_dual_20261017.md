# All mixed three-port charges still leave a budget obstruction

The enlarged rank-13 model over `Q(sqrt(5))` has an exact dual lower bound

\[
\frac{46073614230483840}{3603636833546389}
=12.785310051663\ldots>12
\]

even after adding **every mixed three-port indicator charge** from the
[shared face-incidence lemma](fresh_mixed_face_port_incidence_20261017.md).
Thus these charges, together with the inherited pair, simplex packing and
disjoint-port charges, cannot supply the target budget12 on this saved domain.
This is a limitation of this particular family of weighting rules. It is not
a dimension bound or an existence theorem.

## Exact certificate

The starting point is the exact dual in the
[combined disjoint-port model](fresh_sqrt5_rank13_combined_disjoint_ports_20261017.md).
Its2,397 positive face multipliers comprise815 simplex rows and1,582 quad
rows. Their corners occupy620 forest types. For each of these forests, the
new checker enumerates all286 literal triples of13 ports and groups them
by the complete marked-forest isomorphism code:47,725 orbit columns in total.
The cost of a column is the size of its literal orbit.

At each actual aggregate vertex, a triple contributes1 precisely when the
aggregate occupies at least two of its three port pairs. Different aggregates
have disjoint occupied pairs, so at most one can contribute. Summing over
an orbit gives the claimed vertex cost. Single triangular and quadrilateral
corners alone do not activate this indicator.

The checker reconstructs the actual footprints of every positive face row.
For quad rows, it checks all1,858 source provenances independently and
reproduces their full inherited pair coefficients. It initially takes the
coefficient-wise minimum over provenances, then verifies the stronger fact
that **every positive row has identical complete mixed coefficients across
all its provenances**. Therefore the certificate also applies when each
full star has its own row; a positive merged multiplier can be assigned to
one of its provenances. No first-provenance equality is assumed.

The old face multipliers violate only three new column constraints, at
forests

- `7 A1 + A2 + A4`,
- `3 A1 + A2 + E6 + I2(5)`,
- `3 A1 + A2 + 4 I2(5)`.

Increasing the corresponding vertex multipliers by a total of
`142693/14000000000` repairs all new constraints. The checker independently
reconstructs the old column loads, verifies the old vertex multipliers are
their exact maxima, and checks every new column against the repaired
multipliers. All other forests have zero positive face load. The unchanged
face objective divided by the repaired total vertex mass gives the lower
bound above, compared with the previous `12.7854403641055`.

## Scope and reproduction

This is an exact rational dual certificate, obtained without a new optimizer
run. It covers the saved combined domain and all three-port charges of the
stated kind. Other global incidence constraints, larger port subsets,
different charging rules and missing face families remain separate avenues.
In particular, this does not establish coverage of all compact arithmetic
chambers, all fields, or all higher dimensions.

The [certificate](fresh_mixed_three_port_dual_20261017.json) binds its producer,
the incidence lemma and the source models. Recompute it with:

```sh
python3 scripts/fresh_mixed_three_port_dual_20261017.py --verify
```

Generation and a separate full replay both returned
`PASS_EXACT_ALL_THREE_PORT_MIXED_CHARGE_OBSTRUCTION`. This package is frozen.
