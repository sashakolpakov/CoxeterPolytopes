# Independent cover check at the new restricted face masses

**FINAL/FROZEN.** The [standalone Fraction checker](../scripts/review_sqrt5_restricted_flag_local_cover_20261017.py)
and its [saved certificate](review_sqrt5_restricted_flag_local_cover_20261017.json)
pass generation and fresh replay. They import no producer module, numerical
optimizer or packing enumerator.

The input masses are matched exactly to the frozen `legacy_export_seed`
restriction. All1,864 positive simplex rows are covered, as checked against
the input support's face-kind records. Each row contributes exactly one
corner per simplex vertex. The weighted demand is independently recounted
from the input face obligations, and all3,485 block-orbit marginal demands
are reconstructed from those unchanged masses.

For every one of the593 contributing forests, the checker verifies the
membership and disjointness of every literal block in every saved rational
packing mixture. There are2,830 positive mixture terms and12,278 checked
literal blocks. All marginal cover inequalities hold exactly. The final
vertex-mass vector and its sum are recomputed, giving

\[
 {267632410797552528000000000000000
  \over21640812474569839372569203000501}
 =12.367022315453632\ldots>12.
\]

This proves the upper-cover part of the
[complete restricted-cone result](fresh_sqrt5_restricted_flag_local_cones_20261017.md).
If a nonnegative orbit charge q has disjoint-packing capacity c(q), an
explicit mixture of total mass L covering the required marginals bounds
its load by L*c(q). Raising each forest mass to at least L therefore pays
every such charge simultaneously. This argument needs neither completeness
of a count-vector enumeration nor equality with a local lower bound.

The root also read the producer's new feature reconstruction, literal
orbit enumeration and unchanged-alpha verification of all428,310 inherited
model columns. No source-level gap was identified in that scope. The
independent checker does not duplicate the spherical-forest orbit
classification, the524 local optimality claims or the full matrix replay.
The producer supplies the latter checks separately. The69 exact covering
brackets suffice for the universal-charge argument without claiming an
exact optimum in those local cases.

The result preserves every current flag transfer and uses the new face
masses throughout. It does not reuse a cover from the earlier740-forest
witness. The obstruction concerns simplex-block charges alone, in addition
to the already represented constraints. New joint simplex/prism geometry
may impose further restrictions. No common ambient realization or dimension
theorem follows from these abstract packing mixtures.

```sh
python3 scripts/review_sqrt5_restricted_flag_local_cover_20261017.py --verify
```
