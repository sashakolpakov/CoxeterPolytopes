# Independent review of the escape-cutting wall3 cap

The [separate checker](../scripts/fresh_nonarithmetic_wall3_cap_20261016_review.py)
verifies the new [cap construction](fresh_nonarithmetic_wall3_cap_20261016.md)
without importing its producer. The resulting chamber has **34 actual
facets, 95 finite vertices and 257 finite edges**. The cap cuts the old
escape, but a different exact strict escape proves that the chamber still
has infinite volume.

The reviewer uses its own formal arithmetic for `a+b sqrt(1+sqrt5)` and
reverse-order Gaussian elimination. It shares the frozen octic field and
rigorous base-field sign evaluator with the earlier bridge review. No claim
of another independent implementation of that base field is made.

The rank-nine contact system and the particular/kernel norm identities are
checked against the literal old form. They verify the new norm-two normal.
All 33 contacts pass: eight right, one of label three and 24 strictly dotted.
Acceptance therefore needs no conjectural finite-angle alphabet. Direct
center projection verifies all 34 actual facets and their strict signs.

All 92 inherited finite vertices remain strict at the new wall. The cap has
only nine old finite neighbors, so its finite vertices must occur among the
nine possible eight-neighbor supports. A reverse-order positive-definiteness
test finds exactly three spherical supports. Independently reconstructed
physical vertices have precisely their stated nine-wall zero sets. Together
with the frozen complete 33-wall census, this establishes the complete new
finite-vertex count. The incidence calculation gives 257 finite edges and
341 open finite-corner directions.

Finally, a separate linear solve reconstructs the specified old hyperideal
ray. Adding `H/8192` gives a spacelike point with negative center pairing and
strictly negative pairings with all 34 walls. This certifies infinite volume
without requiring a complete new exterior-ray census.

```sh
python3 scripts/fresh_nonarithmetic_wall3_cap_20261016_review.py
python3 scripts/fresh_nonarithmetic_wall3_cap_20261016_review.py --verify
```

The [saved certificate](fresh_nonarithmetic_wall3_cap_20261016_review.json)
binds the reviewed source and data. Generation has passed; the fresh saved
replay also passes.
