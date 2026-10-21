# Independent review of the bounded first-cap exclusion

The [separate checker](../scripts/fresh_nonarithmetic_two_caps_20261016_review.py)
reconstructs the completeness argument and all sixteen branches in the
[primary report](fresh_nonarithmetic_two_caps_20261016.md), without importing
its producer. Nine branches are obtuse to an old wall; seven have finite
angles strictly between consecutive Coxeter angles. Thus all 128 specified
ordered two-cap recipes fail at their first cap.

For each of the two anchors, direct Gaussian elimination verifies that the
nine wall equations have rank nine. The serialized particular vector and
kernel satisfy the literal form equations, are orthogonal, and have negative
and positive square respectively. These facts prove that the two formal
radical vectors are the complete real norm-two solution set for each of the
four prescribed anchor labels. No singular or continuous case is omitted.

A separately written radical comparison routine checks all 528 old-wall
pairing signs, the center and escape signs, and the seven forbidden angle
intervals. The exceptional label-seven/eight gap uses an exact rational
bracket and the cubic for the label-seven Cartan cycle. This excludes every
integer Coxeter label at the offending pairing, without assuming that its
angle belongs to the original field.

This implementation shares the frozen octic arithmetic, decoder and rigorous
field-sign evaluator of the earlier bridge review. It independently checks
the branch equations, rank, completeness and radical comparisons; it is not
a second independent implementation of the underlying number field.
The result concerns precisely the documented anchors, labels and right
contacts. Other labels, supports and larger cap families remain open.

Generation and fresh saved replay are run with:

```sh
python3 scripts/fresh_nonarithmetic_two_caps_20261016_review.py
python3 scripts/fresh_nonarithmetic_two_caps_20261016_review.py --verify
```

The [saved certificate](fresh_nonarithmetic_two_caps_20261016_review.json)
binds the reviewed data, source code and primary report by their hashes.
