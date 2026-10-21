# Separate arithmetic review of all 82 minimal dotted domains

The [reviewer](../scripts/fresh_sqrt5_global_b1_dotted_parameters_20261016_review.py)
and [saved results](fresh_sqrt5_global_b1_dotted_parameters_20261016_review.json)
reconstruct all 82 classes in the
[parameter census](fresh_sqrt5_global_b1_dotted_parameters_20261016.md), retaining
the same 377 necessary arithmetic seeds and checking 2,578 integral cyclic
products. The reviewer imports neither the parameter producer nor the new
attachment queue.

It uses the frozen independent tuple-Fraction arithmetic, a different minimal
radical gauge for the one-dot matrices, full principal inverses and inertias,
and elementary simple-cycle enumeration. A larger symmetric integer rectangle
is searched, including negative coefficients; its accepted parameter set is
compared with the producer's exact-floor strip. The shared mathematical
embedding-box argument remains an explicit premise.

For the two-dot class it separately checks the smallest parameter 3+phi,
reconstructs the product bound and all twelve possible parameter values,
and checks all 144 ordered pairs in a different diagonal gauge. It finds
the same seventeen surviving pairs. The field-valued scaling introduces no
condition that square roots of the two dotted invariants lie in the field.

All accepted matrices have the required physical and conjugate signatures
and integral cycle products. The singular matrices' relations have alternating
black signs. Four minimal one-dot classes fail conjugate principal restrictions.
The recursive source manifests are checked as part of the review.

Generation and a fresh replay of the saved certificate both pass.

This independently reviews the arithmetic distance domains. It does not
review the new arbitrary-attachment queue, evaluate a weighting rule, prove
global compact realizability, or establish a dimension bound.

```sh
python3 scripts/fresh_sqrt5_global_b1_dotted_parameters_20261016_review.py --verify
```
