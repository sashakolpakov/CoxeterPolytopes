# Independent review of the dimension-24 hub-context weights

The independent exact calculation gives maximum vertex weight
**267,527/5,000 = 53.5054 < 54**, agreeing with the producer. It verifies all
**2,516 face rows and 65,300 common spectator cases**, and every one of the
**157** possible pairs of entire A1 and I2(5) counts at a rank-24 vertex.
The [reviewer](../scripts/fresh_sqrt5_dimension24_i25_cross_20261015_review.py)
imports no optimizer or producer implementation. Its
[certificate](fresh_sqrt5_dimension24_i25_cross_20261015_review.json)
binds the [primary weights](fresh_sqrt5_dimension24_i25_cross_20261015_hub_context.json)
and their complete face inputs.

## What was checked

The separate graph engine reconstructs the literal labelled trees, recognizes
their intrinsic finite types, canonicalizes every marked pair and root, and
checks the full variable ranges. It includes the rank-two aliases A2/I2(3)
and B2/I2(4). All 2,397 original literal-input maps are rechecked, including
their exact dotted parameters and the dimension filter. Every one of the
2,070 singular finite-queue stopping states is accounted for before that
filter; expanded singular states are retained.

For a face, the same external A1 and I2(5) counts are added at all four corners.
Every choice leaving residual rank zero or at least two is checked. Connected
and A1-cross bonuses depend on the entire A1 count and on whether any entire
I2(5) component is present. A1 cross payments divide by the actual A1 count;
I2(5) cross payments divide by twice the actual I2(5) count. Every resulting
face sum is at least one. Nonnegative bonuses preserve the independent
76-marking triangular-face proof and the short-path/H4 payments used in
attachment pruning.

At fixed counts `(a,b)`, each non-A1 destination root has exactly `a` A1/root
pairs. Each destination root outside A1 and I2(5) has exactly `2b` I2(5)/root
pairs. The reviewer checks both normalization identities explicitly. These
two cross families are disjoint; A1-to-I2(5) pairs belong only to the first.
The component sums therefore remain additive. A reverse-type unbounded
dynamic program, distinct from the producer's rank-first recurrence, checks
all residual finite forests. The A1 and I2(5) components themselves are
counted explicitly. Every saved maximizing forest is also checked directly.

Generation and a fresh saved replay pass. A separate source review found no
problem with the normalization, common spectators, variable ranges or forest
enumeration. This verifies the weight rule, not the geometric completeness
of the input branches. That is the role of the separate face-coverage binder.
There is no numerical optimality, formal proof-assistant or priority claim.

```sh
python3 scripts/fresh_sqrt5_dimension24_i25_cross_20261015_review.py
python3 scripts/fresh_sqrt5_dimension24_i25_cross_20261015_review.py --verify
```
