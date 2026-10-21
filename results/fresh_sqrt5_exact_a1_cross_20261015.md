# Exact A1 counts with normalized E7 cross weights

The [partial optimization](../scripts/fresh_sqrt5_exact_a1_cross_20261015.py)
still has exact optimum **73** in dimension 27, above **143/2**. Its
[certificate](fresh_sqrt5_exact_a1_cross_20261015.json) passes generation and
saved replay. No independent review or dimension theorem is claimed.

Connected-pair bonuses now depend on the exact number a of entire A1
components. When a is positive, each marked pair joining an A1 to a rooted
E7 receives its context-dependent bonus divided by a. The total cross
contribution is then additive over E7 roots: each root makes exactly a such
pairs. For each fixed a, an exact unbounded rank dynamic program enumerates
all forests of non-A1 components on the remaining 27−a dimensions.

The calculation imposes 14,128 dimension-feasible cases from 1,343 parameter-
marked face graphs and has 3,729 variables. Exact rational primal and dual
certificates agree. The new dual obstruction uses seven face constraints
and the vertex forest **A1+D5+3D7**, of rank 27 and forced total 73. Thus
refining only the A1 count and the E7 cross term does not suffice. Allowing
normalized A1 crosses with other whole components is a separate next test.

The finite attachment catalogue is not asserted complete here. Some input
families are conservative necessary configurations, without realizability
claims. The optimum is an obstruction only within this specified model.

```sh
python3 scripts/fresh_sqrt5_exact_a1_cross_20261015.py --verify
```
