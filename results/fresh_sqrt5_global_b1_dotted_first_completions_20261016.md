# Exact dotted completions of the underpaid first extensions

The [source](../scripts/fresh_sqrt5_global_b1_dotted_first_completions_20261016.py)
and [certificate](fresh_sqrt5_global_b1_dotted_first_completions_20261016.json)
check all 122 symbolic dotted first extensions underpaid by the
[global proposal](fresh_sqrt5_global_discount_20261016.md). Generation and a
fresh saved replay passed.

Four rows fail a conjugate principal-matrix test. Another 73 have no
field-valued completion with no further touching white roots. On the remaining
45 rows the exact corank-one distance equations give 42 physical dotted
completions: 28 fail cyclic integrality and 14 pass the necessary arithmetic
tests. Every survivor has physical inertia (k-2,1,1), conjugate inertia
(k-1,0,1), an exact mixed-sign root relation, and alternating signs on the
four black roots. All 14 remain underpaid in their saved common spectator
contexts and have been supplied to the repair calculation.

Absence of a zero-further-contact completion does not reject a nonsingular
core that could acquire more touching white roots. These 14 survivors are
necessary arithmetic local data; this screen does not prove their global
compact realizability. It also does not enumerate all dotted parameters for
growing cores or establish attachment closure. Those are separate obligations.

```sh
python3 scripts/fresh_sqrt5_global_b1_dotted_first_completions_20261016.py --verify
```
