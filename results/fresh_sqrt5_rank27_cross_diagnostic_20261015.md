# Exact obstruction for the current rank-27 weight family

The specified connected-pair and unique-A1/E7 weight family has exact optimum
**73** in dimension 27, above the strict face-count threshold **143/2**. This
family therefore does not exclude dimension 27. The result concerns the weight
model and the conservative frozen input catalogue, not existence of a polytope
or an obstruction to other weighting rules.

The [diagnostic](../scripts/fresh_sqrt5_rank27_cross_diagnostic_20261015.py) and
[exact certificate](fresh_sqrt5_rank27_cross_diagnostic_20261015.json) retain
only dimension-27-feasible spectator cases from the frozen
[parameter weight catalogue](fresh_sqrt5_parameter_weight_stress_20261015_cross.md).
They use the same nonnegative baseline and pair orbits. Exact rational primal
and dual certificates agree, and separation checks all allowed finite forests
of total rank 27. Generation and saved replay pass; no independent review of
this new diagnostic is claimed.

The maxima in the zero, one, and at-least-two A1 contexts are respectively
141/2, 73, and 73. The latter two are attained on A1+3E8+I2(10) and 3A1+3E8.
The dual certificate uses nine old E8 constraints and dotted state 107, all
with one external A1, together with the first of these forests. They force
the context-one E8 component total to be at least 24. In dotted state 107,
the core A1+E7 corners acquire a second A1 from the spectator; the current
unique-A1 cross bonus is therefore unavailable there.

This suggests testing weights with finer A1 counts and cross bonuses that
remain available when several A1 components occur. Such a refinement needs
fresh exact face and vertex checks. It is not established by this diagnostic.

```sh
python3 scripts/fresh_sqrt5_rank27_cross_diagnostic_20261015.py --verify
```
