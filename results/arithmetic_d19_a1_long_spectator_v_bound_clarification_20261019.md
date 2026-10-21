# Explicit doubled variable in the D19+A1 long-spectator bound

19 October2026. Paper-only clarification of Section7 of
`arithmetic_d19_a1_two_fork_half_refund_20261019.md`, SHA256
`fcb32a65ac266c5aae7a7345edff0a0f69e52b3a5bb9743221ce79306039a739`.
The pinned source is unchanged. No scientific job was run.

Set **V=2y** throughout that section's norm argument. The printed
phrase “its absolute value is less than9/20” refers to V', not y'.
Here is the explicit substitution, so the factor2 can be checked.

For spectator bond3 the strict conjugate cost is

    19*y'^2/2+y'+3/4<1.

Substituting y'=V'/2 and multiplying by8 gives

    19*V'^2+4V'-2<0,
    (-2-sqrt(42))/19<V'<(-2+sqrt(42))/19.

Thus |V'|<(2+sqrt(42))/19<9/20. The last comparison follows, for
example, from sqrt(42)<13/2 and17/38<9/20.

For spectator bond4 the strict conjugate cost becomes

    19*V'^2+4V'<0,
    -4/19<V'<0,

which has the same weaker absolute bound9/20. The native cycle
makes V integral. In the finite-y branch its physical value is
0<V<2, so its nonzero integer norm in degree d>=2 has absolute
value below

    2*(9/20)^(d-1)<=9/10<1.

No extra native-label exclusion or change to the proof is needed:
the small interval was already for the doubled integral variable.
The independent old-plus-d block forces K=Q(V), since its cost is
a rational polynomial in V and is greater than1 physically.
