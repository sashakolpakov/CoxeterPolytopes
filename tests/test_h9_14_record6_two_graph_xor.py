"""Regression test for the exploratory record-6 parity encoding."""

from itertools import product

import z3

from scripts.explore_h9_14_record6_two_graph import _xor_all as full_xor
from scripts.explore_h9_14_record6_two_graph_core import _xor_all as core_xor


def test_three_input_xor_keeps_the_third_operand() -> None:
    variables = z3.Bools("record6_xor_a record6_xor_b record6_xor_c")
    for helper in (full_xor, core_xor):
        expression = helper(*variables)
        for values in product((False, True), repeat=3):
            substitutions = [
                (variable, z3.BoolVal(value))
                for variable, value in zip(variables, values)
            ]
            observed = z3.is_true(z3.simplify(z3.substitute(expression, substitutions)))
            assert observed is (sum(values) % 2 == 1)
