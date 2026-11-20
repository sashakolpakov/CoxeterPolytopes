#!/usr/bin/env python3
"""A source-bound 26-feature family containing the current rank-58 law.

Importing this module reads no data and imports no scientific worker.
The coordinating runner alone loads the previously verified catalogue and
evaluates features. No forest list, old optimum, face-payment theorem, or
weight-dependent geometric normalization is used.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys
from types import MappingProxyType
from typing import Mapping, Sequence

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "finite-volume-parametric-corner-features-v1"
MIN_DIMENSION = 3
CATALOGUE_DIMENSION = 58
CURRENT_DIMENSION = 58
FEATURE_COUNT = 26
MARGIN_COLUMN = "rho"
MAX_RECEIPT_BYTES = 2 * 1024**2
MAX_SOURCE_BYTES = 16 * 1024**2
Type = tuple[str, int]
Mark = tuple[int, int]
Vector = tuple[Fraction, ...]
ZERO = Fraction(0)
ONE = Fraction(1)

PINS = {
    "scripts/finite57_literal_weights_audit_20261024.py":
        "f6eba7af77af50b821e9e99c589388c1352f9ba5a99052aaeee64393f5ada1e0",
    "scripts/finite_volume_corner_tiling_weights_20261106.py":
        "141d4d13cbd163c298117f132a6e27d39a36e583e12b083913c33707856b7e19",
    "scripts/finite_volume_corner_selected_weights_20261106.py":
        "004193401323f59f6e30887dab934ed4eaf037322ce95a12c30b04de013368ec",
    "scripts/finite_volume_corner_template_catalogue_20261106.py":
        "a3689db2cfbbdeb511297a3665ba2b50460ef06e3dd1c4c8605482253f874a97",
    "scripts/check_finite_volume_corner_template_catalogue_20261106.py":
        "301ba42695a6aa170ce68d6361b979d11282512d1afb8a60687c152f92cb8ec0",
    "scripts/finite_volume_lp_ideal_rows_20261106.py":
        "c6b0b859440360e35478d939607515480ff8018ac6ea38236ea71d37c11d996a",
    "notes/finite_volume_corner_template_catalogue_20261106.md":
        "9fadd98f1c8b1b1717531a573f3dc6d09848b0c40b36198504b91d64dda9075a",
}
RECEIPT_SOURCES = {
    name: PINS[name] for name in (
        "notes/finite_volume_corner_template_catalogue_20261106.md",
        "scripts/finite57_literal_weights_audit_20261024.py",
        "scripts/finite_volume_corner_template_catalogue_20261106.py",
        "scripts/finite_volume_lp_ideal_rows_20261106.py",
    )
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _sha(value):
    require(type(value) is str and re.fullmatch(r"[0-9a-f]{64}", value),
            "A lowercase SHA256 is required")
    return value


def file_digest(path, limit=MAX_SOURCE_BYTES):
    path = Path(path)
    require(path.is_file() and path.stat().st_size <= limit,
            "Missing or oversized bound input")
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def verify_sources():
    for relative, expected in PINS.items():
        require(file_digest(ROOT / relative) == expected,
                "A pinned feature source changed: " + relative)


def _rational(value):
    """Accept exact numeric values only; never silently rationalize a float."""
    if isinstance(value, Fraction):
        return value
    if type(value) is int:
        return Fraction(value)
    require(type(value) is str and len(value) <= 256
            and re.fullmatch(r"-?[0-9]+(?:/[1-9][0-9]*)?", value),
            "Expected a bounded exact integer or fraction")
    return Fraction(value)


def _unique(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "Duplicate receipt JSON key")
        result[key] = value
    return result


def _read_receipt(path, expected_sha):
    path = Path(path).resolve()
    require(path.is_relative_to(ROOT), "The bound receipt must be inside the repository")
    require(path.is_file() and path.stat().st_size <= MAX_RECEIPT_BYTES,
            "Missing or oversized catalogue verification receipt")
    raw = path.read_bytes()
    require(len(raw) <= MAX_RECEIPT_BYTES
            and hashlib.sha256(raw).hexdigest() == _sha(expected_sha),
            "Consumed catalogue receipt differs from its supplied hash")
    result = json.loads(raw, object_pairs_hook=_unique,
                        parse_constant=lambda _: (_ for _ in ()).throw(
                            ValueError("Nonfinite receipt JSON")))
    require(type(result) is dict, "A catalogue verification object is required")
    return path, result


def _load_table_adapter():
    """Load only the frozen definition/data-reader module, at its pinned path."""
    path = ROOT / "scripts/finite_volume_corner_tiling_weights_20261106.py"
    name = "_parametric_corner_table_adapter_20261107"
    if name in sys.modules:
        module = sys.modules[name]
        require(Path(module.__file__).resolve() == path.resolve(),
                "The cached table adapter has an unexpected path")
        return module
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, "Cannot load the pinned table adapter")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        del sys.modules[name]
        raise
    return module


def _count(component):
    return ("count", component[0], component[1])


GUARD_EXPRESSIONS = MappingProxyType({
    "has_a1": ("at_least", ("count", "A", 1), 1),
    "has_small_a": ("at_least", ("add", ("count", "A", 1),
                                 ("count", "A", 2), ("count", "A", 3)), 1),
    "two_rank15": ("at_least", ("add", ("count", "A", 15), ("count", "D", 15)), 2),
    "three_rank15": ("at_least", ("add", ("count", "A", 15), ("count", "D", 15)), 3),
    "a1_and_a10": ("multiply", ("at_least", ("count", "A", 1), 1),
                    ("at_least", ("count", "A", 10), 1)),
    "a1_and_a9": ("multiply", ("at_least", ("count", "A", 1), 1),
                   ("at_least", ("count", "A", 9), 1)),
    "a9_and_a2": ("multiply", ("at_least", ("count", "A", 9), 1),
                   ("at_least", ("count", "A", 2), 1)),
    "a11_or_a9_a2": ("either", ("at_least", ("count", "A", 11), 1),
                      ("multiply", ("at_least", ("count", "A", 9), 1),
                       ("at_least", ("count", "A", 2), 1))),
    "a10_and_two_a1": ("multiply", ("at_least", ("count", "A", 10), 1),
                        ("at_least", ("count", "A", 1), 2)),
})

# This order is the frozen CurrentWeightTables.same_deltas order.
SAME_LAYOUT = (
    (("A", 13), (0, 12), "has_a1", "1"),
    (("A", 14), (0, 13), "has_a1", "1/2"),
    (("A", 15), (0, 13), "has_a1", "1/4"),
    (("A", 15), (1, 14), "has_a1", "1/4"),
    (("D", 15), (1, 13), "has_a1", "1/4"),
    (("D", 15), (1, 14), "has_a1", "1/4"),
    (("A", 16), (1, 12), "has_a1", "5/8"),
    (("A", 16), (3, 14), "has_a1", "5/8"),
    (("D", 16), (0, 14), "has_a1", "1/2"),
    (("D", 16), (0, 15), "has_a1", "1/2"),
    (("D", 14), (0, 12), "has_small_a", "1/2"),
    (("D", 14), (0, 13), "has_small_a", "1/2"),
    (("D", 13), (0, 11), "three_rank15", "-1/2"),
    (("D", 13), (0, 12), "three_rank15", "-1/2"),
    (("D", 28), (0, 12), "two_rank15", "1/2"),
)


@dataclass(frozen=True)
class CrossContext:
    feature_id: str
    left: Type
    left_nodes: tuple[int, ...]
    right: Type
    right_nodes: tuple[int, ...]
    guard: str
    current_coefficient: str

    @property
    def orbit_size(self):
        return len(self.left_nodes) * len(self.right_nodes)

    def aggregate_expression(self):
        if self.left == self.right:
            require(self.left == ("A", 2) and self.orbit_size == 4,
                    "Unsupported equal-type cross feature")
            # 4 * binom(m_A2, 2), not 4*m_A2^2 and not an oriented-pair count.
            raw = ("multiply", 2, _count(self.left), ("add", _count(self.left), -1))
        else:
            raw = ("multiply", self.orbit_size, _count(self.left), _count(self.right))
        return ("multiply", raw, GUARD_EXPRESSIONS[self.guard])


CROSS_CONTEXTS = (
    CrossContext("cross_a4_a12", ("A", 4), (0, 3), ("A", 12), (1, 10), "has_a1", "1/8"),
    CrossContext("cross_a2_a4", ("A", 2), (0, 1), ("A", 4), (0, 3), "a1_and_a10", "1/4"),
    CrossContext("cross_a2_d5", ("A", 2), (0, 1), ("D", 5), (0,), "a1_and_a9", "1/2"),
    CrossContext("cross_a1_d5_first", ("A", 1), (0,), ("D", 5), (0,), "a11_or_a9_a2", "3/8"),
    CrossContext("cross_a1_d5_second", ("A", 1), (0,), ("D", 5), (0,), "a9_and_a2", "5/8"),
    CrossContext("cross_a12_d16", ("A", 12), (0, 11), ("D", 16), (0,), "two_rank15", "1/2"),
    CrossContext("cross_a2_a2_distinct", ("A", 2), (0, 1), ("A", 2), (0, 1),
                 "a10_and_two_a1", "9/32"),
)
FEATURE_IDS = (
    "constant", "finite_baseline", "ideal_baseline", "ordinary_cross",
    *(f"same_{component[0]}{component[1]}_{pair[0]}_{pair[1]}_{guard}"
      for component, pair, guard, _ in SAME_LAYOUT),
    *(item.feature_id for item in CROSS_CONTEXTS),
)
COEFFICIENT_COLUMNS = tuple("theta_" + identifier for identifier in FEATURE_IDS)


def expression_support(expression):
    if type(expression) is int or isinstance(expression, Fraction):
        return frozenset()
    op = expression[0]
    if op == "count":
        return frozenset(((expression[1], expression[2]),))
    children = expression[1:2] if op == "at_least" else expression[1:]
    require(op in ("add", "multiply", "at_least", "either"), "Unknown count-expression operation")
    return frozenset().union(*(expression_support(child) for child in children))


def evaluate_count_expression(expression, counts):
    """Exact expression evaluator; callers validate the count domain first."""
    if type(expression) is int or isinstance(expression, Fraction):
        return Fraction(expression)
    op = expression[0]
    if op == "count":
        return Fraction(counts.get((expression[1], expression[2]), 0))
    if op == "add":
        return sum((evaluate_count_expression(child, counts) for child in expression[1:]), ZERO)
    if op == "multiply":
        value = ONE
        for child in expression[1:]:
            value *= evaluate_count_expression(child, counts)
        return value
    if op == "at_least":
        return Fraction(int(evaluate_count_expression(expression[1], counts) >= expression[2]))
    if op == "either":
        return Fraction(int(any(evaluate_count_expression(child, counts) for child in expression[1:])))
    raise ValueError("Unknown count-expression operation")


@dataclass(frozen=True)
class FeatureTerm:
    feature_index: int
    expression: tuple | int | Fraction

    @property
    def support(self):
        return expression_support(self.expression)


def coefficient_dot(features, coefficients):
    require(len(features) == len(coefficients) == FEATURE_COUNT, "Expected 26 feature coordinates")
    return sum((_rational(a) * _rational(b) for a, b in zip(features, coefficients)), ZERO)


def _json_expression(expression):
    if isinstance(expression, Fraction):
        return {"rational": str(expression)}
    if type(expression) is int:
        return expression
    return [expression[0], *(_json_expression(child) if isinstance(child, (tuple, Fraction))
                             else child for child in expression[1:])]


class FeatureFamily:
    """Native-numbered corner features with authenticated immutable tables."""

    feature_ids = FEATURE_IDS
    coefficient_columns = COEFFICIENT_COLUMNS
    margin_column = MARGIN_COLUMN
    supported_dimension_interval = (MIN_DIMENSION, CATALOGUE_DIMENSION)

    def __init__(self, tables, bindings, source_sha):
        require(len(FEATURE_IDS) == len(set(FEATURE_IDS)) == FEATURE_COUNT,
                "The feature identifiers are not a complete disjoint index list")
        observed = tuple((item.component, item.pair, item.guard, str(item.increment))
                         for item in tables.same_deltas)
        require(observed == SAME_LAYOUT, "The source-bound same-context order or increments changed")
        self.within = MappingProxyType({
            key: MappingProxyType(dict(pairs)) for key, pairs in tables.within.items()
        })
        self.finite_types = frozenset(key[1:] for key in self.within if not key[0])
        self.affine_types = frozenset(key[1:] for key in self.within if key[0])
        self.finite_base_sums = MappingProxyType({
            component: sum(self.within[False, *component].values(), ZERO)
            for component in self.finite_types
        })
        self.affine_base_sums = MappingProxyType({
            component: sum(self.within[True, *component].values(), ZERO)
            for component in self.affine_types
        })
        ordinary, totals = {}, {}
        for left, a, right, b, value in tables.ordinary_cross:
            require(left < right, "The ordinary catalogue must use canonical distinct-type pairs")
            key = (left, a, right, b)
            require(key not in ordinary, "Repeated ordinary marked entry")
            ordinary[key] = value
            totals[left, right] = totals.get((left, right), ZERO) + value
        self.ordinary_pairs = MappingProxyType(ordinary)
        self.ordinary_component_totals = MappingProxyType(totals)
        self.current_coefficients = (
            ZERO, ONE, ONE, ONE,
            *(Fraction(item[3]) for item in SAME_LAYOUT),
            *(Fraction(item.current_coefficient) for item in CROSS_CONTEXTS),
        )
        self.bindings = MappingProxyType(dict(bindings))
        self.source_sha256 = source_sha

    def assert_bindings_unchanged(self):
        verify_sources()
        require(file_digest(__file__) == self.source_sha256, "The feature provider source changed")
        for path, expected in self.bindings.items():
            # The receipt was separately size-checked when consumed. Do not
            # infer an input's role or byte limit from its chosen filename.
            require(file_digest(path, 64 * 1024**2) == expected, "A bound feature input changed")

    def _dimension(self, n):
        require(type(n) is int and MIN_DIMENSION <= n <= CATALOGUE_DIMENSION,
                "This catalogue-backed family supports dimensions 3 through 58 only")

    def _type(self, value, *, affine, n):
        require(isinstance(value, (tuple, list)) and len(value) == 2
                and type(value[0]) is str and type(value[1]) is int,
                "A component must be a (family, Coxeter-rank) pair")
        component = (value[0], value[1])
        domain = self.affine_types if affine else self.finite_types
        require(component in domain and component[1] + int(affine) <= n,
                "Component absent from the supported rank domain")
        return component

    def _counts(self, counts, *, affine, n):
        self._dimension(n)
        require(isinstance(counts, Mapping), "A type-to-multiplicity mapping is required")
        result = {}
        for raw_type, multiplicity in counts.items():
            component = self._type(raw_type, affine=affine, n=n)
            require(type(multiplicity) is int and multiplicity >= 0,
                    "Component multiplicities must be nonnegative integers")
            require(component not in result, "Repeated normalized component type")
            if multiplicity:
                result[component] = multiplicity
        required_rank = n - int(affine)
        require(sum(component[1] * multiplicity for component, multiplicity in result.items()) == required_rank,
                "The complete vertex rank does not equal its ambient rank requirement")
        return result

    def finite_corner(self, forest: Sequence[Type], first: Mark, second: Mark, *, n: int) -> Vector:
        self._dimension(n)
        require(isinstance(forest, (tuple, list)) and 1 <= len(forest) <= n,
                "A complete ordered list of finite component copies is required")
        forest = tuple(self._type(component, affine=False, n=n) for component in forest)
        require(sum(component[1] for component in forest) == n, "Finite corner forest has wrong total rank")
        marks = []
        for mark in (first, second):
            require(isinstance(mark, (tuple, list)) and len(mark) == 2
                    and type(mark[0]) is int and type(mark[1]) is int
                    and 0 <= mark[0] < len(forest) and 0 <= mark[1] < forest[mark[0]][1],
                    "A finite mark lies outside its component copy")
            marks.append((mark[0], mark[1]))
        require(marks[0] != marks[1], "A corner requires two distinct marked nodes")
        (ca, a), (cb, b) = marks
        left, right = forest[ca], forest[cb]
        counts = Counter(forest)
        features = [ZERO] * FEATURE_COUNT
        features[0] = ONE
        if ca == cb:
            pair = tuple(sorted((a, b)))
            features[1] = self.within[False, *left].get(pair, ZERO)
            for index, (component, target, guard, _) in enumerate(SAME_LAYOUT, 4):
                if left == component and pair == target:
                    features[index] = evaluate_count_expression(GUARD_EXPRESSIONS[guard], counts)
            return tuple(features)
        key = (left, a, right, b) if left < right else (right, b, left, a)
        features[3] = self.ordinary_pairs.get(key, ZERO)
        for index, item in enumerate(CROSS_CONTEXTS, 19):
            forward = left == item.left and a in item.left_nodes and right == item.right and b in item.right_nodes
            reverse = right == item.left and b in item.left_nodes and left == item.right and a in item.right_nodes
            if forward or reverse:
                features[index] = evaluate_count_expression(GUARD_EXPRESSIONS[item.guard], counts)
        return tuple(features)

    def affine_corner(self, component: Type, a: int, b: int, *, n: int) -> Vector:
        self._dimension(n)
        component = self._type(component, affine=True, n=n)
        order = component[1] + 1
        require(type(a) is int and type(b) is int and 0 <= a < order and 0 <= b < order and a != b,
                "An affine corner requires two distinct native nodes of its varying factor")
        features = [ZERO] * FEATURE_COUNT
        features[0] = ONE
        features[2] = self.within[True, *component].get(tuple(sorted((a, b))), ZERO)
        return tuple(features)

    def finite_aggregate(self, counts: Mapping[Type, int], *, n: int) -> Vector:
        counts = self._counts(counts, affine=False, n=n)
        features = [ZERO] * FEATURE_COUNT
        features[0] = Fraction(n * (n - 1), 2)
        features[1] = sum((self.finite_base_sums[t] * m for t, m in counts.items()), ZERO)
        features[3] = sum((value * counts.get(left, 0) * counts.get(right, 0)
                           for (left, right), value in self.ordinary_component_totals.items()), ZERO)
        for index, (component, _, guard, _) in enumerate(SAME_LAYOUT, 4):
            features[index] = counts.get(component, 0) * evaluate_count_expression(GUARD_EXPRESSIONS[guard], counts)
        for index, item in enumerate(CROSS_CONTEXTS, 19):
            features[index] = evaluate_count_expression(item.aggregate_expression(), counts)
        return tuple(features)

    def finite_factor_terms(self, *, n: int) -> tuple[FeatureTerm, ...]:
        """Separate local factors; never turn a full baseline sum into one scope.

        Evaluated on a complete rank-n count vector, the terms summed by their
        feature_index equal finite_aggregate. Types whose ranks exceed n are
        identically absent. Context guards retain their declared source scopes.
        """
        self._dimension(n)
        result = [FeatureTerm(0, Fraction(n * (n - 1), 2))]
        for component in sorted(self.finite_types):
            if component[1] <= n and self.finite_base_sums[component]:
                result.append(FeatureTerm(1, ("multiply", self.finite_base_sums[component], _count(component))))
        for (left, right), value in sorted(self.ordinary_component_totals.items()):
            if left[1] + right[1] <= n:
                result.append(FeatureTerm(3, ("multiply", value, _count(left), _count(right))))
        for index, (component, _, guard, _) in enumerate(SAME_LAYOUT, 4):
            if component[1] <= n:
                result.append(FeatureTerm(index, ("multiply", _count(component), GUARD_EXPRESSIONS[guard])))
        for index, item in enumerate(CROSS_CONTEXTS, 19):
            if item.left[1] + item.right[1] <= n:
                result.append(FeatureTerm(index, item.aggregate_expression()))
        return tuple(result)

    def affine_factor_aggregate(self, component: Type, *, n: int) -> Vector:
        self._dimension(n)
        component = self._type(component, affine=True, n=n)
        rank = component[1]
        features = [ZERO] * FEATURE_COUNT
        features[0] = Fraction(rank * (rank + 1), 2)
        features[2] = self.affine_base_sums[component]
        return tuple(features)

    def ideal_link_size(self, factor_counts: Mapping[Type, int], *, n: int) -> int:
        counts = self._counts(factor_counts, affine=True, n=n)
        value = 1
        for component, multiplicity in counts.items():
            value *= (component[1] + 1) ** multiplicity
        return value

    def ideal_normalized_aggregate(self, factor_counts: Mapping[Type, int], *, n: int) -> Vector:
        counts = self._counts(factor_counts, affine=True, n=n)
        features = [ZERO] * FEATURE_COUNT
        features[0] = Fraction(n - 1, 2)
        features[2] = sum((multiplicity * self.affine_base_sums[component] / (component[1] + 1)
                           for component, multiplicity in counts.items()), ZERO)
        return tuple(features)

    def ideal_aggregate(self, factor_counts: Mapping[Type, int], *, n: int) -> Vector:
        normalized = self.ideal_normalized_aggregate(factor_counts, n=n)
        value = self.ideal_link_size(factor_counts, n=n)
        return tuple(value * entry for entry in normalized)

    def current_value(self, features: Vector, *, n: int) -> Fraction:
        require(n == CURRENT_DIMENSION, "CURRENT compatibility is claimed only at dimension 58")
        return coefficient_dot(features, self.current_coefficients)

    def metadata(self):
        specs = [
            {"index": 0, "id": FEATURE_IDS[0], "kind": "constant_corner",
             "finite_aggregate": "n*(n-1)/2", "ideal_normalized_aggregate": "(n-1)/2"},
            {"index": 1, "id": FEATURE_IDS[1], "kind": "finite_literal_same_baseline",
             "finite_aggregate": "sum_T m_T * sum_{a<b} finite_same_weight(T,a,b)"},
            {"index": 2, "id": FEATURE_IDS[2], "kind": "ideal_literal_varying_affine_baseline",
             "ideal_normalized_aggregate": "sum_Gamma m_Gamma * Q_Gamma/(rank(Gamma)+1)"},
            {"index": 3, "id": FEATURE_IDS[3], "kind": "finite_literal_ordinary_cross",
             "finite_aggregate": "sum_{T<U} m_T*m_U * ordinary_cross_sum(T,U)"},
        ]
        for index, (component, pair, guard, _) in enumerate(SAME_LAYOUT, 4):
            specs.append({"index": index, "id": FEATURE_IDS[index], "kind": "same_context_indicator",
                          "component": list(component), "native_pair": list(pair), "guard": guard,
                          "finite_aggregate_expression": _json_expression(
                              ("multiply", _count(component), GUARD_EXPRESSIONS[guard]))})
        for index, item in enumerate(CROSS_CONTEXTS, 19):
            specs.append({"index": index, "id": FEATURE_IDS[index], "kind": "cross_context_indicator",
                          "left": list(item.left), "left_native_nodes": list(item.left_nodes),
                          "right": list(item.right), "right_native_nodes": list(item.right_nodes),
                          "guard": item.guard, "orbit_multiplicity": item.orbit_size,
                          "distinct_component_copies_required": True,
                          "equal_type_uses_unordered_copy_pairs": item.left == item.right,
                          "finite_aggregate_expression": _json_expression(item.aggregate_expression())})
        for spec, column, coefficient in zip(specs, COEFFICIENT_COLUMNS, self.current_coefficients):
            spec["coefficient_column"] = column
            spec["current_coefficient"] = str(coefficient)
            spec["corner_lower_bound"] = "1" if spec["index"] == 0 else "0"
            spec["corner_upper_bound"] = "1"
        return {
            "schema": SCHEMA, "feature_count": FEATURE_COUNT, "features": specs,
            "margin_column": MARGIN_COLUMN, "coefficient_signs": "free; total corner nonnegativity is a separate obligation",
            "native_node_numbering": True, "coherent_vertex_decorations_required": True,
            "supported_ambient_dimensions": {"minimum": MIN_DIMENSION, "maximum": CATALOGUE_DIMENSION},
            "catalogue_ambient_dimension": CATALOGUE_DIMENSION,
            "finite_component_rank_maximum": CATALOGUE_DIMENSION,
            "affine_component_rank_maximum": CATALOGUE_DIMENSION - 1,
            "larger_dimensions_require_new_complete_catalogue_and_verification": True,
            "current_rule_compatibility_dimension": CURRENT_DIMENSION,
            "other_supported_dimensions": "formal feature laws only; no inherited face-payment or final-rule theorem",
            "finite_I2_alias": "all labels at least 5 have the same features; this does not identify Gram matrices",
            "ideal_features_use_only_varying_affine_factor": True,
            "guard_expressions": {name: _json_expression(value) for name, value in GUARD_EXPRESSIONS.items()},
            "source_sha256": self.source_sha256, "sources": dict(PINS), "input_bindings": dict(self.bindings),
            "literal_coefficient_verification_rerun": False, "forest_enumeration_performed": False,
            "face_payment_proved": False, "weight_dependent_label_normalization_used": False,
            "threshold_monotonicity_proved": False, "dimension_bound_proved": False,
        }


def load_feature_family(catalogue_path, catalogue_sha256,
                        catalogue_verification_path, catalogue_verification_sha256) -> FeatureFamily:
    """Authenticate a previously checked catalogue; do not rerun its checker."""
    verify_sources()
    own_sha = file_digest(__file__)
    catalogue_path = Path(catalogue_path).resolve()
    require(catalogue_path.is_relative_to(ROOT), "The catalogue must be inside the repository")
    _sha(catalogue_sha256)
    receipt_path, receipt = _read_receipt(catalogue_verification_path, catalogue_verification_sha256)
    require(receipt.get("schema") == "uniform-corner-template-catalogue-check-v1"
            and receipt.get("status") == "EXACT_CONNECTED_CORNER_CATALOGUE_VERIFIED"
            and receipt.get("dimension") == CATALOGUE_DIMENSION
            and receipt.get("catalogue_sha256") == catalogue_sha256
            and receipt.get("checker_sha256") == PINS["scripts/check_finite_volume_corner_template_catalogue_20261106.py"]
            and receipt.get("producer_sha256") == PINS["scripts/finite_volume_corner_template_catalogue_20261106.py"]
            and receipt.get("sources") == RECEIPT_SOURCES
            and receipt.get("input_bindings", {}).get(str(catalogue_path)) == catalogue_sha256,
            "The complete independently checked catalogue is not authenticated")
    for flag in ("complete_task_prefix_verified", "all_graphs_and_variants_verified",
                 "all_within_marked_pairs_verified", "all_marked_cross_pairs_in_support_verified"):
        require(receipt.get(flag) is True, "The catalogue receipt omits required coverage: " + flag)
    adapter = _load_table_adapter()
    tables = adapter.read_current_weight_tables(catalogue_path, catalogue_sha256)
    bindings = {str(catalogue_path): catalogue_sha256, str(receipt_path): catalogue_verification_sha256}
    result = FeatureFamily(tables, bindings, own_sha)
    result.assert_bindings_unchanged()
    return result
