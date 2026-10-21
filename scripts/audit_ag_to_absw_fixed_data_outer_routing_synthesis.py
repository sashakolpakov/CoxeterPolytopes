#!/usr/bin/env python3
"""Independent hostile audit of the fixed-data/outer-routing synthesis.

The replay deliberately does not import the owner certifier.  It pins every
owner synthesis document, the owner certificate package, and the three
owner/audit theorem pairs; checks each synthesis document separately; and
reconstructs the narrow logical scope that survives the hostile controls.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from fractions import Fraction
from math import isqrt
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/"
    "ag_to_absw_fixed_data_outer_routing_synthesis_independent_hostile_audit.json"
)
SCHEMA = (
    "coxeter.ag_to_absw_fixed_data_outer_routing_synthesis."
    "independent_hostile_audit.v1"
)

# Terminal owner-document/package pins after the hostile integration repair.
DOCUMENT_PINS = {
    "goal": (
        "GOAL.md",
        "0b1736db46f8872ddf3868159c273620642a328ac94bba0488e806ae83f80b4d",
    ),
    "plan": (
        "RESEARCH_PLAN.md",
        "d4533215663c8398fad35de9d3a764b8f6dcf1e6fdea103eda8a5fc34515f983",
    ),
    "completion": (
        "results/canonical_chamber_four_property_completion_audit.md",
        "f1c56a0ade1b37a498b759262b3309c97f96a54ae00900d11d4035a490a05e36",
    ),
    "inverse": (
        "results/inverse_carrier_transfer_ledger.md",
        "ec572e1ba7f0936f3e1ed1eb2c0837780b58b807c2af5009988a498c8e46b2f9",
    ),
    "nikulin": (
        "results/nikulin_to_ag_translation_ledger.md",
        "f016006bbe9424a757b949feb17ad73cd39c9f7303448e1265912159eb858237",
    ),
    "synthesis": (
        "results/ag_to_absw_fixed_data_closure_and_outer_routing_obstruction.md",
        "9c7fdf14ab62629ae1a47d04a6e7ede91a8605e86b936305d2f2b957a6ff77b8",
    ),
}

OWNER_PACKAGE_PINS = {
    "owner_certifier": (
        "scripts/certify_ag_to_absw_fixed_data_outer_routing_synthesis.py",
        "cb78b46f171f6fcc07c99e6139a4198d0f628be628b79e3317471b8179a57487",
    ),
    "owner_archive": (
        "data/enumerations/ag_to_absw_fixed_data_outer_routing_synthesis.json",
        "d0adf354e185f601d2e96ed94ef8225742cb07afa6d6790626f0550932ee5911",
    ),
    "owner_test": (
        "tests/test_ag_to_absw_fixed_data_outer_routing_synthesis.py",
        "ff941226e19d85b50d881669840b15cd973e2b7dba2f3dc315f8239c0813162d",
    ),
}

EVIDENCE_PINS = {
    "finite_pair_owner": (
        "results/projected_k3_maximal_abstract_chamber_finiteness.md",
        "59dd458acf14d18021f6c5583096da13f39c83115fc136d6ba91bcc850ce3603",
    ),
    "finite_pair_audit": (
        "results/projected_k3_maximal_abstract_chamber_finiteness_independent_hostile_audit.md",
        "edc663c0e7f35880123149d81d832812d5f28837c837df16eaa0f1daf725916a",
    ),
    "x3_owner": (
        "results/projected_k3_universality_rank_saturation_no_go.md",
        "7ed5176c0c0a36e036c0e12dfa6acf1ee15d90d5ec7143ca5df59ebc3b7143d2",
    ),
    "x3_audit": (
        "results/projected_k3_universality_rank_saturation_no_go_independent_hostile_audit.md",
        "89296f04bac6b0b4ed68d9b86936a08ede515ba5bf9c847abeaadbc628926071",
    ),
    "h21_owner": (
        "results/borcherds_h21_outer_routing_counterexample.md",
        "fa899977f418d6fa93f7d5270e067c28a648c5e55f2636cb2b362e381534821d",
    ),
    "h21_audit": (
        "results/borcherds_h21_outer_routing_counterexample_independent_hostile_audit.md",
        "8d705a619e9c1ffff3ad69765e90bb1a1a2a5d8713015ecca79382e8a1ccbc86",
    ),
}

ALL_PINS = {**DOCUMENT_PINS, **OWNER_PACKAGE_PINS, **EVIDENCE_PINS}

EVIDENCE_BASENAMES = tuple(
    Path(relative).name for relative, _expected in EVIDENCE_PINS.values()
)

LINK_RE = re.compile(r"!?\[[^\]]*\]\((<?[^)\s>]+>?)")
STRIPPED_TEX_RE = re.compile(
    r"(?<![\\A-Za-z])"
    r"(?:qquad|operatorname|mathcal|mathbb|mathfrak|mathbf|"
    r"subseteq|supseteq|varnothing|infty|mathrm)"
    r"(?![A-Za-z])"
)

# These are the pointwise-overclaim forms found by the hostile integration
# pass.  Correct prose may say that one X3 counterexample refutes a theorem
# universally quantified over all finite-volume inputs; it may not say that
# X3 obstructs each input individually.
DEFECTIVE_X3_PATTERNS = tuple(
    re.compile(pattern, re.IGNORECASE)
    for pattern in (
        r"x3.{0,350}obstructs every exact-module",
        r"x3.{0,500}obstruction applies.{0,180}every finite-volume input",
        r"x3.{0,500}obstruction for every finite-volume exact-module",
        r"x3.{0,500}return is impossible for every finite-volume input",
    )
)

# These four shorthand forms dropped the admissibility, full
# restriction-of-scalars trace, and finite projective order-module gates from
# the field theorem.  A nearby correct statement elsewhere in a long ledger
# is not enough; the local implication itself must retain all three gates.
DEFECTIVE_FIELD_SHORTHANDS = tuple(
    re.compile(pattern, re.IGNORECASE)
    for pattern in (
        r"exact finite-volume full-trace K3 carrier in n>=2",
        r"every exact finite-volume full-trace K3 carrier in dimension at least two",
        r"once an exact finite-volume full-trace K3 return in dimension at least two",
        r"finite field set for every exact finite-volume full-trace K3 carrier in `?n>=2`?",
    )
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def flatten(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def rational_square(value: Fraction) -> bool:
    if value < 0:
        return False
    return (
        isqrt(value.numerator) ** 2 == value.numerator
        and isqrt(value.denominator) ** 2 == value.denominator
    )


def det3(matrix: tuple[tuple[int, int, int], ...]) -> int:
    (a, b, c), (d, e, f), (g, h, i) = matrix
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def validate_pins() -> dict[str, str]:
    actual: dict[str, str] = {}
    for name, (relative, expected) in ALL_PINS.items():
        value = sha256_file(ROOT / relative)
        if value != expected:
            raise AssertionError(
                f"frozen input {name} changed: expected {expected}, got {value}"
            )
        actual[name] = value
    return actual


def local_link_audit(path: Path, text: str) -> dict[str, Any]:
    targets: list[str] = []
    broken: list[str] = []
    escaped_root: list[str] = []
    root = ROOT.resolve()
    for match in LINK_RE.finditer(text):
        raw = match.group(1).strip("<>")
        if raw.startswith(("http://", "https://", "mailto:", "data:", "#")):
            continue
        target = raw.split("#", 1)[0]
        if not target:
            continue
        targets.append(target)
        resolved = (path.parent / target).resolve()
        try:
            resolved.relative_to(root)
        except ValueError:
            escaped_root.append(target)
            continue
        if not resolved.is_file():
            broken.append(target)
    return {
        "local_link_count": len(targets),
        "broken_local_links": sorted(set(broken)),
        "links_escaping_repository": sorted(set(escaped_root)),
        "all_local_links_resolve": not broken and not escaped_root,
    }


def _has(text: str, pattern: str) -> bool:
    return re.search(pattern, text, re.IGNORECASE) is not None


def document_scope_audit(name: str, relative: str) -> dict[str, Any]:
    path = ROOT / relative
    raw = path.read_text(encoding="utf-8")
    text = flatten(raw)
    lower = text.lower()
    defective_x3 = [p.pattern for p in DEFECTIVE_X3_PATTERNS if p.search(text)]
    defective_field = [
        p.pattern for p in DEFECTIVE_FIELD_SHORTHANDS if p.search(text)
    ]
    link_audit = local_link_audit(path, raw)
    stripped = [m.group(0) for m in STRIPPED_TEX_RE.finditer(raw)]

    checks = {
        "all_six_frozen_evidence_reports_linked": all(
            basename in raw for basename in EVIDENCE_BASENAMES
        ),
        "finite_native_pair_condition_present": _has(
            text, r"(?:finite set of|finitely many) native pairs"
        ),
        "nikulin_fixed_data_import_present": (
            "nikulin" in lower and "fixed-data" in lower
        ),
        "reflection_maximality_present_in_count": _has(
            text,
            r"(?:native pairs|fixed-data).{0,450}reflection[- ]maximal"
            r"|reflection[- ]maximal.{0,450}(?:native pairs|fixed-data)",
        ),
        "classification_import_is_explicit": (
            "classification-import" in lower
            or "imports classification" in lower
            or "classification import" in lower
            or "imports nikulin" in lower
        ),
        "ag_only_global_absw_denial_present": (
            "ag-only" in lower
            and "global absw" in lower
            and _has(
                text,
                r"(?:not|neither|no current result|does not).{0,180}"
                r"(?:ag-only|global ABSW)",
            )
        ),
        "x3_universality_scope_present": (
            "x3" in lower
            and "finite-volume" in lower
            and "univers" in lower
        ),
        "x3_one_counterexample_not_pointwise_failure_is_explicit": _has(
            text,
            r"x3.{0,520}(?:one.{0,80}counterexample|counterexample.{0,180}"
            r"universal|attacks all-finite-volume universality)"
            r"|one.{0,80}counterexample.{0,520}x3"
            r"|attacks all-finite-volume universality rather than the maximal-only quantifier",
        ),
        "x3_maximal_only_scope_is_unsettled": (
            "maximal-only" in lower
            and _has(text, r"(?:does not settle|unsettled|rather than)")
        ),
        "x3_has_no_pointwise_all_inputs_overclaim": not defective_x3,
        "h21_absolute_maximality_present": (
            "borcherds" in lower
            and _has(text, r"absolut(?:e|ely).{0,60}reflection[- ]maximal")
        ),
        "h21_current_complex_stack_failure_present": (
            "borcherds" in lower
            and "current complex" in lower
            and _has(text, r"(?:empty|empties|outside|lies outside)")
        ),
        "projected_ambient_hosts_not_needed_for_abstract_count": (
            "projected ambient" in lower
            and _has(
                text,
                r"projected ambient.{0,180}"
                r"(?:unnecessary|not needed|not required|does not require|"
                r"need not be finite|harmless)"
                r"|without.{0,180}projected ambient",
            )
        ),
        "field_claim_retains_exact_full_trace_scope": (
            "field" in lower
            and "full-trace" in lower
            and ("admissib" in lower or "order-module" in lower)
        ),
        "field_scope_has_no_local_shorthand_regression": not defective_field,
        "stripped_tex_scan_clean": not stripped,
        "local_links_resolve": link_audit["all_local_links_resolve"],
    }

    # The canonical completion links the full H21 owner/audit pair and says
    # that the current complex stack fails, but does not repeat the
    # characteristic-two positive construction verbatim.  The other five
    # synthesis documents must preserve it explicitly.
    h21_positive_explicit = (
        _has(text, r"characteristic.{0,20}(?:2|two)")
        and "cremona" in lower
        and _has(text, r"(?:positive|realiz)")
    )
    h21_positive_required = name != "completion"
    checks["h21_characteristic_two_mixed_positive_preserved"] = (
        h21_positive_explicit
        or (
            not h21_positive_required
            and checks["all_six_frozen_evidence_reports_linked"]
        )
    )

    if name == "goal":
        checks["nikulin_dependency_exception_is_not_erased"] = (
            "nikulin and vinberg are henceforth regression corpora, not dependencies"
            not in lower
        )
    if name == "inverse":
        checks["projected_host_finiteness_not_left_as_open_gate"] = (
            "open projected-host gate" not in lower
        )

    return {
        "path": relative,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "defective_x3_patterns": defective_x3,
        "defective_field_scope_patterns": defective_field,
        "h21_characteristic_two_positive_explicit": h21_positive_explicit,
        "stripped_tex_hits": stripped,
        **link_audit,
    }


def validate_evidence_interfaces() -> dict[str, bool]:
    text = {
        name: flatten((ROOT / relative).read_text(encoding="utf-8"))
        for name, (relative, _expected) in EVIDENCE_PINS.items()
    }
    finite_owner = text["finite_pair_owner"]
    finite_audit = text["finite_pair_audit"]
    checks = {
        "finite_pair_is_conditional_and_classification_importing": (
            "Theorem 5.2.1" in finite_owner
            and "reflection maximality" in finite_owner
            and "PASS, conditional and classification-importing" in finite_audit
        ),
        "finite_pair_denies_ambient_host_and_ag_only_claims": (
            "Finite ambient Neron--Severi hosts are not needed" in finite_owner
            and "not an AG-only proof of ABSW" in finite_audit
        ),
        "x3_is_one_counterexample_to_all_finite_volume_universality": (
            "There is no universal `POLYTOPE -> K3 eigen-nef slice` functor for all"
            in text["x3_owner"]
            and "disproves a theorem quantified over every finite-volume"
            in text["x3_owner"]
            and "disproves the statement" in text["x3_audit"]
        ),
        "x3_does_not_settle_maximal_only_scope": (
            "absolute reflection maximality of the `x3` group"
            in text["x3_owner"].lower()
            and "does not address a theorem whose inputs are required to be absolutely"
            in text["x3_audit"]
        ),
        "x3_mixed_divisor_positive_is_preserved": (
            "already has a K3 carrier in a weaker sense" in text["x3_owner"]
            and "does not exclude the known `X3` mixed-divisor K3 carrier"
            in text["x3_audit"]
        ),
        "h21_absolute_maximality_is_audited": (
            "absolutely reflection-maximal" in text["h21_owner"]
            and "absolute-maximality proof survives" in text["h21_audit"]
        ),
        "h21_current_complex_stack_is_empty": (
            "outside the current complex rational" in text["h21_owner"]
            and "direct and exact `RLogComp` complex rational branches are empty"
            in text["h21_audit"]
            and "rank 22 excludes full-trace and projected-root complex-K3 return"
            in text["h21_audit"]
        ),
        "h21_characteristic_two_mixed_carrier_is_positive": (
            "mixed divisor/action realization is a positive fourth-species calibration"
            in text["h21_owner"]
            and "characteristic-two supersingular K3 mixed divisor/Cremona-action carrier"
            in text["h21_audit"]
        ),
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise AssertionError(f"frozen evidence interface failed: {failed}")
    return checks


def validate_owner_archive(
    actual_pins: dict[str, str],
) -> dict[str, bool]:
    archive_path = ROOT / OWNER_PACKAGE_PINS["owner_archive"][0]
    owner = json.loads(archive_path.read_text(encoding="utf-8"))
    expected_document_hashes = {
        relative: actual_pins[name]
        for name, (relative, _expected) in DOCUMENT_PINS.items()
    }
    source_name_map = {
        "finite_pair_owner": "finite_pair_owner",
        "finite_pair_audit": "finite_pair_audit",
        "x3_owner": "x3_owner",
        "x3_audit": "x3_audit",
        "h21_owner": "h21_owner",
        "h21_audit": "h21_audit",
    }
    expected_source_hashes = {
        owner_name: actual_pins[audit_name]
        for owner_name, audit_name in source_name_map.items()
    }
    scope = owner["scope"]
    checks = {
        "owner_schema_is_repaired_v2": (
            owner["schema"]
            == "coxeter.ag_to_absw_fixed_data_outer_routing_synthesis.v2"
        ),
        "owner_all_checks_pass": owner["all_checks_pass"] is True,
        "owner_document_hashes_match_independent_hashes": (
            owner["document_sha256"] == expected_document_hashes
        ),
        "owner_source_hashes_match_independent_hashes": (
            owner["source_sha256"] == expected_source_hashes
        ),
        "owner_denies_projected_ambient_host_requirement": (
            scope["projected_ambient_host_finiteness_needed_for_abstract_count"]
            is False
        ),
        "owner_denies_universal_complex_projected_root_k3": (
            scope["universal_complex_projected_root_k3"] is False
        ),
        "owner_preserves_x3_mixed_divisor_positive": (
            scope["x3_mixed_divisor_k3_positive"] is True
        ),
        "owner_x3_quantifier_is_one_counterexample_not_every_input": (
            scope["x3_no_go_scope"].startswith(
                "one X3 input refutes the universal assertion"
            )
            and "it does not say every input fails" in scope["x3_no_go_scope"]
            and scope["x3_maximal_only_universality_settled"] is False
            and owner["checks"][
                "x3_refutes_all_fv_universality_but_not_each_input_or_maximal_only"
            ]
            is True
        ),
        "owner_field_and_maximality_repair_checks_pass": (
            owner["checks"][
                "field_scope_retains_exact_admissible_full_trace_order_module"
            ]
            is True
            and owner["checks"]["finite_pair_count_retains_reflection_maximality"]
            is True
        ),
        "owner_preserves_h21_absolute_and_complex_empty": (
            scope["h21_absolute_reflection_maximality"] is True
            and scope["h21_current_complex_stack_empty"] is True
        ),
        "owner_preserves_h21_characteristic_two_positive": (
            scope["h21_characteristic_two_mixed_curve_cremona_positive"] is True
        ),
        "owner_denies_global_absw": scope["global_absw_proved"] is False,
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise AssertionError(f"owner archive cross-check failed: {failed}")
    return checks


def hostile_arithmetic_replay() -> dict[str, Any]:
    # X3: one common rational similarity for norm -2 and norm -6 lines
    # would make 3 a square in Q.
    x3_ratio = Fraction(3, 1)

    # The H21 Fano diagonal-point determinant is nonzero in characteristic
    # zero and vanishes modulo two, matching the preserved positive escape.
    fano_matrix = ((1, 1, 0), (1, 0, 1), (0, 1, 1))
    fano_det = det3(fano_matrix)

    # A common similarity cancels from the reflection coefficient.  This is
    # the precise reason native similarity, not ambient Picard scale, is the
    # relevant fixed datum.
    coefficient = Fraction(2 * -7, -6)
    scale_invariant = all(
        Fraction(2 * scale * -7, scale * -6) == coefficient
        for scale in range(1, 65)
    )

    result = {
        "x3_forced_square_class": 3,
        "x3_three_is_rational_square": rational_square(x3_ratio),
        "h21_required_lorentzian_rank": 22,
        "complex_k3_picard_ceiling": 20,
        "h21_complex_k3_rank_gate_passes": 22 <= 20,
        "fano_diagonal_point_determinant": fano_det,
        "fano_configuration_possible_over_characteristic_zero": fano_det == 0,
        "fano_determinant_mod_two": fano_det % 2,
        "native_reflection_coefficient_similarity_invariant": scale_invariant,
    }
    if result != {
        "x3_forced_square_class": 3,
        "x3_three_is_rational_square": False,
        "h21_required_lorentzian_rank": 22,
        "complex_k3_picard_ceiling": 20,
        "h21_complex_k3_rank_gate_passes": False,
        "fano_diagonal_point_determinant": -2,
        "fano_configuration_possible_over_characteristic_zero": False,
        "fano_determinant_mod_two": 0,
        "native_reflection_coefficient_similarity_invariant": True,
    }:
        raise AssertionError("hostile arithmetic replay changed")
    return result


def build_payload() -> dict[str, Any]:
    actual_pins = validate_pins()
    documents = {
        name: document_scope_audit(name, relative)
        for name, (relative, _expected) in DOCUMENT_PINS.items()
    }
    failed_documents = [
        name for name, audit in documents.items() if not audit["all_checks_pass"]
    ]
    if failed_documents:
        failed_checks = {
            name: [
                check
                for check, passed in documents[name]["checks"].items()
                if not passed
            ]
            for name in failed_documents
        }
        raise AssertionError(f"synthesis document scope failures: {failed_checks}")

    evidence = validate_evidence_interfaces()
    owner_archive = validate_owner_archive(actual_pins)
    replay = hostile_arithmetic_replay()

    scope = {
        "finite_pair_theorem": (
            "finite native pair set plus discrete- or arithmetic-reflection "
            "maximality implies finitely many unmarked normalized maximal chambers"
        ),
        "finite_pair_step_imports_nikulin_1982_classification": True,
        "finite_pair_step_is_ag_only": False,
        "global_absw_proved": False,
        "field_finiteness_requires_exact_admissible_full_trace_order_module_k3_return": True,
        "projected_ambient_host_finiteness_is_true_in_general": False,
        "projected_ambient_host_finiteness_needed_for_abstract_count": False,
        "x3_is_one_counterexample_to_all_finite_volume_universality": True,
        "x3_obstructs_every_finite_volume_input": False,
        "x3_proves_a_maximal_only_no_go": False,
        "x3_absolute_reflection_maximality_proved": False,
        "x3_mixed_divisor_complex_k3_positive": True,
        "h21_absolute_reflection_maximality": True,
        "h21_current_complex_stack_empty": True,
        "h21_every_ag_carrier_empty": False,
        "h21_characteristic_two_mixed_curve_cremona_positive": True,
        "delta_246_projected_return_proved_for_one_packet": True,
        "every_exact_projected_complex_k3_return_has_n_at_most_19": True,
        "universal_mixed_divisor_action_compiler_proved": False,
    }
    return {
        "verdict": (
            "PASS_CONDITIONAL_CLASSIFICATION_IMPORTING_FIXED_PAIR_CLOSURE_"
            "CURRENT_COMPLEX_OUTER_STACK_NONEXHAUSTIVE"
        ),
        "input_sha256": actual_pins,
        "document_audits": documents,
        "evidence_interface_checks": evidence,
        "owner_archive_cross_checks": owner_archive,
        "hostile_arithmetic_replay": replay,
        "scope": scope,
        "all_checks_pass": True,
    }


def build_certificate() -> dict[str, Any]:
    payload = build_payload()
    return {
        "schema": SCHEMA,
        "payload": payload,
        "payload_sha256": digest(payload),
    }


def encoded_certificate() -> str:
    return json.dumps(build_certificate(), indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-output", action="store_true")
    args = parser.parse_args()

    encoded = encoded_certificate()
    if args.check_output:
        if not OUTPUT.is_file():
            raise FileNotFoundError(f"missing independent archive: {OUTPUT}")
        if OUTPUT.read_text(encoding="utf-8") != encoded:
            raise AssertionError("independent archive does not match live replay")
    if args.write:
        OUTPUT.write_text(encoded, encoding="utf-8")
    if not args.write and not args.check_output:
        print(encoded, end="")


if __name__ == "__main__":
    main()
