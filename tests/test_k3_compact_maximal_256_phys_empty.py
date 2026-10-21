"""Regression tests for the maximal compact (2,5,6) K3 PHYS no-go."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.certify_k3_compact_maximal_256_phys_empty import (
    ARCHIVE,
    REPORT,
    Q5,
    build_certificate,
    digest,
    encoded_certificate,
    energy,
    golden_energy_shell,
    is_fraction_square,
    is_q5_square,
    raw_gram_certificate,
    square_class_census,
)


CERTIFIER = ROOT / "scripts/certify_k3_compact_maximal_256_phys_empty.py"


def test_exact_q5_arithmetic_and_square_decision() -> None:
    sqrt5 = Q5(Fraction(0), Fraction(1))
    phi = Q5(Fraction(1, 2), Fraction(1, 2))
    assert sqrt5 * sqrt5 == Q5(Fraction(5), Fraction(0))
    assert phi * phi == phi + Q5(Fraction(1), Fraction(0))
    assert phi.trace() == 1
    assert phi.embedding_sign() == 1
    assert phi.embedding_sign(conjugate=True) == -1
    assert is_fraction_square(Fraction(49, 81))
    assert not is_fraction_square(Fraction(5, 9))
    assert is_q5_square(phi * phi)
    assert is_q5_square(Q5(Fraction(3, 2), Fraction(-1, 2)))
    assert not is_q5_square(Q5(Fraction(3), Fraction(0)))


def test_raw_gram_is_admissible_compact_rank_six_packet() -> None:
    gram = raw_gram_certificate()
    assert gram["leading_minors"] == ["-6", "12", "6*phi"]
    assert gram["LDL_pivot_signs_distinguished"] == [-1, -1, 1]
    assert gram["LDL_pivot_signs_conjugate"] == [-1, -1, -1]
    assert gram["distinguished_inertia"] == [1, 2, 0]
    assert gram["conjugate_inertia"] == [0, 3, 0]
    assert gram["coxeter_orders"] == {
        "r1_r2": 2,
        "r2_r3": 5,
        "r3_r1": 6,
    }
    assert gram["order_module_preserved"] is True
    assert gram["cycle_order"] == "Z[phi^2]=Z[phi]=O_k"
    assert gram["compact_angle_sum"] == "1/2+1/5+1/6=13/15<1"
    assert gram["field_degree"] == 2
    assert gram["module_rank_over_field"] == 3
    assert gram["trace_rank"] == 6


def test_codifferent_trace_minus_two_shell_is_complete() -> None:
    shell = golden_energy_shell()
    assert [row["c"] for row in shell] == [-2, -1, 0, 1, 2]
    assert [row["field_norm"] for row in shell] == [
        "1/5",
        "4/5",
        "1",
        "4/5",
        "1/5",
    ]
    assert all(row["numerator_norm_c2_minus_5"] % 3 for row in shell)
    assert all(row["valuation_at_inert_3"] == 0 for row in shell)
    assert energy(0) == Q5(Fraction(-1), Fraction(0))


def test_all_twenty_five_required_square_classes_fail() -> None:
    census = square_class_census()
    assert census["tested_ordered_pairs"] == 25
    assert census["square_count"] == 0
    assert not any(row["is_square"] for row in census["rows"])


def test_inert_prime_parity_and_compiler_verdict() -> None:
    payload = build_certificate()
    coherence = payload["intrinsic_diagonal_PHYS_criterion"]
    assert coherence["invariant_under_facet_line_rescaling"] is True
    assert coherence["full_even_integral_module_also_required"] is True
    assert coherence["raw_maximality_controls_this_intersection"] is False
    assert coherence["triangle_256_pairwise_intersection_empty"] is True
    inert = payload["inert_prime_obstruction"]
    assert inert["values_at_0_1_2_mod_3"] == [2, 2, 1]
    assert inert["irreducible_mod_3"] is True
    assert inert["three_is_inert"] is True
    assert inert["all_energy_values_are_3_adic_units"] is True
    assert inert["facet_norm_ratio_h1_over_h3"] == 3
    assert inert["left_valuation_parity"] == 0
    assert inert["right_valuation_parity"] == 1
    assert inert["contradiction"] is True

    verdict = payload["compiler_verdict"]
    assert verdict["rank_gate_passes"] is True
    assert verdict["trace_rank"] == 6
    assert verdict["K3_complement_rank_if_PHYS_existed"] == 16
    assert verdict["discriminant_length_upper_bound"] == 6
    assert verdict["stable_embedding_inequality"] == "16>=6+2"
    assert verdict["EMB_would_be_automatic_after_PHYS"] is True
    assert verdict["PHYS"] == "EMPTY"
    assert verdict["EMB_reached"] is False
    assert verdict["ROOT_reached"] is False
    assert verdict["K3Comp_nonempty"] is False
    assert verdict["raw_maximality_implies_K3Comp_nonempty"] is False
    assert verdict["projected_root_or_larger_Picard_return_excluded"] is False
    assert verdict["alternative_AG_carrier_excluded"] is False


def test_maximality_input_and_sharp_positive_repair_are_scoped() -> None:
    payload = build_certificate()
    maximality = payload["published_maximality_input"]
    assert maximality["orientation_preserving_triangle_group"] == (
        "Delta(2,5,6)"
    )
    assert maximality["arithmetic"] is True
    assert maximality["finitely_maximal"] is True
    assert maximality["full_triangle_reflection_group_maximal"] is True
    assert "Maximality is cited" in maximality["classification_boundary"]

    repair = payload["sharp_positive_repair"]
    assert repair["PHYS_and_EMB_must_be_assumed_or_proved"] is True
    assert repair["zero_projection_roots_must_be_excluded"] is True
    assert repair["then_raw_maximality_removes_interior_ROOT_cuts"] is True
    assert repair["facet_witnesses_still_required"] is True
    assert repair["period_or_ample_gate_remaining_after_these_conditions"] is False


def test_report_states_the_counterexample_and_does_not_overclaim() -> None:
    text = REPORT.read_text(encoding="utf-8")
    for needle in (
        "raw-maximal compact arithmetic chamber with empty K3 `PHYS`",
        "(\\pi/2,\\pi/5,\\pi/6)",
        "absolutely maximal hyperbolic reflection group",
        "\\mathscr E_{-2}(\\mathcal O_k)",
        "\\bigcap_a{\\cal C}_a\\ne\\varnothing",
        "weakest square-class precondition",
        "v_{\\mathfrak p}(3z^2)=1+2v_{\\mathfrak p}(z)",
        "\\mathrm{PHYS}(P)=\\varnothing",
        "\\operatorname {K3Comp}(P)=\\varnothing",
        "16\\ge\\ell(A_N)+2",
        "`PHYS` is the only lattice-existence",
        "wall-selective native saturation",
        "projected-root return in a",
        "does not exclude",
        "does not make Property 4 automatic",
    ):
        assert needle in text


def test_archive_is_exact_replay() -> None:
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert archived["payload"] == build_certificate()
    assert archived["payload_sha256"] == digest(archived["payload"])
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()


def test_cli_check() -> None:
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check-output"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
