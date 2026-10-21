import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/independent_hostile_audit_silver_lambda22_h14.py"
ARTIFACT = (
    ROOT
    / "data/enumerations/silver_lambda22_h14_independent_hostile_audit.json"
)


def load_audit():
    spec = importlib.util.spec_from_file_location("silver_lambda22_h14_hostile", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_frozen_independent_audit_records_complete_simultaneous_packet() -> None:
    payload = json.loads(ARTIFACT.read_text())
    assert payload["verdict"] == "PASS"
    assert payload["gram"]["smith_invariants"] == [1] * 20 + [2, 6]
    assert payload["gram"]["rank_mod_2"] == 20
    assert payload["gram"]["rank_mod_3"] == 21
    assert payload["pari_count_only"]["method"].endswith(
        "qfminim(Q,b,m=1,default flag=0)"
    )
    packet = payload["H14_packet"]
    assert packet["ordinary_norm10_shell_count"] == 467_596_800
    assert packet["dyadic_S2_norm20_shell_count"] == 1_368_576
    assert packet["total_root_count"] == 468_965_376
    assert packet["q4_difference_squares"] == [10, 8, 6, 4, 0]
    assert payload["hyperideal_witness"]["strict"] is True
    conclusion = payload["conclusion"]
    assert conclusion["H14_prefix_finite_volume"] is False
    assert conclusion["eventual_chamber_status"] == "unresolved"
    assert conclusion["corrected_next_frontier"].startswith("H=36")


def test_standalone_arithmetic_and_modular_kernel_replay() -> None:
    audit = load_audit()
    gram = audit.parse_gram()
    assert gram.det() == 12
    for prime, rank, index in ((2, 20, 2**20), (3, 21, 3**21)):
        basis = audit.kernel_basis_mod_prime(gram, prime)
        assert abs(int(basis.det())) == index
        assert all(int(value) % prime == 0 for value in gram * basis)
        assert 22 - len(gram.nullspace()) == 22  # positive-definite over Q
        assert rank in (20, 21)


def test_packet_pairing_formulas_are_strict_at_proved_dot_bounds() -> None:
    audit = load_audit()
    ordinary = (
        -audit.A * audit.K2.make(6, 4) ** 2
        + audit.K2.make(3, 2) * audit.K2.make(8)
    )
    dyadic = (
        -audit.A * audit.K2.make(8, 6) ** 2
        + audit.K2.make(3, 2) * audit.K2.make(14)
    )
    cross = (
        -audit.A * audit.K2.make(6, 4) * audit.K2.make(8, 6)
        + audit.K2.make(3, 2) * audit.K2.make(13)
    )
    assert ordinary == audit.K2.make(-4, -4)
    assert dyadic == audit.K2.make(-14, -12)
    assert cross == audit.K2.make(-1, -2)
    assert all(value.sign() < 0 for value in (ordinary, dyadic, cross))
