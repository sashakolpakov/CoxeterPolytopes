from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_gram_belyi_h10_lower_bound.py"
ARCHIVE = ROOT / "data/enumerations/gram_belyi_h10_lower_bound.json"
REPORT = ROOT / "results/gram_belyi_coxeter_dessin_theorem.md"


def load_module():
    spec = importlib.util.spec_from_file_location("gram_belyi_certifier", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_archive_is_exact_replay() -> None:
    module = load_module()
    actual = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert actual == module.build_certificate()
    assert actual["all_checks_pass"] is True


def test_symbolic_and_prime_lower_bounds() -> None:
    data = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert data["symbolic_identity"] == {
        "normalized_squared_gram_u": "(m - 1)**2*(m + 1)**2",
        "packet_coordinate": "lambda=u/(u+1)",
        "pair_c_q_2m": "4*(m - 1)*(m + 1)",
        "q_2m": "c + 2*m^2*a + 2*m*(m-1)*b",
        "q_2m_square": "-4",
    }
    rows = data["prime_rows"]
    assert len(rows) == 24
    for row in rows:
        p = row["prime"]
        assert row["path_copies_N"] == 2 * (p + 1)
        assert row["facet_count"] == 2 * p + 14
        assert row["normalized_squared_gram_u"] == (p * (p + 2)) ** 2
        assert row["p_adic_valuation_lambda"] == 2
        assert row["rodriguez_belyi_height_lower_bound"] == p


def test_report_scope_and_bytes() -> None:
    raw = REPORT.read_bytes()
    assert b"\r" not in raw and b"\x08" not in raw and b"\x0b" not in raw
    text = raw.decode("utf-8")
    assert "finite exhaustive marked Belyi-surface stack" in text
    assert "It does **not** yet prove that reflection maximality bounds this degree" in text
    assert "beta_GB(P)<=B" in text
    assert "ABSW" in text and "Prasad" in text

