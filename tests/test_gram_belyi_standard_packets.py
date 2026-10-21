from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_gram_belyi_standard_packets.py"
ARCHIVE = ROOT / "data/enumerations/gram_belyi_standard_packets.json"
REPORT = ROOT / "results/gram_belyi_standard_packet_calibration.md"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "gram_belyi_standard_packet_certifier", SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rows_by_relation() -> dict[str, dict]:
    data = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    return {row["relation"]: row for row in data["rows"]}


def test_archive_is_exact_replay() -> None:
    module = load_module()
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert archived == module.build_certificate()
    assert archived["all_checks_pass"] is True


def test_all_exact_rows_and_orthogonal_exception() -> None:
    rows = rows_by_relation()
    assert set(rows) == {
        "m=2",
        "m=3",
        "m=4",
        "m=5",
        "m=6",
        "parallel",
        "dotted u=2",
    }
    assert rows["m=2"]["packet_lambda"] is None
    assert rows["m=2"]["gram_belyi_contribution"] == 2
    assert rows["m=2"]["auxiliary_packet"]["h_Bel"] == 2
    assert {relation: rows[relation]["h_Bel"] for relation in ("m=3", "m=4", "m=6", "parallel", "dotted u=2")} == {
        "m=3": 5,
        "m=4": 3,
        "m=6": 7,
        "parallel": 2,
        "dotted u=2": 3,
    }
    for relation in ("m=3", "m=4", "m=6", "parallel", "dotted u=2"):
        row = rows[relation]
        assert row["minimum_status"] == "proved_exact"
        assert row["lower_bound"]["bound"] == row["h_Bel"]
        assert row["upper_witness"]["degree"] == row["h_Bel"]
        assert row["lower_bound"]["p_adic_valuation"] != 0


def test_m5_interval_and_factored_composition() -> None:
    row = rows_by_relation()["m=5"]
    upper = 3_099_548_558_797_153_417_819_718_240_423_033_019_647_072
    assert row["minimum_status"] == "certified_interval_exact_minimum_open"
    assert row["proved_interval"] == [29, upper]
    assert row["lower_bound"]["newton_polygon_slopes"] == ["0/1", "1/1"]
    assert row["lower_bound"]["root_valuations"] == ["0/1", "-1/1"]

    witness = row["upper_witness"]
    coordinate = witness["final_negative_coordinate"]
    a, b = coordinate["a"], coordinate["b"]
    assert a == 3**15 * 5**14 * 29**16
    assert b == 2**32 * 11**30
    r = Fraction(-49, 435)
    t = -(r - 1) ** 15 * (15 * r + 1)
    assert t == Fraction(-b, a)
    assert witness["F_degree"] == a + b
    assert witness["total_degree"] == 2 * 16 * (a + b) == upper
    assert row["nonclaim"].startswith("no degree-29 witness")


def test_report_scope_and_bytes() -> None:
    raw = REPORT.read_bytes()
    assert b"\r" not in raw and b"\x08" not in raw and b"\x0b" not in raw
    text = raw.decode("utf-8")
    assert "five exact normalized-polynomial minima" in text
    assert "h_{\\rm Bel}(\\lambda)\\ge 29" in text
    assert "exact `m=5` minimum remains open" in text
    assert "no broad dessin enumeration" in text
