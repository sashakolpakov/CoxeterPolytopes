import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_silver_odd_unimodular_h27_h29_h1000.py"
ARTIFACT = ROOT / "data/enumerations/silver_odd_unimodular_h27_h29_h1000.json"
SNAPSHOT_DIR = ROOT / "data/sources/nebe_sloane_odd_unimodular_h27_h29"

EXPECTED = {
    "T27": {
        "dimension": 27,
        "snapshot_sha256": "c130fe16cfbfa0a5e7ae94ce60ba1a6fe56c213f605a59c1d8df5a0a7f8f7db0",
        "gram_sha256": "6d5975cdb40df7a3a2c066072e18230ca49bb379bba2462ef9fc913c65590bc0",
        "norm_three": 1_640,
        "norm_four": 119_574,
    },
    "U28": {
        "dimension": 28,
        "snapshot_sha256": "ef857a546128dd039fd891163854fd3bad542edb7aaeb898f39c5f4680a552bb",
        "gram_sha256": "a1a0647ae2e7e24d5663869bf431b0d74756d99580ef0532d3abe1e4145713f3",
        "norm_three": 2_240,
        "norm_four": 98_280,
    },
    "dim29odd": {
        "dimension": 29,
        "snapshot_sha256": "e8c4fec6d2fddab8195261f8e8e9fc5eebf5ea532f9de4a4cf3035611b6a1ff8",
        "gram_sha256": "a96355961e7a65911a6154c9c0f5cbc1d9b7a2a481fc617477dc6b1a19d36d9a",
        "norm_three": 1_856,
        "norm_four": 93_786,
    },
}


def parse_gram_independently(path: Path, dimension: int) -> list[list[int]]:
    """Parse without importing or calling the production freezer."""
    source = path.read_text()
    tail = source[source.index('<a NAME="GRAM"><STRONG>GRAM</STRONG></a><br>') :]
    header = re.search(rf"\b{dimension}\s+{dimension}<br>\s*\n", tail)
    assert header is not None
    rows = []
    for line in tail[header.end() :].splitlines()[:dimension]:
        row = [int(value) for value in re.findall(r"-?\d+", line.partition("<br>")[0])]
        assert len(row) == dimension
        rows.append(row)
    assert len(rows) == dimension
    assert rows == [list(row) for row in zip(*rows)]
    return rows


def bareiss_determinant(matrix: list[list[int]]) -> int:
    """Fraction-free exact determinant, independent of PARI and production code."""
    work = [row[:] for row in matrix]
    previous = 1
    for pivot_index in range(len(work) - 1):
        pivot = work[pivot_index][pivot_index]
        assert pivot != 0  # positive-definite hosts have positive leading minors
        for row in range(pivot_index + 1, len(work)):
            for column in range(pivot_index + 1, len(work)):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                assert numerator % previous == 0
                work[row][column] = numerator // previous
        previous = pivot
    return work[-1][-1]


def catalogue_integer(source: str, anchor: str) -> int:
    match = re.search(
        rf'<a NAME="{anchor}"><STRONG>{anchor}</STRONG></a><br>\s*(\d+)<br>',
        source,
    )
    assert match is not None
    return int(match.group(1))


def test_primary_snapshots_and_exact_gram_determinants() -> None:
    payload = json.loads(ARTIFACT.read_text())
    hosts = {host["name"]: host for host in payload["primary_catalogue_hosts"]}
    assert set(hosts) == set(EXPECTED)
    for name, expected in EXPECTED.items():
        snapshot = SNAPSHOT_DIR / f"{name}.html"
        source = snapshot.read_text()
        gram = parse_gram_independently(snapshot, expected["dimension"])
        stream = "\n".join(" ".join(map(str, row)) for row in gram) + "\n"
        assert hashlib.sha256(snapshot.read_bytes()).hexdigest() == expected["snapshot_sha256"]
        assert hashlib.sha256(stream.encode()).hexdigest() == expected["gram_sha256"]
        assert catalogue_integer(source, "DIMENSION") == expected["dimension"]
        assert catalogue_integer(source, "MINIMAL_NORM") == 3
        assert catalogue_integer(source, "KISSING_NUMBER") == expected["norm_three"]
        assert bareiss_determinant(gram) == 1
        assert min(gram[index][index] for index in range(len(gram))) == 3
        host = hosts[name]
        assert host["pari"]["determinant"] == 1
        assert host["pari"]["minimum"] == 3
        assert host["pari"]["norm_three_shell_signed_count"] == expected["norm_three"]
        assert host["first_packet_wall_count"] == expected["norm_four"]


def test_universal_arithmetic_and_h1000_hyperideal_certificate() -> None:
    payload = json.loads(ARTIFACT.read_text())
    theorem = payload["universal_theorem"]
    assert theorem["local_norm_class_count"] == 6
    assert theorem["height_zero_wall_count"] == 0
    assert "S26" not in theorem["local_norm_class_proof"]
    assert "(6,3,4)" in theorem["first_packet_forcing"]
    assert theorem["control_count"] == 134
    assert theorem["distinct_height_count"] == 104
    assert theorem["nonempty_control_count"] == 30
    assert theorem["necessary_binary_type_count"] == 256
    assert theorem["nonempty_height_count"] == 24
    assert theorem["worst_row"] == {
        "height": "455 + 322*sqrt(2)",
        "norm": "2",
        "control": "16 + 11*sqrt(2)",
        "gram": [460, 226, 322],
        "maximum_pairing": "37 + 26*sqrt(2)",
        "margin": "-51/55 + 37/55*sqrt(2)",
    }
    assert theorem["witness"] == "z=(1,(43/55)d), d^2=3"
    assert theorem["witness_lorentz_norm"] == "8572/3025 - 2*sqrt(2)"
    assert theorem["all_wall_margins_strict"] is True
    assert theorem["strict_hyperideal"] is True
    assert payload["conclusion"] == {
        "new_compact_example": False,
        "new_finite_volume_example": False,
        "all_three_prefixes_have_infinite_volume": True,
        "vinberg_terminated": False,
    }


def test_exact_regeneration() -> None:
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--verify"],
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
        text=True,
        capture_output=True,
        timeout=300,
        check=True,
    )
    assert "PASS silver odd-unimodular H27--H29 H1000 family" in completed.stdout
