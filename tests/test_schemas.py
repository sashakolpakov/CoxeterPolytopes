import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_json_artifacts_parse() -> None:
    paths = [
        ROOT / "data/known_examples/index.json",
        ROOT / "data/known_examples/bugaenko_compact_h8.json",
        ROOT / "data/known_examples/bugaenko_vinberg_transcript.json",
        ROOT / "data/candidates/index.json",
        ROOT / "data/candidates/coupled_discriminant11.json",
        ROOT / "data/candidates/coupled_discriminant19.json",
        ROOT / "data/candidates/h19_vinberg_prefix.json",
        ROOT / "data/candidates/shimada_disc11_rm_filter.json",
        ROOT / "data/trace_lattices/bugaenko_l8_trace.json",
        ROOT / "data/trace_lattices/l9_seed_trace.json",
        ROOT / "schemas/golden-element.schema.json",
        ROOT / "schemas/chamber-certificate.schema.json",
    ]
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            assert isinstance(json.load(handle), dict)
