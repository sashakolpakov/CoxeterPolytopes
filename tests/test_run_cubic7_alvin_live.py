from io import BytesIO

from scripts.run_cubic7_alvin_live import CASES, _canonical_qf, _pump, _split_qf


class _Read1Only:
    def __init__(self, chunks: tuple[bytes, ...]):
        self.chunks = list(chunks)

    def read(self, _size: int) -> bytes:  # pragma: no cover - failure path
        raise AssertionError("_pump used blocking read instead of read1")

    def read1(self, _size: int) -> bytes:
        return self.chunks.pop(0)


def test_pump_prefers_available_pipe_bytes_and_preserves_the_stream():
    source = _Read1Only((b"first\n", b"second\n", b""))
    archive = BytesIO()
    terminal = BytesIO()

    _pump(source, archive, terminal)

    assert archive.getvalue() == b"first\nsecond\n"
    assert terminal.getvalue() == archive.getvalue()


def test_high_dimensional_cases_have_exact_hyperbolic_rank_and_diagonal_shape():
    expected = {
        "ramified-7-h22": (22, "[-3,-2,-2]", (2, -1, 0)),
        "inert-2-h22": (22, "[-2,-2,-2]", (2, 0, 0)),
        "ramified-7-h29": (29, "[-3,-2,-2]", (2, -1, 0)),
        "inert-2-h29": (29, "[-2,-2,-2]", (2, 0, 0)),
    }

    for name, (dimension, terminal, canonical_terminal) in expected.items():
        case = CASES[name]
        entries = _split_qf(case["qf"])
        canonical = _canonical_qf(case["qf"])

        assert case["hyperbolic_dimension"] == dimension
        assert case["ambient_rank"] == dimension + 1
        assert len(entries) == len(canonical) == dimension + 1
        assert entries[0] == "[-1,0,0]"
        assert entries[1:-1] == ("1",) * (dimension - 1)
        assert entries[-1] == terminal
        assert canonical[0] == (0, -1, 0)
        assert canonical[1:-1] == ((1, 0, 0),) * (dimension - 1)
        assert canonical[-1] == canonical_terminal


def test_high_dimensional_qf_strings_have_21_and_28_unit_slots():
    assert _split_qf(CASES["ramified-7-h22"]["qf"])[1:-1].count("1") == 21
    assert _split_qf(CASES["inert-2-h22"]["qf"])[1:-1].count("1") == 21
    assert _split_qf(CASES["ramified-7-h29"]["qf"])[1:-1].count("1") == 28
    assert _split_qf(CASES["inert-2-h29"]["qf"])[1:-1].count("1") == 28
