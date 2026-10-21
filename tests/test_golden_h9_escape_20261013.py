"""Independent replay and corruption check for a mixed-norm volume witness."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import json
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=spec_from_file_location('mixed_escape_replay', ROOT/'scripts/verify_golden_h9_mixed_escape_20261013.py')
V=module_from_spec(spec);spec.loader.exec_module(V)

@pytest.mark.parametrize('stem',[
    'golden_h9_full_norms_h300_20261013',
    'golden_h9_odd_a1_full_norms_h300_20261013',
])
def test_independent_mixed_norm_escape(stem,monkeypatch):
    monkeypatch.chdir(ROOT)
    V.verify(ROOT/'results'/f'{stem}_escape.json')


def test_forged_escape_norm_is_rejected(tmp_path,monkeypatch):
    monkeypatch.chdir(ROOT)
    proof=json.loads((ROOT/'results/golden_h9_full_norms_h300_20261013_escape.json').read_text())
    proof['norm']={'a':'1','b':'0'}
    path=tmp_path/'forged.json';path.write_text(json.dumps(proof))
    with pytest.raises(AssertionError):V.verify(path)
