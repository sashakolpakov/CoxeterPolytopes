"""Certificates and adversarial controls for the independent H4-block audit."""
import json
from pathlib import Path
import pytest
from scripts.fresh_h4_blocks_independent_verify import (
    verify,h4_root_calibration,bugaenko_calibration,finite_cone_calibration,simplex_calibration,h9_isometry_calibrations,
)

ROOT=Path(__file__).resolve().parents[1]

@pytest.mark.parametrize('name',['h9_phi','h9_2phi','h9_sqrt5','h12_2phi','h13_2phi'])
def test_saved_escape(name):
    result=verify(ROOT/f'results/fresh_h4_blocks_{name}_h30_20261013.json')
    assert result['status']=='EXACT_SPACELIKE_ESCAPE_PASS'
    assert result['independently_complete_controls']>0

def test_h4_finite_root_closure():
    assert h4_root_calibration()['h4_exact_root_closure']==120

def test_h4_simplex():
    assert simplex_calibration(ROOT/'results/fresh_h4_blocks_independent_h4_phi_h30.json')['coxeter_chain']==[5,3,3,4]

def test_exact_bugaenko_integral_isometry_and_full_cone():
    assert bugaenko_calibration(ROOT)['bugaenko_exact_wall_gram_match']==11
    assert finite_cone_calibration(ROOT/'results/fresh_h4_blocks_h8_2phi_h30_20261013.json',41)['exact_finite_vertices']==41

def test_all_three_h9_transports_are_integral_isometries():
    assert len(h9_isometry_calibrations(ROOT))==3

@pytest.mark.parametrize('forgery',['norm','timelike_escape','missing_control','false_height'])
def test_rejects_forged_certificate(tmp_path,forgery):
    data=json.loads((ROOT/'results/fresh_h4_blocks_h9_phi_h30_20261013.json').read_text())
    if forgery=='norm':data['hyperideal_witness']['norm']={'a':'1','b':'0'}
    elif forgery=='timelike_escape':data['hyperideal_witness']['vector']=data['strict_interior']
    elif forgery=='missing_control':
        i=next(i for i,row in enumerate(data['levels']) if row['roots']==0)
        del data['levels'][i]
    else:data['levels'][0]['height']={'a':'0','b':'0'}
    bad=tmp_path/'bad.json';bad.write_text(json.dumps(data))
    with pytest.raises(AssertionError):verify(bad)
