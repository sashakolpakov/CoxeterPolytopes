import json
from hashlib import sha256
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
ARTIFACT=ROOT/"data/enumerations/silver_lambda22_h738_cutter.json"
CROSS=ROOT/"data/enumerations/silver_lambda22_h738_cross_rows.json"
FRONTIER=ROOT/"data/enumerations/silver_lambda22_h738_replacement_frontier.json"


def payload():return json.loads(ARTIFACT.read_text())


def test_explicit_h738_cutter_arithmetic():
 r=payload()["explicit_cutter"]
 assert r["Gram_d_g"]==[[4,-2],[-2,10]]
 assert r["coefficient_gram"]==[[738,522],[522,370]]
 assert r["norm"]["display"]=="2"
 assert r["height"]["display"]=="738 + 522*sqrt(2)"
 assert r["primitive"] is True and r["reflection_integral"] is True
 assert r["old_ray_pairing"]["display"]=="-9/20 + 7/20*sqrt(2)"
 assert r["old_ray_pairing"]["signs"][0]>0


def test_full_chronological_cache_has_23_real_rows_and_five_positives():
 c=json.loads(CROSS.read_text())
 assert c["gate_compatible_type_count"]==517
 assert c["nonnegative_row_count"]==len(c["rows"])==23
 ext=[e for r in c["rows"] for e in r["extensions"]]
 assert len(ext)==24
 assert c["positive_extension_count"]==5
 assert sum(e["pairing"]["signs"][0]>0 for e in ext)==5
 assert all(e["H14_g_wall_pairing"]["signs"][0]<=0 for e in ext)
 positives=[(r,e) for r in c["rows"] for e in r["extensions"] if e["pairing"]["signs"][0]>0]
 assert [(r["height"]["display"],r["gram"]) for r,e in positives]==[
  ("21 + 15*sqrt(2)",[60,36,45]),
  ("56 + 40*sqrt(2)",[124,52,80]),
  ("274 + 194*sqrt(2)",[580,260,388]),
  ("478 + 338*sqrt(2)",[916,500,676]),
  ("738 + 522*sqrt(2)",[738,370,522]),
 ]
 h478=positives[3][1]
 assert [h478[k] for k in ("d_dot_u","d_dot_v","g_dot_u","g_dot_v","Schur_A","Schur_B","Schur_N")]==[60,44,-18,-10,0,0,0]
 meta=payload()["chronological_acceptance"]["row_cache_artifact"]
 assert meta["path"]=="data/enumerations/silver_lambda22_h738_cross_rows.json"
 assert meta["sha256"]==sha256(CROSS.read_bytes()).hexdigest()


def test_replacement_frontier_is_complete_and_resolved():
 f=json.loads(FRONTIER.read_text())
 assert f["control_count"]==195
 assert f["retained_necessary_type_count"]==521
 assert len(f["relaxed_nonnegative_rows"])==7
 assert all(r["resolution"] for r in f["relaxed_nonnegative_rows"])
 p=payload()["replacement_hyperideal_ray"]
 assert p["scale"]=="32/99"
 assert p["norm"]["display"]=="13897/9801 - 1*sqrt(2)"
 assert p["norm"]["signs"][0]>0
 assert p["H164_cutter_pairing"]["display"]=="-100/99 + 40/99*sqrt(2)"
 assert p["H738_cutter_pairing"]["display"]=="-6/11 + 28/99*sqrt(2)"
 assert p["H164_cutter_pairing"]["signs"][0]<0
 assert p["H738_cutter_pairing"]["signs"][0]<0


def test_scope_is_bounded():
 p=payload()
 assert "bounded prefix" in p["scope"]
 assert "not Vinberg termination" in p["scope"]
 assert p["discovery"] is None
