#!/usr/bin/env python3
"""Producer-free final checkpoint, linked to the full exact walk replay."""
import argparse
import hashlib
import json
from scripts.fresh_h9_vertex_walk_20261013_verify import point_check
from scripts.fresh_h9_vertex_walk_20261013_base import *

def run():
    source=PREFIX+"_final.json";data=read(source)
    audit=read(PREFIX+"_verify.json")
    assert audit["status"]=="PASS_EXACT_TARGETED_VERTEX_WALK"
    for p,digest in audit["source_sha256"].items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==digest
    G,I,h,ps,case=chamber(3)
    assert data["source_geometry"]==PREFIX+"_third_geometry.json"and data["source_case"]==4
    assert data["wall_names"]==list(ps)and data["facet_count"]==len(ps)==73
    assert data["first_escape"]==case["geometry"]["escapes"][0]
    assert strict_escape(tuple(parse(x)for x in data["first_escape"]["certificate"]["y"]),I,h,list(ps.values()))==data["first_escape"]["certificate"]
    assert not case["affine_census"]["affine_circuits"]
    assert len(data["new_cap_references"])==len(data["all_four_actual_finite_vertices"])==4
    for ref,proof in zip(data["new_cap_references"],data["all_four_actual_finite_vertices"],strict=True):
        record=read(ref["source"])["caps"][ref["index"]]
        assert proof["point_field_from_cap"]==ref["name"]
        row=cap_profile(record,ref["name"])
        other=ps[ref["name"]]
        assert key(row["field"])==key(other["field"])
        assert all(row[k]==other[k]for k in("name","c","p","q"))
        point_check(proof,record,ps,I,h)
    return {"status":"PASS_EXACT_FINAL73_FOUR_ACTUAL_FINITE_VERTICES_AND_ESCAPE",
            "source":source,"source_sha256":hashlib.sha256((ROOT/source).read_bytes()).hexdigest(),
            "full_walk_audit":PREFIX+"_verify.json",
            "full_walk_audit_sha256":hashlib.sha256((ROOT/(PREFIX+"_verify.json")).read_bytes()).hexdigest(),
            "facets":73,"actual_finite_vertices":4,"vertex_wall_inequalities":292,
            "future_strict_spacelike_escape":True,
            "scope":"Four actual finite vertices survive simultaneously in the final73 chamber. Full all-pair and affine checks are anchored to the separately replayed walk audit. The strict escape proves infinite volume."}
if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--verify",action="store_true");args=parser.parse_args()
    result=run();out=ROOT/(PREFIX+"_final_verify.json")
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)
