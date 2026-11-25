#!/usr/bin/env python3
"""Check the complete rooted-host/edge pushforwards under the shared watchdog."""
import argparse
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get("COXETER_CHECK_OWNER"): raise RuntimeError("Use the serial scientific watchdog")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True)
    parser.add_argument("--local-lifts",action="store_true")
    args=parser.parse_args()
    out=(ROOT/args.output).resolve()
    if not out.is_relative_to(ROOT/"results") or out.exists(): raise ValueError("Use a fresh results directory")
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import conic,edge_hosts,edge_compatibility,equation_closure
    paths=[ROOT/"research/global_formulation"/(name+".py") for name in
        ("edge_hosts","edge_compatibility","equation_closure","white_profiles","supporting_counts","vertices","moments","polynomial","conic","data")]
    before={str(p.relative_to(ROOT)):conic.digest_bytes(p.read_bytes()) for p in paths}
    out.mkdir()
    result=dict(status="EXACT_COMPLETE_EDGE_SCHEMA_CHECKS_PASSED",rooted_hosts=edge_hosts.self_checks(),
        incidence_pushforwards=edge_compatibility.self_checks(),accepted_bounds_changed=False,
        source_bound_geometric_review_complete=False)
    if args.local_lifts:
        from research.global_formulation import white_profiles
        result["local_level7_triangle"]=white_profiles.local_star_fixture(level=7)
        result["local_level7_quadrilateral"]=white_profiles.quadrilateral_star_fixture(level=7)
    if {str(p.relative_to(ROOT)):conic.digest_bytes(p.read_bytes()) for p in paths}!=before: raise ValueError("A check source changed")
    result["source_bindings"]=before
    conic.write_json(out/"result.json",result)
    print(__import__('json').dumps(conic.plain(result)),flush=True)


if __name__=="__main__": main()
