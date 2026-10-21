#!/usr/bin/env python3
"""Freeze the best73 chamber and all four retained new finite vertices."""
import json
from scripts.fresh_h9_vertex_walk_20261013_base import *
from scripts.fresh_h9_vertex_walk_20261013_caps import finite_vertex

def run():
    G,I,h,selected,case=chamber(3)
    references=[
        {"source":PREFIX+"_caps.json","index":0,"name":"walk0"},
        {"source":PREFIX+"_second_caps.json","index":0,"name":"bridge0"},
        {"source":PREFIX+"_second_caps.json","index":2,"name":"bridge2"},
        {"source":PREFIX+"_third_caps.json","index":0,"name":"step0"}]
    vertices=[]
    for ref in references:
        record=read(ref["source"])["caps"][ref["index"]]
        proof=finite_vertex(selected[ref["name"]],record,I,h,selected)
        vertices.append(proof)
    result={"status":"EXACT_FOUR_REPLACEMENT_73_WALL_CHAMBER_WITH_ESCAPE",
            "source_geometry":PREFIX+"_third_geometry.json","source_case":4,
            "wall_names":list(selected),"facet_count":len(selected),
            "new_cap_references":references,"all_four_actual_finite_vertices":vertices,
            "removed_original_walls":["vertex18","old456","coupled6","vertex22"],
            "first_escape":case["geometry"]["escapes"][0],
            "scope":"Final selected coherent chamber with four optional-wall replacements and four surviving exact actual finite vertices. A strict spacelike escape proves infinite volume; completion beyond this construction is not excluded."}
    (ROOT/(PREFIX+"_final.json")).write_text(json.dumps(result,indent=2)+"\n")
    print({"status":result["status"],"facets":len(selected),"actual_finite_vertices":len(vertices)},flush=True)
if __name__=="__main__":run()
