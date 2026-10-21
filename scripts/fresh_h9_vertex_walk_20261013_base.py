#!/usr/bin/env python3
"""Exact archive-only loaders for the successive vertex-walk chambers."""
from scripts.fresh_h9_gap_completion_20261013_base import *

PREFIX="data/enumerations/fresh_h9_vertex_walk_20261013"

def chamber(stage=0,case_index=None):
    if case_index is None:case_index=4 if stage==3 else 0
    G,I,h,ps,case=selected_coupled_case(0)
    if stage==0:return G,I,h,ps,case
    data=read(PREFIX+"_geometry.json")
    case=data["cases"][case_index if stage==1 else 0]
    bank=read(data["cap_source"])
    for n in case["dropped"]:ps.pop(n)
    for i in case["selected_caps"]:ps[f"walk{i}"]=cap_profile(bank["caps"][i],f"walk{i}")
    assert list(ps)==case["walls"]
    if stage==1:return G,I,h,ps,case
    data=read(PREFIX+"_second_geometry.json")
    case=data["cases"][case_index if stage==2 else 0]
    bank=read(data["cap_source"])
    for n in case["dropped"]:ps.pop(n)
    for i in case["selected_caps"]:ps[f"bridge{i}"]=cap_profile(bank["caps"][i],f"bridge{i}")
    assert list(ps)==case["walls"]
    if stage==2:return G,I,h,ps,case
    if stage==3:
        data=read(PREFIX+"_third_geometry.json");case=data["cases"][case_index]
        bank=read(data["cap_source"])
        for n in case["dropped"]:ps.pop(n)
        for i in case["selected_caps"]:ps[f"step{i}"]=cap_profile(bank["caps"][i],f"step{i}")
        assert list(ps)==case["walls"]
        return G,I,h,ps,case
    raise ValueError(stage)
