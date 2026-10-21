#!/usr/bin/env python3
"""Fraction-only common-escape check for the complete bounded mixed dictionary."""
import json
from pathlib import Path
import fresh_cubic7_dihedral_independent_audit_20261013 as f

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_cubic7_mixed_scalar_h1000_20261013.json'
OUTPUT=ROOT/'results/fresh_cubic7_mixed_scalar_common_escape_20261013.json'

def covector(g,r):
    return tuple(f.sum_field(f.mul(a,b) for a,b in zip(row,r) if a!=f.ZERO and b!=f.ZERO) for row in g)

def run():
    d=json.loads(SOURCE.read_text());base=json.loads((ROOT/d['source']).read_text())
    g=tuple(f.vecdecode(row) for row in base['gram']);old=tuple(f.vecdecode(row) for row in base['all_walls'])
    roots=tuple(f.vecdecode(row['vector']) for row in d['combined_replacement_dictionary'])
    assert len({f.ray_key(r) for r in (*old,*roots)})==len(old)+len(roots)
    normals=tuple(covector(g,r) for r in (*old,*roots));results=[];seen=set();common=None
    for case,row in d['cases'].items():
        v=f.vecdecode(row['geometry']['hyperideal_witness']['vector'])
        if v in seen:continue
        seen.add(v);q=f.dot(v,covector(g,v))
        assert q==f.decode(row['geometry']['hyperideal_witness']['norm']) and f.sign(q)>0 and f.sign(v[0])>0
        products=tuple(f.dot(v,n) for n in normals)
        failed=[i for i,p in enumerate(products) if f.sign(p)>0]
        results.append({'source_case':case,'failed_dictionary_inequalities':failed})
        if not failed:common=(case,v,q,products);break
    out={'source':str(SOURCE.relative_to(ROOT)),'original_walls':len(old),'dictionary_roots':len(roots),
        'tested_saved_escapes':results,'common_escape_found':common is not None}
    if common is not None:
        case,v,q,products=common
        inside=f.vecdecode(base['strict_interior'])
        while any(f.sign(f.dot(inside,n))>=0 for n in normals):inside=(f.ONE,*(f.scale(z,f.F(1,2)) for z in inside[1:]))
        assert f.sign(f.dot(inside,covector(g,inside)))<0 and f.sign(inside[0])>0
        cap=json.loads((ROOT/'results/fresh_cubic7_dihedral_pell_cap_20261013.json').read_text());r=f.vecdecode(cap['root'])
        cap_pair=f.dot(v,covector(g,r));cap_fits=f.sign(cap_pair)<=0
        if cap_fits:
            assert f.sign(f.dot(inside,covector(g,r)))<0
        out.update(status='EXACT_COMMON_SPACELIKE_ESCAPE',source_case=case,vector=f.vecencode(v),norm=f.encoded(q),
            strict_interior=f.vecencode(inside),wall_pairings=f.vecencode(products),pell_cap_pairing=f.encoded(cap_pair),
            pell_cap_also_satisfied=cap_fits,
            scope='Every compatible subset of the complete saved bounded union dictionary retains this spacelike escape. The full union need not be Coxeter. This is not a full-slice or full-host nonreflectivity proof.')
    else:
        out.update(status='NO_SAVED_COMMON_ESCAPE',scope='Only the saved selected-cone escapes were tested. No impossibility or finite-volume conclusion follows for the full dictionary.')
    OUTPUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('vector','strict_interior','wall_pairings','tested_saved_escapes')}))
    return out

if __name__=='__main__':run()
