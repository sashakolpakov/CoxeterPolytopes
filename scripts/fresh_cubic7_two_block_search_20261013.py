#!/usr/bin/env python3
"""Complete bounded two-block cubic-seven search retaining the 73-wall cone."""
import argparse
from itertools import permutations
import json
from pathlib import Path
import time
from fresh_cubic7_dihedral_h9_20261013 import K,T,Z,O,C,D,BETA,decode_matrix,enc,mv,dot,pair,ray_key,finite_angle_counts,geometry
from fresh_cubic7_asymmetric_cuts_20261013 import compatible
from fresh_cubic7_target_alphabet_20261013 import thin_alphabet

ROOT=Path(__file__).resolve().parents[1]

def append_pell_cap(walls,g):
    source=ROOT/'results/fresh_cubic7_dihedral_pell_cap_20261013.json'
    saved=json.loads(source.read_text());r=tuple(map(K.from_json,saved['root']));q=K.from_json(saved['norm'])
    assert q==O and pair(r,g,r)==q and all(z.is_integral for z in r)
    assert all((2*z/q).is_integral for z in mv(g,r))
    cov=tuple(mv(g,w) for w in walls);norms=tuple(dot(w,n) for w,n in zip(walls,cov))
    assert compatible(r,q,cov,norms)
    assert r not in walls
    return [*walls,r],{'source':str(source.relative_to(ROOT)),'vector':enc(r),'norm':q.to_json(),
        'height':saved['height'],'all_previous_walls_coxeter_compatible':True}

def run(height,output):
    started=time.monotonic();source=ROOT/'results/fresh_cubic7_dihedral_h9_2theta_h30_20261013.json'
    basepath=ROOT/'results/fresh_cubic7_asymmetric_cuts_h1000_20261013.json'
    original=json.loads(source.read_text());base=json.loads(basepath.read_text())
    g=decode_matrix(original['gram']);old25=decode_matrix(original['all_walls']);old=decode_matrix(base['walls'])
    assert len(old25)==25 and len(old)==73 and old[:25]==old25
    cov=tuple(mv(g,r) for r in old);norms=tuple(dot(r,n) for r,n in zip(old,cov))
    v=tuple(map(K.from_json,base['geometry']['hyperideal_witness']['vector']))
    center=tuple(map(K.from_json,base['strict_interior']));alpha=2*T
    rows=[];keys={ray_key(r) for r in old};basekeys=set(keys)
    original_rows=[];original_keys={ray_key(r) for r in old25};census=[]
    for q in (K(2),2*BETA):
        times=tuple(x for x in thin_alphabet((height*q/alpha,-q/alpha,-q/alpha))
                    if x.sign()>0 and (2*alpha*x/q).is_integral)
        weights=tuple(w for w in thin_alphabet((2*q*(height+1),2*q,2*q))
                      if w.sign()>=0 and (2*w/q).is_integral)
        blocks=[];buckets={}
        for i,u in enumerate(weights):
            for w in weights[i:]:
                y=((2*u+C*w)/BETA,(C*u+2*w)/BETA)
                if not all(z.is_integral for z in y):continue
                energy=u*y[0]+w*y[1]
                if (q*(height+1)-energy).sign()<0 or any((q-energy).sign(e)<0 for e in (1,2)):continue
                entry=(energy,y,(u,w));blocks.append(entry);buckets.setdefault(energy,[]).append(entry)
        assert any(e==Z for e,_,_ in blocks) and any(u==w and u!=Z for _,_,(u,w) in blocks)
        count={'norm':q.to_json(),'time_controls':enc(times),'weights':enc(weights),'blocks':len(blocks),
            'ordered_norm_solutions':0,'embedded_norm_solutions':0,'compatible_original25':0,
            'compatible_base73':0,'new_roots':0,'old_escape_cutters':0}
        for x in times:
            target=q+alpha*x*x
            for e1,y1,uw1 in blocks:
                for e2,y2,uw2 in buckets.get(target-e1,()):
                    count['ordered_norm_solutions']+=1
                    for first,second in permutations(range(4),2):
                        r=[x]+[Z]*9;r[1+2*first:3+2*first]=y1;r[1+2*second:3+2*second]=y2;r=tuple(r)
                        count['embedded_norm_solutions']+=1
                        assert pair(r,g,r)==q and all(z.is_integral for z in r)
                        assert all((2*z/q).is_integral for z in mv(g,r))
                        if not compatible(r,q,cov[:25],norms[:25]):continue
                        count['compatible_original25']+=1
                        key=ray_key(r);p=pair(r,g,v)
                        record={'vector':enc(r),'norm':q.to_json(),'height':(alpha*x*x/q).to_json(),
                            'previous_escape_pairing':p.to_json(),'support_blocks':[first,second],
                            'block_weights':[enc(uw1),enc(uw2)],'block_energies':enc((e1,e2))}
                        passes73=compatible(r,q,cov[25:],norms[25:])
                        if key not in original_keys:
                            original_keys.add(key);original_rows.append({**record,'already_in_base73':key in basekeys,
                                'passes_retained73_new_wall_gate':passes73})
                        if not passes73:continue
                        count['compatible_base73']+=1
                        if key in keys:continue
                        keys.add(key);count['new_roots']+=1;count['old_escape_cutters']+=int(p.sign()>0)
                        rows.append(record)
        census.append(count)
        print(json.dumps({k:value for k,value in count.items() if k not in ('time_controls','weights')}),flush=True)
    rows.sort(key=lambda r:K.from_json(r['height']))
    walls=list(old);wcov=list(cov);wnorms=list(norms);selected=[];rejected=[]
    for i,row in enumerate(rows):
        r=tuple(map(K.from_json,row['vector']));q=K.from_json(row['norm'])
        if compatible(r,q,wcov,wnorms):
            walls.append(r);wcov.append(mv(g,r));wnorms.append(q);selected.append(i)
        else:rejected.append(i)
    walls,cap=append_pell_cap(walls,g);wcov=tuple(mv(g,r) for r in walls)
    while any(dot(center,n).sign()>=0 for n in wcov):center=(O,*(x/2 for x in center[1:]))
    assert pair(center,g,center).sign()<0
    out={'status':'BOUNDED_TWO_BLOCK_SEARCH_COMPLETE','original_source':str(source.relative_to(ROOT)),
        'base_source':str(basepath.relative_to(ROOT)),'height_bound':height,'census':census,
        'compatible_roots':rows,'selected_indices':selected,'cross_incompatible_indices':rejected,
        'additional_pell_cap':cap,
        'wall_count':len(walls),'walls':[enc(r) for r in walls],'strict_interior':enc(center),
        'finite_angle_counts':finite_angle_counts(walls,g),'geometry':geometry(walls,g,9,inside=center),
        'scope':'All positive-height primitive reflections through the stated height supported on at most two dihedral blocks, scalar zero, and compatible with the retained73-wall cone. Each block has u<=w, including equality and zero, as forced by the original complete zero-height chamber. Ordered embedding counters include duplicates; saved roots are deduplicated. The separate Pell cap lies outside the bounded two-block census. Cross-pair selection is greedy, not an ambient Vinberg completeness claim.',
        'seconds':round(time.monotonic()-started,3)}
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'walls':len(walls),'compatible':len(rows),'selected':len(selected),'geometry':out['geometry']['status'],'seconds':out['seconds']}),flush=True)
    # Reconsider high caps instead of locking in the previous73-wall choice.
    original_rows.sort(key=lambda r:K.from_json(r['height']))
    replacement=list(old25);rcov=list(cov[:25]);rnorm=list(norms[:25]);chosen=[];rejected=[]
    for i,row in enumerate(original_rows):
        r=tuple(map(K.from_json,row['vector']));q=K.from_json(row['norm'])
        if compatible(r,q,rcov,rnorm):replacement.append(r);rcov.append(mv(g,r));rnorm.append(q);chosen.append(i)
        else:rejected.append(i)
    replacement,cap=append_pell_cap(replacement,g);rcenter=tuple(map(K.from_json,original['strict_interior']))
    while any(pair(rcenter,g,r).sign()>=0 for r in replacement):rcenter=(O,*(x/2 for x in rcenter[1:]))
    assert pair(rcenter,g,rcenter).sign()<0
    escaped_new=[i for i,row in enumerate(original_rows) if not row['already_in_base73'] and
        not row['passes_retained73_new_wall_gate'] and K.from_json(row['previous_escape_pairing']).sign()>0]
    rep={'status':'BOUNDED_TWO_BLOCK_REPLACEMENT_COMPLETE','original_source':str(source.relative_to(ROOT)),
        'replaced_base_source':str(basepath.relative_to(ROOT)),'enumeration_source':str(output.resolve().relative_to(ROOT)),
        'height_bound':height,'compatible_roots':original_rows,'selected_indices':chosen,'cross_incompatible_indices':rejected,
        'new_roots_conflicting_with_old73_and_cutting_its_escape':escaped_new,'additional_pell_cap':cap,
        'wall_count':len(replacement),'walls':[enc(r) for r in replacement],'strict_interior':enc(rcenter),
        'finite_angle_counts':finite_angle_counts(replacement,g),'geometry':geometry(replacement,g,9,inside=rcenter),
        'scope':'Complete original25-compatible two-block dictionary at the unchanged height bound, including equal/zero blocks. Height-ordered cross-compatible replacement cone plus the separately archived norm-one Pell cap. No full-height or full-lattice completion claim.'}
    replacement_output=output.with_name(output.name.replace('two_block_h','two_block_replacement_h'))
    replacement_output.write_text(json.dumps(rep,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'replacement_walls':len(replacement),'dictionary':len(original_rows),'selected':len(chosen),
        'conflicting_new_escape_cutters':len(escaped_new),'geometry':rep['geometry']['status']}),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--height',type=int,default=1000);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.height,a.output)
