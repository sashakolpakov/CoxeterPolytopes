#!/usr/bin/env python3
"""Equal dihedral block plus scalar completion, including wall replacement."""
import argparse
import json
from pathlib import Path
import time
from fresh_cubic7_dihedral_h9_20261013 import K,T,Z,O,C,D,BETA,decode_matrix,enc,mv,dot,pair,ray_key,finite_angle_counts,geometry
from fresh_cubic7_asymmetric_cuts_20261013 import compatible
from fresh_cubic7_target_alphabet_20261013 import thin_alphabet

ROOT=Path(__file__).resolve().parents[1]

def finish_cone(old,candidates,g,center,cap=None):
    walls=list(old);cov=list(mv(g,r) for r in walls);norms=list(dot(r,n) for r,n in zip(walls,cov))
    keys={ray_key(r) for r in walls};selected=[];rejected=[]
    for i,row in enumerate(candidates):
        r=tuple(map(K.from_json,row['vector']));q=K.from_json(row['norm'])
        if ray_key(r) in keys:continue
        if compatible(r,q,cov,norms):
            walls.append(r);cov.append(mv(g,r));norms.append(q);keys.add(ray_key(r));selected.append(i)
        else:rejected.append(i)
    cap_status='NOT_PROCESSED_SEPARATELY'
    if cap is not None:
        r=tuple(map(K.from_json,cap['root']));q=K.from_json(cap['norm'])
        assert pair(r,g,r)==q and all((2*z/q).is_integral for z in mv(g,r))
        if compatible(r,q,cov,norms):walls.append(r);cov.append(mv(g,r));norms.append(q);cap_status='APPENDED'
        else:cap_status='REJECTED_BY_EXACT_CROSS_ANGLE_GATE'
    while any(dot(center,n).sign()>=0 for n in cov):center=(O,*(x/2 for x in center[1:]))
    assert pair(center,g,center).sign()<0
    return {'selected_indices':selected,'cross_incompatible_indices':rejected,'wall_count':len(walls),
        'walls':[enc(r) for r in walls],'strict_interior':enc(center),'pell_cap_status':cap_status,
        'finite_angle_counts':finite_angle_counts(walls,g),'geometry':geometry(walls,g,9,inside=center)}

def run(height,output):
    started=time.monotonic();source=ROOT/'results/fresh_cubic7_dihedral_h9_2theta_h30_20261013.json'
    original=json.loads(source.read_text());g=decode_matrix(original['gram']);old=decode_matrix(original['all_walls'])
    ocov=tuple(mv(g,r) for r in old);onorms=tuple(dot(r,n) for r,n in zip(old,ocov));alpha=2*T
    sources={'retained158':'results/fresh_cubic7_two_block_h1000_20261013.json',
             'replacement122':'results/fresh_cubic7_two_block_replacement_h1000_20261013.json'}
    bases={key:json.loads((ROOT/path).read_text()) for key,path in sources.items()}
    bwalls={key:decode_matrix(d['walls']) for key,d in bases.items()}
    bcov={key:tuple(mv(g,r) for r in walls) for key,walls in bwalls.items()}
    bnorms={key:tuple(dot(r,n) for r,n in zip(bwalls[key],normals)) for key,normals in bcov.items()}
    escapes={key:tuple(map(K.from_json,d['geometry']['hyperideal_witness']['vector'])) for key,d in bases.items()}
    keys={ray_key(r) for r in old};rows=[];census=[]
    for q in (O,K(2)):
        times=tuple(x for x in thin_alphabet((height*q/alpha,-q/alpha,-q/alpha)) if x.sign()>0)
        aa=tuple(a for a in thin_alphabet((q*(height+1)/2,q/2,q/2)) if a.sign()>=0)
        zz=tuple(z for z in thin_alphabet((q*(height+1),q,q)) if z.sign()>=0)
        squares={z*z:z for z in zz};count={'norm':q.to_json(),'times':enc(times),'a_alphabet':enc(aa),'z_alphabet':enc(zz),
            'norm_solutions':0,'original25_compatible_new_roots':0,'cut_common_previous_escape':0,
            'compatible_retained158':0,'compatible_replacement122':0}
        for x in times:
            target=q+alpha*x*x
            for a in aa:
                z=squares.get(target-2*a*a)
                if z is None:continue
                count['norm_solutions']+=1
                for block in range(4):
                    r=[x]+[Z]*9;r[1+2*block]=r[2+2*block]=a/D;r[9]=z;r=tuple(r)
                    assert pair(r,g,r)==q and all(t.is_integral for t in r)
                    assert all((2*t/q).is_integral for t in mv(g,r))
                    if not compatible(r,q,ocov,onorms):continue
                    key=ray_key(r)
                    if key in keys:continue
                    keys.add(key);count['original25_compatible_new_roots']+=1
                    flags={name:compatible(r,q,bcov[name],bnorms[name]) for name in bases}
                    pp={name:pair(r,g,v) for name,v in escapes.items()}
                    for name,flag in flags.items():count['compatible_'+name]+=int(flag)
                    assert len(set(pp.values()))==1
                    count['cut_common_previous_escape']+=int(next(iter(pp.values())).sign()>0)
                    rows.append({'vector':enc(r),'norm':q.to_json(),'height':(alpha*x*x/q).to_json(),
                        'parameters':enc((x,a,z)),'block':block,'compatible_with_previous_cones':flags,
                        'previous_escape_pairings':{name:p.to_json() for name,p in pp.items()}})
        census.append(count);print(json.dumps({k:v for k,v in count.items() if k not in ('times','a_alphabet','z_alphabet')}),flush=True)
    rows.sort(key=lambda r:K.from_json(r['height']))
    cases={}
    for name,base in bases.items():
        candidates=[r for r in rows if r['compatible_with_previous_cones'][name]]
        result=finish_cone(bwalls[name],candidates,g,tuple(map(K.from_json,base['strict_interior'])))
        result['candidate_indices']=[i for i,r in enumerate(rows) if r['compatible_with_previous_cones'][name]]
        cases[name]=result
        print(json.dumps({'case':name,'walls':result['wall_count'],'geometry':result['geometry']['status']}),flush=True)
    combined=[];seen={ray_key(r) for r in old}
    for origin,data in (('two_block',bases['replacement122']['compatible_roots']),('equal_block_scalar',rows)):
        for i,row in enumerate(data):
            r=tuple(map(K.from_json,row['vector']));key=ray_key(r)
            if key in seen:continue
            seen.add(key);combined.append({**row,'origin':origin,'origin_index':i})
    combined.sort(key=lambda r:K.from_json(r['height']))
    cap=json.loads((ROOT/'results/fresh_cubic7_dihedral_pell_cap_20261013.json').read_text())
    cases['replacement_from25']=finish_cone(old,combined,g,tuple(map(K.from_json,original['strict_interior'])),cap)
    result=cases['replacement_from25']
    print(json.dumps({'case':'replacement_from25','walls':result['wall_count'],'geometry':result['geometry']['status'],'cap':result['pell_cap_status']}),flush=True)
    out={'status':'COMPLETE_BOUNDED_EQUAL_BLOCK_SCALAR_SEARCH','source':str(source.relative_to(ROOT)),
        'previous_cone_sources':sources,'height_bound':height,'primitive_norm_pool':enc((O,K(2))),
        'census':census,'compatible_roots':rows,'combined_replacement_dictionary':combined,'cases':cases,
        'scope':'Complete original25-compatible equal-block-plus-scalar primitive normals through the given height, all4 block embeddings. Conditional completions of the158/122 cones and a height-ordered replacement from25 using the two-block union. The Pell cap is separate and is retained only if it passes all cross-angle gates. No full-height or full-lattice completeness claim.',
        'seconds':round(time.monotonic()-started,3)}
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'new_roots':len(rows),'seconds':out['seconds']}),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--height',type=int,default=1000);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.height,a.output)
