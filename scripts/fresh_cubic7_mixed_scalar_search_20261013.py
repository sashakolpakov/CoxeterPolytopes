#!/usr/bin/env python3
"""One free block, one equal block and scalar; bounded complete slice search."""
import argparse
from itertools import permutations
import json
from pathlib import Path
import time
from fresh_cubic7_dihedral_h9_20261013 import K,T,Z,O,C,D,BETA,decode_matrix,enc,mv,dot,pair,ray_key
from fresh_cubic7_asymmetric_cuts_20261013 import compatible
from fresh_cubic7_target_alphabet_20261013 import thin_alphabet
from fresh_cubic7_scalar_search_20261013 import finish_cone

ROOT=Path(__file__).resolve().parents[1]

def run(height,output):
    started=time.monotonic();source=ROOT/'results/fresh_cubic7_dihedral_h9_2theta_h30_20261013.json'
    original=json.loads(source.read_text());g=decode_matrix(original['gram']);old=decode_matrix(original['all_walls'])
    ocov=tuple(mv(g,r) for r in old);onorms=tuple(dot(r,n) for r,n in zip(old,ocov));alpha=2*T
    previous='results/fresh_cubic7_scalar_h1000_20261013.json';previous_data=json.loads((ROOT/previous).read_text())
    sources={name:{'file':previous,'case':case} for name,case in (
        ('retained182','retained158'),('retained146','replacement122'),('replacement146','replacement_from25'))}
    bases={name:previous_data['cases'][p['case']] for name,p in sources.items()}
    bwalls={name:decode_matrix(d['walls']) for name,d in bases.items()}
    bcov={name:tuple(mv(g,r) for r in walls) for name,walls in bwalls.items()}
    bnorms={name:tuple(dot(r,n) for r,n in zip(bwalls[name],normals)) for name,normals in bcov.items()}
    escapes={name:tuple(map(K.from_json,d['geometry']['hyperideal_witness']['vector'])) for name,d in bases.items()}
    keys={ray_key(r) for r in old};rows=[];census=[]
    for q in (O,K(2),BETA,2*BETA):
        times=tuple(x for x in thin_alphabet((height*q/alpha,-q/alpha,-q/alpha))
                    if x.sign()>0 and (2*alpha*x/q).is_integral)
        aa=tuple(a for a in thin_alphabet((q*(height+1)/2,q/2,q/2))
                 if a.sign()>=0 and (2*D*a/q).is_integral)
        zz=tuple(z for z in thin_alphabet((q*(height+1),q,q))
                 if z.sign()>=0 and (2*z/q).is_integral)
        weights=tuple(w for w in thin_alphabet((2*q*(height+1),2*q,2*q))
                      if w.sign()>=0 and (2*w/q).is_integral)
        combined={}
        for a in aa:
            for z in zz:combined.setdefault(2*a*a+z*z,[]).append((a,z))
        blocks=[]
        for i,u in enumerate(weights):
            for w in weights[i:]:
                y=((2*u+C*w)/BETA,(C*u+2*w)/BETA)
                if not all(t.is_integral for t in y):continue
                energy=u*y[0]+w*y[1]
                if (q*(height+1)-energy).sign()<0 or any((q-energy).sign(e)<0 for e in (1,2)):continue
                blocks.append((energy,y,(u,w)))
        count={'norm':q.to_json(),'times':enc(times),'weights':enc(weights),'a_alphabet':enc(aa),'z_alphabet':enc(zz),
            'free_blocks':len(blocks),'equal_scalar_pairs':len(aa)*len(zz),'norm_solutions':0,
            'original25_compatible_new_roots':0,'cut_common_previous_escape':0,
            **{'compatible_'+name:0 for name in bases}}
        for x in times:
            target=q+alpha*x*x
            for energy,y,uw in blocks:
                for a,z in combined.get(target-energy,()):
                    count['norm_solutions']+=1
                    for free,equal in permutations(range(4),2):
                        r=[x]+[Z]*9;r[1+2*free:3+2*free]=y;r[1+2*equal:3+2*equal]=(a/D,a/D);r[9]=z;r=tuple(r)
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
                            'parameters':enc((x,a,z)),'free_block':free,'equal_block':equal,'free_weights':enc(uw),
                            'free_energy':energy.to_json(),'compatible_with_previous_cones':flags,
                            'previous_escape_pairings':{name:p.to_json() for name,p in pp.items()}})
        census.append(count)
        print(json.dumps({k:v for k,v in count.items() if k not in ('times','weights','a_alphabet','z_alphabet')}),flush=True)
    rows.sort(key=lambda r:K.from_json(r['height']));cases={}
    for name,base in bases.items():
        candidates=[r for r in rows if r['compatible_with_previous_cones'][name]]
        result=finish_cone(bwalls[name],candidates,g,tuple(map(K.from_json,base['strict_interior'])))
        result['candidate_indices']=[i for i,r in enumerate(rows) if r['compatible_with_previous_cones'][name]]
        cases[name]=result
        print(json.dumps({'case':name,'walls':result['wall_count'],'geometry':result['geometry']['status']}),flush=True)
    union=[];seen={ray_key(r) for r in old}
    for origin,data in (('previous_union',previous_data['combined_replacement_dictionary']),('mixed_scalar',rows)):
        for i,row in enumerate(data):
            r=tuple(map(K.from_json,row['vector']));key=ray_key(r)
            if key in seen:continue
            seen.add(key);union.append({**row,'origin':origin,'origin_index':i})
    union.sort(key=lambda r:K.from_json(r['height']))
    cap=json.loads((ROOT/'results/fresh_cubic7_dihedral_pell_cap_20261013.json').read_text())
    cases['replacement_from25']=finish_cone(old,union,g,tuple(map(K.from_json,original['strict_interior'])),cap)
    r=cases['replacement_from25']
    print(json.dumps({'case':'replacement_from25','walls':r['wall_count'],'geometry':r['geometry']['status'],'cap':r['pell_cap_status']}),flush=True)
    out={'status':'COMPLETE_BOUNDED_MIXED_BLOCK_SCALAR_SEARCH','source':str(source.relative_to(ROOT)),
        'previous_cone_sources':sources,'height_bound':height,'primitive_norm_pool':enc((O,K(2),BETA,2*BETA)),
        'census':census,'compatible_roots':rows,'combined_replacement_dictionary':union,'cases':cases,
        'scope':'Complete original25-compatible primitive normals through the stated height with at most one free dihedral block, one equal-coordinate block, and scalar, all12 ordered placements. Equal/zero free blocks are included. Conditional completions retain each previous cone; union replacement also includes the earlier complete two-block dictionary. No complete ambient height or full-group reflectivity claim.',
        'seconds':round(time.monotonic()-started,3)}
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'dictionary':len(rows),'union_dictionary':len(union),'seconds':out['seconds']}),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--height',type=int,default=1000);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.height,a.output)
