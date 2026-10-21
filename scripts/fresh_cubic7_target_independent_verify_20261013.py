#!/usr/bin/env python3
"""Independent Fraction-only replay of targeted cubic-seven dictionaries.

Imports only the independent standard-library field implementation, never the
producer. Saved candidates and geometry are replayed; the positive-height
candidate enumeration is not independently repeated here.
"""
import argparse
import copy
from functools import lru_cache
import json
from pathlib import Path
import fresh_cubic7_dihedral_independent_audit_20261013 as f

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/fresh_cubic7_dihedral_h9_2theta_h30_20261013.json'

def read(path):return json.loads((ROOT/path).read_text())
def read_case(source):
    return read(source['file'])['cases'][source['case']] if isinstance(source,dict) else read(source)
def walls_of(d):return tuple(f.vecdecode(r) for r in d.get('walls',d.get('all_walls')))
def escape_of(d):return f.vecdecode(d['geometry']['hyperideal_witness']['vector'])

class Auditor:
    def __init__(self):
        self.base=read(BASE);self.g=f.make_gram(f.scale(f.THETA,2),4,True)
        assert tuple(f.vecdecode(r) for r in self.base['gram'])==self.g
        self.original=walls_of(self.base)
        self.allowed={f.ZERO,f.ONE,f.TWO,f.scale(f.ONE,3),f.add(f.TWO,f.THETA),f.add(f.TWO,f.C)}
        self.pool=(f.ONE,f.TWO,f.BETA,f.scale(f.BETA,2))

    @lru_cache(None)
    def cov(self,r):return f.mv(self.g,r)
    @lru_cache(None)
    def norm(self,r):return f.dot(r,self.cov(r))
    @lru_cache(None)
    def angle_ok(self,r,s):
        p=f.dot(r,self.cov(s))
        if f.sign(p)>0:return False
        value=f.divide(f.scale(f.mul(p,p),4),f.mul(self.norm(r),self.norm(s)))
        return f.sign(f.sub(value,f.FOUR))>=0 or value in self.allowed
    def compatible(self,r,walls):return all(self.angle_ok(r,s) for s in walls)

    def root(self,row,height=None):
        r=f.vecdecode(row['vector']);q=self.norm(r)
        assert len(r)==10 and all(f.integral(x) for x in r) and f.sign(r[0])>0
        assert q in self.pool and all(f.sign(q,e)>0 for e in range(3))
        if 'norm' in row:assert f.decode(row['norm'])==q
        assert all(f.integral(f.divide(f.scale(x,2),q)) for x in self.cov(r))
        h=f.divide(f.mul(f.scale(f.THETA,2),f.mul(r[0],r[0])),q)
        assert h==f.decode(row['height'])
        if height is not None:assert f.sign(f.sub(f.scale(f.ONE,height),h))>=0
        assert self.compatible(r,self.original)
        if 'free_block' in row:
            a,b=row['free_block'],row['equal_block'];assert a!=b
            assert r[1+2*b]==r[2+2*b]
            assert all(r[1+2*k]==r[2+2*k]==f.ZERO for k in range(4) if k not in (a,b))
            if 'free_weights' in row:
                y=r[1+2*a:3+2*a];weights=f.mv(f.B,y)
                assert weights==f.vecdecode(row['free_weights'])
                assert f.dot(y,weights)==f.decode(row['free_energy'])
                x,t,z=f.vecdecode(row['parameters'])
                assert r[0]==x and r[-1]==z and r[1+2*b]==f.divide(t,f.D)
                assert f.sign(weights[0])>=0 and f.sign(f.sub(weights[1],weights[0]))>=0
                assert f.sign(t)>=0 and f.sign(z)>=0
            else:
                assert r[1+2*a]!=r[2+2*a] and r[-1]==f.ZERO
        if 'support_blocks' in row:
            blocks=row['support_blocks'];assert len(blocks)==2 and len(set(blocks))==2
            assert all(r[1+2*k]==r[2+2*k]==f.ZERO for k in range(4) if k not in blocks)
            assert r[-1]==f.ZERO
            for k,weights,energy in zip(blocks,row['block_weights'],row['block_energies']):
                y=r[1+2*k:3+2*k];actual=f.mv(f.B,y)
                assert actual==f.vecdecode(weights) and f.dot(y,actual)==f.decode(energy)
                assert f.sign(actual[0])>=0 and f.sign(f.sub(actual[1],actual[0]))>=0
        if 'coefficients' in row:
            aa=f.vecdecode(row['coefficients']);assert len(aa)==4 and q==f.TWO and r[-1]==f.ZERO
            for k,a in enumerate(aa):assert r[1+2*k]==r[2+2*k]==f.divide(a,f.D)
        if 'block' in row:
            k=row['block'];x,a,z=f.vecdecode(row['parameters'])
            assert r[0]==x and r[-1]==z and r[1+2*k]==r[2+2*k]==f.divide(a,f.D)
            assert all(r[1+2*j]==r[2+2*j]==f.ZERO for j in range(4) if j!=k)
            assert q in (f.ONE,f.TWO) and f.sign(a)>=0 and f.sign(z)>=0
        return r

    def geometry(self,d,walls):
        center=f.vecdecode(d['strict_interior'])
        assert f.sign(center[0])>0 and f.sign(self.norm(center))<0
        assert all(f.sign(f.dot(center,self.cov(r)))<0 for r in walls)
        row=d['geometry']['hyperideal_witness'];assert row is not None
        v=f.vecdecode(row['vector']);assert f.sign(v[0])>0
        norm=self.norm(v);assert norm==f.decode(row['norm']) and f.sign(norm)>0
        assert all(f.sign(f.dot(v,self.cov(r)))<=0 for r in walls)
        eps=f.ONE
        while True:
            strict=tuple(f.add(x,f.mul(eps,y)) for x,y in zip(v,center))
            if f.sign(self.norm(strict))>0:break
            eps=f.scale(eps,f.F(1,2))
        assert all(f.sign(f.dot(strict,self.cov(r)))<0 for r in walls)
        return {'strict_spacelike_witness':f.vecencode(strict),'strict_spacelike_norm':f.encoded(self.norm(strict))}

    def cone(self,d,old,roots,extras=()):
        walls=walls_of(d);assert walls[:len(old)]==old and len(walls)==d['wall_count']
        assert all(r in roots or r in extras for r in walls[len(old):])
        assert len({f.ray_key(r) for r in walls})==len(walls)
        if 'selected_indices' in d:
            selected=tuple(roots[i] for i in d['selected_indices'])
            assert walls[len(old):len(old)+len(selected)]==selected
            assert all(r in extras for r in walls[len(old)+len(selected):])
        for i,r in enumerate(walls):
            assert self.norm(r) in self.pool
            assert all(f.integral(f.divide(f.scale(x,2),self.norm(r))) for x in self.cov(r))
            assert all(self.angle_ok(r,s) for s in walls[:i])
        result=self.geometry(d,walls)
        bad=copy.deepcopy(d);bad['geometry']['hyperideal_witness']['norm']=f.encoded(f.ONE)
        try:self.geometry(bad,walls)
        except AssertionError:pass
        else:raise AssertionError('forged escape norm accepted')
        return {'walls':len(walls),'pairwise_angle_checks':len(walls)*(len(walls)-1)//2,
                'forged_norm_rejected':True,**result}

    def single(self,d):
        base=read(d['base_source']) if 'base_source' in d else self.base;old=walls_of(base)
        target=read(d['replaced_base_source']) if 'replaced_base_source' in d else base
        old_escape=escape_of(target)
        rows=d.get('compatible_roots',d.get('compatible_new_roots'));assert rows is not None
        roots=[];cutters=0
        for row in rows:
            r=self.root(row,d['height_bound']);assert self.compatible(r,old)
            p=f.dot(r,self.cov(old_escape));assert p==f.decode(row.get('previous_escape_pairing',row.get('old_escape_pairing')))
            cutters+=f.sign(p)>0;roots.append(r)
        assert len({f.ray_key(r) for r in roots})==len(roots)
        extras=[]
        if 'additional_pell_cap' in d:
            cap=self.root(d['additional_pell_cap']);assert self.norm(cap)==f.ONE
            assert all(x==f.ZERO for x in cap[1:9]);extras.append(cap)
        return {'new_roots':len(roots),'old_escape_cutters':cutters,**self.cone(d,old,roots,extras)}

    def scalar(self,d):
        bases={name:read_case(path) for name,path in d['previous_cone_sources'].items()}
        bwalls={name:walls_of(base) for name,base in bases.items()}
        roots=[];cutters={name:0 for name in bases}
        for row in d['compatible_roots']:
            r=self.root(row,d['height_bound']);roots.append(r)
            for name,base in bases.items():
                assert self.compatible(r,bwalls[name])==row['compatible_with_previous_cones'][name]
                p=f.dot(r,self.cov(escape_of(base)))
                assert p==f.decode(row['previous_escape_pairings'][name]);cutters[name]+=f.sign(p)>0
        assert len({f.ray_key(r) for r in roots})==len(roots)
        combined=d['combined_replacement_dictionary'];union=[]
        if d['status']=='COMPLETE_BOUNDED_MIXED_BLOCK_SCALAR_SEARCH':
            previous=read(next(iter(d['previous_cone_sources'].values()))['file'])
            origins={'previous_union':previous['combined_replacement_dictionary'],'mixed_scalar':d['compatible_roots']}
        else:
            origins={'two_block':bases['replacement122']['compatible_roots'],'equal_block_scalar':d['compatible_roots']}
        def without_origin(row):return {k:v for k,v in row.items() if k not in ('origin','origin_index')}
        for row in combined:
            origin=origins[row['origin']]
            assert without_origin(row)==without_origin(origin[row['origin_index']])
            union.append(self.root(row,d['height_bound']))
        assert len(set(union))==len(union)
        assert set(union)=={f.vecdecode(row['vector']) for rows in origins.values() for row in rows}
        capdata=read('results/fresh_cubic7_dihedral_pell_cap_20261013.json')
        cap=self.root({'vector':capdata['root'],'norm':capdata['norm'],'height':capdata['height']})
        cases={}
        for name,case in d['cases'].items():
            if name=='replacement_from25':
                candidates=union;old=self.original
                if case['pell_cap_status']=='REJECTED_BY_EXACT_CROSS_ANGLE_GATE':
                    assert not self.compatible(cap,walls_of(case))
                else:assert case['pell_cap_status']=='APPENDED'
                extras=[cap] if case['pell_cap_status']=='APPENDED' else []
            else:
                expected=[i for i,row in enumerate(d['compatible_roots']) if row['compatible_with_previous_cones'][name]]
                assert case['candidate_indices']==expected
                candidates=[roots[i] for i in expected];old=bwalls[name];extras=[]
            cases[name]=self.cone(case,old,candidates,extras)
        return {'new_roots':len(roots),'old_escape_cutters':cutters,'cases':cases}

    def verify(self,d):
        if d.get('status')=='EXACT_COMMON_SPACELIKE_ESCAPE':return self.common_escape(d)
        result=self.scalar(d) if 'cases' in d else self.single(d)
        return {'status':'PASS',**result,'scope':'Independent saved-root, reflection, angle, selection, and geometry replay. Positive-height candidate completeness is not independently rerun.'}

    def common_escape(self,d):
        source=read(d['source']);rows=source['combined_replacement_dictionary']
        roots=tuple(self.root(row,source['height_bound']) for row in rows)
        assert d['original_walls']==len(self.original) and d['dictionary_roots']==len(roots)
        walls=self.original+roots;v=f.vecdecode(d['vector'])
        assert f.vecdecode(d['wall_pairings'])==tuple(f.dot(v,self.cov(r)) for r in walls)
        capdata=read('results/fresh_cubic7_dihedral_pell_cap_20261013.json')
        cap=self.root({'vector':capdata['root'],'norm':capdata['norm'],'height':capdata['height']})
        p=f.dot(v,self.cov(cap));assert p==f.decode(d['pell_cap_pairing'])
        assert d['pell_cap_also_satisfied'] and f.sign(p)<=0
        geometry={'strict_interior':d['strict_interior'],'geometry':{'hyperideal_witness':{'vector':d['vector'],'norm':d['norm']}}}
        result=self.geometry(geometry,(*walls,cap))
        return {'status':'PASS','new_roots':len(roots),'walls':len(walls)+1,
                'scope':'Common exact escape and strict interior for the original25 walls, all855 union dictionary normals, and the separate Pell cap. Pairwise Coxeter compatibility of the entire union is neither asserted nor required.',**result}

def verify(d):return Auditor().verify(d)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('inputs',nargs='+',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    auditor=Auditor();out=[]
    for path in a.inputs:
        row=auditor.verify(json.loads(path.read_text()));row['input']=str(path);out.append(row)
        print(json.dumps({'input':str(path),'status':'PASS','roots':row['new_roots'],
            'walls':row.get('walls',{k:v['walls'] for k,v in row.get('cases',{}).items()})}),flush=True)
    a.output.write_text(json.dumps({'status':'PASS','cases':out},indent=2,sort_keys=True)+'\n')
