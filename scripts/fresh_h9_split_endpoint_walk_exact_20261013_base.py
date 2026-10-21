#!/usr/bin/env python3
"""Archive-only exact coordinates for the split-seed edge walk."""
import hashlib
from scripts.fresh_h9_replaced_core_20261013_base import *
from scripts.fresh_h4_blocks_independent_verify import inverse

PREFIX='data/enumerations/fresh_h9_split_endpoint_walk_exact_20261013'
SCOUT='data/enumerations/fresh_h9_split_endpoint_walk_20261013_scout.json'
OPTIONAL_SCOUT='data/enumerations/fresh_h9_split_endpoint_walk_optional_anchor_20261013_scout.json'
NORMALS='data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json'
BASIS=(0,1,2,3,4,5,6,7,8,12)

def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

def source():
    scout=read(SCOUT);snap=scout['source_snapshot'];data=read(NORMALS)
    g16=tuple(tuple(parse(x)for x in row)for row in data['gram'])
    G=tuple(tuple(g16[i][j]for j in BASIS)for i in BASIS);I=inverse(G)
    ps={}
    for i in range(16):
        p=tuple(g16[j][i]for j in BASIS);c=matvec(K,I,p);name=f'seed{i}'
        ps[name]={'name':name,'field':K,'c':c,'p':p,'q':scale(O,2)}
    s=tuple(K.zero for _ in range(10))
    for i in range(16):s=tuple(add(x,y)for x,y in zip(s,ps[f'seed{i}']['c'],strict=True))
    h=tuple(K.neg(x)for x in matvec(K,G,s))
    banks={NORMALS:data}
    for j,(filename,index)in enumerate(snap['refs'],16):
        path='data/enumerations/'+filename
        if path not in banks:banks[path]=read(path)
        raw=banks[path]['caps'][index];F=Quadratic(K,parse(raw['radicand']))
        p16=tuple(F.parse(x)for x in raw['covector'])
        c16=tuple(F.parse(x)for x in raw['normal'])
        assert len(p16)==len(c16)==16
        gg=tuple(tuple(F.lift(x)for x in row)for row in g16)
        assert all(F.sign(F.sub(x,y))==0 for x,y in zip(matvec(F,gg,c16),p16,strict=True))
        p=tuple(p16[i]for i in BASIS)
        c=matvec(F,tuple(tuple(F.lift(x)for x in row)for row in I),p)
        assert F.sign(F.sub(dot(F,c,p),F.rational(2)))==0
        for i in range(16):
            assert F.sign(F.sub(dot(F,c,tuple(F.lift(x)for x in ps[f'seed{i}']['p'])),p16[i]))==0
        assert all(F.sign(x)<=0 for x in p16)
        assert F.sign(dot(F,c,tuple(F.lift(x)for x in h)))>0
        name=f'wall{j}'
        ps[name]={'name':name,'field':F,'c':c,'p':p,'q':scale(O,2),
                  'reference':{'source':path,'index':index}}
    assert len(ps)==72
    return G,I,h,ps,scout,banks

def wall_name(i):return f'seed{i}'if i<16 else f'wall{i}'

def point_pair(row,F,y):
    common,er=embed(F,row['field'])
    return common,dot(common,tuple(er(x)for x in row['c']),
                      tuple(lift_existing(common,F,x)for x in y))

def load_cap(record):return cap_profile(record,f"walk{record['proposal_index']}")
