#!/usr/bin/env python3
"""Separate Fraction/tower replay and common escape for the 107 H20 cap entries."""
from dataclasses import dataclass
from functools import lru_cache
from fractions import Fraction as Q
from pathlib import Path
from collections import Counter
import argparse,hashlib,json,time,copy
from fresh_dimension_bridge_20261013_verify import E,Z,O,parse,vec,dot,mv,scalar,check as check_base

ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_dimension_bridge_angle_repairs_20261013_exact.json'
OUT=ROOT/'results/fresh_dimension_bridge_angle_repairs_20261013_independent.json'


@dataclass(frozen=True)
class Tower:
    base:object=None
    radicand:object=None

    def lift(self,x):return x if self.base is None else(self.base.lift(x),self.base.lift(Z))
    def rational(self,x):return self.lift(scalar(x))
    def parse(self,x):return parse(x)if self.base is None else tuple(self.base.parse(z)for z in x)
    def encode(self,x):return x.enc()if self.base is None else[self.base.encode(z)for z in x]
    def add(self,x,y):return x+y if self.base is None else tuple(self.base.add(a,b)for a,b in zip(x,y,strict=True))
    def neg(self,x):return -x if self.base is None else tuple(self.base.neg(a)for a in x)
    def sub(self,x,y):return self.add(x,self.neg(y))

    @lru_cache(maxsize=250000)
    def mul(self,x,y):
        if self.base is None:return x*y
        f=self.base;a,b=x;c,d=y
        return(f.add(f.mul(a,c),f.mul(self.radicand,f.mul(b,d))),f.add(f.mul(a,d),f.mul(b,c)))

    @lru_cache(maxsize=100000)
    def sign(self,x):
        if self.base is None:return x.sign()
        a,b=x;f=self.base;sa,sb=f.sign(a),f.sign(b);sd=f.sign(self.radicand)
        assert sd>=0
        if not sd or not sb:return sa
        if not sa or sa==sb:return sb
        norm=f.sub(f.mul(a,a),f.mul(self.radicand,f.mul(b,b)))
        return sa*f.sign(norm)

    def dot(self,a,b):
        s=self.lift(Z);zero=s
        for x,y in zip(a,b,strict=True):
            if x!=zero and y!=zero:s=self.add(s,self.mul(x,y))
        return s

    def mv(self,M,v):return tuple(self.dot(row,v)for row in M)


BASE=Tower()


def field(raw):
    if raw['kind']=='sqrt6':return BASE
    assert raw['kind']=='quadratic';base=field(raw['base']);d=base.parse(raw['radicand'])
    assert base.sign(d)>=0
    return Tower(base,d)


def classify(f,p,q,n):
    assert f.sign(p)<=0 and f.sign(q)>0 and f.sign(n)>0
    w=f.mul(f.rational(4),f.mul(p,p));den=f.mul(q,n)
    difference=f.sign(f.sub(w,f.mul(f.rational(4),den)))
    if difference>0:return'dotted'
    if difference==0:return'parallel'
    for label,value in((2,0),(3,1),(4,2),(6,3)):
        if f.sign(f.sub(w,f.mul(f.rational(value),den)))==0:return str(label)
    # Homogeneous polynomial identities and rational isolating intervals avoid
    # all division in possibly dependent radical towers.
    ww=f.mul(w,w);wd=f.mul(w,den);dd=f.mul(den,den)
    for label,s,t,lo,hi in((5,3,1,Q(26,10),Q(27,10)),(8,4,2,Q(34,10),Q(35,10)),
                           (10,5,5,Q(36,10),Q(37,10)),(12,4,1,Q(37,10),Q(38,10))):
        poly=f.add(f.sub(ww,f.mul(f.rational(s),wd)),f.mul(f.rational(t),dd))
        if f.sign(poly):continue
        if f.sign(f.sub(w,f.mul(f.rational(lo),den)))>0 and f.sign(f.sub(f.mul(f.rational(hi),den),w))>0:return str(label)
    raise AssertionError('unverified finite angle')


def rank_field(rows):
    a=[list(row)for row in rows];rank=0
    for j in range(len(a[0])):
        pivot=next((i for i in range(rank,len(a))if a[i][j]!=Z),None)
        if pivot is None:continue
        a[rank],a[pivot]=a[pivot],a[rank];p=a[rank][j];a[rank]=[x/p for x in a[rank]]
        for i in range(rank+1,len(a)):
            p=a[i][j]
            if p!=Z:a[i]=[x-p*y for x,y in zip(a[i],a[rank],strict=True)]
        rank+=1
        if rank==len(a):break
    return rank


def audit():
    started=time.monotonic();data=json.loads((ROOT/SOURCE).read_text());scout=json.loads((ROOT/data['source']).read_text())
    assert hashlib.sha256((ROOT/data['source']).read_bytes()).hexdigest()==data['source_sha256']
    polar=json.loads((ROOT/data['base_source']).read_text())
    assert data['base_source']==scout['source']and data['base_sha256']==scout['source_sha256']
    assert hashlib.sha256((ROOT/data['base_source']).read_bytes()).hexdigest()==data['base_sha256']
    base=json.loads((ROOT/polar['source']).read_text())
    assert hashlib.sha256((ROOT/polar['source']).read_bytes()).hexdigest()==polar['source_sha256']
    base_audit=check_base(base)
    form=tuple(tuple(map(parse,row))for row in base['form'])
    normals=list(map(vec,base['normals']));covectors=[mv(form,r)for r in normals];norms=[dot(r,c)for r,c in zip(normals,covectors,strict=True)]
    center=vec(base['strict_center']);center_cov=mv(form,center);center_norm=dot(center,center_cov)
    for record in polar['cuts']:
        r=vec(record['vector']);c=mv(form,r);q=dot(r,c)
        assert q==parse(record['norm'])and q.sign()>0
        assert dot(r,center_cov)==parse(record['center_pairing'])and dot(r,center_cov).sign()<0
        pairs=[dot(r,p)for p in covectors]
        assert pairs==list(map(parse,record['pairings']))
        assert [classify(BASE,p,q,n)for p,n in zip(pairs,norms,strict=True)]==record['pair_kinds']
        assert len(normals)==record['index']
        normals.append(r);covectors.append(c);norms.append(q)
    assert len(normals)==polar['wall_count']==82
    escape=polar['strict_escape'];v=vec(escape['vector']);vc=mv(form,v)
    vnorm=dot(v,vc);vcenter=dot(v,center_cov)
    assert vnorm==parse(escape['norm'])and vnorm.sign()>0
    assert vcenter==parse(escape['center_pairing'])and vcenter.sign()<0
    old_pairs=[dot(v,c)for c in covectors]
    assert old_pairs==list(map(parse,escape['all_wall_pairings']))and all(p.sign()<0 for p in old_pairs)
    ranks={};rows=[];labels=Counter();scout_ids=[];target_ids=set()
    for i,row in enumerate(data['caps']):
        scout_index=row['scout_index'];raw=scout['candidates'][scout_index];target=polar['rejected_polar_proposals'][raw['target_index']]
        scout_ids.append(scout_index);ids=[j for j in target['active_indices']if j!=raw['released']]+[raw['anchor']]
        assert ids==row['active_indices']and len(ids)==20 and len(set(ids))==20
        assert row['released']==raw['released']and row['anchor']==raw['anchor']and row['anchor_label']==raw['label']
        key=tuple(ids)
        if key not in ranks:ranks[key]=rank_field([covectors[j]for j in ids])
        assert ranks[key]==20
        f=field(row['field']);r=tuple(map(f.parse,row['vector']));q=f.parse(row['norm']);assert f.sign(q)>0
        M=tuple(tuple(f.lift(z)for z in rr)for rr in form);c=f.mv(M,r)
        assert f.sign(f.sub(f.dot(r,c),q))==0
        center_pair=f.dot(r,tuple(map(f.lift,center_cov)))
        assert f.sign(center_pair)<0 and f.sign(f.sub(center_pair,f.parse(row['center_pairing'])))==0
        tv=vec(target['vector']);tc=mv(form,tv)
        assert dot(tv,tc)==parse(target['norm'])and dot(tv,tc).sign()>0
        assert all(dot(tv,p).sign()<=0 for p in covectors)
        target_ids.add(raw['target_index'])
        target_pair=f.dot(r,tuple(map(f.lift,tc)))
        assert f.sign(target_pair)>0 and f.sign(f.sub(target_pair,f.parse(row['target_pairing'])))==0
        pairs=[f.dot(r,tuple(map(f.lift,p)))for p in covectors]
        kinds=[classify(f,p,q,f.lift(n))for p,n in zip(pairs,norms,strict=True)]
        assert kinds==row['pair_kinds']and all(f.sign(pairs[j])==0 for j in ids[:19])
        assert kinds[row['anchor']]==str(row['anchor_label'])
        if row['anchor_label']!=2:assert f.sign(f.add(pairs[row['anchor']],f.rational(1)))==0
        common=f.dot(r,tuple(map(f.lift,vc)))
        assert f.sign(common)<0,('common escape is cut',i)
        rows.append({'cap_index':i,'scout_index':scout_index,'field':row['field'],
                     'common_escape_pairing':f.encode(common),'common_escape_pairing_sign':-1,
                     'exact_old_pair_count':len(kinds),'twenty_constraint_rank':20,
                     'norm_and_center_verified':True,'target_cut_verified':True,'anchor_label':row['anchor_label']})
        labels[str(row['anchor_label'])]+=1
    assert scout_ids==list(range(len(scout['candidates'])))and len(rows)==107
    assert dict(labels)==data['counts_by_anchor_label']and len(ranks)==data['linear_systems']
    # Explicit dependency cases exercise signs without an independence assumption.
    assert Tower(BASE,O).sign((O,-O))==0 and Tower(BASE,Z).sign((-O,scalar(100)))<0
    assert dot(center,mv(form,center)).sign()<0
    forged=list(r);forged[0]=f.add(forged[0],f.rational(1))
    assert f.sign(f.sub(f.dot(tuple(forged),f.mv(M,tuple(forged))),q))!=0
    return {'status':'PASS_INDEPENDENT_107_CAP_BANK_AND_COMMON_STRICT_ESCAPE',
            'source':SOURCE,'source_sha256':hashlib.sha256((ROOT/SOURCE).read_bytes()).hexdigest(),
            'additional_sources_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in(data['source'],data['base_source'],polar['source'])},
            'base_replay':base_audit,'base_walls':82,'cap_entries':107,'new_old_pair_checks':107*82,
            'linear_systems_rank20':len(ranks),'counts_by_anchor_label':dict(labels),
            'distinct_target_rays':sorted(target_ids),'all_targets_feasible_in_base82':True,
            'strict_center_norm':center_norm.enc(),'common_escape':{'vector':escape['vector'],'norm':vnorm.enc(),
                  'center_pairing':vcenter.enc(),'base_inequality_signs':[-1]*82,'cap_inequality_signs':[-1]*107},
            'cap_checks':rows,'total_strict_inequalities':189,'dependent_radical_sign_regressions':True,
            'forged_normal_rejected':True,
            'seconds':round(time.monotonic()-started,3),
            'scope':'Separate Fraction Q(sqrt6) and recursive tower arithmetic replay every cap norm, common center, 19 orthogonalities, prescribed anchor angle, rank20 constraints, all82 old Coxeter pairs, and target cut. The same exact strict future spacelike vector satisfies the entire82+107 halfspace bank, so no subset of these entries can give finite volume while retaining the82 base walls. Entries need not be distinct mirrors; no new-new compatibility, search completeness, or unrestricted no-completion claim.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=audit()
    if args.verify:
        saved=json.loads(OUT.read_text());saved.pop('seconds');result.pop('seconds');assert result==saved
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print({k:result[k]for k in('status','base_walls','cap_entries','new_old_pair_checks','total_strict_inequalities')},flush=True)
