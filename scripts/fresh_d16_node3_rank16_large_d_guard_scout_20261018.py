#!/usr/bin/env python3
"""Exact candidate geometry and small coupon upper; no universal claim alone."""
import argparse
from fractions import Fraction as F
import fcntl
from itertools import combinations
import json
from pathlib import Path
import signal
import fresh_a16_double_large_d_physical_scout_20261018 as S

V=S.V;M=S.M;U=M.U;C=M.C;ROOT=M.ROOT;A=S.A;Q=S.Q;D=S.D;T=S.T
TAG='fresh_d16_node3_rank16_large_d_guard_scout_20261018'
FULL=[498,1231,1294,1308,2020,2447]
PROPER=[504,1019,1031,1266,1285,1537,1614,1624]


def get(index):
    marker=('{"index":'+str(index)+',').encode();c=S.extract(marker,len(marker))
    core=S.extract((json.dumps(c['core_key'],separators=(',',':'))+':{').encode(),1)
    g=A.matrix(core['gram']);left=c['left'][1];right=[c['right_copy_map'][i]for i in c['right'][1]]
    white=sorted(set(range(len(g)))-set(left+right))
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,white+[b]))==[len(white)+1,0,0]for b in left+right)
    actual=Q.inspect_corner(g,[1]*len(g),left+right,*c['preferred_corner'])
    for k in('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):assert actual[k]==c[k]
    return c,core,g,white,left,right


def corners(g,w,c):
    out=[]
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(w)|set(marks));assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
        bs=D.finite_blocks(g,ids);a,b=[next(b for b in bs if m in b['nodes'])for m in marks]
        selected=marks==c['preferred_corner'];coupon=False;images=[]
        if selected:
            assert a is b and a['type']=='D16';images=T.mark_pairs(g,a['nodes'],marks,16,'D');assert images==[(3,13)]
        if a is not b:
            for sm,la,m,n in((a,b,*marks),(b,a,marks[1],marks[0])):
                if(sm['type'],la['type'])==('A2','A13'):
                    images=T.mark_pairs(g,la['nodes'],[n],13,'A')
                    if images==[(3,),(9,)]:
                        assert any(b['type']=='A1'and not set(marks)&set(b['nodes'])for b in bs);coupon=True
        out.append(dict(marks=marks,components=bs,selected_D16_zero=selected,coupon=coupon,images=images,payment=int(coupon)))
    assert sum(r['payment']for r in out)==1
    return out


def geometry():
    sup=V.load(S.L.TAG)['local_proofs']['domains']['D']['complete_selected_support']
    hits=[r for r in sup['rows']if[3,13]in r['marked_images']]
    assert[r['index']for r in hits]==sorted(FULL+PROPER)
    full=[];proper=[]
    for index in FULL:
        c,core,g,w,left,right=get(index);assert len(c['marked_component'])==16
        if c['method']=='radical':
            n=len(g);assert Q.inertia(g)==[n-2,1,1]
            z=list(map(F,core['radical'][0]));assert all(sum(g[i][j]*z[j]for j in range(n))==0 for i in range(n))
            allowed=c['possible_contacts'];assert all(z[i]for i in allowed)
            du,dv=Q.graph_dist(g,left),Q.graph_dist(g,right)
            ss=[dict(nodes=[i,j],length=min(du[i]+dv[j],du[j]+dv[i])+2)for i,j in combinations(allowed,2)if z[i]*z[j]<0]
            assert ss and max(r['length']for r in ss)<=8
            seal=dict(method='radical',radical=z,allowed_contacts=allowed,mixed_sign_shortcuts=ss)
        else:
            gg,J,du,dv,end=S.X.constraints(c,core);assert gg==g
            rem=sorted(set(c['possible_contacts'])-set(c['marked_component']))
            assert rem==c['omitted_black']and sorted(J[i][i]for i in rem)==[4,8]
            seal=dict(method='whole_D16_projection',inverse=J,end_proof=end,remaining=rem)
        full.append(dict(index=index,core_key=c['core_key'],gram=g,white=w,seal=seal,corners=corners(g,w,c)))
    for index in PROPER:
        c,core,g,w,left,right=get(index);gg,J,du,dv,end=S.X.constraints(c,core);assert gg==g
        row=next(r for r in hits if r['index']==index);maps=[]
        for mp in row['maps']:
            if mp['marked_images']!=[3,13]:continue
            image=dict(zip(mp['images'],row['base_nodes']));missing=sorted(set(range(16))-set(image));entries=[]
            for new in missing:
                for old,node in image.items():
                    edge=abs(new-old)==1 and max(new,old)<=13 or set((new,old))in({13,14},{13,15})
                    if edge:entries.append(dict(new=new,old=old,core=node,diagonal=J[node][node]))
            assert entries
            if index in(504,1537):assert all(e['diagonal']==8 for e in entries);outcome='projection_above_two'
            elif index in(1031,1614):
                leaves=[e for e in entries if e['new']in(14,15)]
                assert missing==[0,14,15]and len(leaves)==2 and len({e['core']for e in leaves})==1
                assert all(e['diagonal']==2 for e in leaves);outcome='two_orthogonal_norm_saturated_leaves'
            elif index in(1285,1624):
                assert missing==[0,1,2]and len(entries)==1 and entries[0]['diagonal']==2
                k=entries[0]['core'];negative=[dict(node=j,bound=min(du[k]+dv[j],dv[k]+du[j])+3)for j in c['possible_contacts']if J[k][j]<0]
                assert negative and max(e['bound']for e in negative)<=8;outcome='norm_saturated_root_cannot_continue'
            else:
                assert index in(1019,1266)and missing==[0]and len(entries)==1 and entries[0]['diagonal']==2
                k=entries[0]['core'];n=len(g);completed=[r+[F(-int(i==k))]for i,r in enumerate(g)]
                completed.append([F(-int(i==k))for i in range(n)]+[F(2)])
                radical=[J[i][k]for i in range(n)]+[F(1)];assert Q.inertia(completed)==[n-1,1,1]
                assert all(sum(completed[i][j]*radical[j]for j in range(n+1))==0 for i in range(n+1))
                assert sorted(radical[i]for i in c['omitted_black'])==[-2,2]
                maps.append(dict(images=mp['images'],missing=missing,forced_entries=entries,
                    outcome='forced_single_leaf_completion',completed_gram=completed,completed_white=w+[n],
                    radical=radical,remaining_contacts=c['omitted_black'],corners=corners(completed,w+[n],c)))
                continue
            maps.append(dict(images=mp['images'],missing=missing,forced_entries=entries,outcome=outcome,
                negative_coefficient_shortcuts=negative if index in(1285,1624)else[]))
        proper.append(dict(index=index,gram=g,inverse=J,white=w,selected_nodes=c['marked_component'],
            possible_contacts=c['possible_contacts'],omitted_black=c['omitted_black'],end_proof=end,maps=maps))
    return dict(complete_patterns=394,complete_maps=578,context_maps=sum(len(p['maps'])for p in hits),
        affected_maps=sum(m['marked_images']==[3,13]for p in hits for m in p['maps']),full=full,proper=proper)


def upper():
    cs=(0,0,0,0,1,0,0,1,0);key=next(key for counts,used,key in C.fixed_domain()if counts==cs and key[9]==0)
    rows=[];tabs=[];checks=preds=0
    for k5 in(0,1):
        fixed=C.fixed_components(cs,key)+[('A',5)]*k5;fc=C.upper_value(fixed)
        weights={r:C.item(r,key,k5)+(128 if r==1 else 256 if r in(2,4)else 0)for r in C.ITEMS}
        tab,a,b=U.table(weights);checks+=a;preds+=b;tabs.append(dict(k5=k5,table=tab))
        for s in range(18,26-5*k5):
            n=28-5*k5-s;assert n>=3;assert C.item(s,key,k5)==C.actual_base(('D',s))
            opts=[]
            for a1 in range(1,n-1):
                for a2 in range(1,(n-a1)//2+1):
                    rest=n-a1-2*a2;v=tab['values'][rest]
                    if v is not None:opts.append((fc+C.actual_base(('D',s))+a1*weights[1]+a2*weights[2]+v,a1,a2,rest))
            best=max(opts);rows.append(dict(k5=k5,reserved_D=s,residual=n,upper=F(best[0],64),explicit_counts=best[1:]))
    return dict(classes=2,tables=tabs,rows=rows,Bellman_inequalities=checks,predecessors=preds,
        increasing_family_upper=max(r['upper']for r in rows),scope='All new coupon forests; all negative A16/D16/D17 changes ignored.')


def run():
    gp=geometry();up=upper();assert up['increasing_family_upper']<372
    paths={Path(__file__),Path(S.__file__),Path(M.__file__),Path(V.__file__),S.OVER,
        ROOT/f'results/{S.L.TAG}.json',ROOT/f'results/{S.TAG}.json',ROOT/f'results/{M.TAG}.json'}
    return V.O.enc(dict(status='PASS_D16_NODE3_RANK16_LARGE_D_GUARD_SCOUT',sources={str(p.relative_to(ROOT)):U.digest(p)for p in sorted(paths)},
        geometry=gp,upper=up,rule='D16(3,13) unit becomes zero with separate whole A16 orD16 and D>=18; add A2/A13fourth unit with A1 and the same external witnesses',
        independent_review_complete=False,universal_composition_complete=False,new_dimension_bound=None))


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);M.R.resources();U.memory_guard()
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run();U.memory_guard();M.R.resources(start=False)
    if args.verify:assert V.load(TAG)==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['geometry']['affected_maps'],out['upper']['increasing_family_upper'],flush=True)
