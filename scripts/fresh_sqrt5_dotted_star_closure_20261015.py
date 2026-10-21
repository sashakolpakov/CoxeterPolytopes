#!/usr/bin/env python3
"""Complete attachment queue for fixed dotted-distance rank completions."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import networkx as nx
import fresh_sqrt5_dotted_completion_20261015 as D
import fresh_sqrt5_all_quad_unions_20261015_classes as C

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_dotted_star_closure_20261015'
F,Q,K=D.F,D.Q,D.K
def dec(x):return K(*x)
def matrix(row,ratio):
    g=Q.edges(row);M=Q.tree_scaled(row['nodes'],{p:2 if m==0 else m for p,m in g.items()})
    if M is None:return None,'REJECT_FINITE_CYCLE_FIELD'
    H=nx.Graph();H.add_nodes_from(range(row['nodes']));H.add_edges_from(p for p,m in g.items()if m not in(0,2))
    z=ratio*M[0][0]*M[2][2]
    if nx.has_path(H,0,2):
        x=D.square_root(z)
        if x is None:return None,'REJECT_DOTTED_CYCLE_FIELD'
        M[0][2]=M[2][0]=-x
    else:
        if z.conj().sign()<=0:return None,'REJECT_DOTTED_CONJUGATE_SQUARE'
        component=set(nx.node_connected_component(H,0))
        for i in component:
            for j in component:M[i][j]=z*M[i][j]
        M[0][2]=M[2][0]=-z
    return M,None
def evaluate(row,ratio):
    M,error=matrix(row,ratio);assert error is None
    n=row['nodes'];assert Q.inertia(M)==[n-2,1,1]
    assert Q.inertia([[x.conj()for x in r]for r in M])==[n-1,0,1]
    assert all(c['integral']for c in D.cycles(M))
    domain,dc=F.profiles(row);trials=[];accepted=[]
    for contact in domain:
        ext=F.extension(row,contact);r=dict(contacts=contact)
        if(w:=F.OLD.h4_payments(ext)):r.update(status='PAID_H4',witness=w)
        elif(w:=F.affine(ext)):r.update(status='REJECT_AFFINE',witness=w)
        else:
            G,error=matrix(ext,ratio)
            if error:r['status']=error
            else:
                pi=Q.inertia(G);co=Q.inertia([[x.conj()for x in z]for z in G]);r.update(physical=pi,conjugate=co)
                if pi[1]>1:r['status']='REJECT_TWO_NEGATIVE'
                elif co[1]:r['status']='REJECT_CONJUGATE'
                else:
                    assert pi==[n-1,1,1]and co==[n,0,1]
                    cy=D.cycles(G);bad=next((z for z in cy if not z['integral']),None)
                    if bad:r.update(status='REJECT_CYCLIC_INTEGRALITY',witness=bad)
                    else:r['status']='ACCEPT_UNPAID';accepted.append(contact)
        trials.append(r)
    return dict(status='OPEN_WITH_CONTACTS'if accepted else'SEALED',domain_counts=dc,trials=trials,accepted=accepted)
def run(limit):
    paths=[f'scripts/{TAG}.py',f'scripts/{D.TAG}.py',f'results/{D.TAG}.json','scripts/fresh_sqrt5_all_quad_unions_20261015_classes.py']
    Q.manifest_check(paths);src=json.loads((ROOT/f'results/{D.TAG}.json').read_text())
    states=[];lookup={};seeds=[]
    def add(row,ratio):
        key,order=C.canonical(row);key=(key,tuple(ratio.enc()))
        if key not in lookup:
            lookup[key]=len(states);states.append(dict(index=len(states),row=row,ratio=ratio.enc(),order=order,result=None))
        j=lookup[key];target=states[j];mp=[None]*row['nodes']
        for a,b in zip(order,target['order']):mp[a]=b
        C.check_map(row,target['row'],mp)
        assert {mp[0],mp[2]}=={0,2}
        return dict(index=j,map=mp)
    for ri,r in enumerate(src['records']):
        for ci,c in enumerate(r['result'].get('candidates',[])):
            if c['status']=='PASS_ARITHMETIC_NECESSARY_COMPLETION':seeds.append(dict(record=ri,candidate=ci,**add(r['row'],dec(c['squared_dotted_ratio']))))
    assert len(seeds)==27
    count=Counter();tc=Counter();i=0
    while i<len(states)and i<limit:
        rec=states[i];z=evaluate(rec['row'],dec(rec['ratio']));children=[]
        for c in z['accepted']:children.append(dict(contacts=c,**add(F.extension(rec['row'],c),dec(rec['ratio']))))
        z['children']=children;rec['result']=z;count[z['status']]+=1;tc.update(r['status']for r in z['trials'])
        print('state',i,z['status'],'contacts',len(children),'total',len(states),flush=True);i+=1
    return dict(status='COMPLETE_FIXED_DOTTED_SEED_ATTACHMENT_QUEUE'if i==len(states)else'PARTIAL_FIXED_DOTTED_SEED_ATTACHMENT_QUEUE',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},processed=i,marked_states=len(states),counts=dict(count),trial_counts=dict(tc),seeds=seeds,records=states,scope='Fixed-distance attachment closure of27 necessary zero-contact rank completions only. Complete refers only to this seed family; distance-varying nonsingular intermediate cores are outside it. A sealed row fixes its affected corner components in an unpaid extension, leaving orthogonal spherical spectators. Face payment and actual compact realization remain separate obligations. No dimension bound.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--max-states',type=int,default=1000);args=ap.parse_args();p=ROOT/f'results/{TAG}.json'
    data=json.loads(json.dumps(run(args.max_states)))
    if args.verify:assert data==json.loads(p.read_text())
    else:p.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['processed'],data['marked_states'],data['counts'],data['trial_counts'])
