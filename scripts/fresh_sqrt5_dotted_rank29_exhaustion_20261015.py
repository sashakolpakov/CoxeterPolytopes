#!/usr/bin/env python3
"""Complete dotted attachment queue relevant to dimensions at most29."""
import argparse,json
from collections import Counter,defaultdict
from hashlib import sha256
from pathlib import Path
import fresh_sqrt5_dotted_star_closure_20261015 as S
import fresh_sqrt5_sealed_subcore_closure_20261015_review as J

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_dotted_rank29_exhaustion_20261015'
PARAM='fresh_sqrt5_dotted_parameters_20261015'
D,F,Q,K,C=S.D,S.F,S.Q,S.K,S.C
def evaluate(row,ratio):
    M,error=S.matrix(row,ratio);assert error is None
    n=row['nodes'];pi=Q.inertia(M);co=Q.inertia([[x.conj()for x in r]for r in M]);assert pi[1]==1 and pi[2]in(0,1)and co[1]==0 and co[2]==pi[2]
    assert all(c['integral']for c in D.cycles(M));singular=bool(pi[2])
    assert n<=31
    if n==31:
        # A full star in dimensiond has d+2 nodes and Gram rankd+1.
        # No extension can fitd≤29; a nonsingular31-node core cannot fit either.
        return dict(status='DIMENSION29_SINGULAR_WEIGHT_OBLIGATION'if singular else'REJECT_RANK_EXCEEDS_DIMENSION29',singular=singular,physical=pi,conjugate=co,domain_counts={},trials=[],accepted=[])
    domain,dc=F.profiles(row);trials=[];accepted=[]
    for contact in domain:
        ext=F.extension(row,contact);r=dict(contacts=contact)
        if(w:=F.OLD.h4_payments(ext)):r.update(status='PAID_H4',witness=w)
        else:
            G,error=S.matrix(ext,ratio)
            if error:r['status']=error
            else:
                p=Q.inertia(G);c=Q.inertia([[x.conj()for x in z]for z in G]);r.update(physical=p,conjugate=c)
                if p[1]>1:r['status']='REJECT_TWO_NEGATIVE'
                elif c[1]:r['status']='REJECT_CONJUGATE'
                else:
                    assert p[1]==1 and p[2]in(0,1)and c[2]==p[2]
                    if singular:assert p[2]==1
                    ker=Q.kernel(G)if p[2]else[]
                    if ker and(all(x.sign()>=0 for x in ker[0])or all(x.sign()<=0 for x in ker[0])):
                        r.update(status='REJECT_NONNEGATIVE_ROOT_DEPENDENCE',witness=[x.enc()for x in ker[0]])
                    else:
                        cy=D.cycles(G);bad=next((z for z in cy if not z['integral']),None)
                        if bad:r.update(status='REJECT_CYCLIC_INTEGRALITY',witness=bad)
                        else:r['status']='ACCEPT_UNPAID';accepted.append(contact)
        trials.append(r)
    status=('EXPANDED_SINGULAR'if singular else'EXPANDED_NONSINGULAR')if accepted else('SEALED_SINGULAR_WEIGHT_OBLIGATION'if singular else'REJECT_SEALED_NONSINGULAR_RANK')
    return dict(status=status,singular=singular,physical=pi,conjugate=co,domain_counts=dc,trials=trials,accepted=accepted)
def run(limit):
    paths=[f'scripts/{TAG}.py',f'scripts/{S.TAG}.py',f'scripts/{PARAM}.py',f'results/{PARAM}.json','scripts/fresh_sqrt5_all_quad_unions_20261015_classes.py','scripts/fresh_sqrt5_sealed_subcore_closure_20261015_review.py','results/arithmetic_dimension_bounds_and_effective_high_dimension_search.md']
    Q.manifest_check(paths);src=json.loads((ROOT/f'results/{PARAM}.json').read_text());states=[];lookup=defaultdict(list);seeds=[]
    def add(row,ratio):
        M=J.matrix(row);key=(J.invariant(M),tuple(ratio.enc()))
        for j,N in lookup[key]:
            mp=J.isomorphism(M,N)
            if mp is not None:break
        else:
            j=len(states);mp=list(range(row['nodes']));states.append(dict(index=j,row=row,ratio=ratio.enc(),result=None));lookup[key].append((j,M))
        target=states[j];C.check_map(row,target['row'],mp);assert{mp[0],mp[2]}=={0,2}
        return dict(index=j,map=mp)
    for ri,r in enumerate(src['records']):
        for ci,c in enumerate(r['accepted']):seeds.append(dict(record=ri,accepted=ci,**add(r['row'],S.dec(c['ratio']))))
    assert len(seeds)==32
    count=Counter();tc=Counter();i=0
    while i<len(states)and i<limit:
        rec=states[i];z=evaluate(rec['row'],S.dec(rec['ratio']));children=[]
        for c in z['accepted']:children.append(dict(contacts=c,**add(F.extension(rec['row'],c),S.dec(rec['ratio']))))
        z['children']=children;rec['result']=z;count[z['status']]+=1;tc.update(r['status']for r in z['trials'])
        print('state',i,'nodes',rec['row']['nodes'],z['status'],'contacts',len(children),'total',len(states),flush=True);i+=1
    obligations=[r['index']for r in states if r['result']and r['result']['singular']]
    return dict(status='COMPLETE_DOTTED_ARITHMETIC_ATTACHMENT_QUEUE_THROUGH_DIMENSION29'if i==len(states)else'PARTIAL_DOTTED_ARITHMETIC_ATTACHMENT_QUEUE_THROUGH_DIMENSION29',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},dimension_bound=29,maximum_star_nodes=31,processed=i,marked_states=len(states),counts=dict(count),trial_counts=dict(tc),seeds=seeds,records=states,weight_obligations=obligations,scope='All32 necessary arithmetic parameter seeds from the24 unpaid minimal dotted cores, for dimensions≤29 only. A full quadrilateral star has d+2 nodes; no relevant core exceeds31nodes. A nonsingular31node core also exceeds ambientrank30. Below this dimension boundary, sealed nonsingular stars are excluded by the full-star rank condition. Every singular state needs payment when no further white touches it; otherwise its accepted child is retained. Orthogonal spectators must fit the same dimension budget. Corner weights, actual realization and independent review are separate obligations. This explicit dimension restriction is justified by the inherited arithmetic upper bound29, not an all-size sealing claim. No improved dimension theorem without payment and the remaining finite-label families.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--max-states',type=int,default=1000);args=ap.parse_args();p=ROOT/f'results/{TAG}.json'
    data=json.loads(json.dumps(run(args.max_states)))
    if args.verify:assert data==json.loads(p.read_text())
    else:p.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['processed'],data['marked_states'],data['counts'],data['trial_counts'])
