#!/usr/bin/env python3
"""Complete integral dotted-parameter domains on24 unpaid minimal cores."""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import fresh_sqrt5_dotted_star_closure_20261015 as S

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_dotted_parameters_20261015'
D=S.D;Q=D.Q;K=D.K
def run():
    paths=[f'scripts/{TAG}.py',f'scripts/{S.TAG}.py',f'results/{D.TAG}.json']
    Q.manifest_check(paths);src=json.loads((ROOT/f'results/{D.TAG}.json').read_text());records=[];counts=Counter()
    for r in src['records']:
        if r['contact_index']!=-1:continue
        row=r['row'];s=r['result']['schur'];a,b,c=[K(*s[v])for v in('a','b','c')]
        M=Q.tree_scaled(row['nodes'],{p:2 if m==0 else m for p,m in Q.edges(row).items()})
        bound=8*(a*b+c*c)/(M[0][0]*M[2][2]);upper=0
        while(K(upper)-bound).sign()<0:upper+=1
        assert(K(upper)-bound).sign()>=0 and(K(upper-1)-bound).sign()<0
        trials=[];accepted=[]
        # t=A+B*phi; t>4>t' implies B>0, and t'≥0 implies A≥0.
        # sqrt5>1 gives B<t≤upper, while A<t≤upper.
        for A in range(upper+1):
            for B in range(1,upper+1):
                t=K(A)+B*Q.PHI
                if(t-4).sign()<=0 or(K(upper)-t).sign()<0 or t.conj().sign()<0 or(4-t.conj()).sign()<0:continue
                trial=dict(A=A,B=B,twice_gram_cycle=t.enc());G,error=S.matrix(row,t/4)
                if error:trial['status']=error
                else:
                    pi=Q.inertia(G);co=Q.inertia([[x.conj()for x in z]for z in G]);trial.update(physical=pi,conjugate=co)
                    if pi[1]>1:trial['status']='REJECT_TWO_NEGATIVE'
                    elif co[1]:trial['status']='REJECT_CONJUGATE'
                    else:
                        n=row['nodes'];assert pi[1]==1 and pi[2]in(0,1)and co[2]==pi[2]and co[0]==n-pi[2]
                        cy=D.cycles(G);bad=next((v for v in cy if not v['integral']),None)
                        if bad:trial.update(status='REJECT_CYCLIC_INTEGRALITY',witness=bad)
                        else:
                            trial.update(status='ACCEPT_SINGULAR'if pi[2]else'ACCEPT_NONSINGULAR',gram=D.enc(G),cycles=cy)
                            accepted.append(dict(trial_index=len(trials),ratio=(t/4).enc(),singular=bool(pi[2])))
                trials.append(trial);counts[trial['status']]+=1
        records.append(dict(source_class_index=r['source_class_index'],union_index=r['union_index'],row=row,schur=s,upper_bound_expression=bound.enc(),integer_upper_bound=upper,trials=trials,accepted=accepted))
        print('class',r['source_class_index'],'upper',upper,'integral_parameters',len(trials),'accepted',len(accepted),flush=True)
    assert len(records)==24
    return dict(status='PASS_COMPLETE_INTEGRAL_DOTTED_PARAMETER_DOMAINS',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},counts=dict(counts),records=records,scope='Complete necessary arithmetic dotted-distance domains on all24 unpaid minimal cores. The two-edge cyclic product t=4cosh(distance)^2 is integral,4<t≤ceil(8(ab+c²)/(d0d2)), and0≤t_conjugate≤4. All field integers in this bounded embedding box are enumerated; physical/conjugate inertia and all integral cyclic products are checked. These are necessary principal cores; arbitrary attachment closure and face payment remain separate. No dimension bound.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'])
