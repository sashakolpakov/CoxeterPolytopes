#!/usr/bin/env python3
"""Independent forward-domain/tree-gauge/cached-corner-Schur queue review.

No imports of the fast producer or its tuple-field/white-leaf arithmetic.
"""
import argparse,json,sqlite3,time
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_singular_extension_frontier_20261015 as F
import fresh_sqrt5_sealed_subcore_closure_20261015_review as J

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_nonsingular_rank29_fast_20261015_review'
PRIMARY='fresh_sqrt5_nonsingular_rank29_fast_20261015'
PREFETCH='data/enumerations/fresh_sqrt5_nonsingular_rank29_fast_20261015_batch3473.json'
Q,T=F.Q,F.T
K=Q.K

def total(xs):return sum(xs,K())
def dot(a,b):return total(x*y for x,y in zip(a,b)if x and y)
def multiply(A,v):return [dot(r,v)for r in A]

def pd_inverse(M):
    """Direct Gauss-Jordan, retaining exact pivots at both embeddings."""
    n=len(M);A=[[K(x)for x in row]+[K(i==j)for j in range(n)]for i,row in enumerate(M)]
    for j in range(n):
        d=A[j][j];assert d.sign()>0 and d.conj().sign()>0
        nz=[k for k in range(j,2*n)if A[j][k]]
        for k in nz:A[j][k]=A[j][k]/d
        for i in range(n):
            if i==j or not A[i][j]:continue
            q=A[i][j]
            for k in nz:A[i][k]=A[i][k]-q*A[j][k]
    inv=[r[n:]for r in A]
    assert all(A[i][j]==K(i==j)for i in range(n)for j in range(n))
    # The finite corner is sparse; this literal product is inexpensive.
    assert all(total(M[i][k]*inv[k][j]for k in range(n)if M[i][k])==K(i==j)for i in range(n)for j in range(n))
    return inv

def inverse2(H):
    d=H[0][0]*H[1][1]-H[0][1]*H[1][0];assert d
    return [[H[1][1]/d,-H[0][1]/d],[-H[1][0]/d,H[0][0]/d]]

def bordered_inverse(Ainv,L,Hinv,ids,outside,n):
    """Inverse of [[A,V],[V^t,B]], where L=A^-1 V."""
    X=[[dot(row,[Hinv[k][j]for k in range(len(outside))])for j in range(len(outside))]for row in L]
    out=[[K()for _ in range(n)]for _ in range(n)]
    for a,i in enumerate(ids):
        for b,j in enumerate(ids):out[i][j]=Ainv[a][b]+dot(X[a],L[b])
        for b,j in enumerate(outside):out[i][j]=out[j][i]=-X[a][b]
    for a,i in enumerate(outside):
        for b,j in enumerate(outside):out[i][j]=Hinv[a][b]
    return out

class Core:
    def __init__(self,row):
        self.row=row;self.n=n=row['nodes'];assert n<=31
        assert row['white']==list(range(4,n)) and F.corners(row)and not F.short(row)and not F.OLD.h4_payments(row)
        M=Q.tree_scaled(n,Q.edges(row));assert M is not None;self.M=M
        assert len(Q.components(Q.graph(range(n),Q.edges(row))))==1
        ids=row['white']+[0,1];outside=[2,3]
        Ainv=pd_inverse(Q.restrict(M,ids));V=[[M[i][j]for j in outside]for i in ids]
        L=[[dot(r,[v[j]for v in V])for j in range(2)]for r in Ainv]
        H=[[M[i][j]-dot([v[a]for v in V],[v[b]for v in L])for b,j in enumerate(outside)]for a,i in enumerate(outside)]
        pi=Q.inertia(H);co=Q.inertia([[x.conj()for x in r]for r in H]);pi[0]+=len(ids);co[0]+=len(ids)
        assert pi[1]==1 and pi[2]in(0,1)and co==[n-pi[2],0,pi[2]]
        self.pi,self.co,self.singular=pi,co,bool(pi[2]);self.kernel=[]
        if not self.singular:
            self.ids=list(range(n));self.inv=bordered_inverse(Ainv,L,inverse2(H),ids,outside,n)
        else:
            q=0 if H[0][0]else 1;assert H[q][q]
            y=[-H[0][1],H[0][0]]if q==0 else[H[1][1],-H[0][1]]
            v=[K()for _ in range(n)]
            for i,r in zip(ids,L):v[i]=-dot(r,y)
            for i,x in zip(outside,y):v[i]=x
            assert any(v)and all(not dot(r,v)for r in M);self.kernel=v
            # Delete the other outside node; this is a nonsingular Lorentz ground.
            ground=ids+[outside[q]];self.ids=ground
            B=bordered_inverse(Ainv,[[r[q]]for r in L],[[K(1)/H[q][q]]],list(range(len(ids))),[len(ids)],len(ground))
            self.inv=B
        self.adj={i:[j for j in range(n)if i!=j and M[i][j]]for i in range(n)}

    def extension(self,G):
        n=self.n;M=self.M;kappa=G[0][0]/M[0][0]
        assert kappa.sign()>0 and kappa.conj().sign()>0
        scales={0:K(1)};stack=[0]
        while stack:
            i=stack.pop()
            for j in self.adj[i]:
                if j not in scales:
                    scales[j]=G[i][j]/(kappa*scales[i]*M[i][j]);stack.append(j)
        assert len(scales)==n and all(x.sign()>0 for x in scales.values())
        for i in range(n):
            for j in range(i,n):
                if M[i][j]:assert G[i][j]==kappa*scales[i]*scales[j]*M[i][j]
                else:assert not G[i][j]
        u=[G[i][-1]/(kappa*scales[i])if G[i][-1]else K()for i in range(n)]
        if self.singular and dot(self.kernel,u):return dict(status='REJECT_INHERITED_RELATION')
        ug=[u[i]for i in self.ids];x=multiply(self.inv,ug)
        s=G[n][n]/kappa-dot(ug,x)
        pi=self.pi[:];co=self.co[:]
        pi[{1:0,-1:1,0:2}[s.sign()]]+=1
        co[{1:0,-1:1,0:2}[s.conj().sign()]]+=1
        if pi[1]>1:status='REJECT_TWO_NEGATIVE'
        elif co[1]:status='REJECT_CONJUGATE'
        elif pi[2]>1:status='REJECT_MULTIPLE_LORENTZ_RELATIONS'
        else:
            assert pi[1]==1 and pi[2]==co[2]
            if self.singular:assert pi[2]==1
            if not self.singular and not s:
                v=[-x[i]/scales[i]for i in range(n)]+[K(1)]
                assert all(not dot(r,v)for r in G)
            elif self.singular:v=self.kernel+[K()]
            else:v=[]
            nonnegative=v and(all(z.sign()>=0 for z in v)or all(z.sign()<=0 for z in v))
            status='REJECT_NONNEGATIVE_ROOT_DEPENDENCE'if nonnegative else'ACCEPT_NECESSARY_UNPAID'
        return dict(status=status,physical=pi,conjugate=co)

def patterns(row):
    g=Q.edges(row);out=[]
    for a,b in product((0,2),(1,3)):
        H=T.graph(row['white']+[a,b],g);remaining=set(H);cs=[]
        while remaining:
            ids=T.component(H,min(remaining));remaining.difference_update(ids)
            cs.append({i:{j:m for j,m in H[i].items()if j in ids}for i in ids})
        c=next(h for h in cs if a in h);connected=b in c
        out.append(dict(marked=[a,b],forest=[T.typename(h)for h in cs],a1_count=sum(len(h)==1 for h in cs),
                        marked_component=T.typename(c)if connected else None,pair_orbit=T.keystr(T.tree_key(c,(a,b)))if connected else None))
    return out

def check_record(rec):
    row,old=rec['row'],rec['result'];C=Core(row)
    assert old['physical']==C.pi and old['conjugate']==C.co and old['singular']==C.singular
    if C.singular:
        assert patterns(row)==old['zero_touch_corner_patterns']
        assert old['orthogonal_spectator_rank_budget']=={str(d):d+2-C.n for d in(28,29)if d+2>=C.n}
    if C.n==31:
        status='DIMENSION29_SINGULAR_WEIGHT_OBLIGATION'if C.singular else'REJECT_RANK_EXCEEDS_DIMENSION29'
        assert old['status']==status and not old['trials']and not old['accepted']
        return 0
    domain,_=F.profiles(row);assert list(map(tuple,domain))==list(map(tuple,sorted(domain)))
    trials={tuple(t['contacts']):t for t in old['trials']};assert len(trials)==len(old['trials'])and set(trials)==set(map(tuple,domain))
    accepted=[]
    for c in domain:
        t=trials[tuple(c)];ext=F.extension(row,c)
        if F.OLD.h4_payments(ext):assert t['status']=='PAID_ENTIRE_H4';continue
        G=Q.tree_scaled(C.n+1,Q.edges(ext))
        if G is None:assert t['status']=='REJECT_FINITE_CYCLE_FIELD';continue
        z=C.extension(G)
        for k in('status','physical','conjugate'):assert t.get(k)==z.get(k),(rec['index'],c,k,t.get(k),z.get(k))
        if z['status']=='ACCEPT_NECESSARY_UNPAID':accepted.append(c)
    assert accepted==old['accepted']
    expected=('EXPANDED_SINGULAR'if C.singular else'EXPANDED_NONSINGULAR')if accepted else('SEALED_SINGULAR_WEIGHT_OBLIGATION'if C.singular else'REJECT_SEALED_NONSINGULAR_RANK')
    assert old['status']==expected
    return len(domain)

def record_hash(rec):return sha256(json.dumps(rec,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def run(prechecked=None):
    path=f'results/{PRIMARY}.json';data=json.loads((ROOT/path).read_text())
    assert data['complete_queue']and data['state']['next_core']==data['marked_states']and not data['counts'].get('PENDING',0)
    checkpoint_path='data/enumerations/fresh_sqrt5_nonsingular_rank29_exhaustion_20261015_batch2007.json'
    checkpoint=json.loads((ROOT/checkpoint_path).read_text());assert checkpoint['state']['next_core']==2007
    for rec in checkpoint['records']:
        target=data['records'][rec['index']];assert target['row']==rec['row']
        if rec['result']is None:continue
        a,b=rec['result'],target['result']
        for k in('status','singular','physical','conjugate','accepted','children','zero_touch_corner_patterns','orthogonal_spectator_rank_budget'):
            assert a.get(k)==b.get(k),(rec['index'],k)
        ta={tuple(t['contacts']):t for t in a['trials']};tb={tuple(t['contacts']):t for t in b['trials']};assert set(ta)==set(tb)
        for c in ta:
            for k in('status','physical','conjugate'):assert ta[c].get(k)==tb[c].get(k),(rec['index'],c,k)
    con=sqlite3.connect(f"file:{ROOT/data['database']}?mode=ro",uri=True)
    assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
    assert sha256((ROOT/data['database']).read_bytes()).hexdigest()==data['database_sha256']
    assert json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0])==data['state']
    assert [json.loads(x[0])for x in con.execute('SELECT payload FROM roots ORDER BY literal_index')]==data['root_maps']
    base=json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text())['labelled_union_supports']
    initial=json.loads((ROOT/'results/fresh_sqrt5_sealed_subcore_closure_20261015.json').read_text())
    initial=[r for r in initial['records']if r['status']=='OPEN_NONSINGULAR'];assert len(initial)==len(data['root_maps'])==311
    reached=set();maps=0
    for i,(raw,z)in enumerate(zip(initial,data['root_maps'])):
        assert z['literal_index']==i and all(z[k]==raw[k]for k in('source_union_index','source_class_index','trial_index','contacts'))
        row=F.extension(base[raw['source_union_index']],raw['contacts'])
        J.check_map(J.matrix(row),J.matrix(data['records'][z['target']]['row']),z['map']);reached.add(z['target']);maps+=1
    assert reached==set(range(199))
    counts,tc=Counter(),Counter();obligations=[];profiles=0
    for i,rec in enumerate(data['records']):
        assert rec['index']==i and i in reached and rec['result']is not None
        row,result=con.execute('SELECT row,result FROM cores WHERE id=?',(i,)).fetchone()
        assert json.loads(row)==rec['row']and json.loads(result)==rec['result']
        if prechecked is not None and i in prechecked:
            digest,number=prechecked[i];assert digest==record_hash(rec);profiles+=number
        else:profiles+=check_record(rec)
        old=rec['result'];counts[old['status']]+=1;tc.update(t['status']for t in old['trials'])
        if old['singular']:obligations.append(i)
        children=old['children'];assert [z['contacts']for z in children]==old['accepted']
        for z in children:
            J.check_map(J.matrix(F.extension(rec['row'],z['contacts'])),J.matrix(data['records'][z['index']]['row']),z['map']);reached.add(z['index']);maps+=1
        if i%100==0:print('independent corner-Schur review',i,flush=True)
    assert reached==set(range(data['marked_states']))and obligations==data['weight_obligations']
    assert dict(counts)==data['counts']and dict(tc)==data['trial_counts']
    assert con.execute('SELECT COUNT(*) FROM cores').fetchone()[0]==data['marked_states'];con.close()
    paths=[f'scripts/{TAG}.py',path,checkpoint_path,PREFETCH,f'scripts/{PRIMARY}.py',f'scripts/{F.TAG}.py',f'results/{F.TAG}.json',
           f'scripts/{Q.TAG}.py',f'results/{Q.TAG}.json',f'scripts/{J.TAG}.py',f'results/{J.TAG}.json',
           'scripts/fresh_sqrt5_sealed_weights_20261015.py','results/fresh_sqrt5_sealed_weights_20261015.json']
    Q.manifest_check(paths)
    return dict(status='PASS_INDEPENDENT_COMPLETE_FINITE_ATTACHMENT_QUEUE_DIMENSION29',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                processed=data['state']['next_core'],marked_states=data['marked_states'],pending=0,complete_queue=True,
                counts=dict(counts),trial_counts=dict(tc),weight_obligations=len(obligations),weight_obligation_indices=obligations,
                exact_profiles=profiles,literal_maps=maps,corner_patterns=4*len(obligations),original_engine_checkpoint_processed=2007,
                scope='Every final state and complete one-white domain independently reconstructed using the older separate Fraction/tree-product field and forward finite-Dynkin classifier. Cached finite-corner Gauss-Jordan and2x2/bordered Schur identities replace the fast producer white-leaf elimination. Gauge congruences checked coefficientwise, both embeddings and all root relations/signs checked. No fast producer arithmetic imports. All311 seeds, every child map, all singular zero-touch forests/spectator budgets and terminal database/counters verified. Necessary principal stars only; payment, other seed families and global dimension inference remain separate.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--prefetch',action='store_true');args=ap.parse_args()
    assert not(args.verify and args.prefetch),'Saved replay reconstructs every arithmetic case afresh.'
    checked=None
    if args.prefetch:
        snap=json.loads((ROOT/PREFETCH).read_text());assert snap['state']['next_core']==3473
        checked={}
        for rec in snap['records'][:3473]:
            number=check_record(rec);checked[rec['index']]=(record_hash(rec),number)
            if rec['index']%100==0:print('independent immutable prefetch',rec['index'],flush=True)
        while True:
            try:
                terminal=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
                if terminal['complete_queue']:break
            except(json.JSONDecodeError,FileNotFoundError):pass
            time.sleep(5)
    data=json.loads(json.dumps(run(checked)));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['processed'],data['weight_obligations'])
