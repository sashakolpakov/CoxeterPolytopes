#!/usr/bin/env python3
"""Physical interval rejection of the declared low-S/low-R continuous-N domain.

The numerical packet supplies recipe indices and sheet names, not rejection
proofs. Every excluded branch gets a strict physical native-wall witness;
all remaining branches get complete interval event quadratics.
"""
import argparse,json,math
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import numpy as np
import fresh_oblique_high_s_polar_interval_20261018 as P
import fresh_oblique_high_s_low_r_interval_20261018 as Q
ROOT=P.ROOT;TAG='fresh_oblique_low_s_continuous_n_interval_20261018'
DISCOVERY='fresh_oblique_low_s_continuous_n_events_20261018'
V,add,sub,mul,div,neg,sqrt,dot=P.V,P.add,P.sub,P.mul,P.div,P.neg,P.sqrt,P.I.dot


@lru_cache(None)
def theta(m):
    assert isinstance(m,int)and m>=2
    if m<=6:return P.theta(m)
    def atan(x):
        z=sum((-1)**k*x**(2*k+1)/F(2*k+1)for k in range(100))
        return z,z+x**201/F(201)
    pi=sub(mul(V(16),atan(F(1,5))),mul(V(4),atan(F(1,239))))
    x=div(pi,V(m));xx=mul(x,x);term=value=V(1)
    for k in range(1,60):
        term=div(mul(neg(xx),term),V((2*k-1)*(2*k)));value=add(value,term)
    error=xx[1]*max(abs(q)for q in term)/F(119*120)
    return mul(V(2),(value[0]-error,value[1]+error))


def refined_primary(prior):
    """Use only exact identities of the frozen original Coxeter scaffold.

    Original wall norms are2. An original off-diagonal enclosure isolating
    zero or one finite native label determines that product exactly. Dotted
    products are never rounded to an algebraic value.
    """
    src=prior['primary_projection'];ids=src['walls'];native={int(i):{int(j):V(z)for j,z in r.items()}for i,r in src['physical_Gram'].items()};proofs=[]
    for i in ids:
        for j in ids:
            if j<i:continue
            z=native[i][j];exact=None;label=None
            if i==j:exact=V(2);label='norm2'
            elif z[0]>-1 and z[1]<1 and z[0]<=0<=z[1]:exact=V(0);label=2
            elif -2<z[0]<=z[1]<0:
                m=int(round(math.pi/math.acos(-float(sum(z)/2)/2)))
                if m>=3 and z[0]>-theta(m+1)[0]and z[1]<-theta(m-1)[1]:exact=neg(theta(m));label=m
            if exact is None:continue
            intersection=max(z[0],exact[0]),min(z[1],exact[1]);assert intersection[0]<=intersection[1]
            native[i][j]=native[j][i]=intersection
            proofs.append(dict(pair=[i,j],old=z,identity=label,identity_interval=exact,intersection=intersection))
    wi,pivots=P.I.inverse([[native[a][b]for b in P.I.WHITE]for a in P.I.WHITE]);assert all(p[0]>0 for p in pivots)
    wc={a:[native[a][w]for w in P.I.WHITE]for a in P.OLD}
    g={a:{b:sub(native[a][b],dot(wc[a],wi,wc[b]))for b in P.OLD}for a in P.OLD}
    return g,wc,wi,dict(identities=proofs,physical_Gram=native,white_inverse=wi,white_pivots=pivots,old_white_covectors=wc,projected_Gram=g)


def parent_extension(g,wc,wi,profiles,row):
    wall,frame,pi,labels=(row[k]for k in('wall','frame','profile','labels'));gram=[[g[i][j]for j in frame]for i in frame];inv,det=P.inverse(gram)
    assert det[1]<0
    rhs=[sub(neg(theta(m)),dot(profiles[pi],wi,wc[i]))for i,m in zip(frame,labels)]
    square=sub(V(2),dot(profiles[pi],wi,profiles[pi]));rem=sub(square,dot(rhs,inv,rhs))
    w,sign=row['witness'],row['witness_sign'];cv=[g[i][w]for i in frame];q=sub(g[w][w],dot(cv,inv,cv));assert q[0]>0 and rem[0]>0
    lin=dot(rhs,inv,cv);rad=sqrt(mul(rem,q));fac=mul(V(sign),sqrt(div(rem,q)));out={}
    for j in g:
        cj=[g[i][j]for i in frame];out[j]=rhs[frame.index(j)]if j in frame else add(dot(rhs,inv,cj),mul(fac,sub(g[w][j],dot(cv,inv,cj))))
    proof=dict(wall=wall,frame=frame,profile=pi,labels=labels,frame_Gram=gram,inverse=inv,determinant=det,rhs=rhs,projected_square=square,
        remainder=rem,witness=w,witness_residual=q,branch_witness_products=[sub(lin,rad),add(lin,rad)],witness_sign=sign,new_products=out)
    for j,z in out.items():g[j][wall]=z
    g[wall]={**out,wall:square};wc[wall]=profiles[pi]
    return proof


def extension(g,wc,wi,profiles,numerical,wall,frame,profile,labels,reference):
    """Define a physical root by its witness sheet, independently of its float.

    Numerical proximity only names the sheet. Each recipe is required below
    to use BOTH signs, so even a changed numerical rejection is not inherited.
    """
    gram=[[g[a][b]for b in frame]for a in frame];inv,det=P.inverse(gram)
    assert det[1]<0 and gram[0][0][0]>0
    rhs=[sub(neg(theta(m)),dot(profiles[profile],wi,wc[i]))for i,m in zip(frame,labels)]
    square=sub(V(2),dot(profiles[profile],wi,profiles[profile]));rem=sub(square,dot(rhs,inv,rhs));assert rem[0]>0
    candidates=[]
    for w in g:
        if w in frame:continue
        cv=[g[i][w]for i in frame];q=sub(g[w][w],dot(cv,inv,cv))
        if q[0]>0:candidates.append((q[0],w,q,cv))
    _,w,q,cv=max(candidates);linear=dot(rhs,inv,cv);radical=sqrt(mul(rem,q))
    branches=[sub(linear,radical),add(linear,radical)]
    observed=float(reference@P.H.T.M@numerical[w]);distances=[abs(observed-float(sum(z)/2))for z in branches]
    selected=int(np.argmin(distances));assert distances[selected]<distances[1-selected]/100
    sign=(-1,1)[selected];factor=mul(V(sign),sqrt(div(rem,q)))
    out={}
    for j in g:
        cj=[g[i][j]for i in frame]
        out[j]=rhs[frame.index(j)]if j in frame else add(dot(rhs,inv,cj),mul(factor,sub(g[w][j],dot(cv,inv,cj))))
    error=max(abs(float(sum(z)/2)-float(reference@P.H.T.M@numerical[j]))for j,z in out.items())
    proof=dict(wall=wall,frame=frame,profile=profile,labels=labels,frame_Gram=gram,inverse=inv,determinant=det,rhs=rhs,
        projected_square=square,remainder=rem,witness=w,witness_residual=q,branch_witness_products=branches,witness_sign=sign,
        numerical_witness_product=observed,branch_comparison=distances,all_product_comparison_error=error,new_products=out)
    for j,z in out.items():g[j][wall]=z
    g[wall]={**out,wall:square};wc[wall]=profiles[profile];numerical[wall]=reference
    return proof


def blocker(g,wc,wi,wall,oldkeys):
    for j in oldkeys:
        z=add(g[wall][j],dot(wc[wall],wi,wc[j]))
        if z[0]>0:return dict(wall=j,product=z,kind='positive')
        if z[0]>-1 and z[1]<0:return dict(wall=j,product=z,kind='between_right_and_simple')
        mid=-float(sum(z)/2)
        if 1<mid<2:
            m=int(math.floor(math.pi/math.acos(mid/2)))
            if m>=3 and z[0]>-theta(m+1)[0]and z[1]<-theta(m)[1]:
                return dict(wall=j,product=z,kind='between_consecutive_native_labels',labels=[m,m+1],
                    lower_magnitude=theta(m),upper_magnitude=theta(m+1))
    return None


def event(g,a):
    frame=[10,67,69];inv,det=P.inverse([[g[i][j]for j in frame]for i in frame])
    r0=[V(-a),V(0),V(0)];r1=[V(0),V(0),V(-1)]
    cv=[g[i][64]for i in frame];rv=[g[i][70]for i in frame]
    sig=sub(g[64][64],dot(cv,inv,cv));assert sig[0]>0
    aa=dot(r0,inv,rv);bb=dot(r1,inv,rv);cc=div(sub(g[64][70],dot(cv,inv,rv)),sqrt(sig))
    rho=[sub(V(2),dot(r0,inv,r0)),neg(mul(V(2),dot(r0,inv,r1))),neg(dot(r1,inv,r1))]
    cc2=mul(cc,cc)
    p=[sub(mul(aa,aa),mul(cc2,rho[0])),sub(mul(V(2),mul(aa,bb)),mul(cc2,rho[1])),sub(mul(bb,bb),mul(cc2,rho[2]))]
    assert p[2][0]>0
    disc=sub(mul(p[1],p[1]),mul(V(4),mul(p[2],p[0])));assert disc[0]>0
    roots=[div(sub(neg(p[1]),sqrt(disc)),mul(V(2),p[2])),div(add(neg(p[1]),sqrt(disc)),mul(V(2),p[2]))]
    decisions=[]
    for rt in roots:
        if rt[1]<theta(6)[0]:decision=dict(kind='below_high_domain')
        elif rt[0]>2:decision=dict(kind='above_finite_domain')
        else:
            mid=float(sum(rt)/2);m=int(math.floor(math.pi/math.acos(mid/2)))
            lo,hi=theta(m),theta(m+1)
            assert rt[0]>lo[1]and rt[1]<hi[0],('unresolved_native_event',a,rt,m)
            decision=dict(kind='strict_native_label_gap',labels=[m,m+1],lower_magnitude=lo,upper_magnitude=hi,
                lower_gap=rt[0]-lo[1],upper_gap=hi[0]-rt[1])
        decisions.append(dict(root=rt,**decision))
    return dict(a=a,N_frame=frame,N_inverse=inv,N_determinant=det,old64_residual_square=sig,
        coefficients=[aa,bb,cc],norm_remainder_coefficients=rho,polynomial=p,discriminant=disc,roots=decisions)


def run(saved=None):
    numerical=json.loads((ROOT/f'results/{DISCOVERY}.json').read_text())
    prior=json.loads((ROOT/f'results/{P.TAG}.json').read_text())
    for d in(numerical,prior):
        for path,h in d['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    assert len(numerical['records'])==len(prior['parents'])==2
    oldN=json.loads((ROOT/'results/fresh_oblique_continuous_high_n_high_s_interval_20261018.json').read_text())
    assert oldN['summary']['covering_boxes']==66 and oldN['summary']['unresolved']==0
    assert F(oldN['summary']['minimum_opposite_N_old8'])>0
    assert oldN['N_families']==[[2,2,'m>=6'],[3,2,'m>=6']]and oldN['parents']==prior['parents']
    profiles=[[V(q)for q in row]for row in prior['profiles']]
    g0,wc0,wi,refined=refined_primary(prior)
    oldnum={int(i):np.array(v)for i,v in numerical['numerical_snapshot']['old_normals'].items()}
    records=[];S_count=R_count=S_rej=R_rej=0;eqs=[]
    for nr,pr in zip(numerical['records'],prior['parents']):
        assert nr['index']==pr['index']and nr['instruction_key']==pr['instruction_key']
        g={i:r.copy()for i,r in g0.items()};wc=wc0.copy()
        refined_steps=[parent_extension(g,wc,wi,profiles,step)for step in pr['steps']]
        assert len(pr['steps'])==5 and 70 in g
        del g[70],wc[70]
        for row in g.values():del row[70]
        nums={**oldnum,**{int(i):np.array(v)for i,v in nr['prefix_normals'].items()}}
        retained=[i for i in g if i not in(63,64)]
        row=dict(index=nr['index'],origins=nr['origins'],instruction_key=nr['instruction_key'],refined_parent_steps=refined_steps,S=[]);records.append(row)
        for sq in nr['S_bank']['recipes']:
            assert len(sq['branches'])==2
            signset=[]
            for sb in sq['branches']:
                gs={i:r.copy()for i,r in g.items()};ws=wc.copy();ns={i:nums[i]for i in g}
                ext=extension(gs,ws,wi,profiles,ns,72,[52,68,71],sq['profile'],sq['labels'],np.array(sb['normal']))
                signset.append(ext['witness_sign']);S_count+=1
                block=blocker(gs,ws,wi,72,retained)
                sr=dict(recipe=sq['index'],branch=sb['branch'],profile=sq['profile'],labels=sq['labels'],extension=ext,blocker=block,R=[]);row['S'].append(sr)
                if block:S_rej+=1;continue
                # A branch retained here is only an overdomain candidate. All
                # later strict proofs are independent of its numerical gates.
                assert not sb['blockers']and 'cap'in sb,('unresolved_S',nr['index'],sq['index'],sb['branch'])
                si=sb['cap'];bank=nr['R_banks'][si];assert bank['S']==si;sr['numerical_cap_index']=si
                for rq in bank['recipes']:
                    assert len(rq['branches'])==2
                    rsigns=[]
                    for rb in rq['branches']:
                        gr={i:r.copy()for i,r in gs.items()};wr=ws.copy();nn=ns.copy()
                        rext=extension(gr,wr,wi,profiles,nn,70,[69,71,72],rq['profile'],rq['labels'],np.array(rb['normal']))
                        rsigns.append(rext['witness_sign']);R_count+=1
                        rblock=blocker(gr,wr,wi,70,retained+[72])
                        rr=dict(recipe=rq['index'],branch=rb['branch'],profile=rq['profile'],labels=rq['labels'],extension=rext,blocker=rblock,events=[]);sr['R'].append(rr)
                        if rblock:R_rej+=1;continue
                        assert not rb['blockers'],('unresolved_R',nr['index'],si,rq['index'],rb['branch'])
                        rr['events']=[event(gr,a)for a in(0,1)];eqs.extend(rr['events'])
                    assert sorted(rsigns)==[-1,1]
            assert sorted(signset)==[-1,1]
        print('PHYSICAL_PARENT',row['index'],'S',S_count,'R',R_count,'events',len(eqs),flush=True)
    assert(S_count,R_count,S_rej,R_rej,len(eqs))==(168,2476,92,2416,120)
    roots=[rt for e in eqs for rt in e['roots']];gaps=[rt for rt in roots if rt['kind']=='strict_native_label_gap']
    assert len(roots)==240 and len(gaps)==4 and{tuple(rt['labels'])for rt in gaps}=={(8,9),(10,11)}
    paths={Path(__file__),Path(P.__file__),Path(Q.__file__),Path(P.I.__file__)}
    paths.update(ROOT/p for d in(numerical,prior)for p in d['sources'])
    paths.update(ROOT/f'{folder}/{stem}.{ext}'for stem in(DISCOVERY,P.TAG)for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    result=P.I.enc(dict(status='PASS_PHYSICAL_LOW_S_LOW_R_CONTINUOUS_HIGH_N_EXCLUSION',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        interval_scale=P.I.PREC,parameter_domain=[theta(6)[0],F(2)],refined_primary=refined,profiles=prior['profiles'],parents=prior['parents'],records=records,
        exact_deletion_debit=numerical['exact_deletion_debit'],deleted_walls=numerical['deleted_walls'],topology=numerical['topology'],
        summary=dict(parents=2,parent_origins=4,S_recipes=84,S_branches=S_count,S_strict_rejections=S_rej,S_overdomain_survivors=S_count-S_rej,
            R_recipes=1238,R_branches=R_count,R_strict_rejections=R_rej,R_overdomain_survivors=R_count-R_rej,event_quadratics=len(eqs),event_roots=len(roots),
            outside_domain=sum(rt['kind']!='strict_native_label_gap'for rt in roots),native_label_gaps=len(gaps),
            minimum_native_gap=min(min(rt['lower_gap'],rt['upper_gap'])for rt in gaps),unresolved=0,accepted=0,new_cones=0),
        scope='Two fixed physical parent chains; both high-N families for all finite labels>=6; complete declared low-S cube and RXright RM/RS2..5 finite duties. Every earlier numerical rejection replaced by a strict physical native-wall interval, every residual event eliminated by full quadratic root intervals. Both N sheets covered using frozen opposite-sheet exclusion. Historical scaffold and seven-profile/finite-duty classification inherited. High RM/RS and other parents/incidences remain outside this certificate.'))
    if saved:assert result==saved
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    assert args.verify or not path.exists();out=run(json.loads(path.read_text())if args.verify else None)
    if not args.verify:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
