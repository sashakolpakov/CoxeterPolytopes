#!/usr/bin/env python3
"""Exact coherent five-cap package selected from changed old-face families."""
import argparse,hashlib,itertools,json,math
from fractions import Fraction as F
from collections import Counter
from fresh_h20_targeted_escape_20261013_third import ROOT,BASE,POLAR,EXT,FIRST,SECOND,line,read,sha
from fresh_dimension_bridge_20261013_verify import E,Z,O,parse,vec,dot,mv
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,fdot,fmv,same,label

THIRD='data/enumerations/fresh_h20_targeted_escape_20261013_third.json'
PREFIX='fresh_h20_edge_family_completion_20261013'
OUT='data/enumerations/'+PREFIX+'_exact.json'
CHOICES=[('scout0',11,F(2)),('scout1',11,F(2)),('scout2',15,F(1)),('scout3',5,F(1)),('three_angles',136,F(1))]

def esqrt(x):
    assert x>=0
    for scale in(1,6):
        y=x/scale;a,b=math.isqrt(y.numerator),math.isqrt(y.denominator)
        if a*a==y.numerator and b*b==y.denominator:return E(F(a,b),F(0))if scale==1 else E(F(0),F(a,b))
    raise AssertionError('Square root not in Qsqrt6')

def positive_solve(M,rhs):
    A=[list(row)for row in M];b=list(rhs);n=len(A);piv=[]
    for i in range(n):
        d=A[i][i]
        if d.sign()<=0:return None
        piv.append(d)
        for j in range(i+1,n):
            u=A[j][i]/d;b[j]-=u*b[i]
            for k in range(i+1,n):A[j][k]-=u*A[i][k]
    x=[Z]*n
    for i in range(n-1,-1,-1):x[i]=(b[i]-sum((A[i][j]*x[j]for j in range(i+1,n)),Z))/A[i][i]
    return x,piv

def extra_pair(f,xc,record):
    g=Extension(f,f.lift(parse(record['field']['radicand'])))
    r=tuple((f.lift(parse(x[0])),f.lift(parse(x[1])))for x in record['normal'])
    p=fdot(g,tuple(map(g.lift,xc)),r)
    return g,p,tuple(f.lift(parse(x))for x in record['norm'])

def run(choices=None,previous_package=None):
    choices=CHOICES if choices is None else choices
    base,polar,ext,first,second,third=map(read,(BASE,POLAR,EXT,FIRST,SECOND,THIRD))
    Q=tuple(map(vec,base['form']));h=vec(base['strict_center']);hc=mv(Q,h)
    old=list(map(vec,base['normals']))+[vec(c['vector'])for c in polar['cuts']];cov=[mv(Q,r)for r in old];ns=[dot(r,c)for r,c in zip(old,cov,strict=True)]
    extras=[{'field':c['field'],'normal':c['vector'],'norm':[c['norm'],['0','0']]}for c in ext['caps']]+[first,second,third]
    out=[];sources=[BASE,POLAR,EXT,FIRST,SECOND,THIRD]
    if previous_package is not None:
        extras.extend(read(previous_package)['new_caps']);sources.append(previous_package)
    retained=82+len(extras)
    for name,index,norm in choices:
        sp='data/enumerations/'+PREFIX+'_'+name+'.json';scout=read(sp);source=scout['caps'][index];target_path=scout['target_source'];geometry=read(target_path);sources.extend([sp,target_path])
        S=source['short_orthogonal_support'];anchors=source['anchors'];labels=source['labels'];ids=S+anchors;assert len(ids)==len(set(ids))==20
        assert all(ns[i].b==0 for i in anchors)
        rhs=[Z]*len(S)+[-esqrt(norm*ns[i].a*{0:F(1),2:F(0),3:F(1,4),4:F(1,2),6:F(3,4)}[m])for i,m in zip(anchors,labels,strict=True)]
        u,v=line([cov[i]for i in ids],rhs);aa=dot(v,mv(Q,v));bb=E(F(2))*dot(u,mv(Q,v));cc=dot(u,mv(Q,u))-E(norm)
        dd=bb*bb-E(F(4))*aa*cc;assert aa.sign()!=0 and dd.sign()>0
        f=Extension(B,dd);Qf=tuple(tuple(map(f.lift,row))for row in Q);fhc=tuple(map(f.lift,hc));accepted=[]
        for branch in(-1,1):
            t=(-bb/(E(F(2))*aa),E(F(branch))/(E(F(2))*aa));r=tuple(f.add(f.lift(x),f.mul(t,f.lift(y)))for x,y in zip(u,v,strict=True));rc=fmv(f,Qf,r);q=fdot(f,r,rc)
            assert same(f,q,f.rational(norm));center=fdot(f,r,fhc)
            if f.sign(center)>=0:continue
            pairs=[]
            try:
                for i,(cv,n)in enumerate(zip(cov,ns,strict=True)):
                    p=fdot(f,r,tuple(map(f.lift,cv)));pairs.append({'old':i,'pairing':f.encode(p),'kind':label(f,p,q,f.lift(n))})
                for i,c in enumerate(extras,82):
                    g,p,n=extra_pair(f,rc,c);pairs.append({'old':i,'pairing':g.encode(p),'kind':label(g,p,g.lift(q),n)})
            except AssertionError:continue
            raw=tuple(E(F(str(x)))for x in geometry['whole_bank_geometry']['point']);target=None
            for power in range(12,61):
                y=tuple(x+E(F(1,2**power))*z for x,z in zip(raw,h,strict=True));yc=mv(Q,y);yq=dot(y,yc);yh=dot(y,hc)
                if yq.sign()<=0 or yh.sign()>=0:continue
                ps=[dot(y,c)for c in cov]
                if any(p.sign()>=0 for p in ps):continue
                ep=[]
                for c in extras:
                    g=Extension(B,parse(c['field']['radicand']));rr=tuple(map(g.parse,c['normal']));p=fdot(g,rr,tuple(map(g.lift,yc)))
                    if g.sign(p)>=0:break
                    ep.append(g.encode(p))
                if len(ep)!=len(extras):continue
                cut=fdot(f,r,tuple(map(f.lift,yc)))
                if f.sign(cut)<=0:continue
                target={'vector':[x.enc()for x in y],'norm':yq.enc(),'center_pairing':yh.enc(),'old_wall_pairings':[x.enc()for x in ps]+ep,'new_pairing':f.encode(cut),'perturbation_power':power};break
            if target:accepted.append((branch,r,rc,q,center,pairs,target))
        assert len(accepted)==1,(name,index,len(accepted));branch,r,rc,q,center,pairs,target=accepted[0]
        cusp_pairs=[fdot(f,r,tuple(map(f.lift,mv(Q,vec(c['ray'])))))for c in base['cusp_rows']]
        assert len(cusp_pairs)==107 and all(f.sign(p)<0 for p in cusp_pairs)
        vertex=None
        for additional in itertools.combinations(anchors,19-len(S)):
            J=S+list(additional);incident=[tuple(map(f.lift,old[i]))for i in J]+[r]
            Mf=[[fdot(f,a,fmv(f,Qf,b))for b in incident]for a in incident]
            if not all(x[1]==Z for row in Mf for x in row):continue
            M=[[x[0]for x in row]for row in Mf];b=[dot(old[i],hc)for i in J]+[center[0]]
            # New center can involve sqrtD; solve the two coefficient systems.
            b0=[fdot(f,a,fhc)[0]for a in incident];b1=[fdot(f,a,fhc)[1]for a in incident]
            sol=positive_solve(M,b0)
            if sol is None:continue
            s0,piv=sol;s1,_=positive_solve(M,b1);lam=tuple(zip(s0,s1,strict=True))
            x=tuple(f.sub(f.lift(h[k]),fdot(f,lam,[a[k]for a in incident]))for k in range(21));xc=fmv(f,Qf,x);xq=fdot(f,x,xc)
            assert f.sign(xq)<0 and f.sign(fdot(f,x,fhc))<0
            vp=[]
            for i,cv in enumerate(cov):
                p=fdot(f,x,tuple(map(f.lift,cv)));assert f.sign(p)==(0 if i in J else-1);vp.append(f.encode(p))
            for c in extras:
                g,p,_=extra_pair(f,xc,c);assert g.sign(p)<0;vp.append(g.encode(p))
            assert f.sign(fdot(f,r,xc))==0
            vertex={'vector':list(map(f.encode,x)),'norm':f.encode(xq),'incident_walls':J+[82+len(extras)],'gram':[[e.enc()for e in row]for row in M],'positive_pivots':[e.enc()for e in piv],'old_pairings':vp};break
        row={'source':sp,'source_index':index,'target_source':target_path,'short_orthogonal_support':S,'anchors':anchors,'labels':labels,
             'field':{'kind':'quadratic','base':{'kind':'sqrt6'},'radicand':dd.enc()},'normal':list(map(f.encode,r)),'norm':f.encode(q),'quadratic_branch':branch,
             'linear_part':[e.enc()for e in u],'line_direction':[e.enc()for e in v],'quadratic_coefficients':[e.enc()for e in(aa,bb,cc)],
             'strict_center_pairing':f.encode(center),'old_pair_checks':pairs,'pair_counts':dict(sorted(Counter(p['kind']for p in pairs).items())),
             'strict_target':target,'inherited_cusp_strict_pairings':list(map(f.encode,cusp_pairs)),'actual_finite_vertex':vertex}
        out.append(row);extras.append(row);print({'cap':len(out),'field':row['field'],'pairs':row['pair_counts'],'finite_vertex':vertex is not None},flush=True)
    return {'status':f'EXACT_COHERENT_{retained+len(out)}_WALL_H20_OLD_FACE_EXTENSIONS_VOLUME_SEPARATE','source_sha256':{p:sha(p)for p in sorted(set(sources))},'retained_walls':retained,'new_caps':out,
            'scope':'Five mutually compatible exact caps from distinct18/17-short direction families. All old/new and new/new angles, strict center,107 cusp signs and sequential strict target cuts are exact. Actual finite vertices are supplied where a PD20 Gram is verified. Large numerical bank and finite volume are not asserted by this artifact.'if previous_package is None else'Exact cap package extending the explicitly retained previous package. Norms, all previous pairings, strict center,107 inherited cusp signs and strict target cuts are exact. Actual finite vertices require saved PD20 Grams. Finite volume is checked separately.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert out==read(OUT)
    else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
