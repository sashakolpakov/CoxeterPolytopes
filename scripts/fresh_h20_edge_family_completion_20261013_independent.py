#!/usr/bin/env python3
"""Independent Fraction replay of the coherent93-wall H20 package."""
import argparse,json
from collections import Counter
from copy import deepcopy
from fractions import Fraction as Q
from fresh_h20_targeted_escape_20261013_second_independent import Field,ROOT,BASE,POLAR,CAPS,FIRST,read,sha,add,sub,qm,qs,enc
SECOND='data/enumerations/fresh_h20_targeted_escape_20261013_second.json'
THIRD='data/enumerations/fresh_h20_targeted_escape_20261013_third.json'
SOURCE='data/enumerations/fresh_h20_edge_family_completion_20261013_exact.json'
VOLUME='data/enumerations/fresh_h20_edge_family_completion_20261013_volume.json'
OUT='results/fresh_h20_edge_family_completion_20261013_independent.json'

def qdiv(x,y):
    d=y[0]*y[0]-6*y[1]*y[1];assert d
    z=qm(x,(y[0],-y[1]));return(z[0]/d,z[1]/d)

def check(data):
    retained=data['retained_walls'];assert (retained,len(data['new_caps']))in((88,5),(93,1))
    assert all(sha(p)==s for p,s in data['source_sha256'].items())
    assert all(p in data['source_sha256']for p in(BASE,POLAR,CAPS,FIRST,SECOND,THIRD))
    b,p,c,first,second,third=map(read,(BASE,POLAR,CAPS,FIRST,SECOND,THIRD))
    extras=[{'field':x['field'],'normal':x['vector'],'norm':[x['norm'],['0','0']]}for x in c['caps']]+[first,second,third]
    if retained==93:
        prior='data/enumerations/fresh_h20_edge_family_completion_20261013_exact.json'
        assert data['source_sha256'][prior]==sha(prior)
        extras+=read(prior)['new_caps']
    volume=read(VOLUME);assert volume['source']==SOURCE and volume['source_sha256']==sha(SOURCE)
    results=[]
    for index,row in enumerate(data['new_caps']):
        f=Field(row['field']['radicand']);zero=f.rational(0)
        M=tuple(tuple(f.base(x)for x in r)for r in b['form']);h=tuple(map(f.base,b['strict_center']));hc=f.mv(M,h)
        assert f.sign(f.dot(h,hc))<0
        roots=[tuple(map(f.base,r))for r in b['normals']]+[tuple(map(f.base,r['vector']))for r in p['cuts']]
        cov=[f.mv(M,r)for r in roots];norms=[f.dot(r,v)for r,v in zip(roots,cov,strict=True)];assert len(roots)==82
        r=tuple(map(f.parse,row['normal']));rc=f.mv(M,r);q=f.dot(r,rc)
        assert f.sign(q)>0 and f.sign(sub(q,f.parse(row['norm'])))==0
        center=f.dot(r,hc);assert f.sign(center)<0 and f.sign(sub(center,f.parse(row['strict_center_pairing'])))==0
        parsed=[]
        for previous in extras:
            d=f.base(previous['field']['radicand']);rr=[(f.base(x[0]),f.base(x[1]))for x in previous['normal']]
            parsed.append((d,rr,tuple(f.base(x)for x in previous['norm'])))
        def extra_pair(xc,entry):
            d,rr,n=entry
            return (f.dot(tuple(x[0]for x in rr),xc),f.dot(tuple(x[1]for x in rr),xc))
        def scalar_kind(pp,n):
            assert f.sign(pp)<=0
            if f.sign(pp)==0:return '2'
            square=f.mul(pp,pp);qn=f.mul(q,n);s=f.sign(sub(square,qn))
            if s>0:return 'dotted'
            if s==0:return 'parallel'
            for m,scale in((3,Q(1,4)),(4,Q(1,2)),(6,Q(3,4))):
                if f.sign(sub(square,f.mul(f.rational(scale),qn)))==0:return str(m)
            raise AssertionError('forbidden scalar pair')
        def ext_kind(pp,n,d):
            assert f.extsign(pp,d)<=0
            if f.extsign(pp,d)==0:return '2'
            square=f.extmul(pp,pp,d);qn=tuple(f.mul(q,x)for x in n)
            def delta(scale):return tuple(sub(x,f.mul(f.rational(scale),y))for x,y in zip(square,qn,strict=True))
            s=f.extsign(delta(1),d)
            if s>0:return 'dotted'
            if s==0:return 'parallel'
            for m,scale in((3,Q(1,4)),(4,Q(1,2)),(6,Q(3,4))):
                if f.extsign(delta(scale),d)==0:return str(m)
            raise AssertionError('forbidden extension pair')
        records=[]
        for i,(cv,n)in enumerate(zip(cov,norms,strict=True)):
            pp=f.dot(r,cv);kind=scalar_kind(pp,n);saved=row['old_pair_checks'][i]
            assert saved['old']==i and saved['kind']==kind and f.sign(sub(pp,f.parse(saved['pairing'])))==0
            records.append({'old':i,'kind':kind,'pairing':enc(pp)})
        for i,entry in enumerate(parsed,82):
            d,rr,n=entry;pp=extra_pair(rc,entry);kind=ext_kind(pp,n,d);saved=row['old_pair_checks'][i]
            assert saved['old']==i and saved['kind']==kind
            expected=tuple(map(f.parse,saved['pairing']));assert f.extsign(tuple(sub(x,y)for x,y in zip(pp,expected,strict=True)),d)==0
            records.append({'old':i,'kind':kind,'pairing':[enc(x)for x in pp]})
        assert len(records)==len(row['old_pair_checks'])==retained+index
        counts=dict(sorted(Counter(x['kind']for x in records).items()));assert counts==row['pair_counts']
        assert all(records[i]['kind']=='2'for i in row['short_orthogonal_support'])
        assert all(records[i]['kind']==('parallel'if m==0 else str(m))for i,m in zip(row['anchors'],row['labels'],strict=True))
        cusp=[f.dot(r,f.mv(M,tuple(map(f.base,x['ray']))))for x in b['cusp_rows']]
        assert len(cusp)==107 and all(f.sign(x)<0 for x in cusp)
        assert all(f.sign(sub(x,f.parse(y)))==0 for x,y in zip(cusp,row['inherited_cusp_strict_pairings'],strict=True))
        def signs(x,zeros=()):
            xc=f.mv(M,x)
            for i,cv in enumerate(cov):assert f.sign(f.dot(x,cv))==(0 if i in zeros else-1)
            for entry in parsed:assert f.extsign(extra_pair(xc,entry),entry[0])<0
            return xc
        target=row['strict_target'];w=tuple(map(f.base,target['vector']));wc=signs(w);wq=f.dot(w,wc);wh=f.dot(w,hc);cut=f.dot(r,wc)
        assert f.sign(wq)>0 and f.sign(wh)<0 and f.sign(cut)>0
        assert f.sign(sub(wq,f.base(target['norm'])))==f.sign(sub(wh,f.base(target['center_pairing'])))==f.sign(sub(cut,f.parse(target['new_pairing'])))==0
        vertex=row['actual_finite_vertex'];vr=None
        if vertex:
            J=vertex['incident_walls'];assert len(J)==len(set(J))==20 and J[-1]==retained+index and all(0<=i<82 for i in J[:-1])
            incident=[roots[i]for i in J[:-1]]+[r];G=[[f.dot(a,f.mv(M,b))for b in incident]for a in incident]
            expected=[[f.base(x)for x in line]for line in vertex['gram']]
            assert all(f.sign(sub(x,y))==0 for line,want in zip(G,expected,strict=True)for x,y in zip(line,want,strict=True))
            reduced=[[tuple(map(Q,x))for x in line]for line in vertex['gram']];pivots=[]
            for j in range(20):
                pivot=reduced[j][j];assert qs(pivot)>0;pivots.append(enc(pivot))
                for k in range(j+1,20):
                    for ell in range(j+1,20):reduced[k][ell]=sub(reduced[k][ell],qdiv(qm(reduced[k][j],reduced[j][ell]),pivot))
            x=tuple(map(f.parse,vertex['vector']));xc=signs(x,J[:-1]);xq=f.dot(x,xc)
            assert f.sign(xq)<0 and f.sign(f.dot(x,hc))<0 and f.sign(f.dot(x,rc))==0
            assert f.sign(sub(xq,f.parse(vertex['norm'])))==0
            vr={'incident_walls':J,'norm':enc(xq),'positive_gram_pivots':pivots}
        w=tuple(map(f.base,volume['certificate']['vector']));wc=signs(w);eq=f.dot(w,wc);eh=f.dot(w,hc)
        assert f.sign(eq)>0 and f.sign(eh)<0 and f.sign(f.dot(w,rc))<0
        assert f.sign(sub(eq,f.base(volume['certificate']['norm'])))==f.sign(sub(eh,f.base(volume['certificate']['center_pairing'])))==0
        results.append({'new_wall':retained+index,'norm':enc(q),'pair_counts':counts,'pair_checks':records,'strict_center_pairing':enc(center),'sequential_target_cut_pairing':enc(cut),'preserved_cusps':107,'actual_finite_vertex':vr})
        extras.append(row)
    total=retained+len(results)
    return {'status':f'PASS_INDEPENDENT_{total}_WALL_H20_PACKAGE_VERTICES_TARGETS_AND_ESCAPE','source':SOURCE,'source_sha256':sha(SOURCE),'dependencies':data['source_sha256'],'new_caps':results,
        'new_pair_checks':sum(len(x['pair_checks'])for x in results),'finite_vertices_checked':sum(x['actual_finite_vertex']is not None for x in results),
        'escape':{'source':VOLUME,'source_sha256':sha(VOLUME),'norm':enc(eq),'center_pairing':enc(eh),'strict_inequalities':total},
        'scope':('Separate Fraction arithmetic checks all450 new pairs,5 norms,strict center,535 cusp signs,5 sequential strict target cuts,3 actual finite vertices and all93 strict escape inequalities. Existing88-wall certificate is a frozen dependency. Large numerical banks are not asserted exact.'if retained==88 else'Separate Fraction arithmetic checks all93 new pairs,norm1,strict center,107 cusp signs,a strict target cut,an actual finite vertex and all94 strict escape inequalities. Existing93-wall certificate is a frozen dependency. Large numerical banks are not asserted exact.')}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();data=read(SOURCE);out=check(data)
    for field in('normal','norm','strict_target'):
        bad=deepcopy(data);r=bad['new_caps'][0]
        if field=='normal':r[field][0][0][0]=str(Q(r[field][0][0][0])+1)
        elif field=='norm':r[field][0][0]='0'
        else:r[field]['new_pairing'][0][0]=str(Q(r[field]['new_pairing'][0][0])+1)
        try:check(bad)
        except AssertionError:pass
        else:raise AssertionError('forged '+field+' accepted')
    out['forged_normal_norm_cut_rejected']=True
    if a.verify:assert read(OUT)==out
    else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_pair_checks'],out['finite_vertices_checked'])
