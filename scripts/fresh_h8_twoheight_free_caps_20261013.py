#!/usr/bin/env python3
"""Target an exact escape of the new two-height seed with free-distance caps."""
import heapq,json
from fractions import Fraction as F
from itertools import combinations,product
import numpy as np
from scripts.fresh_h8_twoheight_lift_20261013 import *
from scripts.fresh_h8_release_h9_free_caps_20261013 import esign,emul,esub,eencode,bsign,bmul
from scripts.fresh_h8_alternate_selection_20261013 import pairing
from scripts.fresh_h8_alternate_core_angles_independent_20261013 import finite_label
from scripts.fresh_h8_alternate_coupled_20261013 import escape
OUTPUT=ROOT/'data/enumerations/fresh_h8_twoheight_free_caps_20261013.json'

def run():
    base=json.loads(OUT.read_text());G=matrix(base['seed_gram']);I=matrix(base['seed_inverse']);center=tuple(map(parse,base['strict_center_y']))
    gate=base.get('all_caps_escape') or strict_escape(I,[tuple(map(parse,c['coefficients'])) for c in base['caps']],center)
    assert gate['status']=='EXACT_STRICT_SPACELIKE_ESCAPE';y=tuple(map(parse,gate['y']))
    In=np.array([[real(t) for t in r]for r in I]);ic=np.array([real(t)for t in mv(I,center)]);iy=np.array([real(t)for t in mv(I,y)])
    heaps={(t,u):[]for t in (0,1)for u in range(10)};serial=0;counts=[]
    for typ in (0,1):
        norm=scale(O if typ==0 else Q,2);same=RIGHT if typ==0 else LEFT
        fixed=(O,P,add(O,P),scale(P,2));samevals=tuple(mul(Q,x)for x in fixed)if typ else fixed
        cross=(Q,mul(Q,P),scale(Q,2));alphabets=[samevals if i in same else cross for i in range(10)]
        for u in range(10):
            candidates=0;A=In[u,u];others=[i for i in range(10)if i!=u]
            for k in range(5):
                for supp in combinations(others,k):
                    choices=list(product(*(range(len(alphabets[i]))for i in supp)))
                    ws=np.zeros((len(choices),10))
                    for j,row in enumerate(choices):
                        for i,ii in zip(supp,row):ws[j,i]=real(alphabets[i][ii])
                    B=2*(ws@In[:,u]);C=np.einsum('bi,ij,bj->b',ws,In,ws)-real(norm);D=B*B-4*A*C
                    good=D>1e-10
                    for branch in (-1,1):
                        xs=(-B+branch*np.sqrt(np.maximum(D,0)))/(2*A);ws[:,u]=xs
                        cuts=ws@iy;orient=ws@ic
                        ok=good&(xs>np.sqrt(real(norm)*real(G[u][u]))+1e-7)&(cuts>1e-7)&(orient< -1e-7)
                        for j in np.flatnonzero(ok):
                            candidates+=1;serial+=1;score=float(cuts[j]/np.linalg.norm(ws[j]@In));item=(score,serial,(typ,u,supp,choices[int(j)],branch))
                            h=heaps[typ,u]
                            if len(h)<5:heapq.heappush(h,item)
                            elif score>h[0][0]:heapq.heapreplace(h,item)
                    ws[:,u]=0
            counts.append({'norm_type':typ,'free':u,'numeric_cutters':candidates});print(counts[-1],flush=True)
    caps=[];seen=set()
    for score,_,(typ,u,supp,choice,branch)in sorted([r for h in heaps.values()for r in h],reverse=True):
        norm=scale(O if typ==0 else Q,2);same=RIGHT if typ==0 else LEFT;fixed=(O,P,add(O,P),scale(P,2));samevals=tuple(mul(Q,x)for x in fixed)if typ else fixed;cross=(Q,mul(Q,P),scale(Q,2));f=[Z]*10
        for i,j in zip(supp,choice):f[i]=(samevals if i in same else cross)[j]
        f=tuple(f);A=I[u][u];B=scale(mv(I,f)[u],2);C=sub(pair(f,I,f),norm);D=sub(mul(B,B),scale(mul(A,C),4));assert sign(D)>0
        x=(div(neg(B),scale(A,2)),div(scale(O,branch),scale(A,2)));pa=tuple(neg(x[0]if i==u else t)for i,t in enumerate(f));pb=tuple(neg(x[1])if i==u else Z for i in range(10));ca=mv(I,pa);cb=mv(I,pb)
        if (D,ca,cb)in seen:continue
        seen.add((D,ca,cb));q=(add(dot(pa,ca),mul(D,dot(pb,cb))),scale(dot(pa,cb),2));assert esign(esub(q,(norm,Z)),D)==0
        labels=[]
        for i in range(10):
            assert esign((pa[i],pb[i]),D)<=0
            if i==u:assert esign(esub(emul((pa[i],pb[i]),(pa[i],pb[i]),D),(mul(norm,G[i][i]),Z)),D)>0;labels.append('dotted')
            else:
                rel=angle(pa[i],norm,G[i][i]);assert rel not in(None,'parallel');labels.append(rel)
        assert esign((dot(ca,center),dot(cb,center)),D)>0
        assert esign((dot(ca,y),dot(cb,y)),D)<0
        caps.append({'norm':enc(norm),'free_seed_index':u,'fixed_magnitudes':vec(f),'polynomial_descending':vec([A,B,C]),'discriminant':enc(D),'parameter':eencode(x),
                     'pairings_rational':vec(pa),'pairings_radical':vec(pb),'coefficients_rational':vec(ca),'coefficients_radical':vec(cb),'seed_angles':labels,'target_cut':eencode((dot(ca,y),dot(cb,y))),'score':score})
    nodes=[]
    for c in base['caps']:nodes.append({'q':parse(c['norm']),'d':Z,'a':tuple(map(parse,c['coefficients'])),'b':(Z,)*10,'pa':tuple(map(parse,c['pairings'])),'pb':(Z,)*10})
    for c in caps:nodes.append({'q':parse(c['norm']),'d':parse(c['discriminant']),'a':tuple(map(parse,c['coefficients_rational'])),'b':tuple(map(parse,c['coefficients_radical'])),'pa':tuple(map(parse,c['pairings_rational'])),'pb':tuple(map(parse,c['pairings_radical']))})
    data=json.loads((ROOT/'results/fresh_h8_alternate_core_angles_independent_20261013.json').read_text());polys=[(r['label'],tuple(map(F,r['coefficients_descending'])),F(r['largest_root_interval'][0]),F(r['largest_root_interval'][1]))for r in data['real_cyclotomic_polynomials']]
    cn=np.array([[real(x)+np.sqrt(real(n['d']))*real(z)for x,z in zip(n['a'],n['b'])]for n in nodes]);gn=np.array([[real(x)for x in r]for r in G]);pn=cn@gn@cn.T;norms=np.array([real(n['q'])for n in nodes]);rat=4*pn**2/(norms[:,None]*norms[None,:]);targets=np.array([float((lo+hi)/2)for _,_,lo,hi in polys]);nn=[set()for n in nodes];checks={}
    for i,j in combinations(range(len(nodes)),2):
        if pn[i,j]>1e-6 or (rat[i,j]<=4+1e-7 and min(abs(targets-rat[i,j]))>1e-6):continue
        a,b=nodes[i],nodes[j];p=pairing(a,b);d,e=a['d'],b['d'];s=bsign(p,d,e);label=None
        if s>0:kind='positive';ok=False
        else:
            ratio=tuple(mul(div(scale(O,4),mul(a['q'],b['q'])),z)for z in bmul(p,p,d,e));cmp=bsign((sub(ratio[0],scale(O,4)),*ratio[1:]),d,e)
            if cmp>0:kind='dotted';ok=True
            elif cmp==0:kind='parallel';ok=False
            else:label=finite_label(ratio,d,e,polys);ok=label is not None;kind='finite'if ok else'noncoxeter'
        checks[f'{i},{j}']={'compatible':ok,'kind':kind,'label':label}
        if ok:nn[i].add(j);nn[j].add(i)
    cases=[];seen=set();rng=np.random.default_rng(20261013);nold=len(base['caps'])
    for trial in range(min(40,len(caps))):
        chosen=[nold+trial];pool=set(nn[chosen[0]]);noise=rng.random(len(nodes))
        while pool:
            i=max(pool,key=lambda j:(len(nn[j]&pool),noise[j]));chosen.append(i);pool&=nn[i]
        chosen.sort()
        if tuple(chosen)in seen:continue
        seen.add(tuple(chosen));fixed=[nodes[i]['a']for i in chosen if i<nold];new=[caps[i-nold]for i in chosen if i>=nold]
        gate=escape(I,fixed,new,center,tries=1000);cases.append({'chosen_indices':chosen,'facets':10+len(chosen),'volume_gate':gate});print({'free_caps':len(caps),'case':len(cases),'facets':10+len(chosen),'gate':gate['status']},flush=True)
    allgate=escape(I,[n['a']for n in nodes[:nold]],caps,center,tries=1000)
    out={'schema':'two-height-free-distance-caps-v1','base_source':str(OUT.relative_to(ROOT)),'target_escape':base['all_caps_escape'],'target_y':vec(y),'numeric_counts':counts,'caps':caps,'exact_pair_checks':checks,'cases':cases,'all_dictionary_escape':allgate,'scope':'Bounded one-free-distance proposals and exact selected compatibility/escape checks. Missing graph edges are not exhaustive incompatibility certificates. If all_dictionary_escape is strict, it excludes every subset of this finite bank.'}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print({'new_caps':len(caps),'cases':len(cases)},flush=True)

if __name__=='__main__':run()
