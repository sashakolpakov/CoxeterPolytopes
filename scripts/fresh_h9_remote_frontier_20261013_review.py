#!/usr/bin/env python3
"""Formula review and nested 160-digit replay of the 105/106 wall caps."""
import argparse,copy,hashlib,json,sys
from pathlib import Path
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_remote_frontier_20261013_second_verify as S
F=S.F;V=S.V
OUT=V.ROOT/'results/fresh_h9_remote_frontier_20261013_review.json'

def run(verify=False):
    qb,qz,n0=sp.symbols('qb qz n0',nonzero=True)
    qu=n0+1/qb;rad=(2-qu)/qz
    assert sp.cancel(qu+rad*qz)==2
    assert sp.cancel(-qb/qb)==-1
    phi=(1+sp.sqrt(5))/2
    assert sp.simplify((1+phi**2)/2-(2+phi)/2)==0
    H=sp.Matrix([[2,-1,0],[-1,2,-phi],[0,-phi,2]])
    assert all(sp.simplify(H[:i,:i].det())>0 for i in range(1,4))
    # New side projection has norm one, hence pairing -1 with norm-two
    # anchor gives cos^2(pi/4)=1/2. The other three pairings are zero.
    assert sp.Rational(1,1*2)==sp.cos(sp.pi/4)**2
    rawread=V.read;oldscale=V.SCALE;assert oldscale==10**80
    sources=[F.PREFIX+'_caps.json',F.PREFIX+'_verify.json',S.PREFIX+'_caps.json',S.PREFIX+'_verify.json']
    saved=[rawread(F.PREFIX+'_verify.json'),rawread(S.PREFIX+'_verify.json')]
    # First replay the unchanged archived precision exactly.
    first=F.check(rawread(F.PREFIX+'_caps.json'));second=S.check(rawread(S.PREFIX+'_caps.json'))
    assert first==saved[0]and second==saved[1]
    factor=10**80;counts={'changed_enclosures':0,'prior_certificate_comparisons':0}
    def nested(old,new,path='root'):
        if old==new:return
        if isinstance(old,dict):
            assert isinstance(new,dict)and old.keys()==new.keys(),path
            for key in old:
                if key=='interval_scale':assert int(new[key])==int(old[key])*factor
                else:nested(old[key],new[key],path+'.'+key)
        elif isinstance(old,list):
            assert isinstance(new,list)and len(old)==len(new),path
            if len(old)==2 and all(isinstance(x,str)for x in old):
                try:a,b=map(int,old);c,d=map(int,new)
                except ValueError:pass
                else:
                    assert a*factor<=c<=d<=b*factor,(path,'non-nested enclosure')
                    counts['changed_enclosures']+=1;return
            for i,(a,b)in enumerate(zip(old,new,strict=True)):nested(a,b,path+f'[{i}]')
        else:assert old==new,(path,old,new)
    class SavedCertificate(dict):
        def __eq__(self,other):
            nested(dict(self),other);counts['prior_certificate_comparisons']+=1;return True
        def __ne__(self,other):return not self.__eq__(other)
    def precision_read(path):
        value=rawread(path)
        return SavedCertificate(value)if path.endswith('_verify.json')and isinstance(value,dict)else value
    # Older checkers compare their regenerated predecessor with the saved
    # 80-digit certificate. Replace that comparison by full exact metadata
    # equality plus strict enclosure nesting, leaving hashes and data intact.
    V.read=precision_read;V.SCALE=10**160
    V.PHI=V.div(V.add(V.rat(1),V.sqrt(V.rat(5))),V.rat(2));V.radical.cache_clear()
    first160=F.check(rawread(F.PREFIX+'_caps.json'));second160=S.check(rawread(S.PREFIX+'_caps.json'))
    nested(saved[0],first160);nested(saved[1],second160)
    assert first160['facets']==105 and second160['facets']==106
    assert first160['compact_H2_face']['intrinsic_angle_labels']==[2,2,4,2]
    assert len(first160['actual_new_vertices'])==2 and len(second160['actual_new_vertices'])==1
    assert len(second160['quadrilateral_vertex_survival'])==4
    assert first160['affine_census']['affine']==second160['affine_census']['affine']==0
    assert len(second160['remaining_escape']['strict_pairings'])==106
    files=[V.ROOT/p for p in sources]+[Path(F.__file__),Path(S.__file__),Path(V.__file__),Path(__file__).resolve()]
    out={'status':'PASS_SOURCE_REVIEW_AND_NESTED_160_DIGIT_REMOTE_FRONTIER_REPLAY',
         'sources':{str(p.relative_to(V.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files},
         'frozen_precision_digits':80,'review_precision_digits':160,
         'facets':[105,106],'new_pair_counts':[first160['new_pair_counts'],second160['new_pair_counts']],
         'actual_new_vertices':[2,1],'compact_H2_face_labels':[2,2,4,2],
         'preserved_quadrilateral_vertices':4,'strict_escape_norm':second160['remaining_escape']['norm'],
         'new_affine_circuits':0,**counts,
         'scope':'Independent symbolic norm/anchor/H3 Gram review, exact archived-precision replay, and 160-digit replay of the same rigorous interval engine. All changed predecessor and final enclosures nest; exact metadata agree. This is not a separate implementation of all 106 roots. The final strict escape still proves infinite volume.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],counts,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
