#!/usr/bin/env python3
"""Numerical triage of frozen candidate banks against the exact30-wall chamber."""
import argparse,gzip,json
from hashlib import sha256
from pathlib import Path
import mpmath as mp
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_nine_contact_rescore_20261015'
BASE='results/fresh_h9_nine_contact_extension_20261015.json'
TAGS=['fresh_h9_after_third_vertex_free_scout_20261015',
      'fresh_h9_octic_base_scout_20261015','fresh_h9_octic_germ_scout_20261015']
BASIS=[0,1,2,3,4,5,6,7,8,11]

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def scalar(x):return mp.mpf(x[0][0])+mp.mpf(x[0][1])*mp.sqrt(5)+(mp.mpf(x[1][0])+mp.mpf(x[1][1])*mp.sqrt(5))*mp.sqrt(6+16*mp.sqrt(5)/5)
def smooth(x):return float(format(float(x),'.13g'))

def run():
    mp.mp.dps=100
    d=json.loads((ROOT/BASE).read_text())
    for p,h in d['sources'].items():assert sha(p)==h,p
    Q=mp.matrix([[scalar(x)for x in row]for row in d['form']])
    N=[mp.matrix([scalar(x)for x in row])*mp.sqrt(2/scalar(d['norms'][i]))for i,row in enumerate(d['normals'])]
    M=mp.matrix([[x for x in Q*N[i]]for i in BASIS]);J=M**-1
    assert mp.norm(M*J-mp.eye(10))<mp.mpf('1e-90')
    E=mp.matrix([scalar(x)for x in d['strict_escape']['point']])
    H=mp.matrix([scalar(x)for x in d['center']])
    paths=[BASE,f'scripts/{TAG}.py'];records=[]
    for tag in TAGS:
        meta_path=f'results/{tag}.json';meta=json.loads((ROOT/meta_path).read_text())
        bp=meta['candidate_bank']['path'];assert sha(bp)==meta['candidate_bank']['sha256']
        bank=json.loads(gzip.decompress((ROOT/bp).read_bytes()));paths += [meta_path,bp]
        V=np.array([r['normalized_wall_pairings']for r in bank]);old=V.shape[1]
        assert old in(28,29)and len(bank)==meta['candidate_bank']['count']
        transform=np.array([[float(x)for x in v.T*Q*J]for v in N[old:]+[E,H]])
        new=V[:,BASIS]@transform.T
        finite=-2*np.cos(np.pi/np.arange(2,257))
        ok=new[:,-1]<-1e-8
        for j in range(30-old):
            err=np.min(abs(new[:,j,None]-finite[None,:]),axis=1)
            ok &= (new[:,j]<-2-1e-7)|(err<1e-7)
        cutters=ok&(new[:,-2]>1e-7)
        ids=list(map(int,np.where(cutters)[0]));details=[]
        for i in ids:
            row={k:v for k,v in bank[i].items()if k not in('normalized_wall_pairings','endpoint_other_pairings')}
            details.append(dict(source_index=i,recipe=row,new_pairings=list(map(smooth,new[i]))))
        details.sort(key=lambda r:(-r['new_pairings'][-2],r['source_index']))
        records.append(dict(bank=bp,old_wall_count=old,candidates=len(bank),compatible=int(ok.sum()),
                            escape_cutters=len(ids),maximum_compatible_escape_pairing=smooth(max(new[ok,-2])),
                            accepted_mask_sha256=sha256(np.packbits(ok).tobytes()).hexdigest(),
                            cutter_mask_sha256=sha256(np.packbits(cutters).tobytes()).hexdigest(),cutters=details))
    return dict(status='PASS_BOUNDED_NUMERICAL_30_WALL_RESCORE',sources={p:sha(p)for p in paths},
                physical_precision_digits=100,float_stage='NumPy pairing transform and finite-label gates',
                pairing_basis=BASIS,finite_gate_labels=[2,256],finite_fit_tolerance=1e-7,
                dotted_gap=1e-7,escape_margin=1e-7,center_margin=1e-8,banks=records,
                scope='Only these three frozen finite numerical banks are rescored. Earlier pair gates are inherited, labels are bounded and no interval error certification is supplied. No exact cap or all-label/all-parameter exclusion follows. Exactify surviving recipes before use.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],[(r['candidates'],r['compatible'],r['escape_cutters'])for r in d['banks']])
