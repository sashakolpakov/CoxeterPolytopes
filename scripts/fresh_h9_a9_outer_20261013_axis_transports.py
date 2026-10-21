#!/usr/bin/env python3
"""Exact frame transports involving the two new transverse endpoints."""
import argparse,hashlib,json
from itertools import product
from fresh_h9_a9_outer_20261013_profiles import ROOT,K,Z,O,G,inverse,mv,dot,pair,menc
from fresh_h9_quartic_mixed_20261013_endpoint import square_root
SOURCE=ROOT/'results/fresh_h9_a9_outer_20261013_axis_geometry.json'
OLD=ROOT/'results/fresh_h9_a9_outer_20261013_frontier.json'
OUT=ROOT/'results/fresh_h9_a9_outer_20261013_axis_transports.json'

def paths(gram):
    adj={i:{j for j in range(9)if i!=j and gram[i][j]!=Z}for i in range(9)}
    assert all(gram[i][j]==-O for i in adj for j in adj[i])
    single=[i for i in adj if not adj[i]];assert len(single)==1
    todo=set(adj)-set(single);ps=[]
    while todo:
        start=min(i for i in todo if len(adj[i])==1);path=[start]
        while True:
            nxt=adj[path[-1]]-set(path)
            if not nxt:break
            assert len(nxt)==1;path.append(nxt.pop())
        assert len(path)==4;todo-=set(path);ps.append(path)
    assert len(ps)==2
    return ps,single[0]

def permutations(epath,fpath):
    ep,ei=epath;fp,fi=fpath
    for swap,r0,r1 in product((False,True),repeat=3):
        targets=fp[::-1]if swap else fp;targets=[targets[0][::-1]if r0 else targets[0],targets[1][::-1]if r1 else targets[1]]
        iso={ei:fi}
        for a,b in zip(ep,targets):iso.update(zip(a,b))
        yield [iso[i]for i in range(9)]

def run(verify=False):
    data,fd=[json.loads(p.read_text())for p in(SOURCE,OLD)]
    walls=[tuple(map(K.from_json,r))for r in data['all_walls']];frames=[]
    for row in fd['vertices'][1:]:frames.append(([walls[i]for i in row['incident_indices']],tuple(map(K.from_json,row['point'])),row['name']))
    for row in data['actual_finite_vertices']:
        frames.append(([walls[i]for i in row['certificate']['indices']],tuple(map(K.from_json,row['point'])),f"axis{row['case']}{'J'if row['reversed']else''}"))
    cached=[]
    for E,p,name in frames:
        ec=[mv(G,r)for r in E];eg=tuple(tuple(dot(r,c)for c in ec)for r in E)
        scale=square_root(-20*(K(-1)+K(0,1))/pair(p,G,p));assert scale is not None and scale>0
        cached.append((ec,eg,paths(eg),scale))
    entries=[];accepted=[];I=tuple(tuple(O if i==j else Z for j in range(10))for i in range(10))
    J=tuple(tuple(O if(i==0 and j==0)or(i>0 and j==10-i)else Z for j in range(10))for i in range(10))
    order=[(1,5),(6,5)]+[(i,j)for i in range(10)for j in range(10)if(i,j)not in((1,5),(6,5))]
    for ai,(E,p,name)in enumerate(frames):
        if ai<23:continue # The earlier exact frame censuses cover these.
        ec,eg,epath,pscale=cached[ai];ei=inverse(eg)
        duals=tuple(mv(ei,tuple(c[j]for c in ec))for j in range(10));pcov=tuple(z/pair(p,G,p)for z in mv(G,p))
        for bi,(F,q,label)in enumerate(frames[:ai+1]):
            fc,fg,fpath,qscale=cached[bi];rho=qscale/pscale
            for perm in permutations(epath,fpath):
                f=[F[i]for i in perm];assert tuple(tuple(fg[i][j]for j in perm)for i in perm)==eg
                row={'source_frame':ai,'target_frame':bi,'permutation':perm,'rho':rho.to_json()};hm={};bad=None
                for i,j in order:
                    value=sum((f[k][i]*duals[j][k]for k in range(9)),Z)+rho*q[i]*pcov[j];hm[i,j]=value
                    if not value.is_integral:bad=(i,j);break
                if bad is not None:
                    i,j=bad;row.update(status='NOT_INTEGRAL',nonintegral_entry={'row':i,'column':j,'value':hm[i,j].to_json()})
                else:
                    H=tuple(tuple(hm[i,j]for j in range(10))for i in range(10));hi=inverse(H);assert all(z.is_integral for rr in hi for z in rr)
                    assert all(mv(H,r)==s for r,s in zip(E,f))and mv(H,p)==tuple(rho*z for z in q)
                    cols=tuple(zip(*H));gcols=[mv(G,c)for c in cols]
                    assert tuple(tuple(dot(c,gc)for gc in gcols)for c in cols)==G
                    kind='IDENTITY'if H==I else'KNOWN_REVERSAL'if H==J else'NEW_INTEGRAL_MAP'
                    row.update(status='INTEGRAL',kind=kind,matrix=menc(H),inverse=menc(hi),trace=sum((H[i][i]for i in range(10)),Z).to_json());accepted.append(row)
                    if kind=='NEW_INTEGRAL_MAP':print('NEW INTEGRAL MAP',name,'to',label,'trace',row['trace'],flush=True)
                entries.append(row)
        print('checked frame',ai,name,flush=True)
    out={'status':'COMPLETE_TRANSVERSE_ORIENTED_FRAME_TRANSPORT_CENSUS','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,OLD)},
      'frames':[name for E,p,name in frames],'previous_frames_excluded':23,'transports':entries,'integral_count':len(accepted),
      'new_integral_count':sum(r['kind']=='NEW_INTEGRAL_MAP'for r in accepted),
      'scope':'Every oriented equal-norm simple-frame identification involving a new corner is tested for integrality, with early exact rejection of a single nonintegral matrix entry. Positive accepted maps preserve the complete source/target chamber germs. Any global conclusion from an infinite-order map is a separate certificate.'}
    if verify:assert json.loads(OUT.read_text())==out
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'maps',len(entries),'integral',len(accepted),'new',out['new_integral_count'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
