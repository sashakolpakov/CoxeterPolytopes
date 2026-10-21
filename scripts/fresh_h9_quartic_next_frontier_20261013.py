#!/usr/bin/env python3
"""New-field relative-unit proposals for the four corrected outgoing edges."""
from fresh_h9_quartic_edge_walk_20261013_next import *
import hashlib

SOURCE=ROOT/'results/fresh_h9_quartic_first_endpoint_20261013_completion.json'
FRONTIER=ROOT/'results/fresh_h9_quartic_first_endpoint_20261013_frontier.json'
OUT=ROOT/'results/fresh_h9_quartic_next_frontier_20261013.json'


def run():
    started=time.monotonic();data=json.loads(SOURCE.read_text());front=json.loads(FRONTIER.read_text())
    g,fg=Host(2).g,Host(1).g
    walls=[tuple(map(K.from_json,r))for r in data['all_walls']]
    face=[tuple(map(K.from_json,r))for r in data['forced_H5_face']['all_walls']]
    cv=[mv(g,w)for w in walls];qn=[dot(w,c)for w,c in zip(walls,cv)]
    oldkeys={ray_key(w)for w in walls};rows=[];edge_records=[];found={}
    units={1:(K(-5,-2,8,4),K(-2,0,4,2)),
           3:(K(-8561928901343,-4633682531136,14616126887520,7910190870768),
              K(2736975019280,1481240206632,-4672308615352,-2528635201308))}
    for edge_index,record in enumerate(front['edges'][:4]):
        edge_ids=record['edge_indices'];edge=[face[i]for i in edge_ids];released=record['released_starting_wall']
        eg=tuple(tuple(pair(r,fg,s)for s in edge)for r in edge);ei=inverse(eg)
        coeff=mv(ei,tuple(pair(face[released],fg,e)for e in edge))
        base=tuple(face[released][j]-sum((c*e[j]for c,e in zip(coeff,edge)),Z)for j in range(6))
        base_norm=pair(base,fg,base);start=tuple(map(K.from_json,record['start']));direction=tuple(map(K.from_json,record['direction']))
        basis=nullspace(tuple(mv(fg,e)for e in edge));A=tuple(tuple(pair(r,fg,s)for s in basis)for r in basis)
        aa,bb,cc=A[0][0],A[0][1],A[1][1];delta=bb*bb-aa*cc
        multiplier=3 if edge_index<2 else 1;D=multiplier*ALPHA*(2+T);k=square_root(delta/D);assert k and k>0
        B=tuple(zip(*basis));j2=((-bb/k,-cc/k),(aa/k,bb/k))
        P=matmul(matmul(matmul(B,inverse(A)),basis),fg)
        J=matmul(matmul(matmul(matmul(B,j2),inverse(A)),basis),fg)
        assert matmul(P,P)==P and matmul(J,J)==tuple(tuple(D*z for z in row)for row in P)
        y0,s0=units[multiplier];assert y0*y0-D*s0*s0==1 and y0>1
        maxpower=1 if multiplier==3 else 2;y,s=O,Z;translations={0:tuple(tuple(K(i==j)for j in range(6))for i in range(6))}
        for power in range(1,maxpower+1):
            y,s=y*y0+D*s*s0,y*s0+s*y0
            for orientation in(-1,1):
                H=tuple(tuple(K(i==j)+(y-1)*P[i][j]+orientation*s*J[i][j]for j in range(6))for i in range(6))
                assert matmul(matmul(tuple(zip(*H)),fg),H)==fg
                translations[power*orientation]=H
        profiles=[];candidates=[];compatible=[]
        for q in NORMS:
            alph=[sorted({r for ratio in ALLOWED if(r:=square_root(q*pair(e,fg,e)*ratio/4))is not None and r.is_integral})for e in edge]
            for fixed in product(*alph):
                cf=mv(ei,tuple(-z for z in fixed));r0=tuple(sum((c*e[j]for c,e in zip(cf,edge)),Z)for j in range(6))
                residual=q-pair(r0,fg,r0)
                if any(residual.sign(e)<=0 for e in range(4)):continue
                factor=square_root(residual/base_norm);index=len(profiles)
                profiles.append({'norm':q.to_json(),'fixed_magnitudes':enc(fixed),'fixed_component':enc(r0),
                                 'residual':residual.to_json(),'seed_ratio_square':None if factor is None else factor.to_json()})
                if factor is None:continue
                for power,H in translations.items():
                    b=mv(H,base)
                    for sign in(-1,1):
                        r=tuple(x+sign*factor*y for x,y in zip(r0,b));assert pair(r,fg,r)==q
                        if r[0]<=0 or not all(z.is_integral for z in r):continue
                        if all((z/T).is_integral for z in r)or not all((2*z/q).is_integral for z in mv(fg,r)):continue
                        r9=(r[0],*r[1:5],*(Z,)*4,r[5]);key=ray_key(r9)
                        item={'edge_index':edge_index,'profile_index':index,'unit_power':power,'orthogonal_sign':sign,
                              'vector':enc(r9),'norm':q.to_json(),'height':(ALPHA*r[0]*r[0]/q).to_json(),
                              'start_pairing':pair(start,fg,r).to_json(),'direction_pairing':pair(direction,fg,r).to_json()}
                        candidates.append(item)
                        if key in oldkeys or not coxeter(r9,q,cv,qn):continue
                        sp,dp=pair(start,fg,r),pair(direction,fg,r)
                        if sp>=0 or dp<=0:continue
                        t=-sp/dp;point=tuple(x+t*y for x,y in zip(start,direction))
                        if pair(point,fg,point)>=0:continue
                        if any(pair(point,fg,w)>0 for w in face):continue
                        new=sorted(orbit(r9))
                        if not all(coxeter(a,q,[mv(g,b)],[q])for i,a in enumerate(new)for b in new[i+1:]):continue
                        item.update(edge_parameter=t.to_json(),point=enc(point),point_norm=pair(point,fg,point).to_json(),
                                    new_roots=menc(new))
                        compatible.append(item)
                        if key not in found:found[key]=item;rows.append(item)
        edge_records.append({'edge_indices':edge_ids,'released_wall':released,'edge_gram':menc(eg),
             'relative_discriminant':D.to_json(),'complement_basis':menc(basis),'complement_gram':menc(A),
             'discriminant_square_factor':k.to_json(),'projection':menc(P),'generator':menc(J),
             'seed':enc(base),'seed_norm':base_norm.to_json(),'unit_y':y0.to_json(),'unit_s':s0.to_json(),
             'unit_powers':sorted(translations),'profiles':profiles,'integral_candidates':candidates,
             'compatible_outgoing_proposals':compatible})
        print({'edge':edge_ids,'profiles':len(profiles),'square_seed_ratio':sum(p['seed_ratio_square']is not None for p in profiles),
               'integral_candidates':len(candidates),'outgoing':len(compatible),'seconds':round(time.monotonic()-started,3)},flush=True)
    out={'status':'EXACT_NEW_RELATIVE_FIELD_ENDPOINT_PROPOSALS','source':str(SOURCE.relative_to(ROOT)),
         'frontier':str(FRONTIER.relative_to(ROOT)),'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,FRONTIER)},
         'edges':edge_records,'compatible_dictionary':rows,'seconds':round(time.monotonic()-started,3),
         'scope':'Finite Coxeter coupling profiles at all four new edges, and stated powers of exact relative norm-one units in the two new binary norm fields. Only residual/seed square-ratio classes are used for proposals. No completeness of binary norm orbits, nearest-root claim, or finite-volume claim.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n')
    print({'dictionary':len(rows),'seconds':out['seconds']},flush=True)


if __name__=='__main__':run()
