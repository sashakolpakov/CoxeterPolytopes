#!/usr/bin/env python3
"""Exact discriminant compatibility of literal representatives in a triangle-only dual."""
import argparse, hashlib, json
from collections import Counter
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import fresh_sqrt5_adjacent_discount_witness_20261016 as A
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D
import fresh_sqrt5_rank13_all_selected_nonade_separation_20261017 as E

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_triangle_discriminant_gate_20261017'
INPUT='fresh_sqrt5_rank13_active_triangle_dual_20261017'
ADE_ENDPOINT='fresh_sqrt5_rank13_high_ade_endpoints_20261017'
ADE_TRIANGLE='fresh_sqrt5_rank13_high_ade_triangles_20261017'
K,Z,O=A.K,A.Z,A.O
PHI=K(F(1,2),F(1,2))

def qsquare(x):
    if x<0:return None
    p,q=isqrt(x.numerator),isqrt(x.denominator)
    return F(p,q) if p*p==x.numerator and q*q==x.denominator else None

def ksqrt(x):
    if not x:return Z
    if not x.b:
        r=qsquare(x.a)
        if r is not None:return K(r)
        r=qsquare(x.a/5)
        if r is not None:return K(0,r)
        return None
    norm=qsquare(x.a*x.a-5*x.b*x.b)
    if norm is None:return None
    for s in (norm,-norm):
        a=qsquare((x.a+s)/2)
        if a:
            r=K(a,x.b/(2*a))
            if r*r==x:return r
    return None

def det_positive(M):
    mat=[row[:] for row in M];det=O
    while mat:
        p=mat[0][0];assert p.sign()>0 and p.conj().sign()>0;det*=p
        mat=[[mat[i][j]-mat[i][0]*mat[0][j]/p for j in range(1,len(mat))] for i in range(1,len(mat))]
    return det

def run():
    d=B.read(INPUT);B.dual_bound=D.sparse_dual_bound
    cert=next(c for c in d['certificates'] if c['kind']=='EXACT_RESTRICTED_MODEL_DUAL')
    B.verify_certificate(d,cert)
    alpha={i:F(a) for i,a in cert['face_multipliers']}
    catalogues={};component_cache={};finite={B.C.code(g):name for name,g in B.C.models() if len(g)<=13}
    variable_lookup={(tuple(v['forest']),B.dump(v['marked_forest_code'])):v['index'] for v in d['variables']}
    def read(tag):
        if tag not in catalogues:catalogues[tag]=B.read(tag)
        return catalogues[tag]
    def literal(provenance):
        if provenance['kind']=='old_triangle':
            witness=read(provenance['source'])['triangle_patterns'][provenance['index']]['witness']
            target=read(ADE_ENDPOINT)['targets'][witness['target']]
            occurrence=read(ADE_TRIANGLE)['targets'][witness['target']]['accepted_occurrences'][witness['occurrence']]
            g={i:{} for i in range(14)}
            for i,j in target['edges']:g[i][j]=g[j][i]=3
            port=occurrence['port'];attach=target['ports'][port]['attachments'][occurrence['attachment']]
            for i in attach['contacts']:g[i][13]=g[13][i]=3
            if occurrence['label']!=2:g[port][13]=g[13][port]=occurrence['label']
            return g,[port,occurrence['other'],13]
        target=read(provenance['source'])['targets'][provenance['target']]
        return E.literal(target,target['triangle_occurrences'][provenance['occurrence']])
    def finite_inverse(piece):
        ids=sorted(piece);key=tuple((ids.index(i),ids.index(j),m) for i in ids for j,m in sorted(piece[i].items()) if i<j)
        key=(len(ids),key)
        if key not in component_cache:
            M=[[K(2) if i==j else -{3:O,5:PHI}[piece[i][j]] if j in piece[i] else Z for j in ids] for i in ids]
            component_cache[key]=(det_positive(M),A.inverse(M))
        return ids,component_cache[key]
    records=[];classes=[];counts=Counter()
    for ordinal,i in enumerate(sorted(alpha)):
        face=d['faces'][i];g,black=literal(face['provenance'][0]);parts=B.C.pieces(g)
        corners=[]
        for removed in black:
            marks=[j for j in black if j!=removed]
            end={j:{k:m for k,m in row.items() if k!=removed} for j,row in g.items() if j!=removed}
            pp=B.C.pieces(end);forest=tuple(sorted(finite[B.C.code(p)] for p in pp))
            code=B.dump(sorted(B.C.code(p,marks) for p in pp));corners.append(variable_lookup[forest,code])
        assert sorted(Counter(corners).items())==[tuple(v) for v in face['coefficients']]
        sizes=sorted(len(p) for p in parts);even=all(s%2==0 for s in sizes)
        record=dict(face=i,original_face=face['original_face'],multiplier=str(alpha[i]),provenance=face['provenance'][0],
            black=black,edges=[[j,k,m] for j,row in g.items() for k,m in row.items() if j<k],component_sizes=sizes,
            fixed_discriminant=even,corner_variables=corners)
        if even:
            scale=[O]*14;special=[(j,k,m) for j,row in g.items() for k,m in row.items() if j<k and m not in (3,5)]
            assert len(special)<=1
            if special:
                j,k,m=special[0];assert k==13 and m in (4,6,10)
                cut={a:{b:z for b,z in row.items() if {a,b}!={j,k}} for a,row in g.items()}
                side=next(p for p in B.C.pieces(cut) if 13 in p);assert j not in side
                z={4:K(2),6:K(3),10:K(2)+PHI}[m]
                for a in side:scale[a]=z
            old={j:{k:m for k,m in row.items() if k!=13} for j,row in g.items() if j!=13}
            vector=[]
            for j in range(13):
                m=g[j].get(13,2)
                if m==2:v=Z
                elif m in (3,5):assert scale[j]==scale[13];v=-scale[j]*{3:O,5:PHI}[m]
                else:v=-scale[13]
                vector.append(v)
            determinant=O;schur=2*scale[13]
            for piece in B.C.pieces(old):
                ids,(dd,inv)=finite_inverse(piece);s=scale[ids[0]]
                assert all(scale[j]==s for j in ids) and s.sign()>0 and s.conj().sign()>0
                for _ in ids:determinant*=s
                determinant*=dd
                v=[vector[j] for j in ids]
                schur-=A.dot(v,A.mv(inv,v))/s
            assert schur.sign()<0 and schur.conj().sign()>0
            determinant*=schur
            assert determinant.sign()<0 and determinant.conj().sign()>0
            match=None;root=None
            for ci,cl in enumerate(classes):
                ratio=determinant/K(*cl['representative'])
                root=ksqrt(ratio)
                if root is not None:match=ci;break
            if match is None:
                match=len(classes);classes.append(dict(index=match,representative=determinant.enc(),first_face=i));root=O
            assert root*root==determinant/K(*classes[match]['representative'])
            record.update(discriminant_class=match,canonical_determinant=determinant.enc(),square_ratio_root=root.enc(),
                root_norm_scales=[x.enc() for x in scale],last_schur_pivot=schur.enc())
            counts['fixed']+=1
        else:counts['flexible_odd_component']+=1
        records.append(record)
        if ordinal%250==0:print('DISCRIMINANT_ROWS',ordinal,len(alpha),len(classes),dict(counts),flush=True)
    # Restrict the actual saved dual; this is an exact lower bound, not a new optimum.
    restrictions=[]
    for ci in list(range(len(classes)))+[None]:
        selected={r['face']:alpha[r['face']] for r in records if not r['fixed_discriminant'] or r['discriminant_class']==ci}
        bound,beta,tight,num,den=D.sparse_dual_bound(d,selected)
        restrictions.append(dict(discriminant_class=ci,meaning='any unrepresented squareclass' if ci is None else 'the displayed represented squareclass',
            retained_faces=len(selected),bound=str(bound),objective_numerator=str(num),normalizer=str(den),
            face_multipliers=[[i,str(a)] for i,a in sorted(selected.items())],vertex_multipliers=[[i,str(a)] for i,a in sorted(beta.items())],
            exceeds12=bound>12,exceeds14=bound>14))
    paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{D.TAG}.py',f'scripts/{E.TAG}.py',
        f'results/{INPUT}.json']+[f'results/{tag}.json' for tag in sorted(catalogues)]
    return dict(status='PASS_EXACT_TRIANGLE_DISCRIMINANT_COMPATIBILITY',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        counts=dict(counts,dual_faces=len(records),represented_fixed_discriminant_classes=len(classes)),
        original_bound=cert['bound'],records=records,classes=classes,restricted_duals=restrictions,
        scope='One saved literal representative per positive dual row is classified. Odd connected components retain their arbitrary totally positive scaling freedom. Each class result only restricts the existing exact dual; a failure to exceed the target proves neither payment feasibility nor ambient realization. Other literal representatives of a merged row may add compatibility; no global dimension theorem is asserted.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()))
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],[(x['discriminant_class'],x['bound']) for x in d['restricted_duals']],flush=True)
