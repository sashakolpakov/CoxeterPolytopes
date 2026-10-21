#!/usr/bin/env python3
"""One bounded solve of the preserved model with labelled incidence/packing.

New nonnegative indicator payments are restricted to their exported full-star
classes. All inherited face inequalities, variables and vertex budgets stay.
Positive-support quad rows are split into all their actual provenances.
"""
import argparse,hashlib,json,time
from collections import Counter,defaultdict
from fractions import Fraction as F
from pathlib import Path
import fresh_sqrt5_forced_base_dual_repair_20261017 as R
import fresh_sqrt5_actual_provenance_dual_20261017 as A
import fresh_sqrt5_vertex_refined_base_capacity_20261017 as V
import fresh_sqrt5_vertex_refined_integer_union_20261017 as U

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_labelled_packing_full_model_20261017'
OUT=ROOT/f'results/{TAG}.json'


def build():
    d,patch,shared=R.domain();actual=R.P.read(A.TAG);vertex=R.P.read(V.TAG);packing=R.P.read(U.PACKING)
    oldfaces=len(d['faces']);oldvars=len(d['variables']);assert (oldfaces,oldvars)==(82530,211339)
    variables=[dict(index=v['index'],kind='inherited',vertex=v['vertex'],cost=v['multiplicity'],source_index=v['index']) for v in d['variables']]
    original=[dict(index=f['index'],required=f['required'],kind=f['kind'],coefficients=f['coefficients']) for f in d['faces']]
    faces=[dict(f,coefficients=list(f['coefficients'])) for f in original]
    copies={};pmap={};basefeatures={}
    for row in actual['merged_row_mass_constraints']:
        fi=row['face'];ids=[]
        assert original[fi]['kind']=='canonical_dotted_quad_aggregate'
        for offset,j in enumerate(row['provenance_indices']):
            assert actual['provenance_masses'][j]['merged_face']==fi
            dest=fi if offset==0 else len(faces)
            if offset:faces.append(dict(original[fi],index=dest,coefficients=list(original[fi]['coefficients'])))
            ids.append(dest);pmap[j]=dest
        copies[fi]=ids
    assert len(pmap)==1858 and len(copies)==1582 and len(faces)==82806
    def destinations(fi):return copies.get(fi,[fi])
    blocks=[];zeros=[]
    def add_positive(label,columns,features):
        start=len(variables)
        for c in columns:
            assert c['index']==len(variables)-start
            variables.append(dict(index=len(variables),kind=label,vertex=c['vertex'],cost=c['vertex_cost'],source_index=c['index']))
        for rec in features:
            for fi in destinations(rec['face']):
                for j,m in rec['coefficients']:
                    assert m>0;faces[fi]['coefficients'].append([start+j,m])
        blocks.append(dict(kind=label,start=start,count=len(columns),feature_rows=len(features)))
    add_positive('mixed_three_port',actual['triple_orbit_columns'],actual['positive_face_mixed_features'])
    add_positive('integer_simplex_packing',packing['columns'],packing['positive_simplex_features'])
    def add_signed(kind,source,fi,cap,positive):
        j=len(variables);assert cap>0
        variables.append(dict(index=j,kind=kind,vertex=None,cost=0,source_index=source))
        assert faces[fi]['kind'] in ('maximal_simplex_aggregate','maximal_edge_incidence')
        faces[fi]['coefficients'].append([j,-cap]);terms=[]
        for dest,m in positive:
            assert m>0 and faces[dest]['kind']=='canonical_dotted_quad_aggregate'
            faces[dest]['coefficients'].append([j,m]);terms.append([dest,m])
        zeros.append(dict(variable=j,simplex_face=fi,capacity=cap,positive_coefficients=terms))
    start=len(variables)
    for c in d['transfers']:
        positive=[[dest,m] for fi,m in c['positive_quad_coefficients'] for dest in destinations(fi)]
        add_signed('inherited_shared_transfer',c['index'],c['simplex_face'],c['capacity'],positive)
    blocks.append(dict(kind='inherited_shared_transfer',start=start,count=len(d['transfers'])))
    start=len(variables);sum_positive=defaultdict(Counter);sum_cap=Counter();sum_sf={}
    for c in vertex['labelled_transfer_columns']:
        positive=[[pmap[j],m] for j,m in c['provenance_coefficients']]
        add_signed('vertex_labelled_base_transfer',c['index'],c['simplex_face'],c['capacity'],positive)
        t=c['base_type'];sum_positive[t].update(dict(c['provenance_coefficients']));sum_cap[t]+=c['capacity'];sum_sf[t]=c['simplex_face']
    blocks.append(dict(kind='vertex_labelled_base_transfer',start=start,count=len(vertex['labelled_transfer_columns'])))
    for c in actual['actual_shared_base_constraints']:
        t=c['base_type']
        assert sum_positive[t]==dict(c['positive_provenance_coefficients'])
        assert sum_cap[t]==c['capacity'] and sum_sf[t]==c['simplex_face']
    assert len(variables)==350744 and len(zeros)==17090
    # Exact inherited-submatrix preservation on every original/copy row.
    new_to_old={dest:fi for fi,ds in copies.items() for dest in ds}
    for fi,face in enumerate(faces):
        assert face['index']==fi
        counts=Counter()
        for j,m in face['coefficients']:counts[j]+=m
        face['coefficients']=[[j,m] for j,m in sorted(counts.items()) if m]
        src=new_to_old.get(fi,fi)
        assert face['required']==original[src]['required'] and face['kind']==original[src]['kind']
        assert [[j,m] for j,m in face['coefficients'] if j<oldvars]==original[src]['coefficients']
        assert all(m>0 for j,m in face['coefficients'] if variables[j]['vertex'] is not None)
    budgets=[dict(index=v['index'],forest=v['forest'],coefficients=list(v['variables'])) for v in d['vertices']]
    for var in variables[oldvars:]:
        if var['vertex'] is not None:budgets[var['vertex']]['coefficients'].append([var['index'],var['cost']])
    assert len(budgets)==6558
    for b in budgets:
        assert sum(m for j,m in b['coefficients'] if j<203565)==78
    active=sorted({j for f in faces for j,m in f['coefficients'] if m>0})
    assert all(variables[j]['cost']>=0 for j in active)
    # A column absent from every positive face entry either only spends
    # budget or has nonpositive signed entries; setting it to zero is safe.
    repair_indices=sorted(j for j in active if j<oldvars and (j<203565 or d['variables'][j].get('kind')=='maximal_simplex_disjoint_port_charge'))
    eligible=set(repair_indices)
    assert all(sum(m for j,m in f['coefficients'] if j in eligible)>0 for f in faces)
    h=hashlib.sha256()
    for group in (variables,faces,budgets,zeros):
        for rec in group:h.update(R.B.dump(rec).encode());h.update(b'\n')
    nnz=sum(len(f['coefficients']) for f in faces)
    paths=[f'scripts/{TAG}.py',f'scripts/{R.TAG}.py',f'results/{R.N.TAG}.json',f'results/{R.S.TAG}.json',
        f'results/{R.N.I.LP}.json',f'results/{A.TAG}.json',f'results/{V.TAG}.json',
        f'results/{U.PACKING}.json',f'results/{U.TAG}.json',f'results/{U.TAG}.md']
    prepared=dict(sources={p:R.P.sha(p) for p in paths},dimension=13,target=12,
        counts=dict(faces=len(faces),old_faces=oldfaces,extra_provenance_rows=len(faces)-oldfaces,
            variables=len(variables),old_positive_cost_variables=oldvars,active_variables=len(active),
            zero_only_variables=len(variables)-len(active),vertex_forests=len(budgets),signed_variables=len(zeros),
            face_matrix_nonzeros=nnz,actual_provenances=len(pmap),affected_merged_quad_rows=len(copies)),
        model_sha256=h.hexdigest(),column_blocks=blocks,
        old_quad_to_actual_faces=[[fi,ds] for fi,ds in sorted(copies.items())],
        provenance_to_face=[[j,fi] for j,fi in sorted(pmap.items())],
        redundant_actual_shared_columns_checked=4969,
        redundancy='Each actual unlabelled shared column is exactly the sum of its vertex-labelled columns, checked coefficient by coefficient including the negative simplex capacity. Separate-k capacities follow from shared capacity because k*C(d+1,k)>=d+1.',
        restricted_features='New mixed and integer indicator columns use their complete exported coefficients only on the selected full-star classes, and zero on other classes. All actual variants of every selected merged Q row are expanded. This intentionally underpays unselected faces; it is a valid restricted charge family, not missing coefficient inference.',
        preservation='Every old variable and all82530 old inequalities retain their original coefficients.1582 affected Q inequalities are replaced by1858 actual variants with identical inherited coefficients. All6558 complete78-pair budgets remain.5693 inherited shared transfer columns remain.',
        scope='Preserved partial face/forest domain with new valid class-restricted charges.231 untouched finite forests, fixed discriminant classes and other missing faces remain outside. No complete dimensional exclusion follows even from a budget12 primal.')
    return dict(variables=variables,faces=faces,budgets=budgets,zeros=zeros,active=active,repair_indices=repair_indices),prepared


def verify_primal(d,c):
    den=c['integer_denominator'];x=dict(c['rounded_variable_numerators']);delta=F(c['repair_delta']);eligible=set(d['repair_indices'])
    assert den>0 and delta>=0 and all(n>=0 for n in x.values())
    assert c['repair_variable_indices']==d['repair_indices']
    for f in d['faces']:
        value=F(sum(m*x.get(j,0) for j,m in f['coefficients']),den)
        value+=delta*sum(m for j,m in f['coefficients'] if j in eligible)
        assert value>=f['required'],(f['index'],value)
    budgets=[]
    for b in d['budgets']:
        value=F(sum(m*x.get(j,0) for j,m in b['coefficients']),den)
        value+=delta*sum(m for j,m in b['coefficients'] if j in eligible)
        budgets.append(value)
    assert str(max(budgets))==c['maximum']
    assert [[v,str(a)] for v,a in enumerate(budgets) if a]==c['vertex_budgets']
    return max(budgets)


def recover_primal(d,values):
    eligible=set(d['repair_indices']);best=None
    for den in (10**8,10**10,10**12):
        x={j:max(0,round(float(a)*den)) for j,a in zip(d['active'],values[:-1])}
        x={j:n for j,n in x.items() if n};delta=F(0)
        for f in d['faces']:
            deficit=f['required']*den-sum(m*x.get(j,0) for j,m in f['coefficients'])
            if deficit>0:
                count=sum(m for j,m in f['coefficients'] if j in eligible);assert count>0
                delta=max(delta,F(deficit,den*count))
        budgets=[]
        for b in d['budgets']:
            value=F(sum(m*x.get(j,0) for j,m in b['coefficients']),den)+delta*sum(m for j,m in b['coefficients'] if j in eligible)
            budgets.append(value)
        c=dict(kind='EXACT_PARTIAL_MODEL_PRIMAL',integer_denominator=den,rounded_variable_numerators=sorted(map(list,x.items())),
            repair_delta=str(delta),repair_variable_indices=d['repair_indices'],maximum=str(max(budgets)),
            vertex_budgets=[[v,str(a)] for v,a in enumerate(budgets) if a],meets_target12=max(budgets)<=12)
        verify_primal(d,c)
        if best is None or F(c['maximum'])<F(best['maximum']):best=c
        if c['meets_target12']:break
    return best


def dual_certificate(d,alpha):
    for c in d['zeros']:
        need=sum((m*alpha.get(i,F(0)) for i,m in c['positive_coefficients']),F(0))/c['capacity']
        fi=c['simplex_face'];alpha[fi]=max(alpha.get(fi,F(0)),need)
    alpha={i:a for i,a in alpha.items() if a};loads=defaultdict(F)
    for fi,a in alpha.items():
        for j,m in d['faces'][fi]['coefficients']:loads[j]+=a*m
    mu=defaultdict(F)
    for c in d['variables']:
        j=c['index']
        if c['vertex'] is None:assert loads[j]<=0,(j,loads[j])
        else:mu[c['vertex']]=max(mu[c['vertex']],loads[j]/c['cost'])
    mu={v:a for v,a in mu.items() if a};num=sum((d['faces'][fi]['required']*a for fi,a in alpha.items()),F(0));den=sum(mu.values(),F(0))
    assert den>0
    return dict(kind='EXACT_PARTIAL_MODEL_DUAL',bound=str(num/den),objective_numerator=str(num),normalizer=str(den),
        face_multipliers=[[i,str(a)] for i,a in sorted(alpha.items())],vertex_multipliers=[[v,str(a)] for v,a in sorted(mu.items())],
        excludes_target12=num/den>12)


def verify_dual(d,c):
    rebuilt=dual_certificate(d,{i:F(a) for i,a in c['face_multipliers']});assert rebuilt==c


def optimize(d):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    active=d['active'];index={j:i for i,j in enumerate(active)};T=len(active)
    rr=[];cc=[];vv=[];rhs=[]
    for f in d['faces']:
        row=len(rhs);rhs.append(-f['required'])
        for j,m in f['coefficients']:
            if j in index:rr.append(row);cc.append(index[j]);vv.append(-m)
    for b in d['budgets']:
        row=len(rhs);rhs.append(0)
        for j,m in b['coefficients']:
            if j in index:rr.append(row);cc.append(index[j]);vv.append(m)
        rr.append(row);cc.append(T);vv.append(-1)
    matrix=coo_matrix((np.asarray(vv,dtype=float),(rr,cc)),shape=(len(rhs),T+1)).tocsr()
    objective=np.zeros(T+1);objective[T]=1
    print('ONE_HIGHS_IPM_BEGIN',matrix.shape,matrix.nnz,'seconds300',flush=True)
    start=time.monotonic();result=linprog(objective,A_ub=matrix,b_ub=np.asarray(rhs,dtype=float),bounds=(0,None),method='highs-ipm',options={'time_limit':300})
    solver=dict(status=int(result.status),message=str(result.message),success=bool(result.success),
        objective=None if result.fun is None else float(result.fun),elapsed_seconds=time.monotonic()-start,
        method='highs-ipm',time_limit_seconds=300,iterations=int(getattr(result,'nit',0)))
    print('ONE_HIGHS_IPM_END',solver,flush=True)
    cert=[]
    if result.x is not None:cert.append(recover_primal(d,result.x))
    marg=getattr(getattr(result,'ineqlin',None),'marginals',None)
    if marg is not None:
        den=10**10;alpha={i:F(max(0,round(-float(marg[i])*den)),den) for i in range(len(d['faces']))}
        cert.append(dual_certificate(d,alpha))
    return solver,cert


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--solve',action='store_true');p.add_argument('--verify',action='store_true');a=p.parse_args()
    assert sum((a.prepare,a.solve,a.verify))==1
    d,prepared=build();prepared=json.loads(json.dumps(prepared))
    print('PREPARED_FULL_LABELLED_PACKING_MODEL',prepared['counts'],prepared['model_sha256'],flush=True)
    if a.verify:
        saved=json.loads(OUT.read_text());assert saved['prepared']==prepared
        for c in saved.get('certificates',[]):
            if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL':verify_primal(d,c)
            else:verify_dual(d,c)
        print('PASS_SAVED_FULL_LABELLED_PACKING_MODEL',saved['status'],flush=True)
    elif a.prepare:
        if OUT.exists():assert not json.loads(OUT.read_text()).get('solver_invocations',0)
        OUT.write_text(R.B.dump(dict(status='PREPARED_FULL_LABELLED_PACKING_MODEL',prepared=prepared,solver_invocations=0,solver=None,certificates=[]))+'\n')
    else:
        previous=json.loads(OUT.read_text());assert previous['prepared']==prepared and previous['solver_invocations']==0
        previous.update(status='ONE_SOLVE_STARTED',solver_invocations=1);OUT.write_text(R.B.dump(previous)+'\n')
        solver,certificates=optimize(d)
        status='PASS_EXACT_PARTIAL_PRIMAL_TARGET12' if any(c.get('meets_target12') for c in certificates) else 'PASS_EXACT_PARTIAL_MODEL_CERTIFICATES' if certificates else 'UNRESOLVED_SINGLE_BOUNDED_SOLVE'
        previous.update(status=status,solver=solver,certificates=certificates);OUT.write_text(R.B.dump(previous)+'\n')
        print(status,[(c['kind'],c.get('maximum',c.get('bound'))) for c in certificates],flush=True)
