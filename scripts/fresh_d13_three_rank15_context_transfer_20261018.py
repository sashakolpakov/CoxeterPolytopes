#!/usr/bin/env python3
"""Ambient58-only D13 halves, with two narrowly witnessed finite repairs.

No forest optimization: the frozen complete eight-obstruction certificate
is combined with exact costs of the only two increasing forest families.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,combinations_with_replacement,product
import sympy as sp
import fresh_b13_endpoint_retention_20261018 as Z
import fresh_d17_a16_coupled_ambient58_frontier_20261018 as C
import fresh_affine_d12_leaf_pair_transfer_20261017 as U
import fresh_d12_remaining_selected_half_replacements_20261018 as E

ROOT=Z.ROOT;S=Z.S;T=S.T;A=S.A;Q=S.Q;D=S.D;B=S.B
TAG='fresh_d13_three_rank15_context_transfer_20261018'
PAIRS={(0,11),(0,12)}
STEMS=tuple(dict.fromkeys((Z.TAG,C.TAG,U.TAG,E.TAG,*Z.STEMS,
 'fresh_distance9_supporting_extensions_review_20261017')))
enc=S.enc
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def context_delta(fs):
    co=Counter(fs);k=co['A',15]+co['D',15]
    return (-64*co['D',13]if k>=3 else 0)+(32*co['D',28]+64*co['A',12]*co['D',16]if k>=2 else 0)
def value(fs):return Z.value(fs)+context_delta(fs)

def support(over):
    target=T.target(13,'D');rows=[];hits=[];maps=0
    for c in over['isolation']['cases']:
        if(c.get('component_order_bound')or 0)<13 or len(c['marked_component'])>13:continue
        core=over['isolation']['cores'][c['core_key']]
        mm=T.T.mapped(A.matrix(core['gram']),c['marked_component'],c['preferred_corner'],target)
        aa=[m for m in mm if tuple(m['marked_images'])in PAIRS]
        rows.append(dict(index=c['index'],maps=mm,affected_maps=aa));maps+=len(mm)
        if aa:hits.append((c,core,aa))
    assert len(rows)==801 and maps==862 and len(hits)==33
    assert sum(len(x[2])for x in hits)==76
    return dict(patterns=801,maps=862,affected_contexts=33,affected_maps=76,rows=rows),hits

def basic(c,core):
    g=A.matrix(core['gram']);n=len(g)
    assert core['gauge']==[1]*n
    left=c['left'][1];right=[c['right_copy_map'][i]for i in c['right'][1]]
    white=sorted(set(range(n))-set(left+right))
    assert Q.inertia(g)==core['signature']
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,white+[b]))==[len(white)+1,0,0]for b in left+right)
    actual=Q.inspect_corner(g,core['gauge'],left+right,*c['preferred_corner'])
    for k in('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):assert actual[k]==c[k]
    blocks=D.finite_blocks(g,sorted(set(white)|set(c['preferred_corner'])))
    small=[b for b in blocks if not set(c['preferred_corner'])&set(b['nodes'])]
    return g,white,left,right,small

def selected_proof(c,core,maps):
    g,white,left,right,small=basic(c,core);n=len(g);nodes=c['marked_component']
    common=dict(index=c['index'],core_key=c['core_key'],gram=g,white=white,
        selected_marks=c['preferred_corner'],selected_nodes=nodes,displayed_unmarked_blocks=small)
    if c['method']=='radical':
        assert len(nodes)==13 and core['signature']==[n-2,1,1]
        z=list(map(F,core['radical'][0]));assert all(sum(g[i][j]*z[j]for j in range(n))==0 for i in range(n))
        free=sorted(set(c['possible_contacts'])-set(nodes));assert all(z[i]for i in free)
        du,dv=Q.graph_dist(g,left),Q.graph_dist(g,right);esc=[]
        for i,j in combinations(free,2):
            if z[i]*z[j]<0:
                length=min(du[i]+dv[j],du[j]+dv[i])+2;assert length<=8
                esc.append(dict(nodes=[i,j],length=length))
        assert small and all(len(b['nodes'])<=4 for b in small)
        return common|dict(method='radical_full_core_seal',radical=z,remaining_positions=free,
            mixed_sign_shortcuts=esc,whole_small_factor=small[0],maximum_small_rank=4)
    gg,J,du,dv,end,ww,ll,rr=E.validate(c,core);assert gg==g and ww==white and ll==left and rr==right
    if len(nodes)==12:
        forced=[];target=T.target(13,'D')
        for mp in maps:
            image={j:i for i,j in zip(nodes,mp['images'])};missing=sorted(set(range(13))-set(image));assert len(missing)==1
            contacts=[(image[j],target[missing[0]][j])for j in image if target[missing[0]][j]]
            assert len(contacts)==1 and contacts[0][1]==3
            q=contacts[0][0];assert q in c['possible_contacts']and J[q][q]>2
            forced.append(dict(images=mp['images'],marked_images=mp['marked_images'],missing=missing,
                forced_contact=q,projection_lower=J[q][q]))
        return common|dict(method='genuine_new_fork_leaf_projection_exclusion',end_projection=end,maps=forced)
    assert len(nodes)==13 and small
    free=sorted(set(c['possible_contacts'])-set(nodes));possible={i for i in free if J[i][i]<=2}
    assert all(J[i][i]>0 for i in free)
    common|=dict(end_projection=end,remaining_inverse_diagonals=[[i,J[i][i]]for i in free],possible_contacts=sorted(possible))
    sealed=[b for b in small if not set(b['nodes'])&possible]
    if sealed:
        assert len(sealed[0]['nodes'])<=4
        return common|dict(method='untouched_small_factor',whole_small_factor=sealed[0],maximum_small_rank=len(sealed[0]['nodes']))
    assert len(small)==1
    b=small[0];touch=possible&set(b['nodes']);assert len(touch)==1
    p=next(iter(touch));r=J[p][p];size=len(b['nodes'])
    assert(size,r)in((1,F(4,3)),(2,F(3,2)),(3,F(2)))
    assert all(J[p][i]>0 for i in next(e for e in end['end_sets']if p in e))
    # Each new root touching the small component has exactly one unit
    # core contact: an additional same-end contact has square>2.
    twos=[]
    for e in end['end_sets']:
        for i,j in combinations(sorted(set(e)&set(free)),2):
            lower=J[i][i]+2*J[i][j]+J[j][j];assert lower>2
            twos.append(dict(nodes=[i,j],projection_lower=lower))
    others=sorted(possible-{p});mixed=[]
    for q in others:
        assert r==J[q][q]==F(4,3)and J[p][q]==-F(4,3)
        # The residual Cauchy upper is -4/3+sqrt((2/3)^2)=-2/3.
        upper=J[p][q]+(2-r);assert upper==-F(2,3)<0
        length=min(du[p]+dv[q],du[q]+dv[p])+3;assert length<=8
        mixed.append(dict(nodes=[p,q],native_contact='1',projection_cross=J[p][q],
            maximum_root_product=upper,forced_bond_shortcut=length))
    assert not others or size==1
    # Without the short escape all new roots in the whole small factor
    # meet H only at p. The attached forest's inverse energy is<=1/r.
    maxextra=int(1/(r-1));assert size+maxextra==4
    assert r*F(maxextra+1,maxextra+2)>1
    return common|dict(method='one_contact_inverse_rank_bound',old_small_factor=b,contact=p,
        inverse_diagonal=r,unit_contact_forced=True,two_contact_bounds=twos,
        opposite_contact_shortcuts=mixed,maximum_added_rank=maxextra,maximum_small_rank=4,
        first_forbidden_rank_energy=r*F(maxextra+1,maxextra+2))

def rank15_spectators():
    rows=[];diff=[]
    for kind in('A','D'):
        lab=T.target(15,kind);g=2*sp.eye(15)
        for i,j in combinations(range(15),2):
            if lab[i][j]:assert lab[i][j]==3;g[i,j]=g[j,i]=-1
        J=g.inv();ad=[i for i in range(15)if J[i,i]<2]
        assert ad==([0,1,13,14]if kind=='A'else[0])
        assert all(x>=0 for x in J)
        pairs=[]
        for i,j in combinations(ad,2):
            v=J[i,i]+2*J[i,j]+J[j,j];assert v>=2;pairs.append([i,j,str(v)])
        profiles=[]
        for i in ad:
            for q in(1,sp.sqrt(2)):
                en=sp.simplify(q*q*J[i,i])
                if en<2:profiles.append((i,q,en))
            # The next native square above2 is(3+sqrt5)/2, already excluded.
            assert sp.simplify((3+sp.sqrt(5))*J[i,i]/2)>2
        assert len(profiles)==(6 if kind=='A'else 1)
        for i,q,en in profiles:
            assert en>=sp.Rational(15,16)
            if en==1:assert kind=='D'and i==0 and q==1
        for(i,q,_),(j,t,_)in product(profiles,profiles):
            if i==j:continue
            en=sp.simplify(q*q*J[i,i]+t*t*J[j,j]-2*q*t*J[i,j])
            assert en>=sp.Rational(15,16)and en!=1
            diff.append(dict(kind=kind,left=[i,str(q)],right=[j,str(t)],difference_energy=str(en)))
        rows.append(dict(kind=kind,inverse_diagonals=[str(J[i,i])for i in range(15)],
            admissible_positions=ad,two_contact_bounds=pairs,
            strict_profiles=[dict(node=i,magnitude=str(q),energy=str(en))for i,q,en in profiles]))
    assert 2*F(15,16)>1
    return dict(components=rows,disjoint_same_component_differences=diff,
        minimum_nonzero_component_energy='15/16',two_component_lower='15/8',
        energy_one_only='one unit contact at the long endpoint of exactly one whole D15')

def one_sided_models():
    f,_=U.affine();out=[]
    for attach in(9,10):
        g=2*sp.eye(30);g[:13,:13]=f
        lab=T.target(15,'D')
        for i,j in combinations(range(15),2):
            if lab[i][j]:g[15+i,15+j]=g[15+j,15+i]=-1
        for i,j in((13,attach),(14,12),(13,15)):g[i,j]=g[j,i]=-1
        gg=[[F(x)for x in row]for row in g.tolist()];white=sorted(set(range(30))-{9,11,13,14})
        assert Q.inertia(gg)==[28,1,1]
        assert Q.inertia(A.sub(gg,white))==[26,0,0]
        assert all(Q.inertia(A.sub(gg,white+[x]))==[27,0,0]for x in(9,11,13,14))
        models=[]
        for marks in([9,11],[11,13],[13,14],[14,9]):
            ids=sorted(set(white)|set(marks));sig=Q.inertia(A.sub(gg,ids));assert sig[1]==0
            bs=D.finite_blocks(gg,ids)if sig[2]==0 else[]
            models.append(dict(marks=marks,signature=sig,components=bs))
        quiet=models[-1];b=next(b for b in quiet['components']if set(quiet['marks'])<=set(b['nodes']))
        assert b['type']=='D13'and set(T.mark_pairs(gg,b['nodes'],quiet['marks'],13,'D'))==PAIRS
        pay=models[2];bs=pay['components'];images=[]
        if attach==9:
            assert sorted(b['type']for b in bs)==['A12','D16']
            for b in bs:
                m=next(m for m in pay['marks']if m in b['nodes']);n=len(b['nodes']);im=T.mark_pairs(gg,b['nodes'],[m],n,b['type'][0])
                assert im==([(0,),(11,)]if b['type']=='A12'else[(0,)])
                images.append(dict(type=b['type'],images=im))
        else:
            assert len(bs)==1 and bs[0]['type']=='D28'
            im=T.mark_pairs(gg,bs[0]['nodes'],pay['marks'],28,'D');assert im==[(0,12)]
            images.append(dict(type='D28',images=im))
        out.append(dict(active_attachment=attach,gram=gg,white=white,corners=models,
            quiet_corner=quiet['marks'],new_paying_corner=pay['marks'],paying_images=images,
            quiet_half='1/2',new_half='1/2',total='1',
            additional_whole_spectators='two A15/D15 components, orthogonal to this entire30-root star'))
    return out

def quiet_models(docs):
    rows=[];old=docs[U.TAG]
    for i,p in enumerate(old['nonextreme_zero_cases']):
        if p['extreme_profile']!='mate3':continue
        g=A.matrix(p['gram']);a,b=p['marks'];active=13 if p['nonextreme_side']=='left'else 14
        quiet=14 if active==13 else 13;chosen=[quiet,a]if active==13 else[b,quiet]
        marks=[b,active]if active==13 else[active,a]
        profile=p['nonextreme_profile'];proof=dict(index=i,profile=profile,quiet_corner=chosen)
        if profile=='isolated4':
            ids=sorted(set(p['white'])|set(marks));parts=A.components(g,ids)
            aa,bb=[next(z for z in parts if m in z)for m in marks]
            sm,la=(aa,bb)if len(aa)==1 else(bb,aa);assert len(sm)==1 and len(la)==12
            lm=next(m for m in marks if m in la)
            assert T.mark_pairs(g,la,[lm],12,'D')==[(10,),(11,)]
            proof|=dict(method='distinct_retained_A1_D12_shortleaf_unit',paying_corner=marks,weight='1')
        elif profile=='mate4':
            # Gauge-aware induced target matching identifies the whole immediate B13.
            ids=sorted(set(p['white'])|set(marks));parts=A.components(g,ids);assert len(parts)==1
            import fresh_order13_affine_b_singleton_transfer_20261017 as BB
            target=docs[BB.TAG]['target']
            labels=T.T.labels(g,ids);mm=T.P.embeddings(labels,target);assert mm
            images=sorted({tuple(sorted(m[ids.index(x)]for x in marks))for m in mm})
            proof|=dict(method='independent_current_affine_B13_priority',paying_corner=marks,images=images)
        else:
            assert profile=='center3';other=[13,14];ids=sorted(set(p['white'])|set(other))
            assert set(T.mark_pairs(g,ids,other,13,'D'))==PAIRS
            proof|=dict(method='two_distinct_guarded_D13_halves',paying_corner=other,weights=['1/2','1/2'],total='1')
        rows.append(proof)
    assert len(rows)==24 and Counter(r['profile']for r in rows)=={'isolated4':8,'mate4':8,'center3':8}
    return rows

def long_images(g,white,marks):
    ids=sorted(set(white)|set(marks));parts=A.components(g,ids)
    shared=next((p for p in parts if set(marks)<=set(p)),None)
    return T.mark_pairs(g,shared,marks,13,'D')if shared is not None and len(shared)==13 else[]

def current_interfaces(docs):
    literals=[];later=[];newer=[]
    for st,case,g,w,cs in T.literal_inputs(docs):
        for marks,*_ in cs:
            ims=long_images(g,w,marks);assert not set(ims)&PAIRS
            literals.append(dict(source=st,case=case,marks=marks,images=ims))
    for st,keys in B.FINITE.items():
        for key in keys:
            for path,p,g,w,bs in B.corners(docs[st][key],key):
                ims=long_images(A.matrix(g),w,p['marks']);assert not set(ims)&PAIRS
                later.append(dict(source=st,path=path,marks=p['marks'],images=ims))
    assert len(later)==564
    # Reconstruct the same complete newer physical families as the B13
    # audit, without inferring a D13 statement from that different role.
    collections=[(S.R.TAG,docs[S.R.TAG]['local_proof']['whole_proofs']),
        (S.H.TAG,docs[S.H.TAG]['local_proof']['whole_proofs']),
        (S.S.TAG,docs[S.S.TAG]['local_proof']['whole_proofs']),
        (S.N.TAG,docs[S.N.TAG]['local_proof']['proofs']),
        (S.R.R.TAG,docs[S.R.R.TAG]['local_seals']),
        (S.H.Z.U.TAG,[docs[S.H.Z.U.TAG]['local_proof']]),
        (S.TAG,docs[S.TAG]['local_proof']['old_ideal17_seals']),
        (S.TAG,docs[S.TAG]['local_proof']['A16_seals'])]
    for st,ps in collections:
        for p in ps:
            ph=p.get('old_complete_seal',p)
            for c in p['corners']:
                ims=long_images(A.matrix(ph['gram']),ph['white'],c['marks']);assert not set(ims)&PAIRS
                newer.append(dict(source=st,marks=c['marks'],images=ims))
    variants=[]
    for p in docs[S.L.C.R.S.TAG]['selected_proofs']:
        for v in p.get('variants',[]):
            ims=[long_images(A.matrix(v['gram']),v['white'],c['marks'])for c in v['corners']]
            assert not any(set(im)&PAIRS for im in ims);variants.append(dict(index=p['index'],images=ims))
    assert len(variants)==24
    d12=[]
    for st,n in((B.PART,280),(B.REST,52)):
        ps=[(path,p)for path,p in S.W.walk(docs[st])if F(p['payment'])>0];assert len(ps)==n
        for path,p in ps:
            if any(b['type']=='D13'for b in p['roles']):
                assert not set(map(tuple,p['images']))&PAIRS
                d12.append(dict(source=st,path=path,images=p['images'],payment=p['payment']))
    return dict(literal_corner_rows=literals,literal_configurations=216,later_corners=later,
        newer_corners=newer,D14_variants=variants,D12_positive_domain=332,D12_D13_outputs=d12,
        no_changed_same_component_role_in_these_domains=True,
        universal_unstored_exception='Affine-D12 quiet and arbitrary one-sided outputs, repaired separately')

def universal_roles(docs):
    low=docs['fresh_d12_low_protected_adaptive_interfaces_20261018']
    assert len(low['low'])==66 and all(r['inherited']['component_order_bound']<=12 for r in low['low'])
    assert len(low['protected'])==74 and len(low['adaptive'])==9
    for r in low['protected']+low['adaptive']:
        assert r['distance']==9 and len(r['retained_path'])==10
    # The center/leaf ideal13 theorem's finite D13 mark is the fork center,
    # not a short leaf. Reconstruct that literal universal role directly.
    f,_=U.affine();g=2*sp.eye(14);g[:13,:13]=f;g[12,13]=g[13,12]=-1
    gg=[[F(x)for x in row]for row in g.tolist()];white=sorted(set(range(13))-{0,11});marks=[0,13]
    ids=sorted(set(white)|set(marks));assert Q.inertia(A.sub(gg,ids))==[13,0,0]
    im=T.mark_pairs(gg,ids,marks,13,'D');assert im==[(0,10)]and not set(im)&PAIRS
    # In the coupled-order11 unmixed boundary the finite D_(11+r)
    # output has its mark r places from the long end, so distance9.
    order11=[]
    for r in(1,2,3):
        n=11+r;images=[(r,n-2),(r,n-1)]
        assert n-2-r==9
        if n==13:assert not set(images)&PAIRS
        order11.append(dict(spectator_rank=r,finite_rank=n,images=images,distance=9))
    return dict(low_core_bound=12,low_contexts=66,protected_paths=74,adaptive_paths=9,
        all_protected_adaptive_marked_paths=9,center_leaf_D13_role=im,
        coupled_order11_unmixed_roles=order11,
        cycle_tree_outputs='Finite8 cycles/trees and cycle12-14 give short paths or marked distance9; second-node D_m outputs with marked distance11 have m>=18, not D13.',
        newer_short_outputs='Current depth6 finite6/7/8 escape paths or whole-leaf/A-type outputs never have the D13 long-end/short-leaf distance11 role.',
        older_A16_repair='The eight affine-D11 replacements use D13(1,11)/(1,12), already bound and checked in the frozen A16/D16 packet.',
        current_independent_B13_priority='The frozen B13 endpoint repair uses A1/D12 shortleaf units; the current BC13/B10 finite roles are B-type or small cross-component roles, not a changed D13 pair.',
        inventory_scope='Written complete chosen-output retention inherits the named universal routing and classification theorems; the finite literal domains alone are not asserted to prove exhaustiveness.')

def finite(docs):
    old=docs[C.TAG];assert old['exact_actual_finite_maximum']=='378'
    ob=old['finite']['exact_obstructing_forests'];assert len(ob)==8
    rows=[]
    for rr in ob:
        fs=list(map(tuple,rr));assert sum(n for _,n in fs)==58 and S.value(fs)==378*64
        co=Counter(fs);assert co['B',13]+co['D',13]==1 and co['A',15]+co['D',15]==3
        assert value(fs)==377*64;rows.append(dict(forest=fs,old='378',new='377'))
    inc=[]
    for rank15 in combinations_with_replacement([('A',15),('D',15)],2):
        for fixed in([('D',28)],[('A',12),('D',16)]):
            fs=fixed+list(rank15);assert sum(n for _,n in fs)==58
            cost=F(value(fs),64);expected=F(741,2)if len(fixed)==1 else F(751,2)
            assert cost==expected<377
            inc.append(dict(forest=fs,previous=F(Z.value(fs),64),increment=F(context_delta(fs),64),new=cost))
    assert S.actual_cross(('A',12),('D',16))==0
    return dict(frozen_complete_count=C.TAG,all_old_costly_forests=rows,
        complete_increasing_families=inc,maximum='377',actual_attainers=rows,
        upper_argument='Outside the two exact rank58 increasing families cost is pointwise no larger than the frozen count, except the B13 monotone deletion; all eight old378 forests now equal377.',
        new_DP=False,ideal_weights_unchanged=True,ideal_profiles=2552,cusp_surplus='36',
        curvature='377',strict_global_constant='4205',
        conditional_contradiction='0 >= 4205 + 36*v_i > 0',
        global_dimension_composition='Separate obligation after full mathematical reviews of this local theorem')

def run():
    docs={st:json.loads((ROOT/f'results/{st}.json').read_text())for st in STEMS}
    for st,d in docs.items():
        assert d['status'].startswith('PASS')
        for p,h in d.get('sources',{}).items():assert digest(p)==h,(st,p)
    over=docs[A.OVER];sup,hits=support(over)
    selected=[selected_proof(c,g,ms)for c,g,ms in hits]
    print('SELECTED',Counter(r['method']for r in selected),flush=True)
    spec=rank15_spectators();one=one_sided_models();quiet=quiet_models(docs)
    interfaces=current_interfaces(docs);roles=universal_roles(docs);fin=finite(docs)
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{st}.{ext}'for st in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return dict(status='PASS_AMBIENT58_D13_THREE_RANK15_CONTEXT_LOCAL_TRANSFER_AND_FINITE_CLOSURE',
        sources={p:digest(p)for p in sorted(paths)},ambient_dimension=58,
        changed_pairs=sorted(PAIRS),selected_support=sup,selected_proofs=selected,
        rank15_contact_classification=spec,one_sided_whole_models=one,quiet_repairs=quiet,
        current_interfaces=interfaces,universal_output_inventory=roles,finite=fin,
        scope='Dimension-specific local rule and exact finite-count composition; no new global dimension theorem is claimed by this packet alone.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=ROOT/f'results/{TAG}.json';assert args.verify or not out.exists()
    d=json.loads(json.dumps(enc(run())))
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['finite']['maximum'],flush=True)
