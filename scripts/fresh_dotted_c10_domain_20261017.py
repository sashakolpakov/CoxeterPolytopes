#!/usr/bin/env python3
"""Continuous dotted-C10 endpoint domain; exact native tail certificate."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_release49_62_replacement_20261017 as C
import fresh_old30_connector_native_gate_20261017 as E

G,P,U,I,Q,ROOT=E.G,E.P,E.U,E.I,E.Q,E.ROOT
TAG='fresh_dotted_c10_domain_20261017'
S=[0,1,2,3,5,6]
DELETE=[12,47,52,64,65]
FRAME=[7,30,10]

def data():
    z=C.geometry();ns=z['normals'];ids=sorted(ns);h=z['common_center']
    wg=[[C.literal(z,i,j)for j in S]for i in S]
    wi,wp=P.C53.inverse(wg);P.pd(wg)
    cov={i:[C.literal(z,i,j)for j in S]for i in ids}
    wc={i:U.imv(wi,cov[i])for i in ids}
    def pair(i,j):
        return U.rnd(C.literal(z,i,j)-U.idot(cov[i],wc[j]))
    fg=[[pair(i,j)for j in FRAME]for i in FRAME]
    fi,fp=P.C53.inverse(fg);rhs=[I(),I(),I(-1)]
    vc=U.imv(fi,rhs);alpha=U.rnd(-U.idot(rhs,vc));assert alpha.l>0
    choices=[]
    for seed in ids:
        if seed in S+FRAME:continue
        sc=[pair(seed,j)for j in FRAME];kc=U.imv(fi,sc)
        kq=U.rnd(pair(seed,seed)-U.idot(sc,kc))
        if kq.l>0:choices.append((float((kq.h-kq.l)/kq.l),seed,kc,kq))
    assert choices
    _,seed,kc,kq=min(choices);scale=Q.sqrt_closed(2/kq)
    # The exact sign is fixed below by the old55 product, not by
    # a numerical comparison. The opposite sign is separately tested.
    v={};k={}
    for i in ids:
        col=[pair(j,i)for j in FRAME]
        v[i]=U.rnd(U.idot(vc,col))
        k[i]=U.rnd(scale*(pair(seed,i)-U.idot(kc,col)))
    for i,t in zip(FRAME,rhs):v[i]=t;k[i]=I()
    for i in S:v[i]=k[i]=I()
    if k[55].l>0:scale=-scale;k={i:-q for i,q in k.items()}
    assert k[55].h<0
    hc={i:U.rnd(z['pair'](h,ns[i]))for i in ids}
    swc=U.imv(wi,[hc[i]for i in S])
    hp={i:U.rnd(hc[i]-U.idot(cov[i],swc))for i in ids}
    vh=U.rnd(U.idot(vc,[hp[i]for i in FRAME]))
    kh=U.rnd(scale*(hp[seed]-U.idot(kc,[hp[i]for i in FRAME])))
    projected={i:[U.rnd(a-b)for a,b in zip(ns[i],P.combine(wc[i],[ns[j]for j in S]))]for i in FRAME+[seed]}
    vp=P.combine(vc,[projected[i]for i in FRAME])
    kp=[U.rnd(scale*(a-b))for a,b in zip(projected[seed],P.combine(kc,[projected[i]for i in FRAME]))]
    return dict(z=z,ids=ids,frame_gram=fg,frame_pivots=fp,white_pivots=wp,
                alpha=alpha,seed=seed,kernel_square=kq,kernel_scale=scale,
                V=v,K=k,Vh=vh,Kh=kh,Vphysical=vp,Kphysical=kp,pair=pair)

def run():
    d=data();alpha=d['alpha'];v,k=d['V'],d['K']
    limit=Q.sqrt_closed(alpha/2);start=Q.sqrt_closed(1+2*alpha)
    rows=[];deleted=[];full_tail=[]
    start3=Q.sqrt_closed(1+F(9,2)*alpha)
    for i in d['ids']:
        if i in S+[7,30]:
            assert v[i].l==v[i].h==k[i].l==k[i].h==0
            continue
        a,b=v[i],k[i]
        initial=U.rnd(2*a+b*start)
        derivative=U.rnd(a+I(max(F(0),b.h))*limit)
        row=dict(wall=i,V_product=a.enc(),K_product=b.enc(),value_at_u2=initial.enc(),
                 strict_derivative_upper=derivative.enc())
        assert derivative.h<0,(i,'derivative',derivative.enc())
        value3=U.rnd(3*a+b*start3);assert value3.h<-2
        full_tail.append(dict(wall=i,value_at_u3=value3.enc()))
        if i in DELETE:deleted.append(row);continue
        assert derivative.h<0,(i,'derivative',derivative.enc())
        if i==10:assert initial.l==initial.h==-2
        else:assert initial.h<-2,(i,'initial',initial.enc())
        rows.append(row)
    ch=U.rnd(2*d['Vh']+start*d['Kh'])
    dh=U.rnd(d['Vh']+I(max(F(0),d['Kh'].h))*limit)
    assert ch.h<0 and dh.h<0
    # On the other sheet old55 has product a*u+b*sqrt(1+alpha*u^2/2),
    # with a<0,b>0. A positive asymptotic slope proves it is >0 for
    # every u>0, using sqrt(1+alpha*u^2/2)>u*sqrt(alpha/2).
    oa,ob=v[55],-k[55];slope=U.rnd(oa+ob*limit)
    assert ob.l>0 and slope.l>0
    control=json.loads((ROOT/'results'/f'{C.TAG}.json').read_text())
    escape=G.decode(control['strict_escape']['point'])
    ev=U.rnd(d['z']['pair'](escape,d['Vphysical']))
    ek=U.rnd(d['z']['pair'](escape,d['Kphysical']))
    escape3=U.rnd(3*ev+start3*ek);assert escape3.l>0
    paths=[Path(__file__),Path(C.__file__),Path(E.__file__),Path(Q.__file__),
           Path(G.__file__),Path(C.OLD.__file__),C.PARENT,C.GATE,G.CONTROL,
           ROOT/'results'/f'{C.TAG}.json',ROOT/'results/fresh_nonarithmetic_release49_62_boundary_20261017.json']
    return dict(status='PASS_EXACT_COMPLETE_DOTTED_C10_NATIVE_CURVE',
      sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      retained_control_walls=len(d['ids']),deleted_walls=DELETE,remaining_old_walls=len(d['ids'])-len(DELETE),
      parameter=dict(u='-<C,old10>',domain='u>2',formula='C(u)=u*V+sqrt(1+alpha*u^2/2)*K',
        V_square='-alpha',K_square=2,V_K_product=0,frame=FRAME,
        frame_gram=G.encmat(d['frame_gram']),frame_inverse_pivots=d['frame_pivots'],
        kernel_seed=d['seed'],kernel_square=d['kernel_square'].enc(),kernel_scale=d['kernel_scale'].enc(),
        alpha=alpha.enc(),right_walls=S+[7,30],native_V_contacts={str(i):v[i].enc()for i in d['ids']},
        native_K_contacts={str(i):k[i].enc()for i in d['ids']}),
      physical_recipe=dict(computational_basis=d['z']['basis_ids'],V=P.enc(d['Vphysical']),K=P.enc(d['Kphysical'])),
      retained_dotted_gates=rows,deleted_wall_diagnostic=deleted,
      full_control_tail=dict(parameter_domain='u>=3',deleted_walls=[],dotted_gates=full_tail,
        consequence='For everyu>=3 the same C can be added coherently to all63 old walls, preserving every old finite corner and compact face by acute-center projection. This is a sufficient whole-control tail, not a classification of its possible extra finite contacts for2<u<3.'),
      common_center=dict(V_product=d['Vh'].enc(),K_product=d['Kh'].enc(),
        value_at_u2=ch.enc(),strict_derivative_upper=dh.enc()),
      opposite_sheet_rejection=dict(wall=55,V_product=oa.enc(),K_product=ob.enc(),positive_asymptotic_slope=slope.enc()),
      saved_control_escape=dict(V_product=ev.enc(),K_product=ek.enc(),C3_product=escape3.enc(),
        statement='C(3) strictly cuts the exact saved control escape. This single cut does not establish finite volume or improved complete boundary.'),
      proof='Project off the positive S6 Gram. The invertible three-frame has a positive orthogonal complement, so it is Lorentzian. Its prescribed contacts determine uV. The norm equation gives exactly the two signs of sqrt(1+alpha*u^2/2)K. The chosen sign has negative K.55. For every retained nonright wall, f(u)=a*u+b*sqrt(1+alpha*u^2/2) has fprime bounded above by a+max(0,b)*sqrt(alpha/2)<0, and f(2)<-2 except for old10, where f(2)=-2 exactly. Thus all such contacts are strictly dotted for u>2. The fixed common center satisfies the same negative-start and derivative tests. The other sheet is strictly obtuse to55 for every u>0, since its positive square-root coefficient and positive asymptotic slope bound its product below by a positive multiple ofu. Acute-center projection proves all retained old facets and the new norm2 facet are actual. No new ideal point or finite-volume conclusion follows.',
      scope='Complete continuous native admissibility for norm2 C rightS6+7+30 with C10<-2, retaining precisely the current63 control minus12,47,52,64,65. The chosen branch is admissible for every realu>2; the other branch is impossible. The sufficient subdomainu>=3 retains all63 old walls and their finite/compact faces. This is a one-wall construction domain, not a closed patch, ambient finite-vertex census, or finite-volume chamber.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],len(out['retained_dotted_gates']),flush=True)
