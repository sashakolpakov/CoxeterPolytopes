#!/usr/bin/env python3
"""All-label obstruction for paired top/bottom endcaps retaining current47."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import fresh_nonarithmetic_replace47_close12_20261017 as C
P,U,I,D=C.P,C.U,C.I,C.D
ROOT=C.ROOT;TAG='fresh_nonarithmetic_replace47_b_pair_gate_20261017'


def gap(q):
    if q.l>0:return dict(kind='obtuse',product=q.enc())
    if q.h<-2:return dict(kind='dotted_no_finite_pair',product=q.enc())
    for m in range(2,21):
        a,b=-D.theta(m),-D.theta(m+1)
        if q.h<a.l and b.h<q.l:
            return dict(kind='strict_consecutive_angle_gap',product=q.enc(),labels=[m,m+1],thresholds=[a.enc(),b.enc()])
    return None


def run():
    z=C.geometry();ns=z['normals'];pair=z['pair'];h=z['common_center']
    e7=[U.rnd(a+(b+2*c)/3)for a,b,c in zip(ns[7],ns[5],ns[6])]
    c=C.literal(z,7,52);d=C.literal(z,30,47)
    frame_data=[dict(ids=[7,30,47],roots=[e7,ns[30],ns[47]],
        gram=[[I(F(4,3)),-D.theta(5),I(-1)],[-D.theta(5),I(2),d],[I(-1),d,I(2)]],
        seed=ns[52],seed_square=I(2),seed_products=[c,-D.theta(4),I()],
        profiles=[(2,2,m)for m in [2,3,4]]),
        dict(ids=[30,47,52],roots=[ns[30],ns[47],ns[52]],
        gram=[[I(2),d,-D.theta(4)],[d,I(2),I()],[-D.theta(4),I(),I(2)]],
        seed=e7,seed_square=I(F(4,3)),seed_products=[-D.theta(5),I(-1),c],
        profiles=[(m30,m47,2)for m30,m47 in product([2,3],[2,3,4,5])]+[(2,m47,3)for m47 in [2,3,4,5]])]
    records=[];groups=[];normals=[]
    for group,fr in enumerate(frame_data):
        inv,ip=P.C53.inverse(fr['gram']);seed_coeff=U.imv(inv,fr['seed_products'])
        k=[U.rnd(a-b)for a,b in zip(fr['seed'],P.combine(seed_coeff,fr['roots']))]
        kq=U.rnd(fr['seed_square']-U.idot(fr['seed_products'],seed_coeff));assert kq.l>0
        accepted=[]
        for ms in fr['profiles']:
            contacts=[-D.theta(m)for m in ms];coeff=U.imv(inv,contacts);p=P.combine(coeff,fr['roots'])
            rem=U.rnd(2-U.idot(contacts,coeff));assert rem.l>0
            scale=U.sqrt(rem/kq)
            for sign in [-1,1]:
                n=[U.rnd(a+sign*scale*b)for a,b in zip(p,k)]
                prescribed={**{i:2 for i in C.S},**dict(zip(fr['ids'],ms))}
                raw={i:U.rnd(pair(n,v))for i,v in ns.items()if i not in prescribed}
                rejection=None
                for i,q in raw.items():
                    g=gap(q)
                    if g is not None and g['kind']!='dotted_no_finite_pair':
                        rejection=dict(wall=i,**g);break
                row=dict(index=len(records),group=group,frame=fr['ids'],labels=list(ms),branch=sign,
                    frame_inverse_pivots=ip,kernel_square=kq.enc(),remaining_norm=rem.enc(),positive_scale=scale.enc(),normal=P.enc(n))
                if rejection is not None:row['rejection']=rejection
                else:
                    assert all(q.h<-2 for q in raw.values())
                    assert pair(n,h).h<0
                    accepted.append(len(records));normals.append((len(records),n))
                    row['accepted_native_contacts']=[[i,m]for i,m in sorted(prescribed.items())]
                    row['strict_dotted_products']={str(i):q.enc()for i,q in raw.items()}
                records.append(row)
        groups.append(accepted)
    assert len(records)==30 and len(groups[0])==3 and len(groups[1])==11
    nd=dict(normals);pairs=[]
    for a,b in product(*groups):
        q=U.rnd(pair(nd[a],nd[b]));g=gap(q)
        assert g is not None,(a,b,q.enc())
        pairs.append(dict(indices=[a,b],rejection=g))
    assert len(pairs)==33
    left0=next(i for i in groups[0]if records[i]['labels']==[2,2,2]);a=nd[left0]
    # Finite rank3 diagrams containing an edge m>=6 have that I2(m)
    # as an isolated component. Thus high right47 labels require left0
    # and mutual rightness. The four products then determine B uniquely.
    ad=U.rnd(pair(a,ns[52]))
    tail_gram=[[I(2),d,-D.theta(4),I()],
               [d,I(2),I(),I()],[-D.theta(4),I(),I(2),ad],[I(),I(),ad,I(2)]]
    ti,tp=P.C53.inverse(tail_gram);tails=[]
    for m30 in [2,3]:
        v=D.theta(m30)
        quadratic=ti[1][1];linear=U.rnd(2*v*ti[0][1]);constant=U.rnd(v*v*ti[0][0]-2)
        assert quadratic.l>0 and linear.l>=0
        c6=D.theta(6);value=U.rnd(quadratic*c6*c6+linear*c6+constant)
        derivative=U.rnd(2*quadratic*c6+linear)
        assert value.l>0 and derivative.l>0
        tails.append(dict(contact30=m30,polynomial_ascending=P.enc([constant,linear,quadratic]),
            value_at_sqrt3=value.enc(),derivative_at_sqrt3=derivative.enc(),coefficient_square_positive=True,linear_nonnegative=True))
    return dict(status='PASS_ALL_LABEL_PAIRED_B_ENDPOINT_CAP_OBSTRUCTION',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in[Path(__file__),Path(C.__file__),ROOT/'results'/(C.TAG+'.json'),ROOT/'results/fresh_nonarithmetic_replace47_close12_boundary_20261017.json']},
        support=C.S,records=records,accepted_groups=groups,low_pairs=pairs,left_right_branch_for_tail=left0,
        tail_frame=[30,47,52,'left_all_right_cap'],tail_gram=[P.enc(r)for r in tail_gram],tail_inverse_pivots=tp,high_label_tails=tails,
        profile_completeness='The left cap terminates7/30 and7/47. H4² at7/30 forces rightS6,7,30; the remaining47 label is2,3,4 by positivity with projected7²=4/3. For the right cap here assumed rightS6, endpoints30/52(B2) and47/52(A1²) permit (m30,m47,m52)=(2or3,m,2) for arbitrarym>=2, or(2,m,3) withm2..5. Lowm2..5 give12 profiles. For highm>=6, finite corner47/left/right requires left47=2 and mutualright by finite rank3 classification. The norm2 equation then has no root for c=2cos(pi/m)>=sqrt3: its polynomial is positive there and strictly increasing.',
        scope='Excludes the specified two-cap construction in which one rightS6 cap handles both7-side open edges, one rightS6 cap handles both52-side open edges, and their mutual finite contact supplies the two connecting endpoints on30 and47. No common-center exclusion is used for rejected branches. It does not exclude a right-cap nonzeroA2 profile, other duty assignments, three or more caps, or wall replacements.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['records']),len(out['low_pairs']),len(out['high_label_tails']),flush=True)

if __name__=='__main__':main()
