#!/usr/bin/env python3
"""NEW SOURCE-ONLY preparation: literal-table/aggregate binding, no native imports.

Root must review and launch ONE named batch through the existing bounded
research runner --parent-lock. This does not regenerate certificates or call
an optimizing recurrence. Its receipt proves only its explicitly named bridge.
"""
import argparse
from collections import Counter
import fcntl
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import importlib.util
from itertools import combinations, combinations_with_replacement, product
import json
import os
from pathlib import Path
import signal

ROOT = Path(__file__).resolve().parents[1]
LIB = 'scripts/finite57_literal_weights_audit_20261024.py'
PROJECTOR = 'scripts/check_d13_transfer_batches_20261020.py'
ARRAY = 'scripts/check_finite57_b13_frontier_batches_20261021.py'
COARSE = 'results/fresh_a16_d16_witness_context_ambient58_frontier_20261018.json'
FRONT = 'results/fresh_d17_a16_coupled_ambient58_frontier_20261018.json'
CONDITIONAL = 'results/fresh_d17_a9_a2_context_zero_ambient58_frontier_20261018.json'
FINAL = 'results/fresh_d13_three_rank15_context_transfer_20261018.json'
IDEAL = 'results/fresh_order10_low8_ambient58_frontier_20261018.json'
PINS = {
    LIB:'f6eba7af77af50b821e9e99c589388c1352f9ba5a99052aaeee64393f5ada1e0',
    PROJECTOR:'be8a42d1d3687204af4d559c6bd8a28779dbfac68aeccdd95d41a7172645d555',
    ARRAY:'9052f571ed8642e7fa15d02cf6727e10362e8b1d51ebd3d4a12f9953e551eb15',
    COARSE:'9e5838e5eb9afec9c5b050d4d753e2ec28dc962c788af0fc7e7d5cc144ccc70e',
    FRONT:'0fb359b3b07d6fbabb31adbf611d769f41ebe855075076301975e086055b1223',
    CONDITIONAL:'f80edecac4af2021e39806a194befebcba0ab83c748138a6e40556a7dedc117c',
    FINAL:'2c7d615cd38123bef20cffad68de8ab30baa91826c147a6bb084892590e8ca86',
    IDEAL:'7a01260324b2dff2c81bece2f953f9602e1c38537a389f2da2caba439d5e7a43',
    'manuscripts/finite_volume57.tex':'c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0',
}
BATCHES = ('domination','residual-zero','context-rows','final-forests','affine-small',
    'affine-large-0','affine-large-1','affine-large-2','affine-large-3',
    'coarse-index') + tuple(
    f'{name}-{i}' for name in ('coarse','fibers','conditional') for i in range(4))
W = PROJ = ARR = None
READ = set()


def digest(relative):
    p = (ROOT/relative).resolve()
    assert p.is_relative_to(ROOT) and p.is_file()
    h = sha256()
    with p.open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2),b''): h.update(block)
    return h.hexdigest()


def module(relative,name):
    spec = importlib.util.spec_from_file_location(name,ROOT/relative)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def document(relative,*paths):
    READ.add(relative)
    return PROJ.JsonProjection(relative).load(paths)


class View:
    def __init__(self,relative):
        READ.add(relative)
        self.view = ARR.ArrayView(ROOT/relative)
    def __enter__(self): return self.view
    def __exit__(self,*unused): self.view.close()


def typ(t):
    k,n = t
    return ('I' if k == 'I2*' else k,int(n))


def types():
    return tuple(sorted([('A',r) for r in range(1,59)]+[('B',r) for r in range(2,59)]
        +[('D',r) for r in range(4,59)]+[('E',6),('E',7),('E',8),('F',4),('H',3),('H',4),('I',2)]))


def representative(t):
    if t == ('A',17): return t
    if t in (('A',12),('B',12)): return ('A',12)
    if t == ('D',12): return t
    if t[1] == 6: return ('E',6)
    if t[1] == 9: return ('B',9)
    return ('A' if t[1]<=8 else 'D',t[1])


@lru_cache(None)
def base(t):
    return 64*W.same_component_sum(t,preliminary=True)


@lru_cache(None)
def cross(a,b):
    # Exhaustive host support of the printed ordinary cross table.
    hosts = {('A',r) for r in range(9,18)} | {('B',r) for r in range(9,13)} | {('D',r) for r in range(9,17)}
    if min(a[1],b[1]) > 8 or not ({a,b} & hosts): return F(0)
    return 64*W.ordinary_cross_sum(a,b)


def ordinary(fs):
    c = Counter(fs)
    return sum(base(t)*m for t,m in c.items()) + sum(
        cross(a,b)*ma*mb for (a,ma),(b,mb) in combinations(c.items(),2)) + sum(
        cross(t,t)*m*(m-1)/2 for t,m in c.items())


def value(fs,final=False):
    # The final ordinary B13 restoration is removed separately.
    out = ordinary(fs) + 64*W.context_sum(fs,final=final)
    if final: out -= 64*Counter(fs)[('B',13)]
    return out


BAR_BASE = {('A',5):F(10),('B',9):F(36),('D',10):F(45),('D',11):F(55),
    ('A',12):F(66),('D',12):F(65),('D',13):F(78),('D',14):F(179,2),
    ('D',15):F(201,2),('D',16):F(219,2),('A',17):F(3613,32),('D',17):F(919,8)}
HOSTS = (('A',5),('B',9),('D',10),('D',11),('A',12),('D',12),('D',13),('D',14),('D',15),('D',16),('A',17))
BAR_ROWS = {
 ('A',1):('3/8','1/2','1','5/2','5/2','5/2','137/64','33/16','25/16','3/2','1'),
 ('A',2):('0','1','4','5','0','3','151/32','4','3','0','2'),
 ('A',3):('0','1','4','9/2','0','9/2','9/2','9/2','3','0','3'),
 ('A',4):('0','1','4','8','17/2','6','8','6','4','0','4'),
 ('A',5):('0','5/8','5','19/2','27/4','21/2','15/2','15/2','5','0','5'),
 ('E',6):('0','5/8','0','4','0','0','0','0','0','0','0'),
 ('A',7):('0','5/8','0','0','0','0','0','0','0','0','0'),
 ('A',8):('0','5/8','0','0','0','0','0','0','0','0','0')}


def bar_base(t):
    if t in BAR_BASE: return 64*BAR_BASE[t]
    k,r = t
    if r <= 8: return 32*r*(r-1)
    assert k == 'D' and r>=18
    return 64*(8*r-32-int(r>=20)-max(0,r-20)-max(0,r-21)-max(0,r-22))


def bar_cross(a,b):
    out = []
    for sm,la in ((a,b),(b,a)):
        if sm in BAR_ROWS and la in HOSTS: out.append(64*F(BAR_ROWS[sm][HOSTS.index(la)]))
    assert not out or all(v==out[0] for v in out)
    return out[0] if out else F(0)


def bar_fixed(fs):
    return sum(bar_base(t) for t in fs)+sum(bar_cross(a,b) for a,b in combinations(fs,2))


ITEMS = (1,2,3,4,6,7,8)+tuple(range(17,59))
def item(r,fs):
    t = ('E',6) if r==6 else ('A' if r<=8 else 'D',r)
    return bar_base(t)+sum(bar_cross(t,u) for u in fs)


def fixed(cs,a12,k5):
    out = [('A',5)]*k5
    for r,count in zip(range(9,18),cs,strict=True):
        if r==9: out += [('B',9)]*count
        elif r==12: out += [('A',12)]*a12+[('D',12)]*(count-a12)
        elif r==17: out += [('A',17)]*count
        else: out += [('D',r)]*count
    return out


def domain():
    out = set()
    def visit(r,left,cs):
        if r==18:
            for a12 in range(cs[3]+1):
                for k5 in range(left//5+1): out.add((tuple(cs),k5,a12))
            return
        for count in range(left//r+1): visit(r+1,left-r*count,cs+[count])
    visit(9,58,[])
    assert len(out)==3298
    return out


def domination():
    g = {('A',13):64,('D',14):64,('D',16):64,('A',14):32,('A',15):32,('D',15):32,('A',16):80}
    assert W.finite_cross_weight(('A',1),0,('D',10),7)==0
    assert W.ordinary_cross_sum(('A',1),('D',10))==1
    rows = document(FRONT,('representative_cover','actual_type_base_domination'))['representative_cover']['actual_type_base_domination']
    assert {typ(r['type']) for r in rows}==set(types()) and len(rows)==len(types())
    saved = {typ(r['type']):r for r in rows}
    for t in types():
        literal = base(t)+g.get(t,0)
        upper = bar_base(representative(t))+(22 if t==('A',17) else 0)
        assert literal==saved[t]['maximum_guarded_base']
        assert upper==saved[t]['representative_upper'] and literal<=upper
        for u in types():
            h = 24 if {t,u}=={('A',1),('D',5)} else 32 if {t,u}=={('A',4),('A',12)} else 0
            assert cross(t,u)+h<=bar_cross(representative(t),representative(u))
    # Infinite I2(m),m>=5 is one symbolic payment type, exactly equal to B2.
    assert base(('I',2))==base(('B',2))
    assert all(cross(('I',2),t)==cross(('B',2),t) for t in types())
    return dict(actual_types=len(types()),ordered_pairs=len(types())**2,equality='saved guarded bases',inequality='complete barred cross domination',D10_center_regression='A1 zero')


def residual_zero():
    residual = [t for t in types() if representative(t) in {('E',6)}|{('A',r) for r in (1,2,3,4,7,8)}|{('D',r) for r in range(17,59)} and t not in (('A',1),('A',2))]
    assert not any(t[1]==5 or 9<=t[1]<=16 or t==('A',17) for t in residual)
    for t in residual:
        assert cross(t,('A',1))==cross(t,('A',2))==0
        for u in residual: assert cross(t,u)==0
    return dict(residual_types=residual,separate_fixed_A1_A2_zero=True,ordinary_residual_zero=True,context_exceptions=['A4/A12 with A1','A1/A2/A4/A10 family'])


def context_rows():
    checks = 0
    # Entire marked-node product for each stated cross guard stratum.
    roles = [(("A",4),("A",12),[('A',1)],F(1,2)),
        (("A",2),("A",4),[('A',1),('A',10)],F(1)),
        (("A",2),("D",5),[('A',1),('A',9)],F(1)),
        (("A",1),("D",5),[('A',11)],F(3,8)),
        (("A",1),("D",5),[('A',9),('A',2)],F(1)),
        (("A",1),("D",5),[('A',11),('A',9),('A',2)],F(1)),
        (("A",2),("A",2),[('A',1),('A',1),('A',10)],F(9,8))]
    for a,b,witness,expected in roles:
        # Check every submultiset, with duplicate witness counts retained.
        seen = set()
        for mask in range(1<<len(witness)):
            ws = tuple(sorted(witness[i] for i in range(len(witness)) if mask>>i&1))
            if ws in seen: continue
            seen.add(ws)
            fs = [a,b]+list(ws)
            delta = sum(W.finite_corner_weight(fs,(0,i),(1,j))-W.finite_cross_weight(a,i,b,j) for i in range(a[1]) for j in range(b[1]))
            c=Counter(ws)
            if (a,b)==(('A',4),('A',12)): wanted=F(1,2)*bool(c['A',1])
            elif (a,b)==(('A',2),('A',4)): wanted=F(bool(c['A',1] and c['A',10]))
            elif (a,b)==(('A',2),('D',5)): wanted=F(bool(c['A',1] and c['A',9]))
            elif (a,b)==(('A',1),('D',5)):
                joint=bool(c['A',9] and c['A',2])
                wanted=F(3,8)*bool(c['A',11] or joint)+F(5,8)*joint
            else:
                assert a==b==('A',2)
                wanted=F(9,8)*bool(c['A',1]>=2 and c['A',10])
            assert delta==wanted,(a,b,ws,delta,wanted)
            if ws==tuple(sorted(witness)): assert delta==expected
            checks += 1
    for t,increment in ((('A',13),1),(('A',14),F(1,2)),(('A',15),F(1,2)),(('D',15),F(1,2)),(('A',16),F(5,4)),(('D',16),1)):
        for count in (0,1,2):
            fs=[t]+[('A',1)]*count
            delta=sum(W.finite_corner_weight(fs,(0,a),(0,b))-W.finite_same_weight(t,a,b) for a,b in combinations(range(t[1]),2))
            assert delta==(increment if count else 0)
            checks+=1
    for mask in range(8):
        witness=[('A',i+1) for i in range(3) if mask>>i&1]
        fs=[('D',14)]+witness
        delta=sum(W.finite_corner_weight(fs,(0,a),(0,b))-W.finite_same_weight(('D',14),a,b) for a,b in combinations(range(14),2))
        assert delta==bool(witness)
        checks+=1
    return dict(guard_stratum_products=checks,scope='literal role/Boolean guard checks; universal polynomial counting uses the separately reviewed source argument')


def final_forests():
    ff=document(FINAL,('finite',))['finite']
    old={tuple(sorted(((k,13),)+p)) for k in ('B','D') for p in combinations_with_replacement((('A',15),('D',15)),3)}
    got={tuple(sorted(map(typ,r['forest']))) for r in ff['all_old_costly_forests']}
    assert got==old and len(got)==8
    for fs in old:
        assert value(fs)==378*64 and value(fs,True)==377*64
    positive={tuple(sorted(core+p)) for core in ((('D',28),),(('A',12),('D',16))) for p in combinations_with_replacement((('A',15),('D',15)),2)}
    rows=ff['complete_increasing_families']
    assert {tuple(sorted(map(typ,r['forest']))) for r in rows}==positive and len(rows)==6
    for r in rows:
        fs=list(map(typ,r['forest']))
        assert value(fs,True)==64*F(r['new'])<=64*377
        assert value(fs,True)-value(fs)==64*F(r['increment'])
    return dict(negative_forests=8,positive_forests=6,rank_guard_completeness='58=13+3*15=28+2*15=12+16+2*15')


def affine(batch):
    env={int(r):F(v) for r,v in document(IDEAL,('ideal','current_envelopes'))['ideal']['current_envelopes'].items()}
    orders=range(2,18) if batch=='affine-small' else [m for m in range(18,59) if (m-18)%4==int(batch[-1])]
    rows=[]
    for m in orders:
        kinds=['A']+(['C'] if m>=3 else [])+(['B'] if m>=4 else [])+(['D'] if m>=5 else [])
        kinds += {3:['G2'],5:['F4'],7:['E6'],8:['E7'],9:['E8']}.get(m,[])
        costs={k:W.affine_sum(k,m) for k in kinds}
        assert max(costs.values())/m==env[m-1],(m,costs,env[m-1])
        rows.append(dict(order=m,costs={k:str(v) for k,v in costs.items()},envelope=str(env[m-1])))
    return dict(rows=rows,exceptional='all order<=9 factors have binomial(order,2) cost; classification external',scope='literal costs exactly bind the frozen cusp envelopes; historical per-pair receipt provenance is separate')


def coarse_index():
    with View(COARSE) as view:
        rows=view.array(('finite','expanded_alternatives'))
        got={(tuple(rows[i]['counts']),rows[i]['fixed_A5_count'],rows[i]['fixed_A12_count']) for i in range(len(rows))}
        assert got==domain() and len(rows)==len(got)==3298
        costly={(tuple(rows[i]['counts']),rows[i]['fixed_A5_count'],rows[i]['fixed_A12_count']) for i in range(len(rows))
            if rows[i]['sixtyfourths_upper']+22*rows[i]['counts'][-1]>377*64}
    with View(FRONT) as view:
        classes=view.array(('finite','classes'))
        selected={(tuple(classes[i]['counts']),classes[i]['k5'],classes[i]['a12']) for i in range(len(classes))}
        assert costly==selected and len(classes)==len(costly)==39
    return dict(complete_fixed_choices=3298,complete_costly_classes=39,
        threshold='saved coarse upper +22 per fixed A17 >377*64',
        domain='rank58 counts9..17, every A12 split, every rank5 count')


def coarse(batch):
    shard=int(batch[-1]);checked=inequalities=0;costly=set()
    tabs=document(COARSE,('finite','tables'))['finite']['tables']
    index={(tuple(t['key']),t['fixed_A5_count'],t['exclude2']):t for t in tabs}
    assert len(index)==len(tabs)==6596
    with View(COARSE) as view:
        rows=view.array(('finite','expanded_alternatives'))
        assert len(rows)==3298
        for ri in range(len(rows)):
            if ri%4!=shard: continue
            r=rows[ri];cs=r['counts'];k5=r['fixed_A5_count'];a12=r['fixed_A12_count'];fs=fixed(cs,a12,k5)
            rem=58-sum(t[1] for t in fs);key=tuple(r['key'])
            ws={rank:item(rank,fs) for rank in ITEMS}
            tables=[index[key,k5,flag]['values'] for flag in (False,True)]
            for flag,vs in zip((False,True),tables):
                assert len(vs)==59 and vs[0]==0
                for n in range(1,59):
                    for rank in ITEMS:
                        if rank<=n and not(flag and rank==2):
                            assert vs[n]>=vs[n-rank]+ws[rank];inequalities+=1
            candidates=[tables[0][rem]]
            for reserve in (1,2):
                g0=reserve==2 and bool(cs[1]);g1=bool(cs[0] and k5)
                if rem<reserve or not(g0 or g1) or (reserve==2 and not g0): continue
                for a2 in range((rem-reserve)//2+1):
                    candidates.append(reserve*ws[1]+a2*ws[2]+tables[1][rem-reserve-2*a2]+36*a2*(a2-1)*g0+64*a2*k5*g1)
            bound=bar_fixed(fs)+max(candidates)+22*cs[-1]
            assert bound<=r['sixtyfourths_upper']+22*cs[-1]
            if bound>377*64: costly.add((tuple(cs),k5,a12))
            checked+=1
    return dict(shard=shard,checked_choices=checked,current_bar_transition_inequalities=inequalities,costly_classes_in_shard=len(costly),scope='upper relaxation binding, no maximizing recurrence')


def fibers(batch):
    shard=int(batch[-1]);classes=rows_count=weights_count=0
    residual_set={t for t in types() if t not in (('A',1),('A',2)) and (representative(t)[1] in ITEMS) and representative(t)!=('A',17)}
    all_fibers={r:tuple(t for t in types() if representative(t)==r) for r in set(map(representative,types()))}
    with View(FRONT) as view:
        cls=view.array(('finite','classes'));assert len(cls)==39
        for ci in range(len(cls)):
            if ci%4!=shard: continue
            c=cls[ci];cs=c['counts'];k5=c['k5'];a12=c['a12'];reps=Counter(fixed(cs,a12,k5))
            expected={tuple(sorted(sum((list(x) for x in groups),[]))) for groups in product(*(tuple(combinations_with_replacement(all_fibers[t],n)) for t,n in reps.items()))}
            got={tuple(sorted(map(typ,r['fixed_actual']))) for r in c['rows']}
            assert got==expected and len(got)==len(c['rows'])
            for r in c['rows']:
                fs=list(map(typ,r['fixed_actual']));assert sum(t[1] for t in fs)+c['residual_rank']==58
                seen=set()
                for tab in r['tables']:
                    has1,no4,no3=tab['key'];assert all(type(x) is bool for x in tab['key'])
                    key=tuple(tab['key']);assert key not in seen;seen.add(key)
                    ts=list(map(typ,tab['types']))
                    allowed={t for t in residual_set if not(no4 and t==('A',4)) and not(no3 and t==('A',3))}
                    assert set(ts)==allowed and len(ts)==len(allowed)
                    for t,w in zip(ts,tab['weights'],strict=True):
                        literal=base(t)+sum(cross(t,b) for b in fs)+(32*Counter(fs)[('A',12)] if has1 and t==('A',4) else 0)
                        assert literal==w,(ci,t,literal,w)
                        assert cross(t,('A',1))==cross(t,('A',2))==0
                        weights_count+=1
                attaining=list(map(typ,r['attaining_forest']))
                assert value(attaining)==64*F(r['maximum'])
                rows_count+=1
            classes+=1
    return dict(shard=shard,classes=classes,actual_fixed_fibers=rows_count,exact_saved_item_weights=weights_count,scope='literal equality of actual residual coefficients and full attaining costs; frozen Bellman/flag-coverage proofs remain premises')


def conditional(batch):
    shard=int(batch[-1]);checked=0;expected={(cs,k5,a12) for cs,k5,a12 in domain() if cs[0] and k5>=1}
    assert len(expected)==759
    with View(CONDITIONAL) as view:
        rows=view.array(('increasing_family_upper','rows'));assert len(rows)==759
        got={(tuple(rows[i]['counts']),rows[i]['fixed_A5_count'],rows[i]['key'][9]) for i in range(len(rows))};assert got==expected
        for i in range(len(rows)):
            if i%4!=shard: continue
            r=rows[i];cs=r['counts'];k5=r['fixed_A5_count'];fs=fixed(cs,r['key'][9],k5);rem=r['residual_rank']
            ws={rank:item(rank,fs)+(40*k5 if rank==1 else 0) for rank in ITEMS}
            assert ws=={int(k):v for k,v in r['weights'].items()}
            assert rem==58-sum(t[1] for t in fs)
            v=[t['values'] for t in r['tables']]
            options=[v[0][rem]]
            for reserve in (1,2):
                g0=reserve==2 and bool(cs[1])
                if rem<reserve or (reserve==2 and not g0): continue
                for a2 in range((rem-reserve)//2+1):
                    options.append(reserve*ws[1]+a2*ws[2]+v[1][rem-reserve-2*a2]+64*a2*k5+36*a2*(a2-1)*g0)
            assert bar_fixed(fs)+max(options)+22*cs[-1]<=r['upper']+22*cs[-1]<=64*F(2837,8)
            checked+=1
    return dict(shard=shard,conditional_rows=checked,complete_domain=759,scope='literal barred coefficient equality and conditional-family bound; frozen saved-table upper validity remains an explicit premise')


def parent_lock():
    owner=int(os.environ.get('COXETER_CHECK_OWNER','0'))
    assert owner>1 and os.getppid()==owner,'bounded parent required'
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('r+') as stream:
        try: fcntl.flock(stream,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError: return owner
        else:
            fcntl.flock(stream,fcntl.LOCK_UN)
            raise AssertionError('shared parent lock required')


def main():
    global W,PROJ,ARR
    assert __debug__,'assertions must remain enabled'
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',choices=BATCHES,required=True)
    p.add_argument('--parent-lock',action='store_true',required=True)
    args=p.parse_args();owner=parent_lock();signal.alarm(50)
    for path,pin in PINS.items(): assert digest(path)==pin,('source/input pin',path)
    W=module(LIB,'finite57_literal_weights')
    PROJ=module(PROJECTOR,'finite57_bounded_projection')
    ARR=module(ARRAY,'finite57_bounded_arrays');ARR.PROJ=PROJ
    name=args.batch
    if name=='domination': result=domination()
    elif name=='residual-zero': result=residual_zero()
    elif name=='context-rows': result=context_rows()
    elif name=='final-forests': result=final_forests()
    elif name.startswith('affine'): result=affine(name)
    elif name=='coarse-index': result=coarse_index()
    elif name.startswith('coarse-'): result=coarse(name)
    elif name.startswith('fibers-'): result=fibers(name)
    elif name.startswith('conditional-'): result=conditional(name)
    else: raise AssertionError('unimplemented batch')
    for path,pin in PINS.items(): assert digest(path)==pin,('changed source/input',path)
    assert os.getppid()==owner,'parent exited'
    own=Path(__file__).relative_to(ROOT).as_posix()
    print(json.dumps(dict(status='PASS_DECLARED_LITERAL_BINDING_SCOPE',batch=name,
        source_sha256=digest(own),library_sha256=PINS[LIB],inputs_read=sorted(READ),result=result,
        full_manuscript_verified=False,native_worker_imports=False),sort_keys=True))


if __name__=='__main__': main()
