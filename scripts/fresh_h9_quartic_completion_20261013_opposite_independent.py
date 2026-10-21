#!/usr/bin/env python3
"""Independent reconstruction of the next cap on seed support {6,9}."""
import copy
import json
from fractions import Fraction as Q
from pathlib import Path
from scripts.fresh_h9_quartic_completion_20261013_verify import read, matrix, load_old, node
from scripts.fresh_h9_quartic_completion_20261013_arithmetic import (
    K, Quadratic, O, Z, parse, matvec, pairing, dot)
from scripts.fresh_h9_quartic_completion_20261013_polar_independent import (
    reconstruct as reconstruct_polar, check_archive as check_polar_archive, embed_K)

ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = 'data/enumerations/fresh_h9_quartic_completion_20261013_second_geometry.json'
POLAR = 'data/enumerations/fresh_h9_quartic_completion_20261013_polar.json'
SOURCE = 'data/enumerations/fresh_h9_quartic_completion_20261013_opposite.json'
OUT = ROOT/'results/fresh_h9_quartic_completion_20261013_opposite_independent.json'


def compositum(base, previous):
    """Append a previous real quadratic tower, preserving its chosen roots."""
    if previous is K:
        return base, lambda x: embed_K(base, x)
    parent, embed = compositum(base, previous.base)
    field = Quadratic(parent, embed(previous.d))
    return field, lambda x: (embed(x[0]), embed(x[1]))


def lift_from(field, base, value):
    if field is base:
        return value
    return field.lift(lift_from(field.base, base, value))


def reconstruct():
    polar_result, polar_context = reconstruct_polar()
    check_polar_archive(read(POLAR), polar_result, polar_context)
    polar_field, polar_coefficients = polar_context[:2]
    geometry = read(GEOMETRY)
    case = geometry['cases'][0]
    source = read(geometry['source'])
    base = read(source['base_source'])
    G, I = matrix(base['seed_gram']), matrix(base['seed_gram_inverse'])
    center = tuple(map(parse, base['strict_center_y']))
    assert I[6][9] == I[9][6] == Z
    assert I[6][6] == (Q(4,5), Q(2,5))
    assert I[9][9] == (Q(-9,5), Q(11,10))
    D = K.div(K.sub(K.rational(2), K.mul(K.rational(2), I[6][6])), I[9][9])
    assert D == (Q(12), Q(20)) and K.sign(D) > 0
    E = Quadratic(K, K.rational(2))
    F = Quadratic(E, E.lift(D))
    lift = lambda x: embed_K(F, x)
    pa = [F.zero]*10
    pa[6] = F.neg(F.lift((Z,O)))
    pa[9] = F.neg((E.zero,E.one))
    pa = tuple(pa)
    cf = matvec(F, tuple(tuple(lift(x) for x in row) for row in I), pa)
    assert matvec(F, tuple(tuple(lift(x) for x in row) for row in G), cf) == pa
    assert F.sign(F.sub(dot(F, pa, cf), F.rational(2))) == 0
    assert F.sign(dot(F, cf, tuple(lift(x) for x in center))) > 0
    assert all(F.sign(x) <= 0 for x in pa)
    assert F.mul(pa[6],pa[6]) == F.rational(2)
    assert F.sign(F.sub(F.mul(pa[9],pa[9]),F.rational(4))) > 0

    old = load_old(G)
    previous = []
    for j in case['retained_old_caps']:
        n = old[j]
        field = Quadratic(K,n['d'])
        previous.append(({'kind':'old','index':j},field,tuple(zip(n['a'],n['b'])),n['q']))
    for j in case['selected_new_caps']:
        ref = geometry['cap_references'][j]
        source = read(ref['source'])
        e = Quadratic(K,parse(source['anchor_radicand']))
        row = source['caps'][ref['index']]
        field = Quadratic(e,e.parse(row['discriminant_over_E']))
        previous.append(({'kind':'new','index':j},field,tuple(map(field.parse,row['seed_coefficients'])),K.rational(2)))
    previous.append(({'kind':'polar','index':0},polar_field,polar_coefficients,K.rational(2)))
    rows = []
    for identity, field, coefficients, norm in previous:
        combined, embed = compositum(F,field)
        value = dot(combined,tuple(lift_from(combined,F,x) for x in pa),tuple(map(embed,coefficients)))
        assert combined.sign(value) < 0
        assert combined.sign(combined.sub(combined.mul(value,value),embed_K(combined,K.mul(K.rational(2),norm)))) > 0
        rows.append({**identity,'relation':'dotted'})
    assert len(rows) == 35
    ray = [Z]*10
    ray[6],ray[9] = (Q(2),Q(-1)),O
    ray = tuple(ray)
    assert pairing(K,ray,I,ray) == (Q(1),Q(-1,2))
    assert pairing(K,ray,I,center) == (Q(-7,2),Q(0))
    assert K.sign(pairing(K,ray,I,ray)) > 0
    old_signs = [field.sign(dot(field,coefficients,tuple(embed_K(field,x) for x in ray)))
                 for _,field,coefficients,_ in previous]
    assert all(x >= 0 for x in old_signs)
    assert [previous[j][0] for j,x in enumerate(old_signs) if x == 0] == [
        {'kind':'old','index':67},{'kind':'old','index':310}]
    assert F.sign(dot(F,cf,tuple(lift(x) for x in ray))) < 0
    result = {'status':'PASS_EXACT_CAP','norm':2,'radicand':{'a':'12','b':'20'},
              'seed_relations':['right']*6+['pi/4','right','right','dotted'],
              'all_35_nonseed_pairings':rows,'cuts_exact_opposite_ray':True,
              'strict_common_center':True}
    return result,(F,cf,pa,G,I,center,previous)


def distinct_check(F,pa,G):
    old = load_old(G)
    old += [node(row,G) for row in read('data/enumerations/fresh_h9_joint_completion_20261013.json')['caps']]
    preceding = []
    for n in old:
        field = Quadratic(K,n['d'])
        preceding.append((field,tuple(zip(n['pa'],n['pb'])),n['q']))
    for suffix in ('targeted','adaptive','second'):
        source = read(f'data/enumerations/fresh_h9_quartic_completion_20261013_{suffix}.json')
        e = Quadratic(K,parse(source['anchor_radicand']))
        for row in source['caps']:
            f = Quadratic(e,e.parse(row['discriminant_over_E']))
            preceding.append((f,tuple(map(f.parse,row['seed_pairings'])),K.rational(2)))
    polar = read(POLAR)
    e = Quadratic(K,parse(polar['D']))
    f = Quadratic(e,e.lift(parse(polar['delta_squared'])))
    preceding.append((f,tuple(map(f.parse,polar['seed_pairings'])),K.rational(2)))
    matches = 0
    for field,oldpa,norm in preceding:
        if [j for j,x in enumerate(oldpa) if field.sign(x)] != [6,9]:
            continue
        matches += 1
        combined,embed = compositum(F,field)
        ratios = [combined.mul(combined.mul(embed(x),embed(x)),embed_K(combined,K.div(K.rational(2),norm))) for x in oldpa]
        assert any(combined.sign(combined.sub(lift_from(combined,F,F.mul(x,x)),y))
                   for x,y in zip(pa,ratios,strict=True))
    assert len(preceding) == 3031
    return len(preceding),matches


def check_archive(saved,result,context):
    F,cf,pa,G,I,center,previous = context
    assert saved['polar_source'] == POLAR and saved['facet_count'] == 46
    assert parse(saved['D9']) == (Q(12),Q(20))
    # The independent reconstruction deliberately adjoins radicals in the
    # opposite order to the producer; swap the two formal basis factors.
    oldE = Quadratic(K,(Q(12),Q(20)))
    oldF = Quadratic(oldE,oldE.rational(2))
    def decode(value):
        x = oldF.parse(value)
        return ((x[0][0],x[1][0]),(x[0][1],x[1][1]))
    assert tuple(map(decode,saved['seed_coefficients'])) == cf
    assert tuple(map(decode,saved['seed_pairings'])) == pa
    assert decode(saved['norm']) == dot(F,pa,cf)
    assert decode(saved['center_inequality']) == dot(F,cf,tuple(embed_K(F,x) for x in center))
    assert saved['seed_angles'] == [2,2,2,2,2,2,4,2,2,'dotted']
    assert [r['index'] for r in saved['old_checks']] == [x[0]['index'] for x in previous[:22]]
    assert all(r['kind'] == 'dotted' for key in ('old_checks','new_checks') for r in saved[key])
    assert len(saved['new_checks']) == 13
    source = read(read(GEOMETRY)['source'])
    assert saved['target_source'] == source['target_source']
    targets = read(saved['target_source'])['targets']
    signs = [F.sign(dot(F,cf,tuple(embed_K(F,parse(x)) for x in row['certificate']['y']))) for row in targets]
    assert signs == saved['target_signs']
    assert saved['target_cuts'] == [j for j,x in enumerate(signs) if x < 0]
    assert signs[1] == signs[30] == -1

    volume = saved['volume']
    assert volume['status'] == 'EXACT_STRICT_SPACELIKE_ESCAPE'
    cert = volume['certificate']
    y = tuple(map(parse,cert['y']))
    assert len(y) == 10 and all(K.sign(x) > 0 for x in y)
    norm,cp = pairing(K,y,I,y),pairing(K,y,I,center)
    assert norm == parse(cert['norm']) and K.sign(norm) > 0
    assert cp == parse(cert['center_pairing']) and K.sign(cp) < 0
    ss = [field.sign(dot(field,cc,tuple(embed_K(field,x) for x in y))) for _,field,cc,_ in previous]
    ss.append(F.sign(dot(F,cf,tuple(embed_K(F,x) for x in y))))
    assert ss == cert['nonseed_inequality_signs'] == [1]*36
    return {'strict_future_spacelike_escape':True,'volume_conclusion':'INFINITE_VOLUME',
            'wall_count':46,'saved_target_cuts':saved['target_cuts']}


def main():
    result,context = reconstruct()
    saved = read(SOURCE)
    result.update(check_archive(saved,result,context))
    count,matches = distinct_check(context[0],context[2],context[3])
    result['preceding_entries_checked'] = count
    result['preceding_entries_matching_zero_pattern'] = matches
    assert saved['mirror_check']['preceding_entries'] == count
    assert saved['mirror_check']['matching_zero_pattern'] == matches
    for kind in ('cap','escape'):
        bad = copy.deepcopy(saved)
        if kind == 'cap':
            bad['seed_coefficients'][0][0][0]['a'] = '1'
        else:
            bad['volume']['certificate']['y'] = [{'a':'0','b':'0'}]*10
        try:
            check_archive(bad,result,context)
        except AssertionError:
            result[f'forged_{kind}_rejected'] = True
        else:
            raise AssertionError('Forged certificate accepted')
    result['status'] = 'PASS_EXACT_CAP_AND_INFINITE_VOLUME'
    OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='all_35_nonseed_pairings'}))


if __name__ == '__main__':
    main()
