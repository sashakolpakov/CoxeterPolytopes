#!/usr/bin/env python3
"""Independent finite-profile/orbit enumeration and all constructed cones."""
from fresh_h9_quartic_edge_walk_20261013_independent import *

OUT=ROOT/'results/fresh_h9_quartic_edge_walk_20261013_profiles_independent.json'


def profile_key(row):return dec(row['norm']),vec(row['fixed_magnitudes'])


def replay_profiles(path,next_edge=False):
    data=json.loads(path.read_text());source=json.loads((ROOT/data['source']).read_text())
    base_data=source['cases'][data['source_case']]if next_edge else source
    host=Host(2);face=Host(1);walls=tuple(vec(r)for r in base_data['all_walls'])
    fw=tuple(vec(r)for r in base_data['forced_H5_face']['all_walls']);edge=tuple(fw[i]for i in data['edge_indices'])
    gram=tuple(tuple(face.pair(r,s)for s in edge)for r in edge);assert gram==matrix(data['edge_gram']);pd(gram)
    inv=a.minverse(gram);E=tuple(zip(*edge))
    projection=tuple(tuple(sub(ident(6)[i][j],z)for j,z in enumerate(row))for i,row in enumerate(mm(mm(mm(E,inv),edge),face.g)))
    if next_edge:
        j=matrix(data['generator']);assert projection==matrix(data['plane_projection'])
        seed=vec(data['unit_norm_seed']);old_seed=fw[data.get('seed_wall_index',12)]
        basis=matrix(data['complement_basis']);bg=tuple(tuple(face.pair(r,s)for s in basis)for r in basis)
        assert bg==matrix(data['complement_gram'])
        assert all(face.pair(r,e)==Z for r in basis for e in edge)
        determinant=sub(sq(bg[0][1]),mul(bg[0][0],bg[1][1]));factor=dec(data['discriminant_square_factor'])
        assert sq(factor)==divide(determinant,ALPHA)
    else:
        trans=json.loads((ROOT/'results/fresh_h9_quartic_edge_walk_20261013_halfstep.json').read_text())
        j=matrix(trans['generator']);assert projection==matrix(trans['plane_projection'])
        seed=vec(data['orthogonal_unit_norm_seed']);old_seed=fw[6]
    assert seed==mv(projection,old_seed)and face.pair(seed,seed)==O
    assert mm(j,j)==tuple(tuple(mul(ALPHA,z)for z in row)for row in projection)
    gj=mm(face.g,j);assert all(add(gj[i][k],gj[k][i])==Z for i in range(6)for k in range(6))
    assert mm(j,projection)==j
    powers=data['unit_exponents']if next_edge else data['unit_orbit_exponents']
    assert powers==list(range(-2,3)if next_edge else range(-4,5))
    y0=divide(T,sub(scale(O,2),T));s0=divide(scale(O,2),sub(scale(O,2),T));y,s=O,Z;hs={0:ident(6)}
    for power in range(1,max(powers)+1):
        y,s=add(mul(y,y0),mul(ALPHA,mul(s,s0))),add(mul(y,s0),mul(s,y0))
        for direction in(-1,1):
            hs[power*direction]=tuple(tuple(add(add(O if i==k else Z,mul(sub(y,O),projection[i][k])),scale(mul(s,j[i][k]),direction))for k in range(6))for i in range(6))
    saved={profile_key(row):row for row in data['profiles']};expected=set()
    cov=[mv(host.g,r)for r in walls];norms=[dot(r,c)for r,c in zip(walls,cov)];oldkeys={old.key(r)for r in walls}
    candidates=set();dictionary={}
    for q in POOL:
        alph=[{r for label in ALLOWED if(r:=field_sqrt(scale(mul(mul(q,face.pair(e,e)),label),Q(1,4))))is not None and integral(r)}for e in edge]
        for fixed in product(*alph):
            coefficients=mv(inv,tuple(neg(z)for z in fixed))
            r0=tuple(sum_field(tuple(mul(c,e[k])for c,e in zip(coefficients,edge)))for k in range(6))
            residual=sub(q,face.pair(r0,r0))
            if any(sign(residual,e)<=0 for e in range(4)):continue
            key=q,fixed;expected.add(key);row=saved[key];factor=field_sqrt(residual)
            assert vec(row['fixed_component'])==r0 and dec(row['orthogonal_norm'])==residual
            assert (None if row['square_factor']is None else dec(row['square_factor']))==factor
            if factor is None:continue
            for power in powers:
                v=mv(hs[power],seed)
                for direction in(-1,1):
                    r=tuple(add(x,scale(mul(factor,y),direction))for x,y in zip(r0,v))
                    assert face.pair(r,r)==q
                    if sign(r[0])<=0 or not all(map(integral,r)):continue
                    if all(integral(divide(z,T))for z in r)or not all(integral(divide(scale(z,2),q))for z in mv(face.g,r)):continue
                    r9=(r[0],*r[1:5],*(Z,)*4,r[5]);candidates.add((key,power,direction,r9))
                    ray=old.key(r9)
                    if ray not in oldkeys and compat(r9,q,cov,norms):dictionary[ray]=(key,power,direction,r9)
    assert set(saved)==expected
    candrows=data['integral_candidates']if next_edge else data['integral_primitive_candidates']
    def row_tuple(row):return(profile_key(data['profiles'][row['profile_index']]),row['unit_power'],row['orthogonal_sign'],vec(row['vector']))
    assert candidates==set(map(row_tuple,candrows))
    assert set(dictionary.values())==set(map(row_tuple,data['compatible_dictionary']))
    cases=[]
    for case in data['cases']:
        row=data['compatible_dictionary'][case['dictionary_index']];r=vec(row['vector']);q=dec(row['norm']);new=tuple(vec(v)for v in case['new_roots'])
        assert set(new)==orbit(r)and len(new)==4
        assert tuple(vec(v)for v in case['all_walls'])==walls+new
        assert tuple(vec(v)for v in case['forced_H5_face']['all_walls'])==tuple((*w[:5],w[9])for w in walls+new if all(z==Z for z in w[5:9]))
        assert divide(mul(ALPHA,sq(r[0])),q)==dec(row['height'])
        cut=face.pair((*r[:5],r[9]),vec(base_data['forced_H5_face']['geometry']['vector']))
        assert cut==dec(row['old_H5_escape_pairing'])and sign(cut)>0
        cases.append({'H9':geometry_check(case,2),'H5':geometry_check(case['forced_H5_face'],1),
                      'H9_vertex':check_vertex(case['finite_H9_vertex'],case,host),
                      'H5_vertex':check_vertex(case['forced_H5_face']['finite_vertex'],case['forced_H5_face'],face)})
    return {'source':str(path.relative_to(ROOT)),'finite_profiles':len(expected),
            'square_residual_profiles':sum(row['square_factor']is not None for row in data['profiles']),
            'unit_powers':powers,'primitive_integral_candidates':len(candidates),'compatible_dictionary':len(dictionary),'cases':cases}


def sum_field(values):
    out=Z
    for value in values:out=add(out,value)
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();start=time.monotonic()
    cases=[replay_profiles(ROOT/'results/fresh_h9_quartic_edge_walk_20261013_profiles.json'),
           replay_profiles(ROOT/'results/fresh_h9_quartic_edge_walk_20261013_next.json',True),
           replay_profiles(ROOT/'results/fresh_h9_quartic_edge_walk_20261013_f4.json',True)]
    out={'status':'PASS_INDEPENDENT_FINITE_PROFILE_UNIT_ORBITS','stages':cases,'seconds':round(time.monotonic()-start,3),
         'scope':'Complete finite-label profiles and complete stated bounded unit-orbit grids only. No claim to enumerate all binary norm solutions. Separate Fraction reconstruction, all saved root/angle/strict-escape and actual finite-vertex checks.'}
    if args.verify:
        previous=json.loads(OUT.read_text());previous.pop('seconds');out.pop('seconds');assert out==previous
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)
