#!/usr/bin/env python3
"""Direct rank14 simplex domain: inherited Lorentz cores, all finite spectators."""
import argparse, hashlib, json
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations
from math import comb
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_rank14_core_spectator_domain_20261017'
INPUT = 'fresh_sqrt5_rank13_full_triangle_simplex_separation_20261017'
C = B.C


def sha(path):
    h = hashlib.sha256()
    with (ROOT/path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def refinement(g):
    colors = {i: 0 for i in g}
    for _ in range(len(g)):
        codes = {i: (colors[i], tuple(sorted((m, colors[j]) for j, m in row.items()))) for i, row in g.items()}
        ids = {v: i for i, v in enumerate(sorted(set(codes.values())))}
        colors = {i: ids[codes[i]] for i in g}
    key = (len(g), tuple(sorted(Counter(colors.values()).items())), tuple(sorted(
        (min(colors[i], colors[j]), max(colors[i], colors[j]), m)
        for i, row in g.items() for j, m in row.items() if i < j)))
    return colors, key


def isomorphism(g, h, cg, ch):
    assert len(g) == len(h)
    mapping, used = {}, set()
    def visit():
        if len(mapping) == len(g):
            return dict(mapping)
        choices = []
        for i in g:
            if i in mapping:
                continue
            candidates = [j for j in h if j not in used and cg[i] == ch[j] and
                          all(g[i].get(k, 2) == h[j].get(v, 2) for k, v in mapping.items())]
            if not candidates:
                return None
            choices.append((len(candidates), -len(g[i]), i, candidates))
        _, _, i, candidates = min(choices)
        for j in candidates:
            mapping[i] = j; used.add(j)
            answer = visit()
            if answer is not None:
                return answer
            used.remove(j); del mapping[i]
        return None
    answer = visit()
    if answer is not None:
        assert len(set(answer.values())) == len(g)
        assert all(g[i].get(j, 2) == h[answer[i]].get(answer[j], 2)
                   for i, j in combinations(g, 2))
    return answer


def run():
    prior = B.read(INPUT)
    models = [(name, graph) for name, graph in C.models() if len(graph) <= 14]
    ranks = {name: len(graph) for name, graph in models}
    codes = {name: C.code(graph) for name, graph in models}
    finite = {code: name for name, code in codes.items()}
    def component_name(g):
        if sum(map(len, g.values())) != 2*len(g)-2:
            return None
        return finite.get(C.code(g))
    def endpoint(g, removed):
        end = {i: {j: m for j, m in row.items() if j != removed}
               for i, row in g.items() if i != removed}
        pieces = C.pieces(end)
        names = [component_name(p) for p in pieces]
        return None if None in names else (end, pieces, tuple(sorted(names)))

    @lru_cache(None)
    def forests(rank, first=0):
        if rank == 0:
            return ((),)
        answer = []
        for j in range(first, len(models)):
            name, graph = models[j]
            if len(graph) <= rank:
                answer.extend(tuple(sorted((name,)+tail)) for tail in forests(rank-len(graph), j))
        return tuple(sorted(answer))
    spectator_sets = {n: forests(n) for n in range(15)}
    for n, lists in spectator_sets.items():
        assert len(lists) == len(set(lists))
        assert all(sum(ranks[name] for name in names) == n for names in lists)
    assert len(spectator_sets[13]) == 6789

    cores, graphs, colors, buckets, references = [], [], [], defaultdict(list), []
    for star in prior['stars']:
        g = {i: {} for i in range(14)}
        for i, j, m in star['edges']:
            g[i][j] = g[j][i] = m
        components = C.pieces(g)
        indefinite = [p for p in components if component_name(p) is None]
        assert len(indefinite) == 1
        oldcore = indefinite[0]; ids = sorted(oldcore); reindex = {i: j for j, i in enumerate(ids)}
        core = {reindex[i]: {reindex[j]: m for j, m in row.items()} for i, row in oldcore.items()}
        spectator = sorted(component_name(p) for p in components if p is not oldcore)
        assert None not in spectator
        oldS = [x['omitted'] for x in star['spherical_deletions']]
        assert set(oldS) <= set(oldcore)
        cg, key = refinement(core); ci = None; mp = None
        for candidate in buckets[key]:
            mp = isomorphism(core, graphs[candidate], cg, colors[candidate])
            if mp is not None:
                ci = candidate
                break
        if ci is None:
            ci = len(cores); S = []; deletion = {}
            for removed in core:
                e = endpoint(core, removed)
                if e is not None:
                    S.append(removed); deletion[removed] = e
            S.sort()
            assert sorted(S) == sorted(reindex[i] for i in oldS)
            assert len(S) >= 3
            corners = []
            for removed in S:
                end, pieces, names = deletion[removed]
                for pair in combinations([i for i in S if i != removed], 2):
                    corners.append(dict(vertex_omitted=removed, marked_corner_roots=list(pair),
                                        forest=list(names), marked_forest_code=sorted(C.code(p, pair) for p in pieces)))
            assert len(corners) == 3*comb(len(S), 3)
            cores.append(dict(index=ci, nodes=len(core), edges=[[i,j,m] for i,row in core.items()
                for j,m in sorted(row.items()) if i<j], original_star=star['index'],
                original_core_roots=ids, spherical_deletions=S, simplex_dimension=len(S)-1,
                required=2*comb(len(S),3), corner_occurrences=corners, spectator_rank=15-len(core)))
            graphs.append(core); colors.append(cg); buckets[key].append(ci); mp={i:i for i in core}
        assert sorted(mp[reindex[i]] for i in oldS) == cores[ci]['spherical_deletions']
        references.append(dict(star=star['index'], core=ci, core_roots=ids,
                               core_isomorphism=[mp[i] for i in core], old_spectator_forest=spectator))
        if len(references)%2000 == 0:
            print('RANK14_CORES',len(references),len(cores),flush=True)

    contexts, whole_forests = [], set()
    for core in cores:
        for spectator in spectator_sets[core['spectator_rank']]:
            contexts.append(dict(index=len(contexts), core=core['index'], spectator_forest=list(spectator)))
            for corner in core['corner_occurrences']:
                forest = tuple(sorted(corner['forest']+list(spectator)))
                assert sum(ranks[name] for name in forest) == 14
                whole_forests.add(forest)
    assert len({(r['core'],tuple(r['spectator_forest'])) for r in contexts}) == len(contexts)
    paths = [f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{C.TAG}.py',f'results/{INPUT}.json']
    return dict(status='PASS_EXACT_RANK14_INHERITED_CORE_SPECTATOR_DOMAIN',
        sources={p:sha(p) for p in paths},dimension=14,target=14,pairs_per_vertex=91,
        counts=dict(inherited_stars=len(prior['stars']),connected_lorentz_cores=len(cores),
            distinct_full_stars=len(contexts),corner_forests=len(whole_forests),
            all_finite_rank14_forests=len(spectator_sets[14]),
            cores_by_order=dict(sorted(Counter(c['nodes'] for c in cores).items())),
            full_stars_by_simplex_dimension=dict(sorted(Counter(cores[c['core']]['simplex_dimension'] for c in contexts).items()))),
        finite_spectator_counts={str(n):len(v) for n,v in spectator_sets.items()},
        cores=cores,inherited_star_maps=references,contexts=contexts,
        whole_forests=[list(x) for x in sorted(whole_forests)],
        proof='Each inherited full star has exactly one connected Lorentz component; the other components are spherical. Weighted graph isomorphisms identify cores and all their spherical one-root deletions. Orthogonal finite spectator forests of complementary rank are enumerated completely. Their direct sum preserves physical Lorentz signature, conjugate positivity, cyclic integrality and the compact projected simplex. A one-root deletion is spherical exactly when its removed root is a listed core deletion. A full-star isomorphism preserves the unique Lorentz component, so a core isomorphism class together with a finite spectator multiset identifies the complete full star without duplicates. Corner occurrences retain the literal omitted root and marked pair; the forthcoming vertex model must count all91 pairs.',
        scope='Direct rank14 partial domain with no determinant partition or optimization. Complete only for all finite spectators of the inherited connected cores. New connected order15 cores, other missing connected cores, quadrilateral stars and general faces are not covered. No dimension bound or simultaneous global realization is asserted.')


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data==json.loads(out.read_text())
    else:
        out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],flush=True)
