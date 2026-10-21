#!/usr/bin/env python3
"""Fresh physical158 export from the frozen independent158 reconstruction."""
import argparse,ast,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_e8_exact_cluster_20261017_geometry'
PREVIOUS='scripts/fresh_h20_compact_quad_cap_20261016_review.py'
def run():
 path=ROOT/PREVIOUS;tree=ast.parse(path.read_text())
 fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name=='run')
 cut=max(i for i,n in enumerate(fn.body)if isinstance(n,ast.Assign)and any(isinstance(t,ast.Name)and t.id=='paths'for t in n.targets))
 fn.body=fn.body[:cut]+ast.parse('def f():\n return locals()').body[0].body
 ns={'__name__':'_frozen_independent158','__file__':str(path)}
 exec(compile(ast.fix_missing_locations(tree),str(path),'exec'),ns)
 D=ns['run']();old=D['L'];V=ns['vec'] if 'vec'in ns else D['vec']
 paths=[f'scripts/{TAG}.py',PREVIOUS,'results/fresh_h20_compact_quad_cap_20261016_review.json']
 return dict(status='PASS_REPLAYED_INDEPENDENT_158_PHYSICAL_GEOMETRY_EXPORT',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},form=[[str(x)for x in row]for row in old['Q']],normals=[V(v)for v in D['full']],norms=list(map(str,D['norms'])),center=V(D['H']),S16=sorted(D['S16']),old_escape=V(old['E7']),old_escape_norm=str(old['q7']),old_cusps=[V(v)for v in old['cusp_vectors']+[old['nu2']]],scope='Replays the frozen independent158 physical construction and exports outward rational enclosures of its exact algebraic vectors. The defining exact recipes and parameter isolation remain inherited from that reviewed source chain.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
