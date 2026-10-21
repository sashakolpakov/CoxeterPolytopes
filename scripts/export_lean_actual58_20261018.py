#!/usr/bin/env python3
"""Serialize all1549 actual residual tables from the complete rank58 count.

Identical numeric (rank,charge) items are identified, retaining a metadata
map for every original type. This is exact transcription, no optimization.
"""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_d17_a16_coupled_ambient58_frontier_20261018.json'
PIN='0fb359b3b07d6fbabb31adbf611d769f41ebe855075076301975e086055b1223'
DEST=ROOT/'formal/CoxeterProofs/Generated/Actual58.lean'
META=ROOT/'formal/certificates/actual58_tables_source.json'
def render():
 assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==PIN
 doc=json.loads(SOURCE.read_text())
 for p,h in doc['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 tables=[];metadata=[]
 for ci,c in enumerate(doc['finite']['classes']):
  for ri,r in enumerate(c['rows']):
   for ti,t in enumerate(r['tables']):
    cap=len(t['values'])-1;assert cap==c['residual_rank']
    items=sorted({(typ[1],w)for typ,w in zip(t['types'],t['weights'],strict=True)if typ[1]<=cap})
    assert all(w>=0 for n,w in items)and all(v is None or v>=0 for v in t['values'])
    metadata.append(dict(class_index=ci,fiber_index=ri,table_index=ti,cap=cap,key=t['key'],
       original_types=t['types'],original_weights=t['weights'],allowed_numeric_items=items))
    vals='['+','.join('none'if v is None else f'some {v}'for v in t['values'])+']'
    tables.append(f'⟨{cap},['+','.join(f'⟨{n},{w}⟩'for n,w in items)+f'],{vals}⟩')
 assert len(tables)==1549
 lines=['import CoxeterProofs.OptionalTableBounds','set_option maxRecDepth 100000',
  'set_option maxHeartbeats 0','namespace CoxeterProofs.Generated.Actual58',
  'open CoxeterProofs CoxeterProofs.OptionalTable'];names=[]
 for k in range(0,len(tables),32):
  name=f'chunk{k//32:02}';names.append(name)
  lines += [f'def {name} : List Table := [\n'+',\n'.join(tables[k:k+32])+']',
    f'theorem {name}_checked : {name}.all tableCheck = true := by decide']
 lines += ['def tables : List Table := '+' ++ '.join(names),
  'theorem all_checked : tables.all tableCheck = true := by',
  '  simp only [tables, List.all_append, '+', '.join(n+'_checked'for n in names)+', Bool.and_self]',
  'theorem table_count : tables.length = 1549 := by decide',
  'theorem every_table_upper (t : Table) (member : t ∈ tables) (xs : List ChargedItem)',
  '    (membership : ∀ a ∈ xs, a ∈ t.allowed) (capacity : totalRank xs ≤ t.cap) :',
  '    ∃ z, lookup t.values (totalRank xs) = some z ∧ totalCharge xs ≤ z :=',
  '  forest_upper t.cap t.allowed t.values ((List.all_eq_true.mp all_checked) t member)',
  '    xs membership capacity',
  'end CoxeterProofs.Generated.Actual58','']
 meta=dict(status='EXACT_SERIALIZATION_NOT_COMPLETE_FINITE_PROOF',source=str(SOURCE.relative_to(ROOT)),
  source_sha256=PIN,table_count=1549,original_type_maps=metadata,
  scope='Used-cap numeric residual inequalities including unreachable ranks; fixed-domain and contextual coverage/geometric interpretation separate.')
 return '\n'.join(lines),json.dumps(meta,separators=(',',':'))+'\n'
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();lean,meta=render()
 for path,body in((DEST,lean),(META,meta)):
  if a.write:path.write_text(body)
  else:assert path.read_text()==body,path
 print('PASS exact export of1549 actual rank58 residual tables; full finite domain/context theorem separate')
