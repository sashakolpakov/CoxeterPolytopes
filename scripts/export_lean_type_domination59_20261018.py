#!/usr/bin/env python3
"""Serialize paired actual/representative costs; Lean checks all comparisons.

Pairs preserve the exact saved actual-type order. The right entry is the
saved representative cost pulled back along the saved representative map.
This is transcription of data, not a new optimization or geometric proof.
"""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/certificates/finite_volume58_arithmetic_20261018.json'
EXPECTED='4478785c83ea7da75752806b0cfc3f78e321b42b2bb6a267fe94442ac57dfe73'
DEST=ROOT/'formal/CoxeterProofs/Generated/TypeDomination59.lean'
def pairs(xs):return '['+','.join(f'({a},{b})'for a,b in xs)+']'
def render():
 assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==EXPECTED
 x=json.loads(SOURCE.read_text())
 for p,h in x['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 rep=x['representative_map'];assert len(rep)==len(x['actual_bases'])==180
 bases=list(zip(x['actual_bases'],[x['upper_bases'][r]for r in rep]))
 rows=[list(zip(row,[x['upper_cross'][rep[i]][rep[j]]for j in range(180)]))for i,row in enumerate(x['actual_cross'])]
 assert len(rows)==180 and all(len(row)==180 for row in rows)
 lines=['import CoxeterProofs.ForestDomination','set_option maxRecDepth 100000',
  'set_option maxHeartbeats 0','namespace CoxeterProofs.Generated.TypeDomination59',
  'open CoxeterProofs CoxeterProofs.ForestDomination',
  'def pairCheck (p : Nat × Nat) : Bool := decide (p.1 ≤ p.2)',
  'def basePairs : List (Nat × Nat) := '+pairs(bases),
  'theorem bases_checked : basePairs.all pairCheck = true := by decide']
 names=[]
 for start in range(0,180,10):
  name=f'rows{start//10:02}';names.append(name)
  lines.append('def '+name+' : List (List (Nat × Nat)) := [\n'+',\n'.join(pairs(row)for row in rows[start:start+10])+']')
  lines.append(f'theorem {name}_checked : {name}.all (fun row => row.all pairCheck) = true := by decide')
 lines += ['def crossPairs : List (List (Nat × Nat)) := '+' ++ '.join(names),
  'theorem crosses_checked : crossPairs.all (fun row => row.all pairCheck) = true := by',
  '  simp only [crossPairs, List.all_append, '+', '.join(n+'_checked'for n in names)+', Bool.and_self]',
  'def baseAt (i : Nat) : Nat × Nat := basePairs[i]?.getD (0,0)',
  'def crossAt (i j : Nat) : Nat × Nat := (crossPairs[i]?.getD [])[j]?.getD (0,0)',
  'def actualBase (i : Nat) : Nat := (baseAt i).1',
  'def upperBase (i : Nat) : Nat := (baseAt i).2',
  'def actualCross (i j : Nat) : Nat := (crossAt i j).1',
  'def upperCross (i j : Nat) : Nat := (crossAt i j).2',
  'theorem base_domination (i : Nat) : actualBase i ≤ upperBase i := by',
  '  apply getD_property basePairs (0,0) (fun p => p.1 ≤ p.2)',
  '  · intro p hp; exact of_decide_eq_true ((List.all_eq_true.mp bases_checked) p hp)',
  '  · decide',
  'theorem pair_domination (i j : Nat) : actualCross i j ≤ upperCross i j := by',
  '  have row : (crossPairs[i]?.getD []).all pairCheck = true :=',
  '    getD_property crossPairs [] (fun r => r.all pairCheck = true)',
  '      (List.all_eq_true.mp crosses_checked) (by decide) i',
  '  apply getD_property (crossPairs[i]?.getD []) (0,0) (fun p => p.1 ≤ p.2)',
  '  · intro p hp; exact of_decide_eq_true ((List.all_eq_true.mp row) p hp)',
  '  · decide',
  'theorem every_ordinary_forest_dominated (xs : List Nat) :',
  '    cost actualBase actualCross xs ≤ cost upperBase upperCross xs := by',
  '  apply ordinary_cost_domination',
  '  · intro i _; exact base_domination i',
  '  · intro i _ j _; exact pair_domination i j',
  'end CoxeterProofs.Generated.TypeDomination59','']
 return '\n'.join(lines)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();s=render()
 if a.write:DEST.write_text(s)
 else:assert DEST.read_text()==s
 print('PASS export replay: 180 paired actual/representative bases and32400 cross pairs; contextual domination separate')
