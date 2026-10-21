#!/usr/bin/env python3
"""Separate exact four-anchor extension; frozen301-parent proof unchanged."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
import fresh_h22_contact_bridge_20261013_verify as E
R=E.R;I=E.I;TAG=E.TAG+'_big_verify';OUT='results/'+TAG+'.json';PRIMARY=E.OUT;GRAPH='results/fresh_h22_covering_cells_20261013.json'
def run(verify=False):
 base=R.read(PRIMARY)
 for p,s in base['sources'].items():assert R.sha(p)==s,p
 source,M,walls,cov,nn,h,points=E.setup();types=source['selected_original_types']+[23,24];anchors=[812+2*types.index(k)for k in(0,2,6,22)];assert anchors==[812,816,824,854]
 graph=R.read(GRAPH)['frozen_860_graph']
 for p,s in graph['source_sha256'].items():assert R.sha(p)==s,p
 big=[r for r in graph['components']if len(r['indices'])==305];assert len(big)==1 and big[0]['parent_null_caps']==list(range(301))and big[0]['additional_types']==[0,2,6,22]
 by={(r['vertex'],r['released']):r for r in base['rows']};assert len(by)==88;rows=[];counts=Counter()
 for vi,src in enumerate(source['compact_H2_vertices']):
  inc=src['incident'];A=[]
  for j in inc:
   row=[]
   for k in inc:
    z=R.pair(walls[j],cov[k]);assert set(z)<={1};row.append(z.get(1,F(0)))
   A.append(row)
  R.pd(A);AI=E.inverse(A);assert all(x>=0 for row in AI for x in row);v=R.vm(points[vi],R.sqrtq(F(-1,src['norm'])));assert R.pair(v,R.mv(M,v))==R.C(-1)
  for k,released in enumerate(inc):
   prior=by[vi,released];edge=[j for j in inc if j!=released];assert prior['edge']==edge;dual=E.combine(AI[k],[walls[j]for j in inc]);e=R.vm(dual,R.sqrtq(1/AI[k][k]));assert R.pair(e,R.mv(M,e))==R.C(1)and R.pair(v,R.mv(M,e))=={}
   assert all(R.pair(e,cov[j])==R.pair(v,cov[j])=={}for j in edge)
   minus=E.vsub(v,e);plus=R.vadd(v,e)
   def coefficients(tag):
    if tag=='upper_height':return E.neg(minus[-1]),E.neg(plus[-1])
    assert isinstance(tag,int);return R.pair(minus,cov[tag]),R.pair(plus,cov[tag])
   lower=prior['lower_source'];li=I.rat(1)if lower=='outgoing'else E.cutoff(*coefficients(lower),True);assert E.enc(li)==prior['lower_cutoff'];L=max(I.S,li[0]);rr={'vertex':vi,'released':released,'lower_cutoff':E.enc(li),'anchors':[]}
   if prior.get('conclusion')=='EMPTY_ACUTE_UPPER_CAP_INTERVAL':
    ui=E.cutoff(*coefficients(prior['upper_source']),False);assert L>ui[1];rr['conclusion']='EMPTY_EDGE_FAMILY';counts['empty_edges']+=1
   else:
    counts['nonempty_edges']+=1
    for anchor in anchors:
     a,b=coefficients(anchor);U=E.finite_upper(E.iv(a),E.iv(b),nn[anchor]);assert U[1]<L,(vi,released,anchor,U,L);rr['anchors'].append({'anchor':anchor,'type':types[(anchor-812)//2],'finite_upper_bound':E.enc(U),'strict_lower_minus_upper':str(L-U[1])});counts['empty_anchor_intervals']+=1
   rows.append(rr)
 assert counts=={'empty_edges':10,'nonempty_edges':78,'empty_anchor_intervals':312}
 out={'status':'PASS_EXACT_RIGHT21_DIRECT_BRIDGE_TO_ENTIRE_305_CAP_COMPONENT_EXCLUDED','sources':{p:R.sha(p)for p in(PRIMARY,E.SOURCE,GRAPH,'scripts/'+E.TAG+'_verify.py','scripts/'+TAG+'.py')},'interval_denominator':str(I.S),'counts':dict(counts),'additional_anchor_types':[0,2,6,22],'additional_physical_anchors':anchors,'component_caps':305,'rows':rows,'scope':'Extends the frozen301-parent obstruction by all4other upper caps in the certified305-cap component. Each of312 remaining real finite/parallel contact intervals is empty. Combined, no positive-height cap right to a21-wall edge at any of the four specified vertices can retain860 acutely and directly contact this305-cap component at a finite angle or in parallel. The restriction does not exclude bridges using nonright profiles, other faces/vertices, lower-height caps here, indirect chains of newcaps or retained-wall replacement.'}
 if verify:assert out==R.read(OUT)
 else:(R.P/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],dict(counts),flush=True);return out
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');run(ap.parse_args().verify)
