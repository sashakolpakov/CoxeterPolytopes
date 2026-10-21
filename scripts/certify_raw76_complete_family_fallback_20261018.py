#!/usr/bin/env python3
"""Full-matrix exact fallback after the single proposal timed out.

No optimizer call. The old primal extends by zero; the already proved
family dual repair is checked against every column of the changed matrix.
The original failed proposal is retained byte-for-byte in a separate file.
"""
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_raw76_complete_family_full_model_20261018 as G

ROOT=G.ROOT;OUT=G.OUT
PROPOSAL=ROOT/f'results/{G.TAG}_proposal.json'
TAG='certify_raw76_complete_family_fallback_20261018'

if __name__=='__main__':
    saved_bytes=OUT.read_bytes();saved=json.loads(saved_bytes)
    assert saved['status']=='UNRESOLVED_SINGLE_RAW76_COMPLETE_FAMILY_SOLVE'
    assert saved['solver_invocations']==1 and saved['certificates']==[]
    assert saved['solver']['status']==1 and not saved['solver']['success']
    assert saved['solver']['objective']is None
    assert not PROPOSAL.exists()
    d,prepared=G.build();assert prepared==saved['prepared']
    old=G.M.read(G.B.TAG);family=G.M.read(G.FAMILY)
    primal=next(c for c in old['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    assert all(j<439856 for j,n in primal['rounded_variable_numerators'])
    rep=family['exported']['repair'];normalizer=sum((F(x)for v,x in rep['vertex_multipliers']),F())
    assert str(normalizer)==rep['normalizer']
    dual=dict(kind='EXACT_PARTIAL_MODEL_DUAL',face_multipliers=rep['face_multipliers'],
        vertex_multipliers=rep['vertex_multipliers'],objective_numerator=rep['objective_numerator'],
        normalizer=str(normalizer),bound=rep['bound'],excludes_target12=F(rep['bound'])>12,
        recovery='Frozen family fixed-face-mass monotone repair, independently checked against the entire changed matrix; not a solver optimum.')
    for c in(primal,dual):G.P.verify(d,c)
    lower=F(dual['bound']);upper=F(primal['maximum']);assert 12<lower<=upper
    PROPOSAL.write_bytes(saved_bytes)
    sources={str(PROPOSAL.relative_to(ROOT)):sha256(saved_bytes).hexdigest(),
             f'scripts/{TAG}.py':sha256(Path(__file__).read_bytes()).hexdigest()}
    saved.update(status='PASS_EXACT_RAW76_FAMILY_FALLBACK_INTERVAL_OPTIMUM_UNRESOLVED',
        certificates=[primal,dual],certification_sources=sources,
        certification_method='No second solve. Literal inherited primal with new variables zero; frozen family dual repair; complete exact primal/dual matrix verification.',
        optimization_resolved=False,changed_face_masses=False,
        exact_interval=dict(lower=str(lower),upper=str(upper),gap=str(upper-lower)))
    OUT.write_text(G.P.BASE.dump(saved)+'\n')
    print(saved['status'],saved['exact_interval'],'decimals',float(lower),float(upper),flush=True)
