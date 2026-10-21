#!/usr/bin/env python3
"""One guarded dual proposal, then rational weak duality on every support column."""
import argparse
import fcntl
from fractions import Fraction as F
from hashlib import file_digest
import json
import math
import os
from pathlib import Path
import resource
import selectors
import signal
import subprocess
import sys
import threading
import time
import warnings
from collections import defaultdict
import fresh_compact21_priority_support_scout_20261018 as S

ROOT=S.ROOT
TAG='fresh_compact21_priority_support_exact_dual_20261018'
LIMIT=S.LIMIT
OUTPUT_LIMIT=S.OUTPUT_LIMIT


def verify_dual(data):
    coefficients=defaultdict(F)
    weights={int(i):F(v)for i,v in data['exact_dual_rows']}
    assert len(weights)==len(data['exact_dual_rows'])
    assert all(i>=0 and y<0 for i,y in weights.items())
    value=F()
    count=0
    for i,(row,rhs,_)in enumerate(S.rows()):
        count+=1
        y=weights.get(i,F())
        value+=y*rhs
        for j,a in row.items():
            coefficients[j]+=y*a
    assert count==2805 and all(i<count for i in weights)
    assert all(v<=F(j==93045)for j,v in coefficients.items())
    assert value==F(data['exact_lower_bound'])
    return value


def forest_metadata():
    seed=json.loads((ROOT/S.SEED).read_text())
    target=seed['target_forest_term']
    result={}
    new=0
    for old,term in S.M.objects('terms',4):
        meta=term['metadata']
        if meta['kind']=='forest':
            result[new]={'types':meta['types'],'old_term':old}
            if old==target:
                result[new]['state']='active'
                result[new+1]={'types':meta['types'],'old_term':old,'state':'inactive'}
        new+=2 if old==target else 1
    assert new==2805
    return result


def worker():
    owner=int(os.environ['COXETER_SUPPORT_SCOUT_OWNER'])
    assert os.getppid()==owner
    def owner_watch():
        while True:
            if os.getppid()!=owner:
                os.killpg(os.getpgrp(),signal.SIGKILL)
            time.sleep(.2)
    threading.Thread(target=owner_watch,daemon=True).start()
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    saved=list(S.rows())
    columns=sorted({j for row,_,_ in saved for j in row})
    index={j:i for i,j in enumerate(columns)}
    assert len(saved)==2805 and len(columns)==6121
    rr,cc,vv=[],[],[]
    for i,(row,_,_)in enumerate(saved):
        for j,a in row.items():
            rr.append(i);cc.append(index[j]);vv.append(float(a))
    assert len(vv)==17387
    matrix=coo_matrix((vv,(rr,cc)),shape=(len(saved),len(columns))).tocsr()
    objective=np.zeros(len(columns));objective[-1]=1
    with warnings.catch_warnings(record=True)as caught:
        result=linprog(objective,A_ub=matrix,
                       b_ub=np.asarray([float(rhs)for _,rhs,_ in saved]),
                       bounds=(0,None),method='highs',
                       options={'time_limit':45.0,'threads':1,'disp':False})
    assert result.success, f'Numerical dual proposal failed: {result.message}'
    grid=10**9
    face,forest={},{}
    for i,(_,_,kind)in enumerate(saved):
        mass=max(0.0,-float(result.ineqlin.marginals[i]))
        value=F(math.floor(mass*grid)if kind=='face'else math.ceil(mass*grid),grid)
        if value:
            (face if kind=='face'else forest)[i]=value
    demand,capacity=defaultdict(F),defaultdict(F)
    for i,mass in face.items():
        for j,a in saved[i][0].items():
            assert j!=93045 and a<=0
            demand[j]-=mass*a
    for i,mass in forest.items():
        assert saved[i][0][93045]==-1
        for j,a in saved[i][0].items():
            if j!=93045:
                assert a>=0
                capacity[j]+=mass*a
    scale=min([F(1)]+[capacity[j]/v for j,v in demand.items()if v])
    normalizer=max(F(1),sum(forest.values(),F()))
    weights={i:-v*scale/normalizer for i,v in face.items()if v*scale}
    weights.update({i:-v/normalizer for i,v in forest.items()})
    bound=sum((y*saved[i][1]for i,y in weights.items()),F())
    metadata=forest_metadata()
    top=sorted(({'support_row':i,'positive_mass':str(-y),**metadata[i]}
                for i,y in weights.items()if i in metadata),
               key=lambda r:(-F(r['positive_mass']),r['support_row']))
    target_mass=sum((F(r['positive_mass'])for r in top
                     if sorted(r['types'])==['A1','D5','E7','E8']),F())
    data={'status':'PASS_EXACT_PRIORITY_SUPPORT_DUAL',
          'solver_numeric_objective':float(result.fun),
          'exact_lower_bound':str(bound),'excess_above_40':str(bound-40),
          'proves_support_infeasible_at_40':bound>40,
          'exact_dual_rows':[[i,str(y)]for i,y in sorted(weights.items())],
          'face_mass_shrink':str(scale),'normalizer':str(normalizer),
          'largest_forest_masses':top[:12],
          'A1_D5_E7_E8_forest_mass':str(target_mass),
          'support_rows':2805,'support_columns':6121,'support_nonzeros':17387,
          'solver_warnings':[str(w.message)[:180]for w in caught[:2]],
          'scope':'Exact lower certificate for the stated support successor only; '
                  'an obstruction to this restricted weight family, not to polytopes.',
          'full_model_primal_verified':False,'new_dimension_bound':None}
    if bound<=40:
        data['retained_numeric_dual']=[[i,float(v).hex()]for i,v in enumerate(result.ineqlin.marginals)if v]
    verify_dual(data)
    data['worker_peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin'else 1024)
    assert data['worker_peak_rss_bytes']<LIMIT
    print(json.dumps(data,separators=(',',':')),flush=True)


def guarded_worker(lock_fd):
    before=S.G.resources()
    S.process_tree_rss()
    with Path(S.G.__file__).open('rb')as stream:
        guard_digest=file_digest(stream,'sha256').hexdigest()
    env=os.environ.copy()
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        env[name]='1'
    env['COXETER_SUPPORT_SCOUT_OWNER']=str(os.getpid())
    child=None;output=bytearray();start=time.monotonic();owner=os.getppid();peak=0
    try:
        child=subprocess.Popen([sys.executable,'-B',str(Path(__file__).resolve()),'--worker'],
                               cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                               start_new_session=True,pass_fds=(lock_fd,))
        with selectors.DefaultSelector()as poll:
            poll.register(child.stdout,selectors.EVENT_READ)
            while poll.get_map()or child.poll()is None:
                for key,_ in poll.select(timeout=.1):
                    block=os.read(key.fd,65536)
                    if block:
                        output.extend(block)
                        assert len(output)<OUTPUT_LIMIT
                    else:poll.unregister(key.fileobj)
                peak=max(peak,S.process_tree_rss())
                assert peak<LIMIT and time.monotonic()-start<59
                assert os.getppid()==owner
                S.G.resources(start=False)
        assert child.wait()==0,bytes(output[-4096:]).decode(errors='replace')
        data=json.loads(output)
        verify_dual(data)
        data['guard']={'elapsed_seconds':round(time.monotonic()-start,3),
                       'process_tree_peak_rss_bytes':peak,'rss_limit_bytes':LIMIT,
                       'wall_limit_seconds':59,'output_bytes':len(output),
                       'thread_limit':1,'initial_resources':before}
        data['tooling_provenance']={'guard_module':'scripts/run_lean_batches_20261018.py',
                                   'guard_sha256_at_start':guard_digest,
                                   'role':'Execution provenance; not a mathematical replay dependency'}
        return data
    finally:
        if child is not None:
            try:os.killpg(child.pid,signal.SIGTERM)
            except ProcessLookupError:pass
            try:child.wait(timeout=.3)
            except subprocess.TimeoutExpired:
                try:os.killpg(child.pid,signal.SIGKILL)
                except ProcessLookupError:pass
                child.wait()
            if child.stdout is not None:child.stdout.close()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args()
    if args.worker:
        worker();return
    signal.signal(signal.SIGTERM,S.G.interrupt)
    signal.signal(signal.SIGINT,S.G.interrupt)
    target=ROOT/f'results/{TAG}.json'
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a')as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if args.verify:
            data=json.loads(target.read_text());verify_dual(data)
            for path,digest in data['sources'].items():
                with (ROOT/path).open('rb')as stream:
                    assert file_digest(stream,'sha256').hexdigest()==digest
        else:
            assert not target.exists()
            data=guarded_worker(lock.fileno())
            data['sources']={}
            for path in (f'scripts/{TAG}.py',f'scripts/{S.TAG}.py',
                         'results/fresh_compact21_priority_support_scout_20261018.json',
                         S.M.MODEL,S.SEED):
                with (ROOT/path).open('rb')as stream:
                    data['sources'][path]=file_digest(stream,'sha256').hexdigest()
            encoded=json.dumps(data,indent=2)+'\n'
            assert len(encoded.encode())<OUTPUT_LIMIT
            temporary=target.with_suffix('.json.tmp')
            try:
                with temporary.open('x')as stream:
                    stream.write(encoded);stream.flush();os.fsync(stream.fileno())
                os.replace(temporary,target)
            finally:temporary.unlink(missing_ok=True)
    print(json.dumps({'status':data['status'],'exact_lower_bound':data['exact_lower_bound'],
                      'above_40':data['proves_support_infeasible_at_40'],
                      'dual_terms':len(data['exact_dual_rows']),
                      'A1_D5_E7_E8_mass':data['A1_D5_E7_E8_forest_mass'],
                      'new_dimension_bound':None}))


if __name__=='__main__':main()
