#!/usr/bin/env python3
"""One guarded sparse support LP; exact replay never reruns the optimizer.

This deliberately does NOT verify all full-model face/forest constraints.
"""
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
import fresh_compact21_extreme_e8_priority_saved_dual_20261018 as M
import run_lean_batches_20261018 as G

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_compact21_priority_support_scout_20261018'
SEED='results/fresh_compact21_extreme_e8_priority_saved_dual_20261018.json'
LIMIT=1024**3
OUTPUT_LIMIT=512*1024


def rows():
    seed=json.loads((ROOT/SEED).read_text())
    beta,target=seed['extreme_E8_bonus_column'],seed['target_forest_term']
    for index,term in M.objects('terms',4):
        row={j:F(v)for j,v in term['coefficients']}
        rhs=F(term['rhs'])
        kind=term['metadata']['kind']
        if index==target:
            assert kind=='forest' and row[beta]==2
            yield row,rhs+1,kind
            inactive=row.copy()
            del inactive[beta]
            yield inactive,rhs,kind
        else:
            yield row,rhs,kind


def verify_vector(data):
    seed=json.loads((ROOT/SEED).read_text())
    z=seed['variables']
    values={int(j):F(v)for j,v in data['positive_bonus_values']}
    assert len(values)==len(data['positive_bonus_values'])
    assert all(0<=j<z and v>0 for j,v in values.items())
    values[z]=F(data['upper_coordinate'])
    assert values[z]>=0
    count=0
    for row,rhs,_ in rows():
        assert sum((a*values.get(j,F())for j,a in row.items()),F())<=rhs
        count+=1
    assert count==2805
    return count


def worker():
    owner=int(os.environ['COXETER_SUPPORT_SCOUT_OWNER'])
    assert owner==os.getppid()
    def owner_watch():
        while True:
            if os.getppid()!=owner:
                os.killpg(os.getpgrp(),signal.SIGKILL)
            time.sleep(.2)
    threading.Thread(target=owner_watch,daemon=True).start()
    # Only the guarded child imports numerical packages.
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    saved=list(rows())
    columns=sorted({j for row,_,_ in saved for j in row})
    assert len(saved)==2805 and len(columns)==6121 and columns[-1]==93045
    index={j:i for i,j in enumerate(columns)}
    rr,cc,vv=[],[],[]
    for i,(row,_,_)in enumerate(saved):
        for j,a in row.items():
            rr.append(i);cc.append(index[j]);vv.append(float(a))
    assert len(vv)==17387
    matrix=coo_matrix((vv,(rr,cc)),shape=(len(saved),len(columns))).tocsr()
    objective=np.zeros(len(columns))
    objective[-1]=1
    with warnings.catch_warnings(record=True)as caught:
        result=linprog(objective,A_ub=matrix,
                       b_ub=np.asarray([float(rhs)for _,rhs,_ in saved]),
                       bounds=(0,None),method='highs',
                       options={'time_limit':45.0,'threads':1,'disp':False})
    assert result.success, f'Support solver status {result.status}: {result.message}'
    grid=10**8
    values={j:F(math.ceil(max(0.0,float(result.x[i]))*grid),grid)
            for i,j in enumerate(columns[:-1]) if result.x[i]>0}
    scale=F(1)
    for row,rhs,kind in saved:
        if kind=='face' and rhs<0:
            payment=-sum((a*values.get(j,F())for j,a in row.items()),F())
            assert payment>0
            scale=max(scale,-rhs/payment)
    values={j:v*scale for j,v in values.items()if v}
    upper=F()
    for row,rhs,kind in saved:
        if kind=='forest':
            assert row[93045]==-1
            upper=max(upper,sum((a*values.get(j,F())for j,a in row.items()if j!=93045),F())-rhs)
    data={'status':'PASS_EXACT_PRIMAL_ON_PRIORITY_SUPPORT_ONLY',
          'support_rows':2805,'support_columns':6121,'support_nonzeros':17387,
          'solver_numeric_upper':float(result.fun),'upper_coordinate':str(upper),
          'upper_at_most_40':upper<=40,'face_scaling':str(scale),
          'positive_bonus_values':[[j,str(v)]for j,v in sorted(values.items())],
          'support_exactly_verified':True,'full_model_verified':False,
          'solver_warnings':[str(w.message)[:180]for w in caught[:2]],
          'scope':'Only the old positive-dual support plus the new inactive row. '
                  'Full face separation and all-forest verification remain outstanding.',
          'new_dimension_bound':None}
    verify_vector(data)
    data['worker_peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
    assert data['worker_peak_rss_bytes']<LIMIT
    print(json.dumps(data,separators=(',',':')),flush=True)


def process_tree_rss():
    output=subprocess.check_output(['ps','-axo','pid=,ppid=,rss='],text=True,timeout=2)
    table=[tuple(map(int,line.split()))for line in output.splitlines()if line.strip()]
    descendants={os.getpid()}
    while True:
        new={pid for pid,parent,_ in table if parent in descendants}-descendants
        if not new:
            break
        descendants|=new
    return sum(rss*1024 for pid,_,rss in table if pid in descendants)


def guarded_worker(lock_fd):
    before=G.resources()
    # Fail before spawning if the environment cannot provide the watchdog.
    process_tree_rss()
    env=os.environ.copy()
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        env[name]='1'
    env['COXETER_SUPPORT_SCOUT_OWNER']=str(os.getpid())
    child=None
    output=bytearray()
    start=time.monotonic()
    owner=os.getppid()
    peak=0
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
                        assert len(output)<OUTPUT_LIMIT,'Scout output cap exceeded'
                    else:
                        poll.unregister(key.fileobj)
                peak=max(peak,process_tree_rss())
                assert peak<LIMIT,'Scout process-tree RSS exceeded one GiB'
                assert time.monotonic()-start<59,'Scout exceeded59 seconds'
                assert os.getppid()==owner,'Scout owner exited'
                G.resources(start=False)
        assert child.wait()==0,bytes(output[-4096:]).decode(errors='replace')
        data=json.loads(output)
        verify_vector(data)
        data['guard']={'elapsed_seconds':round(time.monotonic()-start,3),
                       'process_tree_peak_rss_bytes':peak,'rss_limit_bytes':LIMIT,
                       'wall_limit_seconds':59,'output_bytes':len(output),
                       'thread_limit':1,'initial_resources':before}
        return data
    finally:
        if child is not None:
            try:
                os.killpg(child.pid,signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                child.wait(timeout=.3)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(child.pid,signal.SIGKILL)
                except ProcessLookupError:
                    pass
                child.wait()
            if child.stdout is not None:
                child.stdout.close()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args()
    if args.worker:
        worker()
        return
    signal.signal(signal.SIGTERM,G.interrupt)
    signal.signal(signal.SIGINT,G.interrupt)
    target=ROOT/f'results/{TAG}.json'
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a')as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if args.verify:
            data=json.loads(target.read_text())
            verify_vector(data)
            for path,digest in data['sources'].items():
                with (ROOT/path).open('rb')as stream:
                    assert file_digest(stream,'sha256').hexdigest()==digest
        else:
            assert not target.exists()
            data=guarded_worker(lock.fileno())
            data['sources']={}
            for path in (f'scripts/{TAG}.py',M.MODEL,SEED,
                         'scripts/fresh_compact21_extreme_e8_priority_saved_dual_20261018.py',
                         'scripts/run_lean_batches_20261018.py'):
                with (ROOT/path).open('rb')as stream:
                    data['sources'][path]=file_digest(stream,'sha256').hexdigest()
            encoded=json.dumps(data,indent=2)+'\n'
            assert len(encoded.encode())<OUTPUT_LIMIT
            temporary=target.with_suffix('.json.tmp')
            try:
                with temporary.open('x')as stream:
                    stream.write(encoded)
                    stream.flush();os.fsync(stream.fileno())
                os.replace(temporary,target)
            finally:
                temporary.unlink(missing_ok=True)
    print(json.dumps({'status':data['status'],'exact_support_upper':data['upper_coordinate'],
                      'at_most_40':data['upper_at_most_40'],
                      'positive_bonuses':len(data['positive_bonus_values']),
                      'full_model_verified':False,'new_dimension_bound':None}))


if __name__=='__main__':
    main()
