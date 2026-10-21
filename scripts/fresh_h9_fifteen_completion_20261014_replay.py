#!/usr/bin/env python3
"""Read-only fresh arithmetic replay; never uses the resumable audit cache."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fresh_h9_fifteen_completion_20261014_independent import checked,PRIMARY,initialize,replay
def main(cases,jobs):
 d=checked(PRIMARY);rows=d['records'];allids={r['case']for r in rows}
 if cases is not None:
  assert set(cases)<=allids;rows=[r for r in rows if r['case']in set(cases)]
 with ProcessPoolExecutor(max_workers=jobs,initializer=initialize)as pool:
  jobs_done=[pool.submit(replay,r['path'],r['sha256'])for r in rows]
  for k,(ref,job)in enumerate(zip(rows,jobs_done),1):
   out=job.result();assert out['case']==ref['case']and out['verdict']==ref['verdict']and out['states']==ref['nodes']
   if k%25==0 or k==len(rows):print('fresh replay',k,'of',len(rows),flush=True)
 print('PASS_FRESH_READ_ONLY_SEPARATE_ARITHMETIC',len(rows),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--case',type=int,action='append');ap.add_argument('--jobs',type=int,default=2);a=ap.parse_args();assert 1<=a.jobs<=2
 main(a.case,a.jobs)
