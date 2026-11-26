#!/usr/bin/env python3
"""Private SDPA-GMP build replacing repeated graph merges by marked unions.

Only assembly of the symbolic Schur adjacency changes. An edge still means
that two variable coefficient matrices meet the same SDP/SOCP/scalar block.
The arithmetic and the complete solver input stay unchanged.
"""
import argparse
import difflib
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from research.global_formulation import conic,solver_gate

ARCHIVE=Path.home()/'.local/opt/sdpa-gmp/share/sdpa-gmp/source-subset-after-autoreconf.tar.gz'
ARCHIVE_SHA='fec7e13ebf08a52e94a39065cd583cf9d64f64f8d3b804d57adcc46858c134bd'
CACHE=Path('/tmp/coxeter_sdpa_gmp_sparse_graph_20261126')
SOURCE=CACHE/'sdpa-gmp'

FAST=r'''void Chordal::makeGraph(InputData &inputData, int m) {
    // Same union of block-support cliques as upstream; mark duplicates once.
    IVL_init1(adjIVL, IVL_CHUNKED, m);
    std::vector<int> marks(m, -1), row;
    for (int i = 0; i < m; ++i) {
        row.clear();
        const auto append = [&](const int *members, int count) {
            for (int k = 0; k < count; ++k) {
                const int j = members[k];
                if (marks[j] != i) { marks[j] = i; row.push_back(j); }
            }
        };
        const SparseLinearSpace &A = inputData.A[i];
        for (int k = 0; k < A.SDP_sp_nBlock; ++k) {
            const int b = A.SDP_sp_index[k];
            append(inputData.SDP_constraint[b], inputData.SDP_nConstraint[b]);
        }
        for (int k = 0; k < A.SOCP_sp_nBlock; ++k) {
            const int b = A.SOCP_sp_index[k];
            append(inputData.SOCP_constraint[b], inputData.SOCP_nConstraint[b]);
        }
        for (int k = 0; k < A.LP_sp_nBlock; ++k) {
            const int b = A.LP_sp_index[k];
            append(inputData.LP_constraint[b], inputData.LP_nConstraint[b]);
        }
        std::sort(row.begin(), row.end());
        IVL_setList(adjIVL, i, static_cast<int>(row.size()), row.data());
    }
    Graph_init2(graph, 0, m, 0, IVL_tsize(adjIVL), m,
                IVL_tsize(adjIVL), adjIVL, NULL, NULL);
}

'''


def digest(p):return solver_gate.file_digest(p)


def source_text():
    if digest(ARCHIVE)!=ARCHIVE_SHA:raise ValueError('Changed upstream source archive')
    with tarfile.open(ARCHIVE) as f:
        text=f.extractfile('sdpa-gmp/sdpa_chordal.cpp').read().decode()
    begin=text.index('void Chordal::makeGraph(');end=text.index('int Chordal::countNonZero(',begin)
    old=text[begin:end]
    # Preserve every line outside this routine, except the two standard headers.
    changed=text[:begin]+FAST+text[end:]
    changed=changed.replace('#include <sdpa_chordal.h>','#include <sdpa_chordal.h>\n#include <algorithm>\n#include <vector>',1)
    return text,old,changed


def regression(old,original):
    begin=original.index('void Chordal::margeArray(')
    end=original.index('// make aggrigate sparsity pattern',begin)
    merge=original[begin:end]
    pre=r'''
#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
struct IVL { std::vector<std::vector<int> > rows; };
const int IVL_CHUNKED=0;
void IVL_init1(IVL *a,int,int m){a->rows.resize(m);}
void IVL_setList(IVL *a,int i,int n,const int *p){
 if(n)a->rows[i].assign(p,p+n);else a->rows[i].clear();}
int IVL_tsize(IVL *a){int n=0;for(auto &r:a->rows)n+=r.size();return n;}
void Graph_init2(void*,int,int,int,int,int,int,IVL*,void*,void*){}
#define rMessage(x) throw std::runtime_error("reference graph error")
struct SparseLinearSpace {
 int SDP_sp_nBlock,SOCP_sp_nBlock,LP_sp_nBlock;
 int *SDP_sp_index,*SOCP_sp_index,*LP_sp_index;
};
struct InputData {
 int SDP_nBlock,SOCP_nBlock,LP_nBlock;
 int *SDP_nConstraint,*SOCP_nConstraint,*LP_nConstraint;
 int **SDP_constraint,**SOCP_constraint,**LP_constraint;
 SparseLinearSpace *A;
};
namespace ref {struct Chordal {
 IVL store;IVL *adjIVL=&store;void *graph=nullptr;
 void makeGraph(InputData&,int);void margeArray(int,int*,int,int*);
};}
namespace fast {struct Chordal {
 IVL store;IVL *adjIVL=&store;void *graph=nullptr;
 void makeGraph(InputData&,int);
};}
'''
    main=r'''
int main(){
 std::mt19937 generator(20261126);long long checked=0;
 for(int test=0;test<300;++test){
  int m=1+generator()%96;
  std::vector<std::vector<int> > supports[3],incidence[3];
  std::vector<int> counts[3];std::vector<int*> pointers[3];
  for(int g=0;g<3;++g){
   incidence[g].resize(m);int blocks=generator()%90;
   supports[g].resize(blocks);
   for(int b=0;b<blocks;++b){
    for(int j=0;j<m;++j)if(generator()%7==0){
     supports[g][b].push_back(j);incidence[g][j].push_back(b);}
    counts[g].push_back(supports[g][b].size());
    pointers[g].push_back(supports[g][b].data());
   }
  }
  std::vector<SparseLinearSpace> rows(m);
  for(int i=0;i<m;++i)rows[i]={
   (int)incidence[0][i].size(),(int)incidence[1][i].size(),(int)incidence[2][i].size(),
   incidence[0][i].data(),incidence[1][i].data(),incidence[2][i].data()};
  InputData data={(int)supports[0].size(),(int)supports[1].size(),(int)supports[2].size(),
   counts[0].data(),counts[1].data(),counts[2].data(),
   pointers[0].data(),pointers[1].data(),pointers[2].data(),rows.data()};
  ref::Chordal a;fast::Chordal b;a.makeGraph(data,m);b.makeGraph(data,m);
  if(a.store.rows!=b.store.rows)return 2;
  // Also check the defining co-occurrence relation independently.
  for(int i=0;i<m;++i){std::vector<int> expected;
   for(int j=0;j<m;++j){bool meet=false;
    for(int g=0;g<3;++g)for(auto &s:supports[g])
     if(std::binary_search(s.begin(),s.end(),i)&&std::binary_search(s.begin(),s.end(),j))meet=true;
    if(meet)expected.push_back(j);}
   if(expected!=b.store.rows[i])return 3;checked+=expected.size();}
 }
 std::cout<<"300 exact old/new/independent graph cases passed; adjacency entries "<<checked<<"\n";
}
'''
    return pre+merge.replace('Chordal::','ref::Chordal::')+old.replace('Chordal::','ref::Chordal::')+FAST.replace('Chordal::','fast::Chordal::')+main


def build(args):
    original,old,changed=source_text()
    if not SOURCE.exists():
        CACHE.mkdir(exist_ok=True)
        with tarfile.open(ARCHIVE) as f:f.extractall(CACHE,filter='data')
    with tarfile.open(ARCHIVE) as f:
        core={str(Path(q.name).relative_to('sdpa-gmp')):
            hashlib.sha256(f.extractfile(q).read()).hexdigest() for q in f
            if q.isfile() and Path(q.name).suffix in ('.cpp','.h')}
    for name,expected in core.items():
        if name!='sdpa_chordal.cpp' and digest(SOURCE/name)!=expected:
            raise ValueError('Unexpected numerical source change: '+name)
    if (SOURCE/'sdpa_chordal.cpp').read_text() not in (original,changed):
        raise ValueError('Private source has unexpected changes')
    (SOURCE/'sdpa_chordal.cpp').write_text(changed)
    patch=''.join(difflib.unified_diff(original.splitlines(True),changed.splitlines(True),
        fromfile='a/sdpa_chordal.cpp',tofile='b/sdpa_chordal.cpp'))
    (args.output/'graph_assembly.patch').write_text(patch)
    harness=args.output/'graph_regression.cpp';harness.write_text(regression(old,original))
    history=[];started=time.monotonic()
    env=dict(os.environ,MAKEFLAGS='-j1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
        MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    def command(argv,cwd=SOURCE):
        record=dict(command=argv,cwd=str(cwd));history.append(record)
        with (args.output/'build.stdout').open('a') as stream:
            stream.write(json.dumps(record)+'\n');stream.flush()
            done=subprocess.run(argv,cwd=cwd,env=env,stdout=stream,stderr=subprocess.STDOUT,
                timeout=max(1,810-(time.monotonic()-started)),check=False)
        record['returncode']=done.returncode
        if done.returncode:raise RuntimeError('Private build/test command failed')
    report=dict(status='PRIVATE_GMP_GRAPH_BUILD_NOT_COMPLETED',commands=history,
        source_archive_sha256=ARCHIVE_SHA,source_archive=str(ARCHIVE),
        upstream='https://github.com/nakatamaho/sdpa-gmp',
        upstream_commit='ca110db5ea1cc46e811b70dfea9cbb25db74448d',
        only_core_source_changed='sdpa_chordal.cpp',all_numerical_core_sources_unchanged=True,
        graph_patch_sha256=digest(args.output/'graph_assembly.patch'),
        scientific_input_changed=False,installed_solver_untouched=True,
        source_bindings={str(Path(__file__)):dict(sha256=digest(Path(__file__)))})
    try:
        command(['/usr/bin/clang++','-std=c++11','-O2',str(harness),'-o',str(args.output/'graph_regression')])
        command([str(args.output/'graph_regression')])
        report['actual_compiled_old_new_independent_graph_cases_passed']=300
        report['owned_regression_binary_removed_after_use']=dict(sha256=digest(args.output/'graph_regression'),
            bytes=(args.output/'graph_regression').stat().st_size)
        (args.output/'graph_regression').unlink()
        if not (SOURCE/'configure').exists():command(['autoreconf','--force','--install'])
        if not (SOURCE/'Makefile').exists():
            command(['./configure','--with-system-gmp',
                '--with-gmp-includedir=/opt/homebrew/opt/gmp/include',
                '--with-gmp-libdir=/opt/homebrew/opt/gmp/lib',
                'CC=/usr/bin/clang','CXX=/usr/bin/clang++'])
        command(['make','-j1','config'],SOURCE/'external/spooles')
        compatibility=SOURCE/'external/spooles/work/internal/ETree/src/transform.c'
        before=compatibility.read_text();after=before.replace('IVinit(nfront, NULL)','IVinit(nfront, 0)')
        # Whitespace in the upstream file may differ from the installed patch.
        import re
        after,n=re.subn(r'IVinit\(nfront,\s*NULL\)','IVinit(nfront, 0)',before)
        if n not in (0,3):raise ValueError('Unexpected SPOOLES integer compatibility patch')
        compatibility.write_text(after)
        makeinc=compatibility.parents[2]/'Make.inc'
        text=makeinc.read_text()
        text=re.sub(r'^CFLAGS\s*=.*$', 'CFLAGS = -O2 -std=gnu89', text, flags=re.M)
        makeinc.write_text(text)
        command(['make','-j1'])
        for name,expected in core.items():
            if name!='sdpa_chordal.cpp' and digest(SOURCE/name)!=expected:
                raise ValueError('A numerical source changed during compilation: '+name)
        if (SOURCE/'sdpa_chordal.cpp').read_text()!=changed:raise ValueError('Changed graph patch')
        report['all_unchanged_numerical_source_sha256']={n:h for n,h in core.items() if n!='sdpa_chordal.cpp'}
        report['upstream_chordal_source_sha256']=core['sdpa_chordal.cpp']
        report['patched_chordal_source_sha256']=digest(SOURCE/'sdpa_chordal.cpp')
        binary=SOURCE/'sdpa_gmp'
        if not binary.is_file():raise RuntimeError('No private GMP executable')
        receipt=dict(status='PRIVATE_GRAPH_ASSEMBLY_BUILD',binary_sha256=digest(binary),
            upstream=report['upstream'],commit=report['upstream_commit'],version='7.1.3 + marked graph union',
            graph_patch_sha256=report['graph_patch_sha256'])
        conic.write_json(args.output/'installation.json',receipt)
        report.update(status='PRIVATE_GMP_GRAPH_BUILD_AND_COMPILED_EQUIVALENCE_CHECKS_PASSED',
            executable=str(binary),executable_sha256=digest(binary),
            installation=str(args.output/'installation.json'))
    except (RuntimeError,subprocess.TimeoutExpired) as error:
        report.update(status='PRIVATE_GMP_GRAPH_BUILD_BATCH_INCOMPLETE',reason=str(error))
    conic.write_json(args.output/'result.json',report)


def backend(args):
    from scripts import run_balanced_dual_gmp_20261126 as gmp
    report=json.loads((args.build/'result.json').read_text())
    if report['status']!='PRIVATE_GMP_GRAPH_BUILD_AND_COMPILED_EQUIVALENCE_CHECKS_PASSED':
        raise ValueError('No complete private graph-build proof')
    for p,row in report['source_bindings'].items():
        if digest(p)!=row['sha256']:raise ValueError('Changed graph-build source')
    if digest(report['executable'])!=report['executable_sha256']:
        raise ValueError('Changed private GMP executable')
    # Configure only the backend executable/receipt, preserving every input gate.
    gmp.EXECUTABLE=Path(report['executable'])
    gmp.EXECUTABLE_SHA=report['executable_sha256']
    gmp.INSTALLATION=Path(report['installation'])
    if args.action=='check':gmp.check(args)
    else:gmp.solve(args)
    p=args.output/'result.json';result=json.loads(p.read_text())
    result['private_backend_build']=str(args.build/'result.json')
    result['private_backend_launcher_sha256']=digest(Path(__file__))
    p.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['build','check','solve'])
    parser.add_argument('--build',type=Path)
    parser.add_argument('--prepared',type=Path)
    parser.add_argument('--verified',type=Path)
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--iterations',type=int,default=200)
    parser.add_argument('--seconds',type=int,default=780)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the serial scientific watchdog')
    if args.precision not in (128,256,512) or not 1<=args.seconds<=800:
        raise ValueError('Unsupported precision or batch duration')
    args.output=args.output.resolve()
    if args.output.exists() or not args.output.is_relative_to(ROOT/'results'):raise ValueError('Use a fresh results folder')
    args.output.mkdir()
    if args.action=='build':build(args)
    else:backend(args)
    print(json.dumps(dict(status=json.loads((args.output/'result.json').read_text())['status'],output=str(args.output))),flush=True)


if __name__=='__main__':main()
