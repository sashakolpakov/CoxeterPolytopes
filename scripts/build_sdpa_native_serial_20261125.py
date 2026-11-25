#!/usr/bin/env python3
"""Build the isolated standard SDPA command line serially for warm starts.

Compilation uses the scientific watchdog and does not run a numerical test.
It never installs over a system solver or changes retained model inputs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
CACHE=Path('/tmp/coxeter_sdpa_native_20261125')
SOURCE=CACHE/'sdpa-7.3.21'
SDPA_TAR=Path('/tmp/coxeter_sdpa_binary_20261125/sdpa_7.3.21.tar.gz')
MUMPS_TAR=CACHE/'mumps_5.9.0.orig.tar.gz'
SDPA_SHA='300cae3ac04409c2243f6a5b3fd1ee57b2b18db209a21c2472e8d837951d3db6'
MUMPS_SHA='02c6efdb91749ec0f82351d40f3f860547272a1eb1d899126a4265b4d6bcc4ca'
PREFIX=Path.home()/'.local/opt/coxeter-sdpa-7.3.21'


def digest(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the serial scientific watchdog')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    parser.add_argument('--resume-mumps',action='store_true')
    args=parser.parse_args();out=(ROOT/args.output).resolve()
    if not out.is_relative_to(ROOT/'results') or out.exists():raise ValueError('Use a fresh results directory')
    if digest(SDPA_TAR)!=SDPA_SHA or digest(MUMPS_TAR)!=MUMPS_SHA:raise ValueError('Changed upstream archives')
    if digest(SOURCE/'mumps'/MUMPS_TAR.name)!=MUMPS_SHA:raise ValueError('Changed local MUMPS input')
    if 'MUMPS_VER = 5.9.0' not in (SOURCE/'mumps/Makefile').read_text():raise ValueError('Wrong MUMPS source version')
    out.mkdir();started=time.monotonic();history=[]
    sources=[Path(__file__),*sorted(SOURCE.glob('*.cpp')),*sorted(SOURCE.glob('*.h')),SOURCE/'mumps/Makefile']
    bindings={str(p):dict(bytes=p.stat().st_size,sha256=digest(p)) for p in sources}
    (out/'source_bindings.json').write_text(json.dumps(bindings,indent=2,sort_keys=True)+'\n')
    environment=dict(os.environ,MAKEFLAGS='-j1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
        MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    def command(argv,cwd):
        timeout=max(1,int(760-(time.monotonic()-started)))
        entry=dict(command=argv,cwd=str(cwd),timeout_seconds=timeout)
        history.append(entry)
        with (out/'build.stdout').open('a') as stream:
            stream.write(json.dumps(entry)+'\n');stream.flush()
            try:
                result=subprocess.run(argv,cwd=cwd,env=environment,stdout=stream,stderr=subprocess.STDOUT,
                    timeout=timeout,check=False)
                entry['returncode']=result.returncode
                if result.returncode:raise RuntimeError('Compiler/configure command failed')
            except subprocess.TimeoutExpired:
                entry['timed_out']=True;raise
    report=dict(status='NO_NATIVE_BUILD_RESULT',serial_compilation=True,numerical_test_run=False,
        upstream_archives=[dict(path=str(SDPA_TAR),sha256=SDPA_SHA,
            source='https://downloads.sourceforge.net/project/sdpa/sdpa/sdpa_7.3.21.tar.gz'),
            dict(path=str(MUMPS_TAR),sha256=MUMPS_SHA,
            source='https://deb.debian.org/debian/pool/main/m/mumps/mumps_5.9.0.orig.tar.gz')],
        configure_overrides=['MUMPS_VER = 5.9.0'],commands=history)
    try:
        if not (SOURCE/'Makefile').exists():
            command(['./configure','CC=/usr/bin/clang','CXX=/usr/bin/clang++',
                'FC=/opt/homebrew/bin/gfortran','CFLAGS=-O2 -DNDEBUG -fPIC',
                'CXXFLAGS=-O2 -DNDEBUG -fPIC','FCFLAGS=-O2 -fPIC',
                '--with-blas=-L/opt/homebrew/opt/openblas/lib -lopenblas',
                '--with-lapack=-L/opt/homebrew/opt/openblas/lib -lopenblas'],SOURCE)
        if args.resume_mumps and (SOURCE/'mumps/build/Makefile.inc').exists():
            command(['make','-j1','d'],SOURCE/'mumps/build')
        command(['make','-j1'],SOURCE)
        executable=SOURCE/'sdpa'
        if not executable.is_file():raise RuntimeError('No native executable after successful make')
        target=PREFIX/'bin/sdpa'
        if target.exists() and digest(target)!=digest(executable):raise RuntimeError('A different private native solver already exists')
        target.parent.mkdir(parents=True,exist_ok=True)
        if not target.exists():shutil.copy2(executable,target)
        if digest(target)!=digest(executable):raise RuntimeError('The installed executable failed byte verification')
        report.update(status='NATIVE_SDPA_SERIAL_BUILD_COMPLETED_UNTESTED',executable=str(target),
            executable_bytes=target.stat().st_size,executable_sha256=digest(target))
    except subprocess.TimeoutExpired:report['status']='BUILD_BATCH_TIME_LIMIT_NO_SOLVER_RESULT'
    except RuntimeError as error:report.update(status='NATIVE_BUILD_FAILED',reason=str(error))
    if any(digest(Path(p))!=item['sha256'] for p,item in bindings.items()):raise ValueError('An original source changed during compilation')
    report['elapsed_seconds']=time.monotonic()-started
    (out/'result.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
