#!/usr/bin/env python3
'''Render hash-pinned local documentation under the bounded ownership wrapper.

Use run_bounded_research_check_20261018.py --parent-lock. Supply one --sha256
per page, in page order. This is visual inspection, never mathematics.
Browser groups have a parent-loss bridge; profiles/logs are temporary.
Only top documents are hash-pinned, not every browser subresource.
'''
import argparse
from hashlib import file_digest
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import tempfile
import time

import run_lean_batches_signal_safe_20261018 as safe

ROOT = Path(__file__).resolve().parents[1]
CHROME = Path('/Applications/Google Chrome.app/Contents/MacOS/Google Chrome')
BROWSER_BRIDGE = r'''
import os, signal, subprocess, sys, threading, time
owner = int(os.environ['COXETER_BROWSER_OWNER'])
def stop_group(*unused):
    os.killpg(os.getpgrp(), signal.SIGKILL)
def watch_owner():
    while True:
        if os.getppid() != owner:
            stop_group()
        time.sleep(.1)
threading.Thread(target=watch_owner, daemon=True).start()
signal.signal(signal.SIGALRM, stop_group)
signal.alarm(16)
child = subprocess.Popen(sys.argv[1:])
raise SystemExit(child.wait())
'''


def interrupt(signum, unused):
    raise RuntimeError(f'Documentation render interrupted by signal{signum}')


def no_symlinks(path):
    current = ROOT
    for part in path.relative_to(ROOT).parts:
        current /= part
        if current.is_symlink():
            raise ValueError('symlink input/output ancestor is not allowed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pages', nargs='+')
    parser.add_argument('--sha256', action='append', required=True)
    parser.add_argument('--output-name', required=True)
    parser.add_argument('--width', type=int, choices=(480, 960, 1280, 1440), default=1280)
    parser.add_argument('--height', type=int, default=1500)
    args = parser.parse_args()
    owner = os.environ.get('COXETER_CHECK_OWNER','')
    if not owner.isdecimal() or int(owner) != os.getppid():
        parser.error('Run under the bounded ownership wrapper with --parent-lock')
    if not (1 <= len(args.pages) <= 3 and len(args.sha256) == len(args.pages)
            and 480 <= args.height <= 2400
            and re.fullmatch(r'[a-z0-9_]{1,80}', args.output_name)
            and all(re.fullmatch(r'[0-9a-f]{64}', value) for value in args.sha256)):
        parser.error('Invalid page selection, hashes, viewport or output name')
    if not CHROME.is_file():
        parser.error('The configured Chrome executable is unavailable')
    pages = []
    for name,expected in zip(args.pages,args.sha256):
        lexical = Path(os.path.abspath(ROOT/name))
        if not lexical.is_relative_to(ROOT):
            parser.error('Input escapes repository')
        no_symlinks(lexical)
        page = lexical.resolve()
        if not (page.is_file() and 0 < page.stat().st_size < 1024**2
                and (page.is_relative_to(ROOT/'docs/site/_static') and page.suffix == '.svg'
                     or page.is_relative_to(ROOT/'.publication') and page.suffix == '.html')):
            parser.error('Choose an owned local SVG or built HTML below1MiB')
        with page.open('rb') as stream:
            if file_digest(stream,'sha256').hexdigest() != expected:
                parser.error('Input does not match its expected SHA256')
        pages.append((page,expected))
    out = ROOT/'formal/.proof-batches'/args.output_name
    no_symlinks(out)
    out.mkdir(exist_ok=False)
    for signum in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP,signal.SIGALRM):
        signal.signal(signum,interrupt)
    signal.alarm(50)
    rows = []
    try:
        with tempfile.TemporaryDirectory(prefix='coxeter-doc-preview-',dir='/private/tmp') as temp:
            for index,(page,expected) in enumerate(pages,1):
                target = out/f'page{index}.png'
                profile = Path(temp)/f'profile{index}'
                command = [str(CHROME),'--headless','--disable-gpu','--hide-scrollbars',
                           '--no-first-run','--no-default-browser-check',
                           '--disable-background-networking','--disable-extensions',
                           '--disable-component-update','--disable-sync',
                           '--renderer-process-limit=1',
                           '--js-flags=--max-old-space-size=128',
                           '--disable-breakpad','--disable-crash-reporter',
                           '--disable-component-extensions-with-background-pages',
                           '--disk-cache-size=8388608','--media-cache-size=1048576',
                           '--user-data-dir='+str(profile),
                           f'--window-size={args.width},{args.height}',
                           '--virtual-time-budget=2500','--screenshot='+str(target),
                           page.as_uri()]
                with tempfile.TemporaryFile(dir=temp) as log:
                    env = dict(os.environ,COXETER_BROWSER_OWNER=str(os.getpid()))
                    child = subprocess.Popen([sys.executable,'-B','-c',BROWSER_BRIDGE,*command],
                        cwd=ROOT,env=env,start_new_session=True,
                        stdout=log,stderr=subprocess.STDOUT)
                    started = time.monotonic()
                    try:
                        while child.poll() is None:
                            if os.fstat(log.fileno()).st_size > 65536:
                                raise RuntimeError('Browser log exceeds64KiB sampled bound')
                            if time.monotonic()-started > 14:
                                raise RuntimeError('Browser exceeds14-second limit')
                            if os.getppid() != int(owner):
                                raise RuntimeError('Bounded owner exited')
                            time.sleep(.05)
                        size = os.fstat(log.fileno()).st_size
                        if size > 65536:
                            raise RuntimeError('Browser log exceeds64KiB sampled bound')
                        if child.returncode != 0:
                            log.seek(max(0,size-3000))
                            raise RuntimeError(log.read(3000).decode(errors='replace'))
                    finally:
                        safe.stop_group(child,.3)
                if not target.is_file() or not 0 < target.stat().st_size < 8*1024**2:
                    raise RuntimeError('Missing or oversized screenshot')
                with page.open('rb') as stream:
                    if file_digest(stream,'sha256').hexdigest() != expected:
                        raise RuntimeError('Input changed during rendering')
                rows.append({'source':str(page.relative_to(ROOT)),'source_sha256':expected,
                             'png':str(target.relative_to(ROOT)),'bytes':target.stat().st_size})
    finally:
        signal.alarm(0)
    print(json.dumps({'status':'RENDERED_FOR_VISUAL_REVIEW_NOT_PROOF',
        'viewport':[args.width,args.height],'pages':rows,
        'isolated_temporary_profiles_removed':True,
        'mathematical_verification_claimed':False}))


if __name__ == '__main__':
    main()
