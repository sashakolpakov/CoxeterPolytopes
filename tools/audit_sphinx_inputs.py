#!/usr/bin/env python3
"""Audit explicit Sphinx inputs before a documentation build, never mathematics.

Checks local downloads, images, figures, recursive includes, all source pages,
configuration/templates and static assets. Rejects symlinks, including ancestors.
Rejects manifest-listed unuploaded/third-party payloads by path or
exact byte hash, including aliases. This is not analysis of arbitrary Python,
JavaScript, Sphinx extensions, or compressed/modified copies of unknown payloads.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT/'docs/site'
TICK = chr(96)
ROLES = re.compile(r':download:'+TICK+'([^'+TICK+r'\n]+)'+TICK)
DIRECTIVES = re.compile(
    r'^[ \t]*\.\.[ \t]+(?:\|[^|\n]+\|[ \t]+)?'
    r'(image|figure|include|literalinclude)::[ \t]+([^\n]+?)[ \t]*$', re.M)
FORBIDDEN_PARTS = {'.git','.lake','.proof-batches','.publication','__pycache__','node_modules'}
MAX_INPUT, MAX_TEXT, MAX_TOTAL, MAX_FILES = 16*1024**2, 4*1024**2, 128*1024**2, 4096


def canonical_path(path):
    lexical = Path(path)
    if not lexical.is_absolute():
        lexical = ROOT/lexical
    if not lexical.is_relative_to(ROOT):
        raise ValueError('local reference escapes repository')
    current = ROOT
    for part in lexical.relative_to(ROOT).parts:
        current = current.parent if part == '..' else current/part
        if not current.is_relative_to(ROOT):
            raise ValueError('local reference escapes repository')
        if current.is_symlink():
            raise ValueError('symlink input or ancestor is not allowed')
    target = current.resolve()
    if not target.is_relative_to(ROOT):
        raise ValueError('local reference escapes repository')
    return target


def references(text):
    for body in ROLES.findall(text):
        title = re.fullmatch(r'.*<([^<>\n]+)>',body)
        yield 'download', title[1] if title else body
    yield from DIRECTIVES.findall(text)


def local_target(source, raw):
    parsed = urlsplit(raw.strip())
    if parsed.scheme or parsed.netloc:
        return None
    return canonical_path(source.parent/unquote(parsed.path))


def path_forbidden(relative, blocked):
    return relative in blocked or bool(set(Path(relative).parts) & FORBIDDEN_PARTS)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--require-indexed', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if not __debug__:
        parser.error('Do not disable assertions for the synthetic self-test')
    bulk = json.loads((ROOT/'repro/public-data-manifest.json').read_text())
    external = json.loads((ROOT/'repro/external-inputs-manifest.json').read_text())
    rows = [row for part in bulk['parts'] for row in part['files']] + external['files']
    blocked = {row['path'] for row in rows}
    blocked_hashes = {row['sha256'] for row in [*rows,*bulk['parts']]}
    if args.self_test:
        assert path_forbidden(next(iter(blocked)), blocked)
        assert path_forbidden('formal/.proof-batches/private.olean', blocked)
        assert not path_forbidden('docs/site/_static/finite-payment-flow.svg', blocked)
        assert local_target(SITE/'index.rst', 'https://example.org/file.pdf') is None
        assert local_target(SITE/'index.rst', '_static/test%20figure.svg') == SITE/'_static/test figure.svg'
        try:
            local_target(SITE/'index.rst', '../../../../outside')
        except ValueError:
            pass
        else:
            raise AssertionError('escaping path accepted')
        sample = (':download:'+TICK+'x <../../data/private.json>'+TICK+'\n'
                  ':download:'+TICK+'../../data/private.json'+TICK+'\n'
                  '.. figure:: _static/a figure.svg\n'
                  '.. |a| image:: _static/a.svg\n')
        assert len(list(references(sample))) == 4
        assert rows[0]['sha256'] in blocked_hashes  # aliases are checked by digest below
        with tempfile.TemporaryDirectory(prefix='.sphinx-input-check-',dir=ROOT) as temporary:
            fixture = Path(temporary)
            (fixture/'target').mkdir()
            (fixture/'link').symlink_to(fixture/'target',target_is_directory=True)
            try:
                canonical_path(fixture/'link'/'..'/'file')
            except ValueError as error:
                assert 'symlink' in str(error)
            else:
                raise AssertionError('symlink followed by parent traversal accepted')
        print(json.dumps({'publication_input_self_test':'PASS','checks':9,'mathematical_checks':False}))
        return
    indexed = set()
    if args.require_indexed:
        data = subprocess.check_output(['git','ls-files','-z'],cwd=ROOT,timeout=10)
        if len(data) > 8*1024**2:
            raise ValueError('Git index path listing exceeds8MiB bound')
        indexed = set(data.decode().split('\0'))
    refs, issues, targets, sources = [], [], set(), []
    targets.update([ROOT/'repro/public-data-manifest.json',
                    ROOT/'repro/external-inputs-manifest.json'])
    canonical_path(SITE)
    for directory, names, files in os.walk(SITE,followlinks=False):
        directory = Path(directory)
        if directory == SITE:
            names[:] = [name for name in names if name != '_build']
        for name in [*names,*files]:
            node = directory/name
            if node.is_symlink():
                issues.append({'path':str(node.relative_to(ROOT)),'reason':'symlink source/static entry'})
        names[:] = [name for name in names if not (directory/name).is_symlink()]
        for name in files:
            node = directory/name
            if node.is_symlink():
                continue
            target = canonical_path(node)
            if not target.is_file():
                issues.append({'path':str(target.relative_to(ROOT)),'reason':'non-regular site input'})
                continue
            if target.suffix == '.rst':
                sources.append(target)
            if (target.is_relative_to(SITE/'_static')
                    or target.is_relative_to(SITE/'_templates') or target == SITE/'conf.py'):
                targets.add(target)
    pending, visited = sorted(sources), set()
    while pending:
        source = pending.pop()
        if source in visited:
            continue
        visited.add(source)
        if len(visited) > MAX_FILES:
            issues.append({'reason':'too many included text sources'})
            break
        targets.add(source)
        relative = str(source.relative_to(ROOT))
        if path_forbidden(relative,blocked) or source.stat().st_size > MAX_TEXT:
            issues.append({'path':relative,'reason':'blocked or oversized text source'})
            continue
        text = source.read_text()
        for kind,raw in references(text):
            try:
                target = local_target(source,raw)
                if target is None:
                    continue
                relative = str(target.relative_to(ROOT))
                if path_forbidden(relative,blocked):
                    raise ValueError('unuploaded bulk/excluded/cache path')
                if not target.is_file():
                    raise ValueError('missing local file')
                targets.add(target)
                refs.append({'source':str(source.relative_to(ROOT)),'target':relative})
                if kind == 'include':
                    pending.append(target)
            except ValueError as error:
                issues.append({'source':str(source.relative_to(ROOT)),'reference':raw,'reason':str(error)})
    bindings, total = [], 0
    if len(targets) > MAX_FILES:
        issues.append({'reason':'too many documentation input files'})
    for target in sorted(targets)[:MAX_FILES]:
        target = canonical_path(target)
        relative = str(target.relative_to(ROOT))
        if path_forbidden(relative,blocked):
            issues.append({'path':relative,'reason':'forbidden static/input path'})
            continue
        size = target.stat().st_size
        total += size
        if size > MAX_INPUT or total > MAX_TOTAL:
            issues.append({'path':relative,'reason':'documentation input exceeds per-file/aggregate bound'})
            continue
        with target.open('rb') as stream:
            digest = hashlib.file_digest(stream,'sha256').hexdigest()
        if digest in blocked_hashes:
            issues.append({'path':relative,'reason':'exact byte alias of excluded/unuploaded payload'})
        if args.require_indexed and relative not in indexed:
            issues.append({'path':relative,'reason':'not in actual Git index'})
        bindings.append({'path':relative,'bytes':size,'sha256':digest})
    print(json.dumps({'publication_input_audit_only':True,'source_pages':len(sources),
        'text_sources_scanned':len(visited),'explicit_local_references':refs,'input_bindings':bindings,
        'blocked_manifest_paths':len(blocked),'require_indexed':args.require_indexed,
        'issues':issues,'mathematical_checks_launched':False},indent=2))
    if issues:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
