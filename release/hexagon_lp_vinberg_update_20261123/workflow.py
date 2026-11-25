"""Update the existing LP work with its explicitly authorized Vinberg remark.

Called only by the allowlisted publication client under its watchdog.
Original v1 payload descriptors and receipts are never rewritten.
"""
from __future__ import annotations

import argparse
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LEGACY = ROOT/'release/hexagon_submission_schedule_20261123'
PAPER = 'linear_programming_dimension_bounds'
IDENTIFIER = '2611.00962'
SOURCE_SHA = '01b716fcac5e98c7e57282507af9614790a454d2a33affdbc418b9f93eb7a360'
PDF_SHA = '9d3473c36dab9d32ad8258628bf28068c6813912e5c06c0b4d48de9c4dba9bba'
BUILD = 'release/linear_programming_dimension_bounds/build_vinberg_remark01_20261123.json'
PACKAGE = 'release/linear_programming_dimension_bounds/package_vinberg_remark01_20261123.json'
AUTHORIZATION = 'User explicitly requested a side agent to update the existing Hexagon LP submission with the canonical concise Vinberg E7 remark; bounds, title, author, abstract and bibliography remain unchanged.'


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def save_private(path, raw, api):
    if path.exists():
        if path.is_symlink() or path.read_bytes() != raw:
            raise api.Stop('Existing private update evidence differs')
        return
    fd = os.open(path, os.O_WRONLY|os.O_CREAT|os.O_EXCL|os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, 'wb') as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())


def authorized_inputs(api):
    old = read(LEGACY/'approved_payloads.json')
    previous = next(row for row in old['papers'] if row['id'] == PAPER)
    old_receipt_path = LEGACY/(PAPER+'.submission-receipt.local')
    old_receipt = read(old_receipt_path)
    server = old_receipt.get('response', {}).get('receipt', {})
    if old_receipt.get('http_status') != 202 or server.get('identifier') != IDENTIFIER:
        raise api.Stop('Existing LP accepted receipt does not match the target work')
    metadata_path = LEGACY/previous['metadata_template']
    metadata_raw = metadata_path.read_bytes()
    metadata = json.loads(metadata_raw)
    if sha(api.canonical(metadata)) != old_receipt['metadata_sha256']:
        raise api.Stop('Existing LP metadata changed')
    build_raw, package_raw = (ROOT/BUILD).read_bytes(), (ROOT/PACKAGE).read_bytes()
    build, package = json.loads(build_raw), json.loads(package_raw)
    if (build.get('status') != 'PASS_ISOLATED_TEX_BUILD_ONLY' or build.get('pages') != 36
            or build.get('source_sha256') != SOURCE_SHA or build.get('pdf_sha256') != PDF_SHA
            or package.get('status') != 'PASS_LOCAL_THREE_FILE_SOURCE_PACKAGE'
            or package.get('pages') != 36 or package.get('build_receipt', {}).get('sha256') != sha(build_raw)):
        raise api.Stop('Authorized local build/package receipt differs')
    expected = {row['path']: row for row in package['members']}
    rows, files = [], {}
    for prior in previous['upload_files']:
        relative = prior['path']
        record = expected.get(relative)
        raw = (ROOT/relative).read_bytes()
        if not record or len(raw) != record['bytes'] or sha(raw) != record['sha256']:
            raise api.Stop('Canonical update source differs from its exact package receipt')
        if prior['relativePath'] != 'main.tex' and sha(raw) != prior['sha256']:
            raise api.Stop('Bibliography or compiled bibliography changed')
        row = dict(prior, bytes=len(raw), sha256=sha(raw))
        rows.append(row)
        files[prior['relativePath']] = raw
    if files.keys() != {'main.tex', 'references.bib', 'main.bbl'} or sha(files['main.tex']) != SOURCE_SHA:
        raise api.Stop('Unexpected update source inventory')
    pdf_path = ROOT/previous['reference_pdf']['path']
    pdf_raw = pdf_path.read_bytes()
    if sha(pdf_raw) != PDF_SHA or len(pdf_raw) != build['pdf_bytes']:
        raise api.Stop('Canonical reference PDF differs from the authorized edition')
    pattern = rb'(?ms)^\\begin\{remark\}\[Choice of the weight family\]\\label\{fv:vinberg-weight-obstruction\}\n.*?^\\end\{remark\}\n'
    matches = list(re.finditer(pattern, files['main.tex']))
    if len(matches) != 1:
        raise api.Stop('Expected one concise Vinberg remark')
    match = matches[0]
    candidates = [files['main.tex'][:match.start()]+files['main.tex'][match.end():]]
    if files['main.tex'][match.end():].startswith(b'\n'):
        candidates.append(files['main.tex'][:match.start()]+files['main.tex'][match.end()+1:])
    old_main = next(row for row in previous['upload_files'] if row['relativePath'] == 'main.tex')
    restored = [raw for raw in candidates if len(raw) == old_main['bytes'] and sha(raw) == old_main['sha256']]
    if len(restored) != 1:
        raise api.Stop('Removing the sole remark does not restore the accepted v1 source bytes')
    entry = dict(previous, upload_files=rows,
        reference_pdf=dict(previous['reference_pdf'], bytes=len(pdf_raw), sha256=PDF_SHA))
    inputs = dict(entry=entry, metadata=metadata, metadata_raw=metadata_raw,
        files=files, restored_v1_main=restored[0], old_receipt=old_receipt,
        old_receipt_raw=old_receipt_path.read_bytes(), prior=previous,
        build_raw=build_raw, package_raw=package_raw)
    return inputs


def prepare(api):
    inputs = authorized_inputs(api)
    old = inputs['prior']
    snapshot = HERE/'payload_history.local'
    snapshot.mkdir(mode=0o700, exist_ok=True)
    save_private(snapshot/'v1.main.tex', inputs['restored_v1_main'], api)
    save_private(snapshot/'v1.references.bib', inputs['files']['references.bib'], api)
    save_private(snapshot/'v1.main.bbl', inputs['files']['main.bbl'], api)
    save_private(snapshot/'v1.accepted-receipt.json', inputs['old_receipt_raw'], api)
    for name, raw in inputs['files'].items():
        save_private(snapshot/('updated.'+name), raw, api)
    save_private(snapshot/'metadata.json', inputs['metadata_raw'], api)
    document = dict(schema='hexagon-lp-vinberg-remark-update-authorization-v1',
        authorization=AUTHORIZATION, target_identifier=IDENTIFIER,
        prior_version_id=inputs['old_receipt']['response']['receipt']['versionId'],
        prior_approved_payloads_sha256=sha((LEGACY/'approved_payloads.json').read_bytes()),
        prior_accepted_receipt_sha256=sha(inputs['old_receipt_raw']),
        metadata_sha256=sha(api.canonical(inputs['metadata'])),
        sole_change='The concise remark labeled fv:vinberg-weight-obstruction',
        removing_remark_restores_v1_source_sha256=sha(inputs['restored_v1_main']),
        old_source_inventory=old['upload_files'], updated_source_inventory=inputs['entry']['upload_files'],
        updated_reference_pdf=inputs['entry']['reference_pdf'], upload_names=['main.tex', 'references.bib'],
        build_receipt=dict(path=BUILD, sha256=sha(inputs['build_raw'])),
        package_receipt=dict(path=PACKAGE, sha256=sha(inputs['package_raw'])),
        title_author_abstract_bibliography_unchanged=True, retained_bounds=[57,17,15],
        original_payload_bytes_preserved=True, new_paper_creation_authorized=False,
        local_mathematical_or_TeX_execution_authorized=False)
    target = HERE/'authorization.json'
    raw = json.dumps(document, indent=2)+'\n'
    if target.exists() and target.read_text() != raw:
        raise api.Stop('Existing public update authorization differs')
    if not target.exists():
        target.write_text(raw, encoding='utf-8')
    print(json.dumps(dict(prepared=True, target_identifier=IDENTIFIER,
        source_sha256=SOURCE_SHA, pdf_sha256=PDF_SHA, pages=36,
        exact_v1_source_restored=True, new_paper_creation=False)), flush=True)
    return inputs


def context(api):
    inputs = authorized_inputs(api)
    if not (HERE/'authorization.json').exists():
        raise api.Stop('Prepare and retain the explicit update authorization first')
    # Use the reviewed upload/preview/commit predicates with separate update
    # history. The global unresolved markers and private mutation lock stay
    # shared across workflows; credential and repository roots stay fixed.
    api.HERE = HERE
    api.RECEIPT_PATH = HERE/'execution.local'
    def preparation(paper):
        if paper != PAPER:
            raise api.Stop('This update is restricted to the existing LP paper')
        fresh = authorized_inputs(api)
        return fresh['entry'], fresh['metadata'], fresh['files'], []
    api.preparation = preparation
    def bind(draft, paper, metadata):
        control = read(HERE/'control.local')
        if (paper != PAPER or control.get('draft_id') != draft
                or control.get('target_identifier') != IDENTIFIER
                or control.get('metadata_sha256') != sha(api.canonical(metadata))):
            raise api.Stop('Correction draft is not bound to the existing LP identifier')
    api.bind_draft = bind
    return inputs


def records(api):
    return read(api.RECEIPT_PATH) if api.RECEIPT_PATH.exists() else []


def exchange(api, inputs, action, method, path, body=None, filename=None):
    rows = records(api)
    answered = {row.get('request_id') for row in rows if row.get('state') == 'HTTP_RESPONSE_RECEIVED'}
    if method != 'GET' and any(row.get('state') == 'REQUEST_PENDING_OR_RESULT_UNKNOWN'
            and row.get('request_id') not in answered for row in rows):
        raise api.Stop('An update mutation has an unknown result; reconcile it first')
    args = argparse.Namespace(paper=PAPER, file=filename)
    return api.exchange(args, inputs['entry'], inputs['metadata'], api.credential(), action,
        method, path, None if body is None else api.canonical(body))


def target_status(api, inputs):
    draft = inputs['old_receipt']['draft_id']
    code, value = exchange(api, inputs, 'original-status', 'GET', '/api/v1/submissions/'+draft)
    if code != 200 or value.get('id') != draft:
        raise api.Stop('Could not verify the existing LP submission status')
    print(json.dumps(dict(target_identifier=IDENTIFIER, current_status=value.get('status'),
        reported_version=value.get('versionId'), update_api_supported=True)), flush=True)
    return value


def open_correction(api, inputs):
    control_path = HERE/'control.local'
    if control_path.exists():
        return read(control_path)
    original = target_status(api, inputs)
    parent_id = inputs['old_receipt']['draft_id']
    parent_base = '/api/v1/submissions/'+parent_id
    if original.get('status') in {'under_review', 'processing', 'submitted'}:
        code, _ = exchange(api, inputs, 'pause-for-corrections', 'POST', parent_base+'/corrections')
        if code != 200:
            raise api.Stop('The service refused to pause this package for corrections; no new work created')
    elif original.get('status') not in {'published', 'changes_requested', 'awaiting_changes'}:
        raise api.Stop('Existing package is not in a documented correction/revision state')
    code, value = exchange(api, inputs, 'open-revision', 'POST', parent_base+'/revise')
    if code != 200:
        raise api.Stop('The existing package has no editable correction/revision draft')
    candidate = value.get('draftId')
    if not isinstance(candidate, str) and isinstance(value.get('draft'), dict):
        candidate = value['draft'].get('id')
    draft = api.component(candidate, 'correction draft identifier')
    if draft == parent_id:
        raise api.Stop('The service did not return a distinct editable correction draft')
    control = dict(target_identifier=IDENTIFIER, original_draft_id=parent_id, draft_id=draft,
        metadata_sha256=sha(api.canonical(inputs['metadata'])),
        files_sha256=sha(api.canonical(inputs['entry']['upload_files'])),
        opened_utc=dt.datetime.now(dt.timezone.utc).isoformat())
    api.write_private_json_once(control_path, control)
    return control


def upload(api, inputs, control, operation_lock):
    draft = control['draft_id']
    base = '/api/v1/submissions/'+draft
    # Revision metadata may initially have blank version comments. Restore the
    # already-authorized metadata exactly, including all personal declarations.
    code, _ = exchange(api, inputs, 'update-metadata', 'PATCH', base, inputs['metadata'])
    if code != 200:
        raise api.Stop('Could not preserve the approved LP metadata on the correction draft')
    args = argparse.Namespace(paper=PAPER, action='draft-status', file=None,
        draft_id=draft, upload_id=None, operation_lock=operation_lock)
    if api.execute_args(args) != 0:
        raise api.Stop('Correction metadata verification failed')
    api.bind_verified_metadata(draft, PAPER, inputs['metadata'], inputs['entry'])
    latest = [row for row in records(api) if row.get('action') == 'draft-status'
              and row.get('http_status') == 200][-1]['response']
    # Revise clones the prior package into a distinct editable draft. Remove
    # only its two staged copies, leaving the committed v1 objects/history
    # untouched, then send the two exact authorized source files anew.
    copies = latest.get('uploads') or []
    if copies:
        names = [row.get('relativePath') for row in copies]
        if len(names) != len(set(names)) or not set(names) <= api.HEXAGON_SOURCE_FILES:
            raise api.Stop('Correction draft contains an unexpected staged source')
        prior = {row['relativePath']: row for row in inputs['prior']['upload_files']}
        updated = {row['relativePath']: row for row in inputs['entry']['upload_files']}
        for row in copies:
            filename = row['relativePath']
            first, second = prior[filename], updated[filename]
            if (row.get('filename') != filename or row.get('kind') != first['kind']
                    or row.get('contentType') != first['contentType']
                    or row.get('sizeBytes') not in {first['bytes'], second['bytes']}
                    or row.get('sha256') not in {None, first['sha256'], second['sha256']}):
                raise api.Stop('A staged copy differs from the authorized prior/update inventory')
        for row in copies:
            filename = row['relativePath']
            upload_id = api.component(row.get('id'), 'staged correction upload identifier')
            code, _ = exchange(api, inputs, 'remove-staged-copy', 'DELETE',
                base+'/uploads/'+upload_id, filename=filename)
            if code != 200:
                raise api.Stop('The service refused staged-copy replacement with this credential; owner-browser action may be required')
        if api.execute_args(args) != 0:
            raise api.Stop('Could not verify the cleared correction staging area')
        cleared = [row for row in records(api) if row.get('action') == 'draft-status'
                   and row.get('http_status') == 200][-1]['response']
        if cleared.get('uploads'):
            raise api.Stop('Correction staging still contains files after replacement preparation')
    for filename in sorted(api.HEXAGON_SOURCE_FILES):
        args.file = filename
        args.action = 'create-upload'
        if api.execute_args(args) != 0:
            raise api.Stop('Correction source registration failed')
        row = [row for row in records(api) if row.get('action') == 'create-upload'
               and row.get('http_status') in {200,201} and row.get('manifest_file') == filename][-1]
        args.upload_id = api.component(row['response'].get('uploadId'), 'update upload identifier')
        for action in ('put-part', 'complete-upload'):
            args.action = action
            if api.execute_args(args) != 0:
                raise api.Stop('Correction source upload did not finish')
    args.action, args.file, args.upload_id = 'draft-status', None, None
    if api.execute_args(args) != 0:
        raise api.Stop('Correction source verification failed')
    return args


def finish(api, inputs, control, operation_lock, science_lock):
    draft = control['draft_id']
    args = argparse.Namespace(paper=PAPER, action='draft-status', file=None,
        draft_id=draft, upload_id=None, operation_lock=operation_lock, scientific_lock=science_lock)
    for action in ('draft-status', 'source-check', 'preview-status'):
        args.action = action
        if api.execute_args(args) != 0:
            raise api.Stop('Correction source checks failed before preview')
    if not (HERE/(PAPER+'.preview-resolution.local')).exists():
        args.action = 'preview'
        if api.execute_args(args) != 0:
            raise api.Stop('Correction preview did not complete')
    if not (HERE/(PAPER+'.preview-artifacts.local')).exists():
        args.action = 'preview-artifacts'
        if api.execute_args(args) != 0:
            raise api.Stop('Correction preview artifacts were not retrieved')
    if not (HERE/(PAPER+'.preview-content-review.local')).exists():
        args.action = 'inspect-preview'
        if api.execute_args(args) != 0:
            raise api.Stop('Correction preview comparison failed')
    review = read(HERE/(PAPER+'.preview-content-review.local'))
    if review.get('content_review_passed') is not True or review.get('reference_pages') != 36:
        raise api.Stop('Correction preview does not match the authorized local edition')
    today = dt.datetime.now(dt.timezone.utc).date().isoformat()
    total_commits = sum(row.get('action') == 'commit' and row.get('http_status') == 202
        and row.get('utc_time', '').startswith(today)
        for row in read(LEGACY/'execution.local')+records(api))
    if total_commits >= api.DAILY_COMMIT_LIMIT:
        raise api.Stop('The user-authorized daily allowance is exhausted')
    args.action = 'commit'
    if api.execute_args(args) != 0:
        raise api.Stop('Correction commit was not accepted')
    return report(api, inputs, control)


def report(api, inputs, control):
    receipt_path = HERE/(PAPER+'.submission-receipt.local')
    receipt = read(receipt_path)
    server = receipt.get('response', {}).get('receipt', {})
    if (receipt.get('http_status') != 202 or server.get('identifier') != IDENTIFIER
            or re.fullmatch(re.escape(IDENTIFIER)+r'v[1-9][0-9]*', server.get('versionId', '')) is None):
        raise api.Stop('Accepted correction is not the existing LP work; inspect its receipt')
    code, status = exchange(api, inputs, 'updated-status', 'GET',
                           '/api/v1/submissions/'+control['draft_id'])
    if code != 200:
        raise api.Stop('Accepted correction status could not be refreshed')
    review = read(HERE/(PAPER+'.preview-content-review.local'))
    prior = read(HERE/'authorization.json')
    if (sha((LEGACY/(PAPER+'.submission-receipt.local')).read_bytes())
            != prior['prior_accepted_receipt_sha256']
            or sha((LEGACY/'approved_payloads.json').read_bytes())
            != prior['prior_approved_payloads_sha256']):
        raise api.Stop('Original accepted receipt or approved payload history changed')
    document = dict(schema='hexagon-lp-vinberg-remark-correction-completion-v1',
        target_identifier=IDENTIFIER, accepted_version_id=server['versionId'],
        current_status=status.get('status'), publication_confirmed=status.get('status') == 'published',
        commit_http_status=202, accepted_utc=receipt['accepted_utc'],
        same_existing_work=True, new_paper_created=False,
        authorization_sha256=sha((HERE/'authorization.json').read_bytes()),
        accepted_update_receipt_sha256=sha(receipt_path.read_bytes()),
        old_receipt_and_payload_history_unchanged=True, original_v1_source_bytes_preserved=True,
        source_sha256=SOURCE_SHA, reference_pdf_sha256=PDF_SHA,
        preview={key: review.get(key) for key in ('reference_pages', 'preview_pages',
            'normalized_full_text_equal', 'normalized_bibliography_equal', 'bibtex_executed',
            'final_pass_has_undefined_references', 'content_review_passed', 'pdf_sha256')},
        workflow_source_sha256=sha(Path(__file__).read_bytes()),
        client_source_sha256=sha((LEGACY/'client.py').read_bytes()),
        global_unresolved_markers_clear=not api.PREVIEW_MARKER.exists() and not api.COMMIT_MARKER.exists(),
        no_local_TeX_or_scientific_jobs=True, retained_bounds=[57,17,15],
        title_author_abstract_bibliography_unchanged=True, bulk_upload_performed=False)
    target = HERE/'completion.json'
    if not target.exists():
        target.write_text(json.dumps(document, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: document[key] for key in ('target_identifier', 'accepted_version_id',
        'current_status', 'publication_confirmed', 'same_existing_work', 'new_paper_created',
        'global_unresolved_markers_clear')}), flush=True)
    return 0


def run(api, step):
    if os.environ.get('COXETER_PUBLICATION_IO_GUARD') != '1':
        raise api.Stop('Use the publication I/O watchdog')
    if step == 'prepare':
        prepare(api)
        return 0
    inputs = context(api)
    if step == 'status':
        target_status(api, inputs)
        return 0
    if step not in {'run', 'finish'}:
        raise api.Stop('Unknown update stage')
    if api.PREVIEW_MARKER.exists() or api.COMMIT_MARKER.exists():
        raise api.Stop('A global remote operation is unresolved; reconcile it before continuing')
    rows = records(api)
    answered = {row.get('request_id') for row in rows if row.get('state') == 'HTTP_RESPONSE_RECEIVED'}
    if any(row.get('state') == 'REQUEST_PENDING_OR_RESULT_UNKNOWN'
           and row.get('request_id') not in answered for row in rows):
        raise api.Stop('An earlier update request has an unknown result')
    fd = os.open(LEGACY/'submission.lock.local', os.O_RDWR|os.O_CREAT|os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, 'a') as operation_lock:
        try:
            fcntl.flock(operation_lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise api.Stop('Another publication mutation is running') from None
        with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as science_lock:
            try:
                fcntl.flock(science_lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:
                raise api.Stop('The shared scientific/typesetting lane is occupied') from None
            if (HERE/(PAPER+'.submission-receipt.local')).exists():
                return report(api, inputs, read(HERE/'control.local'))
            if step == 'run':
                control = open_correction(api, inputs)
                upload(api, inputs, control, operation_lock)
            else:
                control = read(HERE/'control.local')
            return finish(api, inputs, control, operation_lock, science_lock)
