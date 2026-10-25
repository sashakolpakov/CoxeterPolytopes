#!/usr/bin/env python3
"""Preserve each complete unchanged large-D result with lossless encoding.

The parent completed its assertions but exceeded its public output cap. This
wrapper calls the same pinned function, checks a complete byte round trip,
and emits a bounded compressed result. It changes no domain or assertion.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
import base64
from hashlib import sha256
import json
from pathlib import Path
import types
import zlib

ROOT = Path(__file__).resolve().parents[1]
PARENT = 'scripts/check_finite57_current_universal_receivers_20261025.py'
PARENT_SHA = '29ddbcf65cae2a7c7ff7703a1b4a81bad3ca1fdb0ade4bdbffda26fb9c8a1631'
MAX_DECODED = 8 * 1024**2
MAX_OUTPUT = 500 * 1024
BATCHES = tuple(f'large-d-{i}' for i in range(4))


def need(condition, message):
    if not condition:
        raise ValueError(message)


def raw(name):
    path = (ROOT/name).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'existing repository source')
    with path.open('rb') as stream:
        data = stream.read(4*1024**2+1)
    need(0 < len(data) <= 4*1024**2, 'bounded exact source')
    return data


def encode_result(result):
    decoded = json.dumps(result, sort_keys=True, separators=(',', ':')).encode('ascii')
    need(0 < len(decoded) <= MAX_DECODED, 'bounded COMPLETE decoded large-D result')
    compressed = zlib.compress(decoded, level=6)
    need(zlib.decompress(compressed) == decoded, 'byte-exact complete result round trip')
    return dict(encoding='zlib-base64-json-v1',
        decoded_bytes=len(decoded), decoded_sha256=sha256(decoded).hexdigest(),
        compressed_bytes=len(compressed),
        data=base64.b64encode(compressed).decode('ascii'))


def decode_result(packet):
    """Strict bounded inverse, also used by the separately reviewed consumer."""
    need(type(packet) is dict and set(packet) == {
        'encoding', 'decoded_bytes', 'decoded_sha256', 'compressed_bytes', 'data'},
        'exact compressed-result schema')
    need(packet['encoding'] == 'zlib-base64-json-v1' and
         type(packet['decoded_bytes']) is int and 0 < packet['decoded_bytes'] <= MAX_DECODED and
         type(packet['compressed_bytes']) is int and 0 < packet['compressed_bytes'] <= MAX_OUTPUT and
         type(packet['data']) is str and len(packet['data']) <= MAX_OUTPUT,
         'bounded compressed-result lengths')
    compressed = base64.b64decode(packet['data'], validate=True)
    need(len(compressed) == packet['compressed_bytes'], 'exact compressed length')
    decoder = zlib.decompressobj()
    decoded = decoder.decompress(compressed, packet['decoded_bytes']+1)
    need(len(decoded) == packet['decoded_bytes'] and decoder.eof and
         not decoder.unused_data and not decoder.unconsumed_tail and
         sha256(decoded).hexdigest() == packet['decoded_sha256'],
         'complete bounded decoded bytes, digest and end of stream')
    # The final consumer uses its own duplicate-safe exact JSON loader on
    # these bytes. This function never discards or silently truncates data.
    return decoded


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', choices=BATCHES, required=True)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    args = parser.parse_args()
    own = str(Path(__file__).resolve().relative_to(ROOT))
    before, parent_bytes = raw(own), raw(PARENT)
    need(sha256(parent_bytes).hexdigest() == PARENT_SHA, 'exact unchanged universal parent')
    u = types.ModuleType('finite57_large_d_unchanged_parent')
    u.__file__ = str(ROOT/PARENT)
    exec(compile(parent_bytes, u.__file__, 'exec'), u.__dict__)
    a, r, c, extra, ledger = u.prepare()
    need(PARENT not in c.PINS or c.PINS[PARENT] == PARENT_SHA, 'parent pin collision')
    c.PINS[PARENT] = PARENT_SHA
    c.READ.add(PARENT)
    result = u.large(a, r, c, int(args.batch[-1]))
    encoded = encode_result(result)
    need(json.loads(decode_result(encoded)) == result, 'complete result object round trip')
    for name in c.READ:
        need(c.digest(name) == c.PINS[name], 'consumed large-D source/input changed')
    need(raw(own) == before and raw(PARENT) == parent_bytes and
         raw(u.LEDGER) == ledger and sha256(raw(u.LOADER)).hexdigest() == u.LOADER_SHA,
         'worker/parent/shared loader/ledger changed')
    need(raw(a.EXTRA) == extra and sha256(raw(a.ROUTE)).hexdigest() == a.ROUTE_SHA and
         sha256(raw(r.COLLECTOR)).hexdigest() == r.COLLECTOR_SHA and
         sha256(c.raw(c.MANIFEST)).hexdigest() == c.MANIFEST_SHA,
         'shared frozen preparation chain changed')
    payload = dict(status='PASS_COMPLETE_LARGE_D_RECEIVER_RESULT_LOSSLESSLY_ENCODED',
        batch=args.batch, source_sha256=sha256(before).hexdigest(),
        parent_source_sha256=PARENT_SHA, result_encoding=encoded,
        consumed_inputs={p: c.PINS[p] for p in sorted(c.READ)}, fuel_used=c.FUEL[0],
        scope='The unchanged complete large-D shard result, losslessly encoded with both byte and object round-trip checks. All parent geometry, classification and numerical inputs retain their declared premise scopes.',
        independent_S8_verified=False, current_protected74_verified=False,
        current_terminal_bridge_complete=False, global_finite57_verified=False)
    output = json.dumps(payload, sort_keys=True)
    need(len(output.encode())+1 <= MAX_OUTPUT, 'complete encoded output cap; never truncate')
    print(output)


if __name__ == '__main__':
    main()
