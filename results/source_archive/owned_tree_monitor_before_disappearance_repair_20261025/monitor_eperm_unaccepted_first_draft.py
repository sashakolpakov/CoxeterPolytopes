#!/usr/bin/env python3
"""Track and clean only an observed owned proof process tree.

This is an additional owner/cleanup bridge, not a proof or a replacement
resource runner. Actual children remain the unchanged bounded runners.
Birth identities are checked before every group signal. Polling is not an
atomic pidfd facility: unobserved adversarial fork/reparent races and the
kernel check-to-signal interval are explicit limitations on macOS.
"""
import os
import sys

if not __debug__ or sys.flags.optimize:
    raise SystemExit('ASSERTIONS_REQUIRED: owned proof monitoring needs normal Python')
if os.environ.get('PYTHONOPTIMIZE', '') not in ('', '0'):
    raise SystemExit('ASSERTIONS_REQUIRED: remove PYTHONOPTIMIZE before verification')

import argparse
import ctypes
from dataclasses import dataclass
import errno
import hashlib
import json
from pathlib import Path
import signal
import struct
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
ALLOWED_RUNNERS = {
    ROOT / 'scripts/run_bounded_research_check_20261018.py',
    ROOT / 'scripts/build_illustrated_ra_manuscript_20261018.py',
    ROOT / 'scripts/build_compact_arithmetic_nine_field19_20261020.py',
}

CONTROLLERS = (
    'scripts/run_proof_batch_20261024.py',
    'scripts/audit_existing_lean_module_20261024.py',
    'scripts/verify_manuscript_audited_20261024.py',
    'scripts/monitor_owned_proof_tree_20261024.py',
    'scripts/run_bounded_research_check_20261018.py',
    'scripts/run_lean_batches_signal_safe_20261018.py',
    'scripts/run_lean_batches_20261018.py',
    'scripts/build_illustrated_ra_manuscript_20261018.py',
    'scripts/build_compact_arithmetic_nine_field19_20261020.py',
    'scripts/verify_manuscript.py',
)


def execution_sources():
    """Bind controller code only, not the Python/OS or mathematical inputs."""
    result = {}
    for name in CONTROLLERS:
        path = (ROOT / name).resolve()
        if path.parent != ROOT / 'scripts' or not path.is_file():
            raise RuntimeError('CONTROLLER_SOURCE_MISSING_OR_EXTERNAL: ' + name)
        before = path.stat()
        with path.open('rb') as stream:
            digest = hashlib.file_digest(stream, 'sha256').hexdigest()
        after = path.stat()
        identity = lambda stat: (stat.st_dev, stat.st_ino, stat.st_size,
                                 stat.st_mtime_ns, stat.st_ctime_ns)
        if identity(before) != identity(after):
            raise RuntimeError('CONTROLLER_SOURCE_CHANGED_WHILE_HASHING: ' + name)
        result[name] = digest
    return result


@dataclass(frozen=True)
class Process:
    pid: int
    ppid: int
    pgid: int
    born: tuple
    zombie: bool


class ProcessReader:
    """Read process metadata; denied monitoring fails before any proof child."""
    def __init__(self):
        if sys.platform == 'darwin':
            # sys/proc_info.h: PROC_PIDTBSDINFO=3, proc_bsdinfo is136 bytes;
            # pid/ppid/pgid are at12/16/100, start sec/usec at120/128.
            self.library = ctypes.CDLL('/usr/lib/libproc.dylib', use_errno=True)
            self.proc = self.library.proc_pidinfo
            self.proc.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_uint64,
                                  ctypes.c_void_p, ctypes.c_int]
            self.proc.restype = ctypes.c_int
        elif not sys.platform.startswith('linux'):
            raise RuntimeError('OWNED_TREE_UNSUPPORTED: need Darwin proc_pidinfo or Linux /proc')

    def process(self, pid):
        if pid <= 0:
            return None
        if sys.platform == 'darwin':
            buffer = ctypes.create_string_buffer(136)
            ctypes.set_errno(0)
            length = self.proc(pid, 3, 0, ctypes.byref(buffer), 136)
            error = ctypes.get_errno()
            if length == 0 and error in (errno.ESRCH, errno.ENOENT):
                return None
            if length == 0 and error == errno.EPERM:
                # Darwin can report EPERM for a descendant which disappeared
                # after the process-table sample. Accept absence only when a
                # fresh, successful ps query confirms that exact PID is gone.
                # A live, denied, malformed or indeterminate query still fails.
                absent = subprocess.run(['ps', '-p', str(pid), '-o', 'pid='],
                                        capture_output=True, text=True, timeout=3)
                if (absent.returncode == 1 and not absent.stdout.strip()
                        and not absent.stderr.strip()):
                    return None
            if length != 136:
                raise RuntimeError('OWNED_TREE_IDENTITY_DENIED_OR_ABI: pid=' + str(pid)
                                   + ' errno=' + str(error))
            raw = buffer.raw
            actual, parent = struct.unpack_from('=II', raw, 12)
            group = struct.unpack_from('=I', raw, 100)[0]
            born = struct.unpack_from('=QQ', raw, 120)
            state = struct.unpack_from('=I', raw, 4)[0]
            if actual != pid or not born[0]:
                raise RuntimeError('OWNED_TREE_INVALID_IDENTITY: pid=' + str(pid))
            return Process(pid, parent, group, born, state == 5)
        try:
            raw = Path('/proc', str(pid), 'stat').read_text()
        except (FileNotFoundError, ProcessLookupError):
            return None
        tail = raw.rsplit(')', 1)[1].split()
        # tail starts at field3; Linux starttime is field22.
        return Process(pid, int(tail[1]), int(tail[2]), (int(tail[19]),), tail[0] == 'Z')

    def table(self):
        result = subprocess.run(['ps', '-axo', 'pid=,ppid=,pgid='],
                                check=True, capture_output=True, text=True, timeout=3)
        if len(result.stdout) > 2 * 1024**2:
            raise RuntimeError('OWNED_TREE_PROCESS_TABLE_TOO_LARGE')
        rows = {}
        for line in result.stdout.splitlines():
            if not line.strip():
                continue
            values = line.split()
            if len(values) != 3:
                raise RuntimeError('OWNED_TREE_MALFORMED_PROCESS_TABLE')
            pid, parent, group = map(int, values)
            rows[pid] = (parent, group)
        return rows

    def probe(self):
        self.table()
        own = self.process(os.getpid())
        if own is None or own.zombie:
            raise RuntimeError('OWNED_TREE_MONITORING_UNAVAILABLE')


class OwnedTree:
    """Track births reached through live owned parents or owned process groups.

    A group is trusted only while at least one already observed birth remains
    in it and its numerical group leader has not been replaced. A recycled
    PID never inherits an old ownership record. Unknown groups are not killed.
    """
    def __init__(self, child, reader=None):
        self.child = child
        self.reader = reader or ProcessReader()
        first = self.reader.process(child.pid)
        if first is None:
            raise RuntimeError('OWNED_TREE_ROOT_DISAPPEARED_BEFORE_TRACKING')
        self.known = {first.pid: first}
        self.depth = {first.pid: 0}
        self.group_leaders = {first.pgid: first.born} if first.pid == first.pgid else {}
        self.signals = []
        self.skipped_reused = set()

    def _live(self):
        live = {}
        for pid, recorded in list(self.known.items()):
            current = self.reader.process(pid)
            if current is not None and current.born != recorded.born:
                self.skipped_reused.add(pid)
                continue
            if current is not None and not current.zombie:
                live[pid] = current
        return live

    def _trusted_groups(self, live):
        groups = set()
        for info in live.values():
            group = info.pgid
            if group <= 1 or group == os.getpgrp():
                continue
            leader = self.reader.process(group)
            recorded = self.group_leaders.get(group)
            if leader is not None:
                if recorded is None or leader.born != recorded:
                    continue
            groups.add(group)
        return groups

    def sample(self):
        rows = self.reader.table()
        live = self._live()
        while True:
            additions = {}
            groups = self._trusted_groups(live)
            for pid, (parent, group) in rows.items():
                if pid in live or (parent not in live and group not in groups):
                    continue
                current = self.reader.process(pid)
                if current is None or current.zombie:
                    continue
                old = self.known.get(pid)
                if old is not None and old.born != current.born:
                    self.skipped_reused.add(pid)
                    continue
                parent_owned = False
                if current.ppid in live:
                    ancestor = self.reader.process(current.ppid)
                    parent_owned = (ancestor is not None
                                    and ancestor.born == live[current.ppid].born)
                if not parent_owned and current.pgid not in groups:
                    continue
                additions[pid] = current
                self.depth[pid] = self.depth.get(current.ppid, 0) + 1
                if pid == current.pgid:
                    previous = self.group_leaders.get(pid)
                    if previous is not None and previous != current.born:
                        self.skipped_reused.add(pid)
                        additions.pop(pid)
                        continue
                    self.group_leaders[pid] = current.born
            if not additions:
                break
            self.known.update(additions)
            live.update(additions)
        return live

    def _signal(self, group, signum, live):
        if group <= 1 or group == os.getpgrp():
            raise RuntimeError('OWNED_TREE_REFUSED_UNSAFE_GROUP: ' + str(group))
        members = [info for info in live.values() if info.pgid == group]
        if not members:
            return
        # Revalidate the birth and group immediately before the signal.
        valid = []
        for member in members:
            current = self.reader.process(member.pid)
            if (current is not None and not current.zombie
                    and current.born == member.born and current.pgid == group):
                valid.append(current)
        if not valid or group not in self._trusted_groups({i.pid: i for i in valid}):
            return
        try:
            os.killpg(group, signum)
            self.signals.append({'pgid': group, 'signal': int(signum),
                                 'verified_members': [i.pid for i in valid]})
        except ProcessLookupError:
            pass

    def _groups_by_depth(self, live):
        groups = self._trusted_groups(live)
        return sorted(groups, key=lambda g: max(self.depth.get(p, 0)
                      for p, i in live.items() if i.pgid == g), reverse=True)

    def stop(self, grace_seconds=2):
        # Inner compiler groups are signalled before their Python monitors.
        live = self.sample()
        for group in self._groups_by_depth(live):
            self._signal(group, signal.SIGTERM, live)
        deadline = time.monotonic() + grace_seconds
        while live and time.monotonic() < deadline:
            time.sleep(.05)
            live = self.sample()
        for _ in range(3):
            if not live:
                break
            for group in self._groups_by_depth(live):
                self._signal(group, signal.SIGKILL, live)
            time.sleep(.05)
            live = self.sample()
        try:
            self.child.wait(timeout=1)
        except subprocess.TimeoutExpired:
            pass
        if live:
            raise RuntimeError('OWNED_TREE_CLEANUP_INCOMPLETE: verified live pids '
                               + ','.join(map(str, sorted(live))))
        return {'observed_processes': len(self.known),
                'signals': self.signals,
                'reused_ids_skipped': sorted(self.skipped_reused),
                'verified_live_owned_processes_remaining': 0}


GATE = """import os,sys
fd=int(sys.argv[1]); ready=os.read(fd,1); os.close(fd)
if ready!=b'1': raise SystemExit(125)
os.execv(sys.argv[2],sys.argv[2:])
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if (len(command) < 4 or Path(command[0]).resolve() != Path(sys.executable).resolve()
            or command[1:3] != ['-E', '-B']
            or Path(command[3]).resolve() not in ALLOWED_RUNNERS):
        parser.error('Only unchanged bounded Python research/TeX runners are allowed')
    reader = ProcessReader()
    reader.probe()  # No child exists before monitoring access is established.
    controller_sources = execution_sources()
    owner = reader.process(os.getppid())
    if owner is None:
        raise RuntimeError('OWNED_TREE_OWNER_MISSING')
    stop_reason = []
    def interrupted(signum, frame):
        stop_reason.append('signal' + str(signum))
    for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(signum, interrupted)
    read_fd, write_fd = os.pipe()
    child, tree, cleanup = None, None, None
    try:
        environment = dict(os.environ)
        environment.pop('PYTHONOPTIMIZE', None)
        child = subprocess.Popen([sys.executable, '-E', '-B', '-c', GATE,
                                  str(read_fd), *command], cwd=ROOT,
                                 env=environment, start_new_session=True,
                                 pass_fds=(read_fd,))
        os.close(read_fd)
        read_fd = None
        tree = OwnedTree(child, reader)
        tree.sample()
        current = reader.process(owner.pid)
        if (stop_reason or os.getppid() != owner.pid or current is None
                or current.born != owner.born):
            raise RuntimeError('OWNED_TREE_CANCELLED_BEFORE_RELEASE')
        started = time.monotonic()
        os.write(write_fd, b'1')
        os.close(write_fd)
        write_fd = None
        while True:
            tree.sample()
            current = reader.process(owner.pid)
            if (os.getppid() != owner.pid or current is None
                    or current.born != owner.born):
                raise RuntimeError('OWNED_TREE_OWNER_EXITED')
            if stop_reason:
                raise RuntimeError('OWNED_TREE_INTERRUPTED: ' + stop_reason[0])
            if time.monotonic() - started >= 59:
                raise RuntimeError('OWNED_TREE_BRIDGE_TIME_LIMIT: 59seconds')
            code = child.poll()
            if code is not None:
                break
            time.sleep(.05)
        cleanup = tree.stop()
        if execution_sources() != controller_sources:
            raise RuntimeError('CONTROLLER_SOURCE_CHANGED_DURING_JOB')
        print(json.dumps({'owned_tree_bridge_only': True,
                          'child_exit_code': code, 'cleanup': cleanup,
                          'controller_source_sha256': controller_sources,
                          'assertions_enabled': __debug__,
                          'python_optimize_flag': sys.flags.optimize,
                          'child_pythonoptimize': environment.get('PYTHONOPTIMIZE')}),
              flush=True)
        if code == 0 and cleanup['signals']:
            raise RuntimeError('OWNED_TREE_STRAGGLERS_AFTER_SUCCESS: cleaned, not a proof PASS')
        return code if code >= 0 else 128 - code
    finally:
        if read_fd is not None:
            os.close(read_fd)
        if write_fd is not None:
            os.close(write_fd)
        if tree is not None and cleanup is None:
            tree.stop()
        elif child is not None and tree is None:
            # This direct child has not been reaped and is still behind GATE.
            # It has not executed a runner or a mathematical worker.
            child.terminate()
            try:
                child.wait(timeout=1)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (RuntimeError, ValueError, OSError, subprocess.SubprocessError) as error:
        print(str(error), file=sys.stderr)
        sys.exit(1)
