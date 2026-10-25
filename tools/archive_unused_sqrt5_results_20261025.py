#!/usr/bin/env python3
"""Byte-preserving streamed archive of 43 unused parked sqrt5 JSON results.

This is neither publication packing nor a scientific/Lean/TeX launcher. It
does not change the40GiB/6GiB proof admission gate. It retains the original
until a streamed gzip round trip matches its fixed content identity.
"""
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
GIB = 1024**3
FILES = {
    'fresh_sqrt5_adjacent_discount_frontier_20261016':
        (54736542, '83947110599f43d4883a35e34833ac14eb0c23de921841ae2bd168335d6dcc07'),
    'fresh_sqrt5_all_quad_contacts_20261015':
        (39193205, 'c7a81ba0eb5689b8b0139d4a78e80b44c81a3bc4f30c539722ffb4ddb3f12feb'),
    'fresh_sqrt5_all_quad_unions_20261015':
        (30417018, '670539bc2fa64195bd497c31efa0b2ff13b6c018880c6b08e888a88c5d38143b'),
    'fresh_sqrt5_dimension21_a4_exceptional_20261016_cap2':
        (144777262, '9645cac9d7a62ce7e56e4692df4d96bbe7b879b14b1244d6dfd3d3f2af306da7'),
    'fresh_sqrt5_dimension21_a4_pair_exceptional_20261016_cap2':
        (144915640, '3660011ee75126289224ce5cf3305dc1d67d036a8f2d752f3a509ba35c731821'),
    'fresh_sqrt5_dimension21_a4_source_orbits_20261016_cap2':
        (27186536, 'de9cfd720adbb9976f6e6c3d0e39594ac0dd3d938477ebf451f3f5ab8c7b773d'),
    'fresh_sqrt5_dimension21_e7e8_context_20261016_cap2':
        (32854301, '994ae57b9502bfe62c176a652fe6ab85a31895d9c18d3dbd649fb67bf0449579'),
    'fresh_sqrt5_dimension21_e7e8_exact_20261016_cap2':
        (55236816, '56ffa14579f372d5508797cadc40cc736cdd27ef985b63774c8ddd0f3c63273b'),
    'fresh_sqrt5_global_b1_dotted_closure_20261016':
        (118658678, 'ca4ec6dff5ac96b883efff02e798c28fee5002d09a694c09f4a5f24e55b033a6'),
    'fresh_sqrt5_global_b1_dotted_fast_20261017':
        (514679613, '4b4168d24a3c998262dd0ff729a5332d97d1eddc570dbc21e0f310f3001944aa'),
    'fresh_sqrt5_global_b1_finite_closure_20261016':
        (826786125, '40b483cd2469c769a440ccc6976df51b0d1db6b2318cee14e283e0a4bea6c2c2'),
    'fresh_sqrt5_h3_hub_coverage_20261016_bound21':
        (26828124, 'dcf4be8a774f5f1961df2b3d882bd79a4d9a22c0e236259c26829882a603c9f1'),
    'fresh_sqrt5_integer_simplex_port_packing_20261017':
        (132103430, '5d25a38ff2f57d50b6081642c8e8996065ac336190e679236608599449d72c75'),
    'fresh_sqrt5_joint_five_port_diagnostic_20261017':
        (49285924, 'e42997ac1e8804adf8ead5a8118f656dcb5e23391b82e66af7fb97cda692bb1e'),
    'fresh_sqrt5_joint_four_port_diagnostic_20261017':
        (29508943, 'e4f7f81c29b12833c14e4485cb534f93d435cef12b306385277626b1044f1fe4'),
    'fresh_sqrt5_labelled_coface_base_packing_20261017':
        (51455513, '59231726adc0a68d3ef1a80bf0a06b8b512087ee3298748b6d64a4461b795817'),
    'fresh_sqrt5_nonsingular_rank29_fast_20261015':
        (68491759, 'f6f14c2b226221a0a507c3c2df6179c0aa1d3ec5bba91136b2db46eeab94d22f'),
    'fresh_sqrt5_prism_simplex_base_incidence_20261017':
        (30646244, '25ac36e08ef7a9ee33e26a6cadc2a28b8fc1ae66c2310d7fa255a6749c076275'),
    'fresh_sqrt5_rank13_active_triangle_dual_20261017':
        (44477706, '116a36c24e7c42be8c55db08d7b578b74ceb92c0e60aa6709ba112bc38ed8a28'),
    'fresh_sqrt5_rank13_all_selected_nonade_triangles_20261017':
        (585686433, '06ee37d5251737d81877ac1c284ca01b99d3579131ed4748e1f983d3ba11acb7'),
    'fresh_sqrt5_rank13_combined_disjoint_ports_20261017':
        (92008813, 'c352738277c6fb46ead52ed621ad3219e208b0e7594f5c6eb74cb23027fa7c74'),
    'fresh_sqrt5_rank13_discriminant_simplex_lp_20261017':
        (64149818, 'adfaf9fceb47ceffcd8a5d4f409c95e62355377a9841a9fcbb3616572fb8e569'),
    'fresh_sqrt5_rank13_dotted_b1_spectators_20261017':
        (169327409, '02f0ce96487423702ca4315f4914c5221e8881b3931149bbb7068e39ad2fd211'),
    'fresh_sqrt5_rank13_full_flexible_simplex_lp_20261017':
        (149368369, '74e659892cf8958f431256c247456f3fd6cf8126f3d4bb04d0afa844934715af'),
    'fresh_sqrt5_rank13_full_triangle_simplex_separation_20261017':
        (181750201, 'e21a488a3254689e93b695bb2e33d1b03e5bfa48f8354c2bac3f9252f3719272'),
    'fresh_sqrt5_rank13_high_ade_triangles_20261017':
        (37391288, '5a35ac8c4f2d2d257e785da4d96113a27c6983fc81f418dae39924580bf76e8b'),
    'fresh_sqrt5_rank13_irreducible_joint_quads_20261017':
        (80008894, '2d9c382f6bbdbc3943308380379d8e7387c162540c898d9141082dd60689a22a'),
    'fresh_sqrt5_rank13_joint_corner_all_nonade_repair_20261017':
        (302717179, '1dfde0ea6d3b6ac665210953f2624fd892a1b2a3721b31d0ff01f1d548c2ceac'),
    'fresh_sqrt5_rank13_joint_corner_double_h4_repair_20261017':
        (165971106, 'ec36267dd66c0a6aa689b6be530a6ca31f717956670235e48ed39c25f16c5881'),
    'fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017':
        (49493871, '16bcf91d7aab9cb177861c3c2fd20a48b0c03e6af5d27ccd5592dd25a33c2f81'),
    'fresh_sqrt5_rank13_joint_corner_repair1_20261017':
        (49495459, '86ccb18fbc821aef08345c20561c2cd76d3510a469db57001262d4baf5f8bf0a'),
    'fresh_sqrt5_rank13_joint_corner_spectator_repair_20261017':
        (162682064, '1ca9c193ffe7523dea442460faa3d9acf745e7918558f651e00083270c0d0ce1'),
    'fresh_sqrt5_rank13_packing_quad_separation_20261017':
        (38723945, '8f9ca7f7a163db92bdbd5b9586d0ee97c907596aded99c28aeadff60054a1734'),
    'fresh_sqrt5_rank13_simplex_aggregate_lp_20261017':
        (46862547, 'b9b6375a4d175510e7d0c3e72ee5f76929e527f6c056487b123b03ba623c1726'),
    'fresh_sqrt5_rank13_simplex_packing_repair_20261017':
        (149388919, 'e8aa86cf6f59e854e2248d6c2667ad4bb07077bf775bca569cc95e484ee8679e'),
    'fresh_sqrt5_rank13_simplex_quad_combined_20261017':
        (86131211, '54e52e3a435d93e92856b8590f21ec3b889043f44344a3d9dd8d009e3fa9a79f'),
    'fresh_sqrt5_rank13_small_d_joint_quads_20261017':
        (189057836, 'aa85f386933150932bb3d7a7ff476dd0466248d0714df825c7200ff1efd62f72'),
    'fresh_sqrt5_rank13_triangle_all_representatives_20261017':
        (29160898, 'dc9d93e0b0d0e69b44aa9f25c1e72b16d9ee7244df8832477309329249059a32'),
    'fresh_sqrt5_rank13_triangle_simplex_flags_20261017':
        (43871477, '36ca67cc0264a1594e9c188affe43c174d88f5dd6e49995337e470f0ec32e8b2'),
    'fresh_sqrt5_rank13_triangle_third_costs_20261017':
        (26310570, 'ede5249677825a8c7e37b3cb6d8c39124be46d11b0ca2cb1d40ed12b7feae711'),
    'fresh_sqrt5_rank14_core_spectator_domain_20261017':
        (36053278, 'ef7906f4729114f95d9b73805a5a752222f74c1983f1abba712d1952a5414331'),
    'fresh_sqrt5_rank14_core_spectator_lp_20261017':
        (960803140, '677f5fe33bda068b1e2354c24cb7a0305a6de9ae0cf01c6127fb5dccefa12b28'),
    'fresh_sqrt5_rank14_triangle_padding_barrier_20261017':
        (48348777, '9238b1c40a4bead8017f43454fcdaf29c0b18d8dd00cb4f550088d23eaddbcf5'),
}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def identity(path):
    s = path.stat()
    return s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns


def main():
    started = time.monotonic()
    peak = 0
    last = 0
    owner = os.getppid()
    group = os.getpgrp()

    def guard(force=False):
        nonlocal peak, last
        require(os.getppid() == owner, 'Storage owner exited')
        require(time.monotonic()-started < 59, 'Storage-maintenance59s limit')
        now = time.monotonic()
        if not force and last and now-last < .2:
            return
        last = now
        disk = shutil.disk_usage(ROOT).free
        require(disk >= 20*GIB, 'Storage maintenance requires20GiB free')
        vm = subprocess.run(['vm_stat'], capture_output=True, text=True,
                            check=True, timeout=3).stdout
        page = int(re.search(r'page size of (\d+) bytes', vm).group(1))
        ram = page*sum(int(re.search(re.escape(k)+r':\s*(\d+)', vm).group(1))
                       for k in ('Pages free', 'Pages inactive', 'Pages speculative'))
        require(ram >= 3*GIB, 'Storage maintenance requires3GiB RAM')
        table = subprocess.run(['ps', '-axo', 'pgid=,rss='], capture_output=True,
                               text=True, check=True, timeout=3).stdout
        rss = sum(int(row[1])*1024 for line in table.splitlines()
                  if len(row := line.split()) == 2 and int(row[0]) == group)
        peak = max(peak, rss)
        require(rss < 256*1024**2, 'Storage-maintenance group RSS limit')
        return {'free_disk_bytes': disk, 'available_ram_bytes': ram}

    initial = guard(True)  # Monitoring access is mandatory before touching files.
    with (ROOT/'results/finite57_terminal_adapter_inputs_small_affine_20261025.json').open('rb') as inp:
        manifest = inp.read(65537)
    require(hashlib.sha256(manifest).hexdigest() ==
            '8a0973d7b7ef08231e199bee9f04f45c6a8dc01ba66ca8de8acb2d8e513ee7b7',
            'Current terminal manifest changed')
    require(not any(x['path'] == 'results/'+name+'.json'
                    for x in json.loads(manifest)['inputs'] for name in FILES),
            'A current terminal input must not be archived')
    archive = ROOT/'formal/.proof-batches/parked_sqrt5_json_archive_20261025'
    receipt = ROOT/'results/parked_sqrt5_json_streamed_archive_20261025.json'
    require(not archive.exists() and not receipt.exists(), 'New maintenance targets required')
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        sources = []
        for name, (size, digest) in FILES.items():
            p = ROOT/f'results/{name}.json'
            require(p.is_file() and not p.is_symlink() and p.stat().st_size == size,
                    'Exact unused historical-result path/size required')
            git = subprocess.run(['git', '-c', 'core.packedGitWindowSize=8m',
                '-c', 'core.packedGitLimit=16m', '-c', 'core.deltaBaseCacheLimit=4m',
                'ls-files', '--error-unmatch', '--', str(p.relative_to(ROOT))],
                cwd=ROOT, capture_output=True, text=True, timeout=3)
            require(git.returncode == 1 and not git.stdout,
                    'Only this untracked historical result may be archived')
            sources.append((name, p, size, digest, identity(p)))
        archive.mkdir()
        records = []
        verified_archives = []
        for name, source, size, expected, before in sources:
            target = archive/(name+'.json.gz')
            original_hash = hashlib.sha256()
            count = 0
            with source.open('rb') as inp, target.open('xb') as out:
                with gzip.GzipFile(filename='', mode='wb', fileobj=out,
                                   compresslevel=1, mtime=0) as packed:
                    for block in iter(lambda: inp.read(1024**2), b''):
                        guard()
                        original_hash.update(block)
                        count += len(block)
                        packed.write(block)
            require(count == size and original_hash.hexdigest() == expected,
                    'Original object content identity mismatch; original retained')
            packed_identity = identity(target)
            restored_hash = hashlib.sha256()
            count = 0
            with gzip.open(target, 'rb') as inp:
                for block in iter(lambda: inp.read(1024**2), b''):
                    guard()
                    restored_hash.update(block)
                    count += len(block)
            require(count == size and restored_hash.hexdigest() == expected,
                    'Gzip round trip failed; original retained')
            require(identity(target) == packed_identity, 'Archive changed during round trip')
            compressed_hash = hashlib.sha256()
            with target.open('rb') as inp:
                for block in iter(lambda: inp.read(1024**2), b''):
                    guard()
                    compressed_hash.update(block)
            require(identity(target) == packed_identity, 'Archive changed during hashing')
            verified_archives.append((target, packed_identity))
            require(identity(source) == before, 'Original object changed; retained')
            record = {'original': str(source.relative_to(ROOT)),
                      'archive': str(target.relative_to(ROOT)), 'bytes': size,
                      'sha256': expected, 'archive_bytes': target.stat().st_size,
                      'archive_sha256': compressed_hash.hexdigest(),
                      'roundtrip_verified': True}
            records.append(record)
        # Verify all archive streams before removing any original.
        # This immutable byte manifest also survives interruption during removal.
        (archive/'manifest.json').write_text(json.dumps(records, indent=2)+'\n')
        for name, source, size, expected, before in sources:
            guard(True)
            require(identity(source) == before, 'Object changed before removal')
        with (ROOT/'results/finite57_terminal_adapter_inputs_small_affine_20261025.json').open('rb') as inp:
            require(inp.read(65537) == manifest, 'Terminal manifest changed before removal')
        for target, verified in verified_archives:
            require(identity(target) == verified, 'Archive changed before original removal')
        for name, source, size, expected, before in sources:
            source.unlink()
        final = guard(True)
        result = {'status': 'BYTE_PRESERVING_UNUSED_HISTORICAL_JSON_ARCHIVE',
                  'records': records, 'initial': initial, 'final': final,
                  'peak_group_rss_bytes': peak,
                  'seconds': round(time.monotonic()-started, 3),
                  'original_raw_path_archived': True,
                  'scientific_byte_contents_preserved': True,
                  'current_terminal_manifest_110_inputs_unchanged': True,
                  'scientific_or_lean_execution': False,
                  'proof_admission_policy_changed': False,
                  'restore': 'Stream gzip decompression to its original results path; verify the recorded SHA256 before use. The current110-input terminal manifest does not use this historical result.',
                  'monitor_precision': 'Sampled256MiB group RSS and20GiB/3GiB storage floors every0.2s; each1MiB block checks owner/deadline. No mathematical admission policy changed.'}
        with receipt.open('x') as out:
            json.dump(result, out, indent=2)
            out.write('\n')
        print(json.dumps(result))


if __name__ == '__main__':
    main()
