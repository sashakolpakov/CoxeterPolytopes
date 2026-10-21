"""Build and verify the solver-independent H9/14 record-23 exclusion."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record23 import (
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
    audit_record23_structure,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
BACKEND_SOURCE = ROOT / "data/enumerations/h9_14_facet_completion.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record23_certificate.json"
A0_SOURCE = ROOT / "scripts/certify_h9_14_record23_a0.cpp"
EXACT_POSITIVE_SOURCE = ROOT / "scripts/certify_h9_14_record23_positive.cpp"
RELAXED_POSITIVE_SOURCE = ROOT / "scripts/independent_record23_positive_a.cpp"
INDEPENDENT_A1_SOURCE = ROOT / "scripts/independent_record23_a1_audit.cpp"

EXPECTED_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_BACKEND_SOURCE_SHA256 = (
    "2d329439e5baeffb574932ff7bbab2dd78daf706b2917a0674bdc9f01a2f5123"
)
EXPECTED_GENBG_SHA256 = (
    "638224266e86d68d1cd720f18e8b7904d3571849f3da8c6130b1d0e58f515f83"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "7b13a6b215c4cf6c2c871361c984f150f044a7107bce5b5fb9f484096360da5c"
)

A0_GENBG_ARGUMENTS = (
    "-q", "-Y1", "-d7:0", "-D7:7", "7", "12", "49:49"
)
A1_GENBG_ARGUMENTS = (
    "-q", "-Y1", "-d7:0", "-D8:6", "6", "12", "43:43"
)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_digest(manifest: Mapping[str, object]) -> str:
    critical = dict(manifest)
    critical.pop("certificate_sha256", None)
    return sha256(
        json.dumps(critical, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _source_record(source: Mapping[str, object]) -> Mapping[str, object]:
    records = source.get("records")
    if not isinstance(records, list):
        raise ValueError("nerve-support source lacks records")
    survivors = [
        record
        for record in records
        if isinstance(record, Mapping) and record.get("survives_second_layer")
    ]
    record = survivors[TARGET_RECORD_INDEX]
    if record.get("graph6") != TARGET_GRAPH6:
        raise ValueError("survivor index 23 has changed")
    return record


def _backend_corroboration() -> dict[str, object]:
    backend = json.loads(BACKEND_SOURCE.read_text(encoding="utf-8"))
    matching = [
        record
        for record in backend["records"]
        if record.get("index") == TARGET_RECORD_INDEX
    ]
    if len(matching) != 1 or matching[0]["solver"]["status"] != "unsat":
        raise ValueError("record-23 backend corroboration changed")
    solver = matching[0]["solver"]
    return {
        "path": "data/enumerations/h9_14_facet_completion.json",
        "sha256": _file_sha256(BACKEND_SOURCE),
        "status": solver["status"],
        "backend": solver["backend"],
        "backend_version": solver["backend_version"],
        "assertion_count": solver["assertion_count"],
        "trusted_for_proof": False,
    }


def _parse_histogram(text: str) -> dict[str, int]:
    result: dict[str, int] = {}
    for item in text.rstrip(",").split(","):
        if not item:
            continue
        key, value = item.split(":")
        result[key] = int(value)
    return result


def _parse_stream_census(output: str) -> dict[str, object]:
    record: dict[str, object] = {
        "raw_output_sha256": sha256(output.encode()).hexdigest(),
        "facet_count_histogram": {},
    }
    histogram = record["facet_count_histogram"]
    if not isinstance(histogram, dict):
        raise AssertionError("unreachable histogram type")
    vector_fields = {
        "maximizing_left_vertices",
        "maximizing_complement_columns",
    }
    for line in output.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "histogram":
            histogram[parts[1]] = int(parts[2])
        elif parts[0] in vector_fields:
            record[parts[0]] = [int(value) for value in parts[1:]]
        elif len(parts) == 2:
            record[parts[0]] = int(parts[1])
        else:
            raise ValueError(f"unparsed stream-census line: {line}")
    return record


def _parse_independent_a1(output: str) -> dict[str, object]:
    record: dict[str, object] = {
        "raw_output_sha256": sha256(output.encode()).hexdigest(),
        "facet_count_histogram": {},
    }
    histogram = record["facet_count_histogram"]
    if not isinstance(histogram, dict):
        raise AssertionError("unreachable histogram type")
    for line in output.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "histogram":
            histogram[parts[1]] = int(parts[2])
        elif parts[0] == "best_B":
            record[parts[0]] = [int(value) for value in parts[1:]]
        elif len(parts) == 2:
            record[parts[0]] = int(parts[1])
        else:
            raise ValueError(f"unparsed independent-a1 line: {line}")
    return record


def _parse_relaxed_positive(output: str) -> dict[str, object]:
    active_types: dict[int, dict[str, object]] = {}
    later: dict[int, dict[str, object]] = {}
    for line in output.splitlines():
        if not line:
            continue
        active_match = re.fullmatch(
            r"a=1 active_size=(\d+) isolate_family_size=(\d+) "
            r"isolate_orbits=(\d+)",
            line,
        )
        if active_match:
            active_size, family_size, count = map(int, active_match.groups())
            entry = active_types.setdefault(active_size, {"active_size": active_size})
            entry.setdefault("isolate_orbits", {})[str(family_size)] = count
            continue
        active_hist_match = re.fullmatch(
            r"a=1 active_size=(\d+) valid_families=(\d+) "
            r"target_96=(\d+) histogram=(.*)",
            line,
        )
        if active_hist_match:
            active_size, valid, target = map(int, active_hist_match.groups()[:3])
            entry = active_types.setdefault(active_size, {"active_size": active_size})
            entry.update(
                valid_family_count=valid,
                target_facet_count_occurrences=target,
                facet_count_histogram=_parse_histogram(active_hist_match.group(4)),
            )
            continue
        maximum_match = re.fullmatch(
            r"a=1 active_size=(\d+) maximum=(\d+) q=(\d+) "
            r"maximum_family_masks=(.*)",
            line,
        )
        if maximum_match:
            active_size, maximum, q = map(int, maximum_match.groups()[:3])
            entry = active_types.setdefault(active_size, {"active_size": active_size})
            entry.update(
                maximum_facet_count=maximum,
                missing_four_representative_mask=q,
                maximum_family_masks=[
                    int(value)
                    for value in maximum_match.group(4).rstrip(",").split(",")
                    if value
                ],
            )
            continue
        orbit_match = re.fullmatch(
            r"a=(\d+) a4_orbits=(\d+) "
            r"shared_isolate_triple_pair_histogram=(.*)",
            line,
        )
        if orbit_match:
            a, count = map(int, orbit_match.groups()[:2])
            later.setdefault(a, {}).update(
                missing_four_count=a,
                a4_orbit_count=count,
                shared_isolate_triple_pair_histogram=_parse_histogram(
                    orbit_match.group(3)
                ),
            )
            continue
        later_match = re.fullmatch(
            r"a=(\d+) valid_families=(\d+) target_(\d+)=(\d+) "
            r"histogram=(.*)",
            line,
        )
        if later_match:
            a, valid, target_value, target_count = map(
                int, later_match.groups()[:4]
            )
            later.setdefault(a, {}).update(
                valid_family_count=valid,
                required_facet_count=target_value,
                target_facet_count_occurrences=target_count,
                facet_count_histogram=_parse_histogram(later_match.group(5)),
            )
            continue
        raise ValueError(f"unparsed relaxed-positive line: {line}")
    return {
        "raw_output_sha256": sha256(output.encode()).hexdigest(),
        "a1_active_types": [active_types[key] for key in sorted(active_types)],
        "a1_maximum_facet_count": max(
            int(entry["maximum_facet_count"]) for entry in active_types.values()
        ),
        "a1_target_facet_count_occurrences": sum(
            int(entry["target_facet_count_occurrences"])
            for entry in active_types.values()
        ),
        "later_branches": [later[key] for key in sorted(later)],
    }


def _compile(source: Path, output: Path) -> None:
    subprocess.run(
        ["clang++", "-O3", "-std=c++20", str(source), "-o", str(output)],
        check=True,
    )


def _run_genbg_filter(arguments: Sequence[str], executable: Path) -> str:
    producer = subprocess.Popen(
        ["genbg", *arguments], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if producer.stdout is None:
        raise AssertionError("genbg stdout pipe was not created")
    consumer = subprocess.run(
        [str(executable)], stdin=producer.stdout, capture_output=True, text=True
    )
    producer.stdout.close()
    stderr = producer.stderr.read().decode() if producer.stderr else ""
    producer_status = producer.wait()
    if producer_status or consumer.returncode:
        raise RuntimeError(
            f"census pipeline failed: genbg={producer_status}, "
            f"filter={consumer.returncode}, stderr={stderr}{consumer.stderr}"
        )
    return consumer.stdout


def recompute_censuses() -> dict[str, object]:
    genbg = shutil.which("genbg")
    if genbg is None:
        raise RuntimeError("recomputation requires nauty genbg")
    if _file_sha256(Path(genbg)) != EXPECTED_GENBG_SHA256:
        raise RuntimeError("recomputation requires the frozen genbg executable")
    with tempfile.TemporaryDirectory(prefix="record23-census-") as directory:
        temporary = Path(directory)
        a0_binary = temporary / "a0"
        exact_positive_binary = temporary / "positive"
        relaxed_binary = temporary / "relaxed"
        independent_a1_binary = temporary / "independent-a1"
        _compile(A0_SOURCE, a0_binary)
        _compile(EXACT_POSITIVE_SOURCE, exact_positive_binary)
        _compile(RELAXED_POSITIVE_SOURCE, relaxed_binary)
        _compile(INDEPENDENT_A1_SOURCE, independent_a1_binary)
        a0_output = _run_genbg_filter(A0_GENBG_ARGUMENTS, a0_binary)

        producer = subprocess.Popen(
            ["genbg", *A1_GENBG_ARGUMENTS],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if producer.stdout is None:
            raise AssertionError("genbg stdout pipe was not created")
        consumer = subprocess.run(
            [str(exact_positive_binary), "1", "5", "96"],
            stdin=producer.stdout,
            capture_output=True,
            text=True,
        )
        producer.stdout.close()
        stderr = producer.stderr.read().decode() if producer.stderr else ""
        producer_status = producer.wait()
        if producer_status or consumer.returncode:
            raise RuntimeError(
                f"a1 census failed: genbg={producer_status}, "
                f"filter={consumer.returncode}, stderr={stderr}{consumer.stderr}"
            )
        relaxed = subprocess.run(
            [str(relaxed_binary)], check=True, capture_output=True, text=True
        ).stdout
        independent_a1_output = _run_genbg_filter(
            A1_GENBG_ARGUMENTS, independent_a1_binary
        )
    return {
        "a0_exact_genbg": _parse_stream_census(a0_output),
        "a1_exact_genbg": _parse_stream_census(consumer.stdout),
        "a1_independent_replay": _parse_independent_a1(independent_a1_output),
        "positive_a_relaxed_independent": _parse_relaxed_positive(relaxed),
    }


def _enumerator_provenance() -> dict[str, object]:
    return {
        "generator": {
            "name": "nauty genbg",
            "path_at_freeze": "/opt/homebrew/bin/genbg",
            "sha256_at_freeze": EXPECTED_GENBG_SHA256,
            "role": (
                "isomorph-free generation of bicoloured incidence graphs; "
                "all command-line properties are rechecked by the filters"
            ),
        },
        "sources": [
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": _file_sha256(path),
            }
            for path in (
                A0_SOURCE,
                EXACT_POSITIVE_SOURCE,
                RELAXED_POSITIVE_SOURCE,
                INDEPENDENT_A1_SOURCE,
            )
        ],
        "commands": {
            "a0": (
                "genbg -q -Y1 -d7:0 -D7:7 7 12 49:49 | "
                "certify_h9_14_record23_a0"
            ),
            "a1_exact": (
                "genbg -q -Y1 -d7:0 -D8:6 6 12 43:43 | "
                "certify_h9_14_record23_positive 1 5 96"
            ),
            "positive_relaxed": "independent_record23_positive_a",
        },
    }


def build(*, recompute: bool = True) -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    structure = audit_record23_structure(
        TARGET_GRAPH6, record["dehn_sommerville_schedule"]
    ).to_record()
    if not recompute:
        raise ValueError("building the record-23 archive requires recomputation")
    censuses = recompute_censuses()
    a0 = censuses["a0_exact_genbg"]
    a1 = censuses["a1_exact_genbg"]
    a1_replay = censuses["a1_independent_replay"]
    relaxed = censuses["positive_a_relaxed_independent"]
    if (
        a0["generated_graph_count"] != 2_807_714
        or a0["distinct_row_graph_count"] != 2_765_889
        or a0["locally_valid_left_choice_count"] != 11_456_748
        or a0["target_facet_count_occurrences"] != 0
        or a0["maximum_facet_count"] != 94
    ):
        raise AssertionError("a=0 census did not reproduce")
    if (
        a1["generated_graph_count"] != 75_890
        or a1["distinct_row_graph_count"] != 75_016
        or a1["schedule_compatible_graph_count"] != 73_919
        or a1["blocker_valid_left_choice_count"] != 86_241
        or a1["target_facet_count_occurrences"] != 0
        or a1["maximum_facet_count"] != 90
    ):
        raise AssertionError("a=1 exact census did not reproduce")
    if (
        a1_replay["generated"] != 75_890
        or a1_replay["structural"] != 73_919
        or a1_replay["left_choices"] != 86_241
        or a1_replay["target"] != 0
        or a1_replay["maximum"] != 90
        or a1_replay["facet_count_histogram"]
        != a1["facet_count_histogram"]
    ):
        raise AssertionError("independent a=1 replay did not reproduce")
    later = {item["missing_four_count"]: item for item in relaxed["later_branches"]}
    if (
        relaxed["a1_target_facet_count_occurrences"] != 0
        or relaxed["a1_maximum_facet_count"] != 92
        or later[2]["a4_orbit_count"] != 17
        or later[2]["target_facet_count_occurrences"] != 0
        or later[3]["a4_orbit_count"] != 144
        or later[3]["target_facet_count_occurrences"] != 0
    ):
        raise AssertionError("positive-a relaxed census did not reproduce")

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record23_solver_independent_exclusion",
        "claim": (
            "H9/14 survivor record 23 cannot be the nerve of a compact "
            "hyperbolic Coxeter 9-polytope."
        ),
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_unsat",
        "proof_backend": (
            "finite integer and isomorph-free incidence enumeration; "
            "no SAT/SMT backend"
        ),
        "source": {
            "path": "data/enumerations/h9_14_nerve_support.json",
            "sha256": _file_sha256(SOURCE),
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
        },
        "proof_dependencies": [
            "every compact Coxeter missing face is Lanner of order at most five",
            "strict generalized lower-bound inequalities in nonsimplex links",
            "the g-theorem Macaulay inequalities for vertex, edge, and triple links",
            "Dehn--Sommerville face-vector schedules",
            (
                "maximal faces are complementary to inclusion-minimal "
                "transversals of the complete minimal-nonface clutter"
            ),
            "McKay's isomorph-free bicoloured graph generation algorithm (genbg)",
        ],
        "structure": structure,
        "blocker_formula": {
            "usable_vertices": 12,
            "minimal_graph_covers": [[0, 1, 2], [12, 13]],
            "facet_count": "tau_3(all A4/A5) + tau_2(A5 disjoint from {0,1,2})",
            "purity_gate": (
                "tau(all A4/A5)>=3 and tau(A5 disjoint from left)>=2; "
                "then every five-transversal is minimal"
            ),
        },
        "censuses": censuses,
        "enumerator_provenance": _enumerator_provenance(),
        "backend_corroboration": _backend_corroboration(),
        "updated_theorem_partition": {
            "independently_certified_impossible": 117,
            "certified_explicit_witness": 2,
            "theorem_undecided": 23,
            "total_graph_shadows": 142,
        },
        "nonclaims": [
            "This excludes one graph shadow, not all compact H9/14 nerves.",
            "The archived Z3 UNSAT result is corroboration and is not trusted.",
            "The bicoloured graph generator is an enumerator, not a SAT solver.",
            "Expanded left-colour choices may repeat isomorphic coloured states; this is safe for absence.",
            "This is not a general non-existence theorem in dimension 9.",
        ],
    }
    manifest["certificate_sha256"] = _canonical_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object], *, recompute: bool = False) -> dict[str, object]:
    if manifest.get("schema_version") != "1.0.0" or manifest.get("kind") != (
        "h9_14_record23_solver_independent_exclusion"
    ):
        raise ValueError("unexpected record-23 certificate metadata")
    if manifest.get("record_index") != TARGET_RECORD_INDEX or manifest.get(
        "graph6"
    ) != TARGET_GRAPH6:
        raise ValueError("record-23 target mismatch")
    if manifest.get("theorem_status") != "certified_unsat":
        raise ValueError("record-23 theorem status was not promoted")
    if _file_sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise ValueError("nerve-support source is not frozen")
    if _file_sha256(BACKEND_SOURCE) != EXPECTED_BACKEND_SOURCE_SHA256:
        raise ValueError("backend corroboration source is not frozen")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    expected_structure = audit_record23_structure(
        TARGET_GRAPH6, record["dehn_sommerville_schedule"]
    ).to_record()
    if manifest.get("structure") != expected_structure:
        raise ValueError("record-23 structural proof does not verify")
    if manifest.get("backend_corroboration") != _backend_corroboration():
        raise ValueError("record-23 backend corroboration mismatch")
    backend = manifest["backend_corroboration"]
    if not isinstance(backend, Mapping) or backend.get("trusted_for_proof") is not False:
        raise ValueError("backend corroboration was incorrectly promoted")
    provenance = manifest.get("enumerator_provenance")
    if provenance != _enumerator_provenance():
        raise ValueError("record-23 enumerator source provenance mismatch")
    digest = _canonical_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-23 certificate digest mismatch")
    if digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("record-23 certificate is not frozen")
    if recompute:
        expected_censuses = recompute_censuses()
        if manifest.get("censuses") != expected_censuses:
            raise ValueError("record-23 incidence censuses do not reproduce")
    return dict(expected_structure)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args(argv)
    if args.verify_only:
        manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
        proof = verify(manifest, recompute=args.recompute)
    else:
        manifest = build(recompute=True)
        OUTPUT.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        proof = verify(manifest, recompute=False)
    print(
        json.dumps(
            {
                "output": str(OUTPUT),
                "certificate_sha256": manifest["certificate_sha256"],
                "verified": proof["verified"],
                "recomputed": bool(args.recompute or not args.verify_only),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
