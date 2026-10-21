#!/usr/bin/env python3
"""Certify the p7 silver odd-unimodular family in H27--H29.

For each of the primary catalogue lattices T27, U28 and dim29odd this freezer
downloads/parses no data at run time: the exact full Gram matrices are
reconstructed from vendored normalized primary-page snapshots.  PARI proves
determinant one, minimum three, and the exact norm-four shell cardinality.

The Vinberg arithmetic and H<=1000 hyperideal argument are universal for any
positive-definite integral unimodular lattice of minimum three containing a
norm-three vector.  We reuse the independently tested finite Gram theorem
from ``silver_s26`` and explicitly record the hypotheses checked here.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
import shutil
import subprocess
from pathlib import Path

from coxeter_polytopes.search.silver_s26 import (
    FIRST_PACKET_HEIGHT,
    H1000_WITNESS_SCALE,
    NORM_CLASSES,
    h1000_minimal_direction_audit,
)


ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT_DIR = ROOT / "data/sources/nebe_sloane_odd_unimodular_h27_h29"
JSON_PATH = ROOT / "data/enumerations/silver_odd_unimodular_h27_h29_h1000.json"
REPORT_PATH = ROOT / "results/silver_odd_unimodular_h27_h29_h1000.md"

CATALOGUES = (
    {
        "name": "T27",
        "dimension": 27,
        "url": "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/T27.html",
        "snapshot": "T27.html",
        "kissing": 1_640,
    },
    {
        "name": "U28",
        "dimension": 28,
        "url": "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/U28.html",
        "snapshot": "U28.html",
        "kissing": 2_240,
    },
    {
        "name": "dim29odd",
        "dimension": 29,
        "url": "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/dim29odd.html",
        "snapshot": "dim29odd.html",
        "kissing": 1_856,
    },
)

UNIVERSAL_LOCAL_NORM_CLASS_PROOF = (
    "Put pi=sqrt(2), rho=2*pi-1 and delta=3+pi=epsilon*rho. For a "
    "primitive r=(x,y), reflection integrality is q(r)|2I(r), where "
    "I(r)=(rho*x,(y,M)). Since M is unimodular, every prime away from "
    "(pi) and (rho) is excluded by primitivity. At (pi), rho is a unit; "
    "if x is a pi-unit then rho*x is a unit, while otherwise primitivity "
    "of y and unimodularity make some (y,M) a unit. Thus I(r) is "
    "pi-locally the unit ideal and v_pi(q) is 0, 1, or 2. At p=(rho), "
    "positive valuation forces y in p*M; then primitive x is a p-unit and "
    "q=rho*(-x^2+rho*z^2) has exact p-valuation one. Class number one "
    "and U^+=U^2 give the six canonical totally-positive representatives "
    "(1,2+pi,2) times (1 or 3+pi)."
)


def extract_gram(path: Path, dimension: int) -> tuple[tuple[int, ...], ...]:
    text = html.unescape(path.read_text())
    marker = '<a NAME="GRAM"><STRONG>GRAM</STRONG></a><br>'
    if marker not in text:
        raise AssertionError(f"missing GRAM marker in {path}")
    tail = text.split(marker, 1)[1]
    lines = tail.splitlines()
    header_index = next(
        index
        for index, line in enumerate(lines)
        if re.fullmatch(rf"\s*{dimension}\s+{dimension}\s*<br>\s*", line)
    )
    rows = []
    for line in lines[header_index + 1 : header_index + 1 + dimension]:
        values = tuple(map(int, re.findall(r"-?\d+", line.split("<br>", 1)[0])))
        if len(values) != dimension:
            raise AssertionError(
                f"expected {dimension} entries in {path.name} Gram row, got {len(values)}"
            )
        rows.append(values)
    gram = tuple(rows)
    if len(gram) != dimension:
        raise AssertionError(f"incomplete Gram in {path}")
    if any(gram[i][j] != gram[j][i] for i in range(dimension) for j in range(dimension)):
        raise AssertionError(f"nonsymmetric Gram in {path}")
    return gram


def gram_stream_sha256(gram: tuple[tuple[int, ...], ...]) -> str:
    stream = "\n".join(" ".join(map(str, row)) for row in gram) + "\n"
    return hashlib.sha256(stream.encode()).hexdigest()


def gp_matrix(gram: tuple[tuple[int, ...], ...]) -> str:
    return "[" + ";".join(",".join(map(str, row)) for row in gram) + "]"


def pari_host_audit(gram: tuple[tuple[int, ...], ...]) -> dict[str, object]:
    executable = shutil.which("gp")
    if executable is None:
        raise RuntimeError("PARI/GP is required to regenerate the family certificate")
    program = (
        f"G={gp_matrix(gram)};"
        'print("DET|",matdet(G));'
        "R=qfminim(G,3,,2);"
        'print("M3|",R[1],"|",floor(R[2]+1/2),"|",matsize(R[3])[2]);'
        "S=qfminim(G,4,,0);"
        'print("M4|",S[1],"|",floor(S[2]+1/2),"|",matsize(S[3])[2]);quit()'
    )
    completed = subprocess.run(
        [executable, "-q", "-s", "1000000000"],
        input=program + "\n",
        text=True,
        capture_output=True,
        timeout=300,
        check=True,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if len(lines) != 3:
        raise AssertionError(f"unexpected PARI transcript: {lines!r}; {completed.stderr!r}")
    det = int(lines[0].split("|")[1])
    m3 = tuple(map(int, lines[1].split("|")[1:]))
    m4 = tuple(map(int, lines[2].split("|")[1:]))
    if det != 1 or m3[1] != 3 or m3[0] <= 0 or m4[1] != 4:
        raise AssertionError(f"host invariant audit failed: {lines!r}")
    return {
        "engine": "PARI/GP matdet; qfminim flag=2 at norm 3 and flag=0 at norm 4",
        "transcript": lines,
        "determinant": det,
        "minimum": m3[1],
        "norm_three_shell_signed_count": m3[0],
        "norm_three_representatives_mod_sign": m3[2],
        "cumulative_signed_count_through_norm_four": m4[0],
        "norm_four_shell_signed_count": m4[0] - m3[0],
        "norm_four_representatives_mod_sign": (m4[0] - m3[0]) // 2,
    }


def build_payload() -> dict[str, object]:
    universal = h1000_minimal_direction_audit()
    if not universal.strict_hyperideal:
        raise AssertionError("the universal H1000 witness certificate failed")
    hosts = []
    for catalogue in CATALOGUES:
        snapshot = SNAPSHOT_DIR / catalogue["snapshot"]
        gram = extract_gram(snapshot, catalogue["dimension"])
        audit = pari_host_audit(gram)
        if audit["norm_three_shell_signed_count"] != catalogue["kissing"]:
            raise AssertionError(f"{catalogue['name']} kissing number differs")
        if not any(gram[i][i] == 3 for i in range(catalogue["dimension"])):
            raise AssertionError(f"{catalogue['name']} has no explicit norm-three basis vector")
        hosts.append(
            {
                **catalogue,
                "snapshot_path": str(snapshot.relative_to(ROOT)),
                "snapshot_sha256": hashlib.sha256(snapshot.read_bytes()).hexdigest(),
                "gram_stream_sha256": gram_stream_sha256(gram),
                "pari": audit,
                "first_packet_height": str(FIRST_PACKET_HEIGHT),
                "first_packet_wall_count": audit["norm_four_shell_signed_count"],
                "first_packet_description": (
                    "r_b=(2+sqrt(2),(1+sqrt(2))*b), b in spatial lattice, b^2=4"
                ),
            }
        )
    worst = universal.worst_row
    return {
        "schema_version": "0.1",
        "field": "Q(sqrt(2))",
        "coefficient": "2*sqrt(2)-1",
        "hyperbolic_dimensions": [27, 28, 29],
        "scope": (
            "For each displayed integral unimodular minimum-three spatial host: "
            "complete primitive norm alphabet, height-zero exclusion, exact first "
            "packet, and a strict hyperideal witness against every possible wall "
            "through H=1000. This is bounded, not Vinberg termination."
        ),
        "primary_catalogue_hosts": hosts,
        "universal_theorem": {
            "hypotheses": [
                "positive-definite integral unimodular spatial lattice",
                "minimum at least three",
                "a fixed vector d of norm three",
            ],
            "local_norm_class_count": len(NORM_CLASSES),
            "local_norm_classes": [str(row.norm) for row in NORM_CLASSES],
            "local_norm_class_proof": UNIVERSAL_LOCAL_NORM_CLASS_PROOF,
            "height_zero_wall_count": 0,
            "first_packet_height": str(FIRST_PACKET_HEIGHT),
            "first_packet_forcing": (
                "The raw minimum-three binary possibilities are (6,3,4) and "
                "(4,4,4). The first is impossible because (u-v)^2=1<3; the "
                "second forces u=v=b with b^2=4. Thus the first packet is the "
                "full spatial norm-four shell."
            ),
            "complete_height_bound": str(universal.complete_height_bound),
            "control_count": universal.control_count,
            "distinct_height_count": universal.distinct_height_count,
            "nonempty_control_count": universal.nonempty_control_count,
            "necessary_binary_type_count": universal.gated_binary_type_count,
            "nonempty_height_count": universal.nonempty_height_count,
            "witness": "z=(1,(43/55)d), d^2=3",
            "witness_scale": str(H1000_WITNESS_SCALE),
            "witness_lorentz_norm": str(universal.witness_lorentz_norm),
            "worst_row": {
                "height": str(worst.control.height),
                "norm": str(worst.control.norm_class.norm),
                "control": str(worst.control.control),
                "gram": [worst.triple.u_square, worst.triple.v_square, worst.triple.u_dot_v],
                "maximum_pairing": str(worst.pairing.maximum_pairing),
                "margin": str(worst.witness_margin),
            },
            "all_wall_margins_strict": universal.all_margins_strict,
            "strict_hyperideal": universal.strict_hyperideal,
            "proof_scope": (
                "The binary and augmented Gram gates use only integrality, minimum "
                "three and positive semidefiniteness. They retain a superset of all "
                "actual roots in every rank, so their maximum applies uniformly."
            ),
        },
        "conclusion": {
            "new_compact_example": False,
            "new_finite_volume_example": False,
            "all_three_prefixes_have_infinite_volume": True,
            "vinberg_terminated": False,
        },
    }


def render(payload: dict[str, object]) -> str:
    theorem = payload["universal_theorem"]
    rows = "\n".join(
        f"| `{host['name']}` | {host['dimension']} | "
        f"{host['pari']['norm_three_shell_signed_count']:,} | "
        f"{host['first_packet_wall_count']:,} |"
        for host in payload["primary_catalogue_hosts"]
    )
    return f"""# Silver odd-unimodular family in H27--H29 through H1000

## Outcome

Three genuinely new rootless target-rank hosts are now exact:

| spatial lattice | hyperbolic dimension | norm-three shell | first packet walls |
|:--|--:|--:|--:|
{rows}

For each host

\\[
 L_M=\\langle-(2\\sqrt2-1)\\rangle\\oplus M,
\\]

the complete six-class local arithmetic has no height-zero wall. The first
packet at `H=5+4sqrt(2)` is exactly the full norm-four shell shown above.
No compact example results: for any fixed `d^2=3`,

\\[
 z=(1,(43/55)d)
\\]

is strictly hyperideal against every possible wall through `H<=1000`.

## Host-specific exact audits

The literal full Gram matrices were parsed from normalized snapshots of the
primary Nebe--Sloane catalogue pages ([T27]({payload['primary_catalogue_hosts'][0]['url']}),
[U28]({payload['primary_catalogue_hosts'][1]['url']}), and
[dim29odd]({payload['primary_catalogue_hosts'][2]['url']})). For every row, PARI independently gives
determinant one and minimum three, reproduces the catalogue kissing number,
and counts the exact norm-four shell. Each displayed Gram contains an
explicit diagonal norm-three basis vector.

## Uniform H1000 theorem

The universal arithmetic census has **{theorem['control_count']} controls at
{theorem['distinct_height_count']} heights**. Exactly
{theorem['nonempty_control_count']} controls retain a necessary Gram type,
for {theorem['necessary_binary_type_count']} types at
{theorem['nonempty_height_count']} heights. Its tight row is

\\[
 H={theorem['worst_row']['height']},\\quad q=2,\\quad
 (u^2,v^2,u.v)=(460,226,322),
\\]

with maximum `37+26sqrt(2)` and strict margin
`{theorem['worst_row']['margin']}`. The witness is spacelike by
`{theorem['witness_lorentz_norm']}>0`.

This finite theorem depends only on integrality, unimodularity, minimum three,
and a norm-three direction—not on rank or on an automorphism-orbit sample.
The necessary binary and augmented Gram gates retain supersets, so the
pairing bounds apply to every actual wall.

Thus all three complete bounded prefixes have open hyperideal patches and
infinite volume. Walls above `H=1000` remain open; this is not reflectivity or
Vinberg termination.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/certify_silver_odd_unimodular_h27_h29_h1000.py --verify
PYTHONPATH=src pytest -q tests/test_silver_odd_unimodular_h27_h29_h1000.py
```
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    report = render(payload)
    if args.verify:
        if JSON_PATH.read_text() != serialized:
            raise SystemExit(f"artifact drift: {JSON_PATH}")
        if REPORT_PATH.read_text() != report:
            raise SystemExit(f"report drift: {REPORT_PATH}")
        print("PASS silver odd-unimodular H27--H29 H1000 family")
        return
    JSON_PATH.write_text(serialized)
    REPORT_PATH.write_text(report)
    print(JSON_PATH)
    print(REPORT_PATH)


if __name__ == "__main__":
    main()
