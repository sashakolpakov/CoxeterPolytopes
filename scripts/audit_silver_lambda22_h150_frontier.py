#!/usr/bin/env python3
"""Freeze the exact silver Lambda22 witness frontier through real height 150.

This is a continuation of ``audit_silver_lambda22_h14.py``.  It does not
enumerate every spatial shell.  Instead it enumerates every reflective
control, every necessary binary coefficient Gram type, and every PSD
extension by the fixed minimal direction ``d=e_1``.  Only seven types after
the frozen H14 packet can possibly meet the ray ``z=(1,d/3)``.  Exact gate,
minimum, residue, and PARI shell arguments exclude all seven.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from hashlib import sha256
from math import lcm
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_lambda22 import (
    FIRST_PACKET_HEIGHT,
    GRAM_FILE_SHA256,
    LAMBDA22_COEFFICIENT,
    SOURCE_URL,
    enumerate_lambda22_controls,
    gated_relaxed_witness_bound,
    lambda22_gram,
    lambda22_pairing_radical_basis,
    witness_norm,
)
from coxeter_polytopes.search.silver_leech import candidate_leech_gram_triples
from coxeter_polytopes.search.silver_leech_manifest import (
    has_short_binary_combination,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_lambda22_h150_frontier.json"
REPORT_PATH = ROOT / "results/silver_lambda22_h150_frontier.md"
H14_JSON_PATH = ROOT / "data/enumerations/silver_lambda22_h14.json"
HEIGHT_BOUND = Silver(150)
MINIMUM_BY_PAIRING_MODULUS = {1: 4, 2: 12, 3: 24, 6: 96}


def value_json(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(0), value.sign(1)],
    }


def gate_data(level: object) -> tuple[int, int, int, int]:
    """Return exact pairing moduli/minima for ``u,v`` in ``y=u+sqrt(2)v``."""

    p_modulus = q_modulus = 1
    if level.dyadic_valuation == 3:
        p_modulus = 2
    elif level.dyadic_valuation == 4:
        p_modulus = q_modulus = 2
    if level.ternary_gate:
        p_modulus = lcm(p_modulus, 3)
        q_modulus = lcm(q_modulus, 3)
    return (
        p_modulus,
        q_modulus,
        MINIMUM_BY_PAIRING_MODULUS[p_modulus],
        MINIMUM_BY_PAIRING_MODULUS[q_modulus],
    )


def relaxed_frontier_scan() -> dict[str, object]:
    """Return a complete necessary-superset scan through ``H<=150``.

    If ``u`` pairs divisibly by ``m`` and ``v`` pairs divisibly by ``n``,
    then ``u.v`` is divisible by ``lcm(m,n)``.  Using a gcd here would be a
    genuine bug; in particular, the apparent H72 tangent has dot 525 and is
    excluded because its required moduli are 6 and 3.
    """

    controls = enumerate_lambda22_controls(HEIGHT_BOUND)
    retained = []
    dangerous = []
    retained_levels: set[tuple[str, str, str]] = set()
    h72_cross_gate_rejections = []
    for level in controls:
        p_modulus, q_modulus, u_minimum, v_minimum = gate_data(level)
        cross_modulus = lcm(p_modulus, q_modulus)
        for item in candidate_leech_gram_triples(level.target):
            if has_short_binary_combination(item):
                continue
            u_square = item.u_square
            v_square = item.v_square
            u_dot_v = item.u_dot_v
            if (
                (u_square and u_square < u_minimum)
                or (v_square and v_square < v_minimum)
                or u_square % p_modulus
                or v_square % q_modulus
            ):
                continue
            if u_dot_v % cross_modulus:
                if level.height == Silver(72, 51):
                    h72_cross_gate_rejections.append(
                        [u_square, v_square, u_dot_v]
                    )
                continue
            try:
                d_dot_u, d_dot_v, spatial_maximum = gated_relaxed_witness_bound(
                    u_square,
                    v_square,
                    u_dot_v,
                    p_modulus,
                    q_modulus,
                )
            except AssertionError:
                continue
            threshold = 3 * LAMBDA22_COEFFICIENT * level.control
            margin = spatial_maximum - threshold
            row = {
                "height": value_json(level.height),
                "norm": value_json(level.norm),
                "control": value_json(level.control),
                "target": value_json(level.target),
                "gram": [u_square, v_square, u_dot_v],
                "pairing_gate": {
                    "u_modulus": p_modulus,
                    "v_modulus": q_modulus,
                    "cross_modulus": cross_modulus,
                    "u_minimum": u_minimum,
                    "v_minimum": v_minimum,
                },
                "relaxed_maximum": {
                    "d_dot_u": d_dot_u,
                    "d_dot_v": d_dot_v,
                    "spatial": value_json(spatial_maximum),
                    "three_times_lorentz_threshold": value_json(threshold),
                    "margin": value_json(margin),
                },
            }
            retained.append(row)
            retained_levels.add(
                (str(level.height), str(level.norm), str(level.control))
            )
            if level.height > FIRST_PACKET_HEIGHT and margin.sign(0) >= 0:
                dangerous.append(row)

    expected_dangers = [
        ("21 + 15*sqrt(2)", [60, 36, 45]),
        ("36 + 51/2*sqrt(2)", [354, 198, 264]),
        ("41 + 29*sqrt(2)", [88, 40, 58]),
        ("42 + 30*sqrt(2)", [228, 144, 180]),
        ("56 + 40*sqrt(2)", [124, 52, 80]),
        ("56 + 40*sqrt(2)", [100, 64, 80]),
        ("70 + 99/2*sqrt(2)", [230, 126, 170]),
    ]
    actual_dangers = [
        (row["height"]["display"], row["gram"]) for row in dangerous
    ]
    if actual_dangers != expected_dangers:
        raise AssertionError(f"the H150 dangerous frontier changed: {actual_dangers!r}")
    if len(controls) != 58 or len(retained_levels) != 18 or len(retained) != 72:
        raise AssertionError(
            "the H150 necessary-superset census changed: "
            f"{len(controls)}, {len(retained_levels)}, {len(retained)}"
        )
    expected_h72 = [
        [780, 354, 525],
        [768, 360, 525],
        [756, 366, 525],
        [744, 372, 525],
        [732, 378, 525],
        [720, 384, 525],
        [708, 390, 525],
    ]
    if h72_cross_gate_rejections != expected_h72:
        raise AssertionError("the H72 cross-gate rejection list changed")
    return {
        "control_count": len(controls),
        "retained_level_count": len(retained_levels),
        "retained_necessary_type_count": len(retained),
        "retained_types": retained,
        "post_H14_nonnegative_type_count": len(dangerous),
        "post_H14_nonnegative_types": dangerous,
        "H72_cross_gate_rejections": h72_cross_gate_rejections,
    }


def radical_norm_residues_mod8() -> dict[str, object]:
    """Enumerate the four classes of ``S2/2Lambda22`` and their norm residues."""

    gram = lambda22_gram()
    basis = lambda22_pairing_radical_basis(2)
    generators = {
        tuple(basis[i][j] % 2 for i in range(22)) for j in range(22)
    }
    classes = {(0,) * 22}
    for generator in generators:
        classes |= {
            tuple((left[i] + generator[i]) % 2 for i in range(22))
            for left in tuple(classes)
        }
    rows = []
    for vector in sorted(classes):
        norm = sum(
            vector[i] * gram[i][j] * vector[j]
            for i in range(22)
            for j in range(22)
        )
        rows.append({"representative": list(vector), "norm_mod_8": norm % 8})
    residues = [row["norm_mod_8"] for row in rows]
    if len(rows) != 4 or sorted(residues) != [0, 4, 4, 4]:
        raise AssertionError(f"the S2 radical norm residues changed: {residues}")
    return {"class_count": 4, "classes": rows, "residue_multiset": sorted(residues)}


def _pari_matrix(matrix: object) -> str:
    return "[" + ";".join(",".join(map(str, row)) for row in matrix) + "]"


def pari_direction_shell_audit() -> dict[str, object]:
    """Enumerate the only short gated shells used by the frontier exclusions."""

    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP is required to regenerate the H150 frontier")
    gram_text = _pari_matrix(lambda22_gram())
    outputs = {}
    for modulus, bound, norms in ((2, 16, (12, 16)), (6, 96, (96,))):
        basis_text = _pari_matrix(lambda22_pairing_radical_basis(modulus))
        lines = [
            f"G={gram_text}",
            f"B={basis_text}",
            "H=B~*G*B",
            "d=vector(22,i,i==1)~",
            "ell=d~*G*B",
            f"R=qfminim(H,{bound},,2)",
            "M=R[3]",
            'print("signed_count=",R[1])',
            'print("maximum=",floor(R[2]+1/2))',
            'print("representatives=",matsize(M)[2])',
        ]
        for norm in norms:
            lines.extend(
                [
                    "P=List()",
                    (
                        "for(j=1,matsize(M)[2],"
                        f"if(M[,j]~*H*M[,j]=={norm},"
                        "listput(P,abs(ell*M[,j]))))"
                    ),
                    f'print("norm_{norm}_representatives=",#P)',
                    f'print("norm_{norm}_abs_pairings=",Set(Vec(P)))',
                ]
            )
        lines.append("quit()")
        process = subprocess.run(
            [gp, "-q", "-s", "512000000"],
            input=";".join(lines) + "\n",
            text=True,
            capture_output=True,
            timeout=180,
            check=True,
        )
        parsed: dict[str, object] = {}
        for raw in process.stdout.splitlines():
            key, value = raw.strip().split("=", 1)
            if key.endswith("abs_pairings"):
                parsed[key] = json.loads(value)
            else:
                parsed[key] = int(value)
        outputs[str(modulus)] = parsed

    expected = {
        "2": {
            "signed_count": 58_344,
            "maximum": 16,
            "representatives": 29_172,
            "norm_12_representatives": 4_224,
            "norm_12_abs_pairings": [0, 2, 4],
            "norm_16_representatives": 24_948,
            "norm_16_abs_pairings": [0, 2, 4, 8],
        },
        "6": {
            "signed_count": 1_782,
            "maximum": 96,
            "representatives": 891,
            "norm_96_representatives": 891,
            "norm_96_abs_pairings": [0, 6, 12],
        },
    }
    if outputs != expected:
        raise AssertionError(f"the PARI direction-shell census changed: {outputs!r}")
    return {
        "algorithm": "PARI qfminim(B^T G B,bound,,flag=2), d=e1",
        "pairing_radical_shells": outputs,
    }


def exact_exclusions(scan: dict[str, object]) -> list[dict[str, object]]:
    """Return the seven exact resolutions in frontier order."""

    dangers = scan["post_H14_nonnegative_types"]
    if len(dangers) != 7:
        raise AssertionError("expected seven relaxed frontier threats")
    return [
        {
            "height": dangers[0]["height"],
            "gram": dangers[0]["gram"],
            "resolution": (
                "The complete q=6 type list has u,v in S3 and nonzero u-v "
                "of norm 4,6,...,16, always below min(S3)=24. The dangerous "
                "type itself has norm 60+36-2*45=6."
            ),
        },
        {
            "height": dangers[1]["height"],
            "gram": dangers[1]["gram"],
            "resolution": (
                "All 17 possible types have U+2V=750 and u.v=264, with "
                "u in S6 and v in S3. Thus w=u-v in S3 has w^2=222-V, "
                "while s=u-2v in S6 has s^2=2V-306. For V<=198, "
                "0<s^2<=90<96; for V>=200, 0<w^2<=22<24."
            ),
        },
        {
            "height": dangers[2]["height"],
            "gram": dangers[2]["gram"],
            "resolution": (
                "Tangency requires w=u-v in S2 with w^2=12 and d.w=6. "
                "The complete S2 norm-12 shell has |d.w| in {0,2,4}."
            ),
        },
        {
            "height": dangers[3]["height"],
            "gram": dangers[3]["gram"],
            "resolution": (
                "Here u,v lie in S6. The nonzero vector u-v lies in S6 and "
                "has norm 228+144-2*180=12<96. The complete q=12 packet is "
                "likewise killed by S6 minimum; its difference norms are 6..48."
            ),
        },
        {
            "height": dangers[4]["height"],
            "gram": dangers[4]["gram"],
            "resolution": (
                "For the relaxed maximizer, s=u-2v lies in S2 with s^2=12 "
                "and d.s=-6, excluded by the complete S2 norm-12 pairing spectrum."
            ),
        },
        {
            "height": dangers[5]["height"],
            "gram": dangers[5]["gram"],
            "resolution": (
                "The nonzero vector u-v lies in S2 and has norm "
                "100+64-2*80=4<min(S2)=12."
            ),
        },
        {
            "height": dangers[6]["height"],
            "gram": dangers[6]["gram"],
            "resolution": (
                "The entire q=4+2sqrt(2) packet requires u in S2, but its six "
                "possible U values 250,246,242,238,234,230 are 2 or 6 mod 8. "
                "The four classes of S2/2Lambda22 have norm residues 0,4,4,4."
            ),
        },
    ]


def build_payload() -> dict[str, object]:
    scan = relaxed_frontier_scan()
    shells = pari_direction_shell_audit()
    radical = radical_norm_residues_mod8()
    d = (1,) + (0,) * 21
    h14_sha = sha256(H14_JSON_PATH.read_bytes()).hexdigest()
    return {
        "schema_version": "0.1.0",
        "seed_id": "silver-H22-Lambda22-a-sqrt2-1",
        "scope": (
            "Complete exact necessary-superset witness audit through real Vinberg "
            "height H<=150, continuing the frozen H14 prefix. It proves that the "
            "fixed hyperideal ray survives this cutoff, not global nonreflectivity."
        ),
        "field": "Q(sqrt(2))",
        "hyperbolic_dimension": 22,
        "lattice": "<-(sqrt(2)-1)> + Lambda22",
        "lambda22_source": SOURCE_URL,
        "gram_file_sha256": GRAM_FILE_SHA256,
        "base_H14_artifact": {
            "path": str(H14_JSON_PATH.relative_to(ROOT)),
            "sha256": h14_sha,
        },
        "height_bound": value_json(HEIGHT_BOUND),
        "fixed_direction": "d=e1, d^2=4",
        "hyperideal_witness": {
            "vector": "z=(1,d/3)",
            "norm": value_json(witness_norm(d)),
        },
        "gate_rule": (
            "For y=u+sqrt(2)v, if Gu is divisible by m and Gv by n, then "
            "u.v is divisible by lcm(m,n)."
        ),
        "necessary_superset_scan": scan,
        "S2_radical_norm_residues_mod8": radical,
        "pari_direction_shell_audit": shells,
        "exact_resolutions": exact_exclusions(scan),
        "redundant_H72_check": {
            "primary_resolution": (
                "All 22 necessary binary types have u.v=525, but u in S6 and "
                "v in S3 force 6|u.v; hence the entire level is empty."
            ),
            "rejected_types": scan["H72_cross_gate_rejections"],
            "shell_cross_check": (
                "The former tangent parametrization would require an S6 norm-96 "
                "vector pairing 18 with d; the complete shell spectrum is only "
                "{0,6,12}."
            ),
        },
        "conclusion": (
            "Every possible post-H14 root through H<=150 pairs strictly negatively "
            "with z, so the exact prefix through this cutoff still has infinite volume."
        ),
        "discovery": None,
        "next_target": (
            "At H=82+58sqrt(2) the apparent q=4 tangent again reduces to the "
            "excluded S2 norm-12 pairing 6. The next unresolved relaxed cutter is "
            "H=103+73sqrt(2), q=4, type (208,104,146)."
        ),
    }


def report(payload: dict[str, object]) -> str:
    scan = payload["necessary_superset_scan"]
    return f"""# Silver Lambda22 strict hyperideal frontier through H<=150

For the exact arithmetic host

`<-(sqrt(2)-1)> + Lambda22`

in hyperbolic dimension 22, the fixed ray `z=(1,e1/3)` remains strictly
hyperideal after every possible reflective root of real Vinberg height at
most 150.  Its norm is `13/9-sqrt(2)>0`.  Thus this **bounded prefix** still
has infinite volume; this is not a global nonreflectivity statement.

## Complete necessary-superset scan

The exact scan has {scan['control_count']} control levels.  After spatial
minimum and pairing-gate filters, {scan['retained_necessary_type_count']}
binary Gram types at {scan['retained_level_count']} levels remain.  Exact
PSD extension by `d=e1`, with all gate divisibilities imposed, leaves only
seven post-H14 types whose relaxed wall pairing is nonnegative.

They are excluded as follows:

1. `H=21+15sqrt(2)`: `u-v` is a nonzero S3 vector of norm 16, below 24.
2. `H=36+51sqrt(2)/2`: for all 17 types, either `u-v in S3` has norm at
   most 22 or `u-2v in S6` has norm at most 90.
3. `H=41+29sqrt(2)`: tangency needs an S2 norm-12 vector pairing 6 with
   `e1`; the complete shell has absolute pairings only `0,2,4`.
4. `H=42+30sqrt(2)`: the dangerous S6 difference has norm 12; in fact the
   whole q=12 packet has positive difference norm between 6 and 48.
5. `H=56+40sqrt(2)`, type `(124,52,80)`: `u-2v` would be an S2 norm-12
   vector pairing -6 with `e1`, excluded by the same shell census.
6. `H=56+40sqrt(2)`, type `(100,64,80)`: `u-v in S2` has norm 4, below 12.
7. `H=70+99sqrt(2)/2`: `u in S2`, but every possible `u^2` is 2 or 6
   modulo 8, whereas the four S2 radical classes have residues `0,4,4,4`.

An apparent H72 tangent disappears once the cross gate is imposed correctly:
`u in S6` and `v in S3` force `6 | u.v`, whereas every necessary type has
`u.v=525`.  As a redundant check, all 1,782 signed minimal S6 vectors have
absolute pairing `0,6,12` with `e1`, never the formerly required 18.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h150_frontier.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h150_frontier.py
```

Source: {SOURCE_URL}
"""


def main() -> None:
    payload = build_payload()
    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    REPORT_PATH.write_text(report(payload))
    print(JSON_PATH)
    print(REPORT_PATH)


if __name__ == "__main__":
    main()
