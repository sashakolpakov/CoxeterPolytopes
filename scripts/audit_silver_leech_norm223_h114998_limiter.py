#!/usr/bin/env python3
"""Freeze/verify realization and acceptance of the norm-223 H114998 limiter."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from coxeter_polytopes.exact.silver import Silver  # noqa: E402
from coxeter_polytopes.search.silver_leech_norm223 import (  # noqa: E402
    norm223_leech_prefix_audit,
)
from coxeter_polytopes.search.silver_leech_norm223_limiter import (  # noqa: E402
    norm223_limiting_root_audit,
)


JSON_PATH = ROOT / "data/enumerations/silver_leech_norm223_h114998_limiter.json"
REPORT_PATH = ROOT / "results/silver_leech_norm223_h114998_limiter.md"
BASE_ARTIFACT = ROOT / "data/enumerations/silver_leech_norm223_h250000.json"
SOURCE_PATHS = (
    "src/coxeter_polytopes/exact/silver.py",
    "src/coxeter_polytopes/search/cubic7_leech.py",
    "src/coxeter_polytopes/search/cubic7_leech_continuation.py",
    "src/coxeter_polytopes/search/silver_leech.py",
    "src/coxeter_polytopes/search/silver_leech_manifest.py",
    "src/coxeter_polytopes/search/silver_leech_norm223.py",
    "src/coxeter_polytopes/search/silver_leech_norm223_limiter.py",
)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def encode(value: object) -> object:
    if isinstance(value, Silver):
        return value.to_json()
    if isinstance(value, Fraction):
        return str(value)
    if is_dataclass(value):
        return encode(asdict(value))
    if isinstance(value, dict):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def build_payload() -> dict[str, object]:
    if not BASE_ARTIFACT.exists():
        raise AssertionError("the frozen H250000 necessary-prefix dependency is missing")
    prefix = norm223_leech_prefix_audit()
    audit = norm223_limiting_root_audit()
    if not audit.accepted:
        raise AssertionError("the limiting root is not Vinberg-accepted")
    return {
        "schema": "coxeter.silver_leech_norm223_limiting_root.v1",
        "field": "Q(sqrt(2))",
        "host": "<-(-13+14*sqrt(2))> + Lambda_24",
        "hyperbolic_dimension": 24,
        "source_sha256": {
            path: file_sha256(ROOT / path) for path in SOURCE_PATHS
        },
        "dependency": {
            "path": str(BASE_ARTIFACT.relative_to(ROOT)),
            "sha256": file_sha256(BASE_ARTIFACT),
            "scope": "complete necessary arithmetic prefix through H=250000",
        },
        "audit": audit,
        "surviving_prefix_witness": {
            "coordinate": prefix.witness_coordinate,
            "minimum_necessary_wall_margin": prefix.minimum_wall_margin,
            "lorentzian_excess": prefix.witness_hyperideal_excess,
            "strict_hyperideal": True,
        },
        "conclusion": {
            "necessary_superset_limiter_resolved": True,
            "explicit_primitive_reflective_root": True,
            "chronologically_vinberg_accepted": True,
            "accepted_root_supports_unshrunk_limiting_ray": True,
            "strict_shrunken_ray_survives_complete_H250000_prefix": True,
            "finite_volume_through_prefix": False,
            "compact_example": False,
            "vinberg_termination": False,
            "global_nonreflectivity": False,
        },
        "scope": (
            "The explicit H=114998+81316*sqrt(2) root is realized and accepted "
            "against the complete earlier necessary superset and its complete "
            "same-height cohort. The surrounding H<=250000 computation remains "
            "a bounded necessary-prefix audit; no termination claim is made."
        ),
    }


def render(payload: dict[str, object]) -> str:
    audit = payload["audit"]
    witness = payload["surviving_prefix_witness"]
    first, second = audit["maximizing_types"]
    return f"""# Silver Leech norm-223: the H114998 limiter is an accepted root

The complete necessary audit through `H=250000` had a unique limiting
component at

`H=114998+81316*sqrt(2), q=2, x=130+92*sqrt(2)`.

Its 216 binary Gram types have exact norm-four-direction maximum
`679+479*sqrt(2)`. Exactly two types attain it:

- `{tuple(first['triple'].values())}`, with augmented determinant
  `{first['augmented_determinant']}`;
- `{tuple(second['triple'].values())}`, with augmented determinant
  `{second['augmented_determinant']}`.

Both maxima uniquely have `(d.u,d.v)=(679,479)`. In the first type, putting
`e=u-170d` and `f=v-120d` gives `(e^2,f^2,e.f)=(10,4,6)`, so
`(e-f)^2=2`; this type is impossible in the Leech lattice. In the second,
`(e^2,f^2,e.f)=(6,6,6)`, hence `e=f`.

## Explicit realization

In standard MOG coordinates scaled by `sqrt(8)`, take

```text
d = {tuple(audit['explicit_d_scaled'])}
e = {tuple(audit['explicit_e_scaled'])}
u = 170*d+e
v = 120*d+e
```

Then `d^2=4`, `e^2=6`, `d.e=-1`, and
`(u^2,v^2,u.v)=(115266,57366,81316)`. Thus

`r=(130+92*sqrt(2), u+sqrt(2)*v)`

has square two and realizes the limiting direction pairing. The vector is
primitive: divisibility by the only possible nonunit `sqrt(2)` would force
`u in 2*Lambda`, whereas `d.u=679` is odd. Norm two makes reflection
integrality automatic.

## Chronological acceptance

There are `{audit['earlier_component_count']}` earlier nonempty necessary
components. Writing the spatial part of `r` as

`((679+479*sqrt(2))/4)d + (1+sqrt(2))(e+d/4)`

gives `(e+d/4)^2=23/4`. Exact Schur/Cauchy bounds are strictly negative for
every earlier component. The tight rejection-square certificate is control
ordinal `{audit['tight_earlier_control_ordinal']}` and equals
`57+41*sqrt(2)>0`. The simultaneous norm-one component is empty. The only
possible positive pairings with a distinct same-height norm-two root are
`1` and `sqrt(2)`; their spatial difference targets `2` and
`4-2*sqrt(2)` have empty Leech Gram censuses. Hence the displayed root is
Vinberg-accepted.

## Scope

This upgrades the old necessary-superset threat to a genuine accepted wall,
but it does not produce a finite-volume example. The archived strict shrink
still survives every necessary wall through `H=250000`, with exact Lorentz
excess `{witness['lorentzian_excess']['a']} + ({witness['lorentzian_excess']['b']})*sqrt(2)>0`.
This is a bounded result, not Vinberg termination or global nonreflectivity.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = encode(build_payload())
    artifact = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    report = render(payload).encode()
    if args.verify:
        if JSON_PATH.read_bytes() != artifact or REPORT_PATH.read_bytes() != report:
            raise SystemExit("frozen norm-223 limiter certificate differs")
        print("PASS", sha256(artifact).hexdigest())
        return
    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_bytes(artifact)
    REPORT_PATH.write_bytes(report)
    print(sha256(artifact).hexdigest())


if __name__ == "__main__":
    main()
