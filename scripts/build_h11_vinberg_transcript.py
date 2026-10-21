"""Print the exact bounded Vinberg transcript for the H11 field lattice."""

from __future__ import annotations

import argparse
import json

from coxeter_polytopes.search.h11_vinberg import build_h11_vinberg_transcript


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--identity-bound", type=int, default=20)
    args = parser.parse_args()
    print(
        json.dumps(
            build_h11_vinberg_transcript(args.identity_bound),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
