"""Print the bounded exact S19 mixed-root Vinberg transcript."""

from __future__ import annotations

import argparse
import json

from coxeter_polytopes.search.s19_vinberg import build_s19_vinberg_transcript


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--control-bound", type=int, default=3)
    parser.add_argument("--attempt-stopping", action="store_true")
    args = parser.parse_args()
    print(
        json.dumps(
            build_s19_vinberg_transcript(
                args.control_bound, attempt_stopping=args.attempt_stopping
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
