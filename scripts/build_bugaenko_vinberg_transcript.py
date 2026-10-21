"""Print Bugaenko's exact Vinberg wall-enumeration transcript."""

from __future__ import annotations

import json

from coxeter_polytopes.search.bugaenko_vinberg import (
    build_bugaenko_vinberg_transcript,
)


if __name__ == "__main__":
    print(json.dumps(build_bugaenko_vinberg_transcript(), indent=2, sort_keys=True))
