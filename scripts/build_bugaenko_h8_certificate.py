"""Print the exact geometry certificate for Bugaenko's compact H8 polytope."""

from __future__ import annotations

import json

from coxeter_polytopes.known_examples.bugaenko import bugaenko_h8_certificate


if __name__ == "__main__":
    print(json.dumps(bugaenko_h8_certificate(), indent=2, sort_keys=True))
