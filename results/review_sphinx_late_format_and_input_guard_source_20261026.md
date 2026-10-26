# Independent source review: final Sphinx format and input guard

Date: 26 October 2026. Reviewer: `/root/finite57_completion_audit`.

**Source-clear for the eleven heading repairs and the final guard below.**
This review used source and existing receipts only. The reviewer launched no
Sphinx, publication, scientific, Lean, TeX, PDF, rendering or test process and
imported no repository worker. Only this review note was written for this task.

## Eleven heading repairs: exact inverse mapping

The first strict build's preserved stdout reports eleven warnings, all
`Title underline too short`, and exit 1 with no watchdog stop reason.
`results/sphinx_short_underline_repair_20261026.json` is 3,774 bytes,
SHA256 `fe4c3b48545a58ac5918c948a6b5ddb06af08eb3e5d24940dc58ba45c508c739`.

For every recorded edit, the actual current line equals the recorded after-line.
Replacing only those eleven adornment lines by their before-lines restores each
full file's before-SHA exactly. Each before-SHA also equals that path's actual
binding in `results/sphinx_inputs_actual_index_guard_20261026.stdout`.
There is no body, formula, link, literal, numerical claim, quantifier or theorem
change in this repair.

| File | Changed lines | Restored before SHA256 | Current after SHA256 |
|---|---|---|---|
| `docs/site/methods.rst` | 234 | `0313362c046a4ca6d45702a92051c8d03efbabc47b6cead3ef4b31c4b6d61a6c` | `9b80bf39fee1a35cf814fdced8209eaac690e07574998b8010d05cec89302d6d` |
| `docs/site/right_angled.rst` | 22 | `3322501fbcf56821bc279bb48f4d914a7050e199c5747a774efa69be0015e1c4` | `dcd7d3e5256525229c560691c0b10f2fa6f872d6618ea496621a04a7c1202828` |
| `docs/site/sources.rst` | 73, 93, 111 | `2fca826d051b6b6b76adb85f313a05beb8d9b1e7ecca705a83b02d5499988ed7` | `42aebd1eb33adaea64ad5f25dac2a4587c30623bab6cc32f29ea6b3d5b6a9746` |
| `docs/site/verification.rst` | 131, 171 | `5a8a41be14bf718c444ddeb372b62d3f2b9540c9e13ad777fce79404bc368f32` | `a47fc0c6257a62808497e0b71d6808593cc28f8cb898ee289246f2cd1ee63177` |
| `docs/site/volume_problems.rst` | 33, 304, 333, 364 | `251369fe86d33b2cd20b18d182cb79c11090fd6bb807e0c5cf95b62d60b108b1` | `223238152885753af41722abae0ab029c22c4e02a4803172ae48c3c427204edd` |

The original guard receipt has 14 RST pages, 190 explicit local references,
145 input bindings and 973 denied manifest paths, with an empty issue list and
`require_indexed: true`. The only changed bindings at the first comparison were
these five pages; their exact inverse mapping preserves the prior source
review and input-provenance evidence. A later receipt can bind their final
adornment bytes directly without suggesting a new mathematical replay.

## Hardened publication-only guard

Read the complete current `tools/audit_sphinx_inputs.py`: 9,307 bytes,
SHA256 `1cc08d0dbb13c65b1895af5eb53d2c5d67bbc54f370e88f39f081e44d13c8923`.

The source binds all current RST sources, both input manifests, configuration,
templates and static assets. It examines titled and untitled downloads, images,
figures, includes and literal includes, and recursively scans included text.
The visited set prevents include cycles from looping. It checks current
source/input membership in the actual Git index when requested. This is index
membership, not a substitute for the separately required staged-blob audit.

Blocked identity checks cover the 936 scientific-file paths/hashes, the
37 third-party-file paths/hashes and the 33 prepared archive hashes. The source
rejects named private/cache components, missing inputs, repository escapes and
symlink entries. Text scans have a 4-MiB per-source bound; input hashing has
16-MiB per-file and 128-MiB aggregate bounds, with a 4,096-file bound. An issue
causes refusal rather than a truncated successful audit. The documented
exclusions for arbitrary Python/JavaScript/extensions and modified or unknown
compressed payloads remain explicit.

An intermediate canonicalizer normalized away raw `link/../name` components
before inspecting symlink ancestors. I reported that concrete case to root;
it did not occur among the existing 145 bound inputs. The accepted successor
retains the raw `Path` components, checks each step before any final resolve,
rejects a traversed symlink immediately, and rejects an intermediate repository
escape even if a later component could return. Normal current `../..` source
references remain inside the repository and retain their intended targets.
The owned temporary symlink-followed-by-parent-traversal fixture exercises the
specific correction. The final source above implements this logic correctly.

Existing execution records retain their exact boundaries:

* `results/sphinx_input_guard_self_test_20261026.stdout` reports eight actual
  synthetic checks, PASS, no mathematical checks and watchdog exit 0.
* `results/sphinx_input_guard_raw_ancestor_regression_20261026.stdout`, SHA256
  `c236a34d043c2ef535346b962454cdd81b97812b029b7e276143d080c570f331`,
  reports nine actual synthetic checks, PASS, no mathematical checks and
  watchdog exit 0 after the added regression.
* `results/sphinx_inputs_after_heading_repairs_actual_guard_20261026.stdout`,
  SHA256 `ce28a169a325cfd4cb4a0ca9a3618388dafbf8adfa8dd1a955a156769ac9ac2c`,
  binds all five repaired pages and has empty issues, 14 pages, 190 references,
  145 inputs and 973 denied paths. It pins the preceding 8,566-byte guard,
  SHA256 `cca58406146cb112268ee7b7cb3dd4d273b7713f709a774ab3c5e760260990b1`.
  At this review's source comparison its sole current mismatch is the guard
  successor itself. A final actual indexed-input receipt must bind the
  9,307-byte successor before being described as the current guard PASS.

This source clearance supplies no new execution credit. Strict Sphinx, the
current indexed-input audit and actual staged-tree audit retain their own
receipts; source publication, deployment and bulk availability stay separate.

## Two small chapter deltas also reviewed

The methods sentence now explicitly requires `b_i - mu*t_i >= 0` at every
required anchor, consistent with its existing individually funded reference
argument. The verification chapter's expanded guard description matches the
reviewed final implementation, including recursive source coverage, source and
template bindings, archive hashes and the stated finite caps. Its six-figure
visual link names the actual sampled page/PDF record rather than transferring
old coverage. Neither change strengthens a mathematical conclusion.

The earlier complete mathematical reviews of these chapters and the other
native pages remain applicable at the exact editorial mappings above.
