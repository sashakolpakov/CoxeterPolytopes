# Documentation publication tooling: source review

26 October 2026. Reviewer: /root/compact_manuscript_completion_audit.
**REQUEST CHANGES to the original two helper sources; proposed repairs cleared at source scope only.**

I read both complete helper sources, the publication-watchdog allowlist change,
the documentation CI change, current Sphinx configuration, the manifest schemas,
and the existing process-tree/watchdog/group-stop implementations they use.
The Git diff was read with small mapped-pack/window/delta caches. No tool under
review was launched or imported. No self-test, browser, renderer, Sphinx, Lean,
TeX or scientific worker ran. This is not runtime or interruption-test evidence.
Only private proposal files and this new review note were written.

## Original findings requiring repair

1. **Static/source symlink boundary.** Explicit role targets are resolved, but
   the original static walk adds unresolved paths and follows file symlinks.
   A static symlink to an unknown private file can evade the forbidden-path
   check and be copied by Sphinx. Directory symlinks are not reliably covered.
   A source-page symlink can also be read without the corresponding path guard.
2. **Incomplete explicit-input traversal.** The original download regex only
   recognizes labelled roles with angle brackets. Ordinary implicit download
   roles are missed. Included text outside the initial site-page list is hashed
   but not recursively scanned for further downloads/includes. Image substitutions
   and directive arguments containing spaces also need recognition.
3. **Source/index bindings.** The original actual-index requirement covers
   referenced/static targets but not the RST source pages themselves. Source
   pages, configuration, templates and decision manifests should be bound along
   with the assets; new untracked source pages must not silently enter the build.
4. **Detached browser ownership.** Chrome starts a new session. The original
   renderer's default SIGALRM can terminate Python before its finally block, and
   abrupt outer-owner loss kills the renderer group without necessarily killing
   that separate Chrome group. A child-group ownership watcher is required.
5. **Captured-output bound.** The original communicate call collects all browser
   output in memory before testing its length. The stated 64KiB check therefore
   does not bound that capture. The original renderer also records an observed
   input hash rather than requiring a caller-pinned source version.

These are publication/process-boundary findings, not defects in any mathematical
proof or a reason to rerun successful proof jobs.

## Proposed repairs

The unapplied complete patch is
/private/tmp/documentation_publication_tools_hardening_20261026.patch,
SHA256 b9f92f06ae34160ae0b96b78fd6377da462179d349b4c54bc90bbec60bef032b.

| Helper | Original SHA256 | Proposed SHA256 |
|---|---|---|
| tools/audit_sphinx_inputs.py | a64cd028f6ba1eea235676e99645acf96ab2a267b8a94ca66e0a1e9dcd58fbad | cca58406146cb112268ee7b7cb3dd4d273b7713f709a774ab3c5e760260990b1 |
| scripts/render_documentation_pages_20261026.py | 791860cac1725c0ebcaa99a9a47b9e766b8d7ea40524e10e1f0b0b661cebe659 | ae95bab9d11db3610c6fb2dce1b7dd8fc0a1736b34d688c5d2c87e460cf9428a |

The scanner proposal rejects symlink files and ancestors before reading them,
prunes symlink directories, rejects non-regular site entries and canonicalizes
all assets. It supports implicit/labelled downloads, image substitutions and
space-containing directive arguments. Included text is scanned through a visited
queue; literalinclude remains literal. All site RST pages, current configuration,
templates, static assets and both manifests are hash/index-bound. Known bulk
archive digests join the original file/external-input digest denylist. Hashing
remains streamed; per-file 16MiB, text 4MiB, aggregate 128MiB and 4096-file limits
bound the intended scan. The Git path listing has a ten-second timeout and an
8MiB post-read bound; it reads the index, not Git objects. The outer publication
watchdog remains the live RSS/time/resource guard. The synthetic parser check is
expanded in source but was not executed by this reviewer.

The renderer proposal accepts at most three owned local SVG/HTML inputs, requires
one expected SHA256 per page, rejects symlink ancestors and checks each source
again afterwards. A bounded-parent ownership environment is required before any
browser launch; the coordinator must still invoke the existing wrapper with
--parent-lock. SIGTERM/INT/HUP/ALRM raise into cleanup. The separate browser group
contains a small watcher process which kills only its own new group if its parent
is lost, and has its own 16-second alarm. Ordinary timeout/error cleanup uses the
existing signal-safe group stopper. No unrelated browser/profile is targeted.

Browser output goes to a temporary file with a sampled 64KiB size bound and only
a 3000-byte error tail is read into memory. The parent 14-second per-page and
50-second total alarms remain below the outer 59-second budget. A fresh temporary
profile is used for each page; 8MiB disk-cache and 1MiB media-cache settings are
browser hints, not enforced allocation caps. Normal/handled-error profile/log
cleanup is owned by the context managers. An uncatchable SIGKILL can still leave
owned temporary files for coordinator cleanup; it must not leave an indefinitely
running browser group. Screenshots remain coordinator-owned until inspection
and cleanup. This ownership reasoning has source review only, not a live signal
or owner-loss exercise.

The renderer CLI intentionally adds repeated --sha256 arguments in page order.
Root should supply hashes from the selected source/build snapshot and use the
existing monitored wrapper with --parent-lock. No new scientific runner, parallel
proof execution, heap claim or weakened admission/abort policy is introduced.

## Watchdog and CI deltas

The actual watchdog diff adds only tools/audit_sphinx_inputs.py to the existing
publication-I/O allowlist. It does not admit the browser renderer, mathematical
checks or arbitrary commands. The watchdog source SHA256 at review is
49b4d010172766458a9d32bbfd412d669cc1a9838bfb9dd33f7ac2aee6badb62.

The actual CI diff adds only the input audit with --require-indexed before the
strict Sphinx step. Checkout does not persist credentials; build permissions
remain contents:read, pip caching remains limited to documentation dependencies,
and the artifact path remains the generated site alone. Pull requests cannot
upload/deploy Pages. Main-branch deployment still depends on a successful build;
the later actual-index audit remains separate. No proof, bulk restoration or
scientific worker was added. CI source SHA256 at review is
490a8355571750ecec9e3a4a8345c1d6bbd60941bbf204edba5f84d0368071b5.

These two narrow integration deltas are acceptable with the repaired input guard.
They do not constitute a new successful CI run, deployment or bulk upload.

## Explicit limits and verdict

The scanner is a source-input gate for the stated local RST/static/template
surface. It does not analyze arbitrary Python/Sphinx extension behavior, raw HTML,
JavaScript, remote resources or unknown transformed/compressed private data. Its
known-path/exact-hash checks must not be called a complete data-loss-prevention
proof. The current configuration is ordinary declared Sphinx input configuration.

The browser hashes the selected top documents, not every subresource. Its
background-networking flag is not a network-isolation guarantee. A screenshot
is visual evidence, not mathematical verification or exhaustive all-page review.

The proposed patch addresses the identified source defects without changing any
proof source or successful evidence scope. Root acceptance and actual bounded
parser/path checks, strict Sphinx build and small browser batches remain required.
No runtime PASS is supplied by this note.
