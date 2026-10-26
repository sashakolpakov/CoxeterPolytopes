# Source review: browser memory flag delta, 26 October 2026

Reviewer: independent manuscript/documentation source reviewer. Decision: the
five-flag delta is accepted for the existing bounded publication renderer.
This is a source-only review; the reviewer launched no browser, scientific
worker, build, rendering job or test.

## Bound inputs

| Item | SHA-256 |
| --- | --- |
| Previously reviewed renderer, 6,994 bytes | `ae95bab9d11db3610c6fb2dce1b7dd8fc0a1736b34d688c5d2c87e460cf9428a` |
| Current `scripts/render_documentation_pages_20261026.py`, 7,276 bytes | `d15dcf5faa90fc25f38fd3b3eeccaf157c851ff2b62b3950947bd93ce9f59be1` |
| Reviewed four-line patch, 869 bytes | `b7c7f41cb9763c5b90e7c6a44bd489edebd406b70229968d1b2dfa7df4e882cc` |

The reviewer read and hashed the current source and retained predecessor.
Removing exactly the four added command-list lines restores the predecessor
bytes. The complete delta adds these five Chrome arguments:

- `--renderer-process-limit=1`
- `--js-flags=--max-old-space-size=128`
- `--disable-breakpad`
- `--disable-crash-reporter`
- `--disable-component-extensions-with-background-pages`

## Findings and scope

The delta attempts to reduce browser startup memory and unnecessary background
activity. It does not add `--no-sandbox` or `--single-process`. The existing
percentage resource admission/abort policy, 1 GiB process-tree RSS stop,
process-tree and ownership monitoring, bridge cleanup, timeouts, bounded log
handling, input path/hash checks and fresh per-page profiles are unchanged.
The root's one-input-per-batch execution plan requires no helper change.

The renderer process setting is not a hard total-process cap, and the V8 setting
limits old space rather than total browser RSS. Crash-handler flags may have
platform/version-dependent effects. Acceptance therefore makes no claim that
these flags alone enforce a 1 GiB memory limit; the existing outer monitoring
remains responsible for that stop condition. These distinctions follow the
[Chromium renderer switch definition](https://chromium.googlesource.com/chromium/src/+/master/content/public/common/content_switches.cc)
and [V8 old-space flag definition](https://chromium.googlesource.com/v8/v8/+/e6902daebfce639b4e48388f5cec2dcd478fa7e0/src/flags/flag-definitions.h).

The coordinating root separately reported successful one-input batches for one
direct SVG, two montages, desktop HTML and mobile HTML, each below the existing
1 GiB stop threshold with zero owned processes remaining. That is root-reported
execution evidence, not execution performed or independently replayed by this
reviewer. This note certifies only the exact source delta and preservation of
the previously reviewed safeguards. The separate input-guard fixture revision
is outside this delta review.

This note supplements
[the publication tooling review](review_documentation_publication_tools_20261026.md).
No shared tool source was edited by this review.
