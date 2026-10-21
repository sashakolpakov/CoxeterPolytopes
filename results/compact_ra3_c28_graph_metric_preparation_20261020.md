# Fixed C28 graph/algebra checker: source-only preparation

20 October 2026. No worker import or execution. Root alone schedules one
bounded generation and a separate fresh-process replay after code review.
The frozen C26 packet, primary manuscript and central status are unchanged.

Worker: `scripts/check_compact_ra3_c28_graph_metric_20261020.py`.
Certificate: `results/compact_ra3_c28_graph_metric_20261020.json`.
Written proof SHA256:
`076864a5d836456e56ac377dffbcaa581a4738ec977dbf58a7905580d9838afc`.

The stdlib-only worker constructs the two specified 16-vertex dual
complexes, checking all 28 triangles, 42 edges with two incidences,
connected cyclic links, Euler characteristic two, twelve degree-five
and four degree-six vertices. The degree-six induced graphs distinguish
the two types. It verifies all 24 index permutations of the tetrahedral
type, its distinguished-edge endpoint incidences, and the second type's
four explicit permutations, face orbits, supplied cyclic links and
two-regular-neighbor propagation. This is a check of two fixed complexes,
not an enumeration of fullerenes.

Exact Fraction arithmetic and a self-contained Q(sqrt5) class check the
tetrahedral squared metric identities and edge value 4/3<phi; the second
type's nondegenerate basis, normal equations, eliminated quartic and its
two roots, rejection of the root below one, H Gram relation/rank and the
one-dimensional unit-normal fiber behind the duplicate obstruction.
The final volume comparison is checked only as conditional scalar algebra.

Written premises remain explicit: the primary BGM2012v2 Table 2 count
two at 28 vertices; actual compact RA3 metric uniqueness and symmetry
realization; polygon capacity and its full L7 equality pattern; actual
normal/edge-length and pentagon closure lemmas; and the inherited RA4
capacity, dodecahedral-facet existence, W gap and volume formula. The
worker neither retrieves the primary paper nor proves those geometric
premises or existence/compactness from a formal Gram matrix. C26 sources
and its saved certificate are hash-bound inputs, never imported or rerun.

Controls follow the frozen C26 parent-lock contract: COXETER_CHECK_OWNER
must equal the actual parent PID, child alarm 45s, wrapper watchdog 59s,
one thread, 1GiB process-tree RSS cap, inherited lock/owner monitoring and
the existing disk/RAM gates. No subprocess or inner lock is created by
this checker. Certificate cap 128KiB; exclusive --write refuses overwrite.
--verify independently recomputes all fixed checks and compares the whole
saved JSON. All collected proof/control paths and the worker's own bytes
are stream-hashed before and again at completion. The bounded saved read
is checked unchanged before exit. Optional --verify-code checks only
current source pins against the saved certificate, with explicit
PROVENANCE_ONLY status; it is not a proof replay and needs no separate run.

Root-only planned commands (not run):

    python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock --ledger results/compact_ra3_c28_graph_metric_generation_20261020.json -- scripts/check_compact_ra3_c28_graph_metric_20261020.py --write

    python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock --ledger results/compact_ra3_c28_graph_metric_replay_20261020.json -- scripts/check_compact_ra3_c28_graph_metric_20261020.py --verify

On any mismatch, resource gate, cap or failed identity, stop without
expanding the domain or replacing evidence. This preparation is not an
execution ledger; no C28 certificate has been generated at handoff.
