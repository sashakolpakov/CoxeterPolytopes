# Bounded proof execution: verified status

The current research priority is lowering dimension bounds and finding
examples. Lean follows settled mathematics as a final check.

Latest20 October batch: every subset-size obligation for the new mixed
right-angled ceiling10 passed written review and exact generation/replay.
Only afterward `RightAngledMixed10` compiled in1.47seconds and passed a
separate1.21second eight-theorem axiom audit, one thread/2048MiB.
Its aggregate real transport and geometric interfaces remain explicit.
The four small exact sources, including the independently retained k9
interval proof, each completed in under0.5seconds and56MB process-tree
RSS. All actual jobs used the shared lock and normal resource gates.
Disk remained about181GiB; the below20GiB cleanup trigger did not apply.
No unused build or raster outputs were created by this scalar batch.

The [runner](../scripts/run_lean_batches_20261018.py) uses one module per
batch, one thread and a2048MiB Lean memory limit, with an exclusive lock.
It verifies source, dependency-object and compiler/configuration hashes
before publishing an object. A512KiB output tail bounds retained logs.
Existing compiled dependencies without a new ledger are explicitly labeled
legacy cached objects; they are not relabeled as freshly capped builds.

The resource gates require40GiB free disk and6GiB available memory to
start a proof. A running job stops below8GiB disk or3GiB available memory.
Below20GiB disk, the coordinator must immediately start or reassign an
agent to cleanup before further jobs. The runner emits a small request
marker; it cannot itself call the agent orchestration API. No cleanup
daemon or additional proof farm is left running.

Synthetic small-process tests passed for source/object/dependency changes,
bounded output, stdout closing before process exit, timeout cleanup,
exclusive locking, all-cached cleanup and resource request/floor behavior.
The source files pass Python syntax checks. The aggregate checker routes
Lean through this runner and cleans TeX intermediates after each PDF.

A real final check then passed in ce11a3/16d309:

```sh
python3 -B scripts/run_lean_batches_20261018.py --audit-axioms
```

The root module compiled in1.90 seconds and its separate axiom audit took
0.25 seconds. Both used the memory cap and one thread. All23 dependency
objects were reused from their prior successful kernel builds. The audit
listed only propext, Quot.sound and Classical.choice. Necessary objects
and small hash ledgers were retained; temporary work and unused build
outputs were removed. No proof process remained after completion.

This validates execution of the current partial formal project. It does
not formalize new geometric results or turn the full dimension57/compact
arithmetic21 theorems into end-to-end Lean proofs. The pending generated
certificate expansions must still be split into smaller modules before
fresh compilation; mathematical research takes precedence over that work.

The separate three-page octahedral-block manuscript compiled with no
final-pass warnings in eeadd1; unused TeX files were removed afterward.

## Signal-safe cleanup follow-up

A later sandboxed numerical scout exposed a real cleanup edge case:
`killpg(pid,0)` can raise PermissionError even when signaling our owned
process group with SIGTERM is permitted. The earlier runner is preserved
byte-identically because frozen execution records pin its hash. Future
Lean batches use the thin
[signal-safe wrapper](../scripts/run_lean_batches_signal_safe_20261018.py).
It retains all earlier caps, locks, resource gates and ledgers, waits on
the actual leader, and sends a final SIGKILL to any remaining members of
the owned group without a signal-zero probe. The aggregate entry point
now calls this wrapper too; the complete aggregate has not been rerun.

Real small-process integration checks passed in281c49: an exited leader
with a SIGTERM-ignoring descendant, and a directly ignoring leader, were
both reaped/terminated. A duplicated output pipe reached EOF only after
the descendant was killed. The test rejected any attempted signal-zero
probe. Both real SIGTERM and SIGKILL worked in the sandbox. No Lean
compilation or new geometric verification was performed for this follow-up.

## New arithmetic20 and illustrated manuscript batches

The new compact arithmetic20 result has exact rational generation, a fresh
independent face/forest check and the final13+24-branch composition audit.
The root independent check took4.121seconds with85.7MB peak process-tree
RSS; the streamed source audit took0.79seconds with74.6MB. It bound107
manifests and228 files while streaming502MB rather than loading their
large arrays. A separate solver-free producer replay also passes.

Two illustrated right-angled manuscripts compiled and were frozen with
their mathematical reviews: nine pages/five figures for the restricted
octahedral-chain proof, and seven pages/two figures for compact facet
capacity. Both final builds had zero warnings, used less than86MB peak
process-tree RSS and removed their own temporary directories.

The new compact20 and mixed-opposition TeX sources have mathematical
translation reviews. Their next PDF builds have stopped at the start
gate, before launching TeX, while available memory was5.2--5.5GiB. A
rounded display of6.0GiB can still be below the exact threshold and was
also rejected. No gate was relaxed. Disk remained above194GiB, so the
user's below20GiB cleanup-agent trigger did not apply.

The [common research-check wrapper](../scripts/run_bounded_research_check_20261018.py)
supplies the live process-tree watchdog for direct checker entry points
that only have their own alarms/periodic checks. Its first execution and
the new compact20 Lean source are queued behind the same start gate.
The Lean source is not yet claimed as kernel checked. The result's full
geometry remains outside existing formal coverage.

## Consolidated ideal-four package

The shared checker wrapper now supports legacy non-locking entrypoints
with an inherited parent-held lock, bounded output, source-hash checks,
live59-second/1GiB process-tree limits, parent-death handling and reliable
owned-group cleanup. All42 catalogue replay/lock modes were read. The
new per-manuscript driver passed its real resource-deferral path and then
the ideal-four core integration run: four exact checks, cached old gap
logic, fresh capped `ManuscriptIdeal4` compilation and actual axiom prints.
The final audit lists only standard foundational axioms. Named published,
script and geometric proof parameters remain in the assumption register;
they are not listed by Lean's axiom printer.

The same source builds to16 pages/five figures with zero warnings in
1.914seconds and84.3MB peak RSS. Temporary files were removed. The next
compact20 PDF attempt stopped before TeX at5.4GiB available RAM and
193.4GiB free disk. No additional proof process was launched after that
deferral. The below20GiB cleanup trigger did not apply.

## Subsequent bounded research and main-paper verification

The compact20 paper now builds to8 pages with zero warnings,82.2MB peak
RSS and automatic temporary-file removal. Its
[ten-job core run](../results/compact20_core_verification_20261019.json)
passes three exact replays, six declared Lean module checks and a fresh
temporary PDF build. The newly compiled compact20 and finite57
conditional composition theorems pass actual printed-axiom audits with
only propext and Quot.sound, no sorryAx. Old imports retain their legacy
cache labels; geometric premises remain explicit.

The A15 monotone refund, D17 physical-face capacity and Figure7 cusp model
each have successful separate generation and replay ledgers. Their new
checker jobs took at most0.421seconds and92MB process-tree memory. The
initial A15 JSON tuple/list replay defect was fixed in a distinct v2
certificate with the original evidence preserved. A D17 attempt also
stopped below6GiB RAM before child launch and resumed only after recovery.
No gates were relaxed, no scientific jobs ran concurrently, and no
optimizer or giant-data reload was needed. Disk remained above192GiB;
the cleanup-agent trigger did not apply.

The new `RightAngledCuspCapacity` module subsequently passed a single-module
compilation and separate printed-axiom audit. Its final warning-free source
checks the59/180 and240/4980 consequences and updated mixed inequalities,
conditionally on the listed geometry. Final compilation/audit took0.45/0.43s,
with192.0GiB disk and7.77/7.76GiB available RAM at their respective gates.
Only standard foundational axioms were printed; no `sorryAx`. The runner
removed its own intermediate outputs. This is not a new dimension ceiling.

The opposed-fork arithmetic model checker then passed generation and fresh
separate replay under the shared lock and live59s/1GiB process-tree guard.
Each run took about0.91s and under56MB, at191.97GiB disk and7.75–7.76GiB
available RAM. Its four exact stars are local models, not new complete
polytopes. The small certificate binds the three written model/capacity
reports; those source bytes must remain unchanged. Both runs exited0 and
no process remains from these jobs.

After the new positive cusp-defect proof passed both written reviews,
`RightAngledCuspDefect` compiled in0.59s and passed its separate axiom
audit in0.43s. The capacity predecessor had a capped cached dependency;
only one new module ran, one thread,2048MiB. Start gates recorded190.89GiB
disk and7.31/7.29GiB available RAM. No warning or `sorryAx`; standard
foundational axioms only. Temporary outputs were removed. These are
conditional62/189 and252/5504 count implications, not dimension exclusions.

The next reviewed facet-propagation proof then received its own one-module
batch:1.02s compilation and0.34s axiom audit, one thread,2048MiB, at190.24GiB
disk and7.93/7.86GiB RAM. Its two dependencies were capped cached modules;
no warning or `sorryAx`. The module checks62/206 and269/5875 implications,
conditionally on explicit geometric count interfaces. Intermediates were
removed. Dimension ceilings remain unchanged.

The26-wall H9 local-algebra checker subsequently passed generation and
separate replay in0.229/0.225s, under58MB process-tree RSS. Start gates
recorded183.05GiB disk and6.38/6.69GiB available RAM. All325 pairs and
the listed actual-facet/vertex witnesses check; two spacelike initial
edge endpoints prove this chamber still has infinite volume. A following
bounded short-pencil search and replay each took0.258s, under57MB, with
186.67GiB disk and6.84/6.66GiB RAM. Its537 admissible controls gave no
field-square solution. A subsequent paper congruence argument excludes
that whole pencil at every height; the bounded search alone did not.
Both source families used the wrapper's shared lock,59s/1GiB caps and
parent-death protection. No build intermediates or concurrent jobs were
created, and the below20GiB cleanup trigger did not apply. See the
[root execution review](../results/review_h19_twenty_six_wall_execution_root_20261019.md).

The34-profile H9 search next passed generation/replay in0.650/0.654s,
under58MB, with186.65GiB disk and7.51/7.81GiB RAM. Its two survivors
cannot coexist; the smaller gives a27-wall chamber. The fixed27-wall
generation subsequently passed in0.234s, under57MB, at187.49GiB disk
and6.19GiB RAM. It checks351 pairs,27 actual facets and eight finite
initial edges; the ninth still escapes. The first generation attempt
and subsequent separate replay attempts stopped before child launch
at5.8--5.9GiB available RAM. After recovery the separate replay also
passed in0.234s, under58MB, at187.49GiB disk and6.19GiB RAM.
No gate was relaxed and no temporary build outputs were created.
The below20GiB disk trigger does not apply.

The18-profile search then passed generation/replay in0.404/0.392s,
under57MB, at187.53GiB disk and6.20/6.29GiB available RAM. An earlier
launch deferred at5.7GiB. Its one compatible survivor gives a28th wall
and cuts the last known initial escape. The search is bounded in root
height; a full finite-volume conclusion is not supplied. The fixed28-wall
and further edge-walking diagnostic is still in preparation. These runs
were serial and created only small certificates/ledgers, no build outputs.

The fixed28-wall check and bounded walk subsequently passed generation
and replay in0.544/0.546s, under58MB, at187.47GiB disk and6.19/6.18GiB
available RAM. It verifies378 pairs,28 actual facets and all nine finite
initial edges, then finds a feasible spacelike endpoint after19 finite
vertices were discovered. This proves the chamber still has infinite
volume; the walk does not claim a complete census. The saved certificate
is34,290bytes. The shared lock,45s child alarm and59s/1GiB wrapper caps
were retained; no build intermediates or concurrent jobs were created.

The next D16/D13 saved physical-output audit is prepared as a54-file,
106,205,768-byte sequential read with four large sections streamed. Its
source was fully reviewed before launch. Several attempts stopped at the
prelaunch resource gate, the latest at187.4GiB disk and5.9GiB available
RAM. No scientific child, new certificate or execution ledger was created
by those deferrals. The source remains queued for one guarded generation
and a separate replay when memory recovers; no limits were relaxed.

After memory recovery, that D16/D13 audit passed generation and a separate
replay in0.506/0.472s, under91MB process-tree RSS. Start gates recorded
187.31/187.32GiB disk and6.48/6.27GiB available RAM. All54 inputs and
the declared full saved interface checked, with zero unknown targets.
The certificate is86,534bytes; both compact ledgers are under1KB. Both
jobs exited0 and released the lock. They created no build intermediates.
The following eleven-short-profile H9 search attempt then stopped before
child launch at5.9GiB RAM. No H9 search certificate or ledger was created
by that deferral. The cleanup trigger still does not apply.

The eleven-short-profile search subsequently passed generation and
separate replay in0.412/0.404s, below58MB process-tree RSS. Both started
with187.29GiB disk and6.12/6.14GiB available RAM. Its11803 exact square
tests produce two individually compatible mirrors; no all-height census
or finite-volume chamber is asserted. The fixed29-wall check is ready,
but its first launch stopped before the scientific child at5.6GiB RAM.
The59s/1GiB caps and shared lock were retained. No intermediates were
created and no cleanup threshold was reached.

The fixed29-wall checker passed separate generation/replay in0.721/0.766s
under58MB process-tree RSS, with normal gates and no intermediates. Its
new feasible spacelike endpoint proves infinite volume. The30-wall
checker then passed generation in0.765s under58MB; its first separate
replays deferred below6GiB RAM and have not launched a scientific child.
The31-wall source is reviewed/prepared; no31-wall execution yet.

The A2/D13-short saved-interface audit passed generation and separate
replay in0.508/0.477s under92MB, all54 inputs and zero unknown targets.
The current ideal5/6 manuscript passed its final bounded12-page build
in1.073s under87MB, zero warnings and temporary files removed. Four
figure pages were inspected, a Figure4 spacing change rechecked, and
all disposable previews removed.

The latest conditional Lean module RightAngledElevenVertexCounts compiled
and passed its separate printed-axiom audit in0.45/0.40s, one thread and
2048MiB, after7.43/7.44GiB RAM gates. No warnings or sorryAx occur in the
successful run. The resolved failed elaboration's temporary log was
removed after its precise correction was documented.

A serial manuscript replay subsequently passed all three exact scripts
and the cached current Lean/axiom check, then deferred only the check-only
PDF job at the memory gate. That aggregate ledger correctly records
DEFERRED_RESOURCE_GATE, not a complete pass. The identical final source
and canonical PDF already have their separate successful bounded build
ledger. No limit was relaxed and no cleanup threshold was reached.

The next continuation performed paper/source work only. Separate attempts
to replay the30-wall certificate, generate31 walls, run the fixed cubic
arithmetic producer, and build the compactRA4 successor all stopped at
the prelaunch memory gate. Latest compact build attempt:184.7GiB disk,
4.8GiB available RAM. No scientific child, successor PDF, corresponding
new certificate or success ledger was created. Existing evidence was
preserved. No proof process is left running by these deferrals.

Meanwhile the compactRA4 W>=16 manuscript has full source review; its
W>=15 source was archived in byte-verified gzip. Three arithmetic dotted
rows now have a uniform all-field written norm obstruction, and the
fixed cubic row has a complete hand proof. Their producer/helper remain
unexecuted support. Disk remains far above the20GiB cleanup trigger;
the former cleanup agent is doing mathematical work. Lean stays behind
settled research and is not used to occupy the blocked execution queue.

Available memory subsequently recovered above the6GiB launch gate. The
coordinator ran five serial research batches:30-wall replay;31-wall
generation and separate replay;D16/D14 repair generation and separate
replay. All passed under93MB process-tree RSS; every elapsed time was
under1.5seconds. The54-file audit has zero unknown outputs. No resource
threshold or lock was bypassed and no background worker was started.

Between those research batches, the reviewed compactRA4 source20e5b4a3
built successfully as a13-page PDF with four vector figures and zero
warnings. The2.745-second build peaked below85MB and removed its TeX
intermediates. Source-bound ledgers retain precise caps and initial
resources. Latest successful batch started with182.34GiB free disk and
6.13GiB available RAM. Disk is far above the20GiB cleanup trigger; the
former cleanup agent remains assigned to mathematical research.

The compactRA4 strengthening was subsequently consolidated in the SAME
paper: W=0 or W>=17 and degree-two=>W>=20. Both earlier source and PDF
were byte-verified gzip-archived before the unfrozen PDF was replaced.
Two complete source reviews pass. Its final14-page build took1.310s,
peaked below86MB,had zero warnings,and removed intermediates. A separate
four-figure render passed the normal gates; two byte-identical pages
reuse the preceding visual review and two revised pages were inspected.
All eight disposable previews from the two visual batches were removed,
reclaiming about2MB. No Lean job was launched for these new proofs.

The generic figure utility initially refused an uppercase output name;
the corrected lowercase invocation passed. The successor PDF builder
initially refused replacement without its explicit unfrozen flag;
the rerun used that flag after the verified archive. Neither initial
refusal started a renderer or TeX process, and neither bypassed a guard.

After the independent third-agent full-packet review, the settled
eight-field19 final logic was checked in a single serial module,
`ArithmeticEightField19`. Compilation took0.44s and the separate
printed-axiom audit0.35s, with one thread and a2048MiB Lean cap.
Both passed; launch resources were181.74GiB disk and7.02/6.93GiB
available RAM. Only propext and Quot.sound were printed; no sorryAx.
The dependency report retains the older CompactHVector cache's missing
capped-build ledger; no unrestricted dependency rebuild was launched.
Actual geometry and coverage remain explicit premises. Research and
complete written review preceded this final verification step.

The next settled result was mixed right-angled ceiling11. Root read the
complete new proof, inherited lower-count geometry and209-line Fraction
checker; a separate reviewer independently rederived the whole packet.
Bounded coefficient generation/replay passed0.086/0.084s,peak below56MB,
181.60GiB disk and6.52/6.49GiB available RAM. The certificate is3333B.
Only after those reviews and checks, `RightAngledMixed11` compiled in
0.56s and passed its separate six-theorem axiom audit in0.34s.
Launch resources were181.58GiB disk,7.47/7.49GiB RAM. It imports
only Std, has no legacy project dependency gap, and prints only
propext/Quot.sound with no sorryAx. All jobs were serial and used
normal resource gates. There are no disposable TeX or raster outputs
from these two new proofs; final source/certificate/kernel records remain.

The focused mixed11 source subsequently passed two complete manuscript
reads. Its final11-page/three-figure PDF has zero warnings; the last
build took1.649s and peaked below87MB. One earlier attempt was stopped
before TeX launch at5.8GiB available RAM; later jobs passed the6GiB
gate without changes. All actual figure pages were visually checked,
layout overlaps fixed,and all nine raster previews removed (1,941,885B).
One read-only process inspection identified another workspace's helpers;
no unrelated process was stopped. Final disk remained above181GiB.

The same compact arithmetic paper now has the scoped19 composition
appendix. A dedicated thin wrapper checks both TeX files and all48
declared source pins before and after the existing guarded builder.
The final18-page PDF has zero warnings. Prior eight-page source/PDF
and intermediate source-ledger evidence remain byte-preserved. The
common verification driver now routes its temporary PDF check through
that included-source wrapper; metadata-only plan checks passed for both
compact20 and mixed11. No optimizer, background farm or unrestricted
Lean rebuild was launched. Temporary build files are removed per batch.
