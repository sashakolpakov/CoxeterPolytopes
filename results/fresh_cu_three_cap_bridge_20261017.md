# Three-cap disk with a moving A companion

The bounded scout finds no compatible left cap, so no three-cap closure is
constructed. It tests one continuously moving companion family rather than
reusing the eight numerical A clusters of the preceding four-cycle search.
The exact63-wall/199-finite/10-open-end infinite-volume control is unchanged.

Take the zero A2-profile moving family with native labels A7=2,A30=3,AC=2.
The analytic moving-frame formula from the
[event search](fresh_cu_companion_contact_events_20261017.md) gives a
parallel A47 threshold u approximately3.28277304328. On the sampled side
u above that threshold, A47 is dotted. All63 original walls, C(u), A and
the original strict center are retained.

The proposed boundary is7–47–52–30–A–7. R owns the two boundary edges at7,
L the two at30, and D the remaining47–52 edge. The nine triangles are

    7/47/R, 7/R/A, A/R/L, A/L/30, 30/L/52,
    52/L/D, 52/D/47, 47/D/R, R/D/L.

They form a disk with five boundary vertices, three interior vertices,
16 edges and9 triangles. The cap duties are distinct; this is not an attempt
to add two endpoints to one old edge. A hypothetical completed disk would
still require actual full-cone incidence and all native contacts.

R uses frame7/47/A and L uses30/52/A. The complete symbolic finite-corner
classification gives24 R and23 L recipes. High labels are sampled only at
6,8,12, producing60 R and46 L norm branches per parameter. Forty-eight
logarithmically spaced values between threshold+0.0001 and32 are tested.
The plan was to bracket a finite R–L contact, solve its parameter continuously,
then solve D and test every remaining finite corner and the complete cone.

Among2,880 R branch tests,1,001 pass all native contacts and the strict
center. Every one of the2,208 L branches fails a native gate;244 are blocked
only by old47 and339 by old7 together with old47. No branch is blocked only
by the new C wall. Thus there are no R–L brackets and no D stage in this run.
There are no unresolved accepted roots or cones; individual R bank entries
are necessary native-admissible prefixes, not claimed complete chambers.

This failure suggests changing L's old47 interaction or its endpoint duties,
rather than increasing the terminal D labels. It is not an exact obstruction:
other A families, parameter values between the samples or above32, isolated
native contact events, tangential mutual-contact roots and unsampled high
labels remain open. The numerical norm solver also discards near-singular
or near-tangent frames; its return values are not a complete root-coverage
certificate. No native-compatible R–L pair exists in the sampled
domain, so the planned continuous root solver is never invoked.

The old deletion debit is zero. The source contains full native/common-center
and complete target-cone/independent-PD checks for any event or third cap that
would survive; those later stages produce no accepted output in this run.
No numerical wall is promoted, and no full ambient finite-volume criterion
is asserted.

[Source](../scripts/fresh_cu_three_cap_bridge_20261017.py) and
[saved data](fresh_cu_three_cap_bridge_20261017.json) contain the literal disk,
recipe lists, parameter grid and native blocker counts. **FINAL/FROZEN:** primary generation
and a separate fresh saved replay both terminated exit0/PASS. A separate
literal topology check confirms V8/E16/F9, the five boundary edges, interior
edge multiplicities, interval links at boundary vertices and cycle links at
R/L/D. That review does not cover numerical geometry or parameter coverage.
Root's full narrow source/report read passed the distinct duties, frame
contacts, native offsets, common center, bounded grid/bracket logic and
conditional later-stage cone checks. It did not repeat the numerical root
or cone calculations or assert complete continuous coverage.

    python3 scripts/fresh_cu_three_cap_bridge_20261017.py --verify
