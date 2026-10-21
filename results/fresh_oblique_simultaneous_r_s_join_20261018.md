# Simultaneous R/S low banks do not close the M corner

Removing both prior R and S before testing their replacements gives
**80 native R candidates and142 native S candidates** on the four fixed
N geometries. Every one of their2,844 pairings is checked:2,838 are dotted
and six are nonnative. No pair supplies a finite M/R/S corner. There is
no completed disk or new finite-volume example.

The [producer](../scripts/fresh_oblique_simultaneous_r_s_join_20261018.py),
[saved result](fresh_oblique_simultaneous_r_s_join_20261018.json) and
[standalone saved checker](../scripts/check_oblique_simultaneous_r_s_saved_20261018.py)
preserve the frozen [joint N/R pilot](fresh_oblique_joint_n_r_pilot_20261018.md).
They test a simultaneous two-wall join, so the old S no longer filters
the R bank and the old R no longer filters the S bank.

## Fixed geometry and independent candidate domains

The four prefixes are exactly the previously selected forward indices
140,133,166,348, with N/X labels3,4,5,5. The complete eight-step chains
are reconstructed, including both former R normals and the former S
needed as recipe anchors. Both R70 and S72 are then removed. J67,L68,
X69,M71,N73 and all retained original walls remain fixed.

There are52 projected original normals and six fixed original whites,
hence58 original facets. Five fixed new normals give63 fixed facets.
Each candidate receives57 explicit product tests (52 originals and five
new normals), plus the six original-white contacts through its profile.
The other moving wall is the only intentionally deferred product.

R is solved in the physical X/M/N frame with finite duties X/M/R and
X/N/R. There is no M/N edge assumption. Independently, S is solved in
52/L/M with finite duties52/L/S and L/M/S; there is no52/M edge.
Each bank contains every seven-profile label triple in{2,3,4,5} allowed
by those two actual finite duties. The inherited list is independently
compared with the entire duty-filtered cube.

|Prefix|R recipes|S recipes|Native R|Native S|Pairs|
|---|---:|---:|---:|---:|---:|
|140|50|42|0|35|0|
|133|38|42|36|35|1260|
|166|26|42|22|36|792|
|348|26|42|22|36|792|

The308 recipes give600 real norm branches. Eight R recipes have no real
root in the numerical quadratic calculation; all other recipes have both
real branches checked. Norm/frame residuals, profile offsets, native
spacings and strict common-center signs precede bank admission. Identical
rounded candidates within a bank may be deduplicated; no different frozen
prefix is merged by numerical proximity. The saved branch records retain
the entire declared domain and its digests.

The candidate normal calculations remain numerical. Symbolic high families
are recorded but unsearched. This is a complete low-cube two-wall join on
four specified geometries, not a search over all1,100 otherwise-native
N branches or a continuous parameter domain.

## Complete pair test and diagnostic cones

For each prefix, every Cartesian R-bank/S-bank pair is tested with its
full physical product, including the new/new white-profile offset.
An actual finite M/R/S corner requires a native finite R/S bond and a
positive-definite corner Gram. None passes even the finite-bond gate.

The six nonnative products occur in the last two prefixes, three in each:

    −1.1771803392, −1.3053058122, −0.9719724788.

The first two lie strictly between the simple and double native bond
magnitudes; the last lies in(−1,0). Their closeness to a native value is
not accepted as an angle relation. Every other pair is dotted.

Full closure candidates would next receive all mutual native/center
tests, all12 prescribed corners, complete target-cone/PD agreement, and
the whole-face ray sets of every recorded wall. That accepted-join stage
is empty. Separately, one closest native pair per nonempty Cartesian bank
is scored as an explicitly marked diagnostic. Thus exactly three joined
cones, not all2,838 dotted pairs, are reconstructed and compared with PD
incidence. Every diagnostic has24 finite rays and14 open ends.

Their whole compact target faces are exactly8,11,J,L,X, with the complete
J,L,X ray sets matching their disk links. M,R,S,N remain open. The sole
failed prescribed triangle is M/R/S. All three diagnostics share the
exterior supports

    (7,47,55), (10,33,N), (33,55,59), (33,55,N), (33,59,66),
    (47,52,N), (47,55,N), (52,N,S), (58,60,66), (59,60,66),
    (M,R,S), (N,R,S).

The corresponding one-open/two-open edge ledger gives14 ends. Four
deduplicated boundary ledgers cover the original control, the two
removed-wall prefix types and these diagnostics. Before adding R/S,
the first prefix has21 finite target rays and13 ends; the other three
have20 and14. The original63-wall control has26 and10.

## Unchanged topology, exact debit and scope

The same12-triangle disk has boundary8–10–N–R–S–52–8 and
V/E/F=10/21/12. Its interior J,L,X,M links are cycles and its boundary
links are paths. The source and saved checker verify this literally.
The prior topology-only peer review applies to this unchanged disk.

Deletion30/44/63/64/65 still loses19 of the original199 finite-support
certificates,12 compact two-face certificates and one compact three-face
certificate, including three finite supports on8. Both programs
reconstruct this exact old debit. It does not guarantee survival of all
remaining ambient finite vertices after a new numerical cut. Target-face
finite-ray counts and ambient H9 finite-support counts remain distinct.

Joint low-bank construction can bring R/S nearer a finite native bond
than the previous fixed-S test, but gives no native finite join and no
reduction in the diagnostic open-end count. The result does not exclude
high-label N motion with a simultaneously altered R/X relation; the
forward high-N obstruction required frozen nonright R/X. That different
mechanism is not executed here and would need a separate packet.

## Verification

Generation84457/chunked2718, separate fresh replay5067/chunk7086e1,
and standalone saved checker chunk20bcb9 terminate exit0/PASS.
The root agent completed a full producer/report/helper mathematical and
scope read (554afb/a1bd10): all chains remove R/S before the gates, both
independent low cubes and all2,844 pairs are covered, the three diagnostics
are separate from finite joins, and the debit, unchanged disk and stated
numerical limits agree. This internal review did not execute the programs,
solve normals, independently prove no-real statuses or reconstruct cones.
The packet is final and frozen after its direct-binding, AST and link audit.

The standalone checker uses NumPy and the standard library, with no
project imports. It verifies both removed-wall chains, complete low-cube
domains, all600 branch values and both bank digests at every prefix,
all2,844 pair values/digests, empty finite join set, the deterministic
three diagnostic choices, whole-face and stored cone/PD agreement,
the disk links and exact original-support debit. It checks stored corner
positive definiteness numerically. It does not solve norm equations,
independently prove no-real statuses, reconstruct cones, repeat the
historical scaffold/classifier or enclose new normals by intervals.
No ambient volume certificate is claimed.

    python3 scripts/fresh_oblique_simultaneous_r_s_join_20261018.py --verify
    python3 scripts/check_oblique_simultaneous_r_s_saved_20261018.py
