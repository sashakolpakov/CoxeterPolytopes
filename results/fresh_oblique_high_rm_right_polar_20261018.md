# High R–M closure fails the necessary right-polar test on the declared N/S domain

The remaining high R–M hypothesis gives no native completion on the20
high-N geometries and their complete low-label S banks. Exactly192 S
branches have S–M right. For each, both possible profile0 R normals
orthogonal to X,N,S are reconstructed. Every branch has

    1.0236068047 <= |R·M| <=1.4572876741 <sqrt(3).

Thus none can have a finite R–M bond of label at least6. All192 branches
with a negative common-center product also fail native contacts with
L68 and M71. No new cone or whole-face closure stage is reached, and
no finite-volume example is produced.

The [producer](../scripts/fresh_oblique_high_rm_right_polar_20261018.py),
[certificate](fresh_oblique_high_rm_right_polar_20261018.json) and
[standalone saved checker](../scripts/check_oblique_high_rm_right_polar_saved_20261018.py)
bind the frozen [high-N/right-R packet](fresh_oblique_high_n_right_r_join_20261018.md).
The result has no imposed upper cutoff on R–M labels within its declared
N/S domain. Its new geometric inputs and conclusions remain numerical;
it is not an exact interval exclusion for actual reconstructed normals.

## Why the polar equations are necessary

The inherited high N–X finite-corner argument forces R–X and R–N right.
Now suppose R–M has a finite native label at least6. The physical norm2
principal matrix on M,R,S has determinant

    8−2(a²+b²+c²)−2abc,

with a²≥3 for R–M and b,c the nonnegative S–M/S–R magnitudes. Either
nonright native bond has magnitude at least1, making the determinant
at most0. A finite M/R/S corner is strictly positive definite, so both
S–M and S–R must be right. The same principal argument with a fixed
white root shows R cannot have any nonzero profile contact: M has
profile0, and a high R–M bond together with any such native contact
already violates positive definiteness. Hence R has profile0.

These necessary conditions are therefore

    R perpendicular to X,N,S; S perpendicular to M; R profile0.

They apply before considering any other original wall. The literal disk
remains the inherited12-triangle disk, boundary8–10–N–R–S–52–8,
V/E/F=10/21/12. No new adjacency is assumed and no old wall is released.
In particular the finite X/M/R, X/N/R and M/R/S corners are the genuine
duties underlying the argument.

## Complete declared S domain and reconstruction

The frozen high-N packet has two distinct J/L/X/M parent instruction
classes retaining all four selected origins. On each it considers
Nprofile0, labels[3,2,m] and[2,2,m] in frame10/J/X for
m=6,7,8,10,12. The present producer reconstructs those20 accepted N
normals from their complete parent chains. Earlier R and S are recipe
ancestors only and are removed from the fixed set.

For each N the frozen complete seven-profile2..5 S bank in frame52/L/M
already tested the genuine52/L/S and L/M/S corners, every retained
native product and the full common center. Restricting that bank to
S–M=2 leaves exactly192 of712 native S normals. The saved subset indices
are checked against the full frozen bank; no favorable S branch is
selected within this necessary subset. Every S is reconstructed from
its saved norm branch before any R calculation.

For each parent, each [3,2,m] family contributes ten S branches.
The [2,2,m] families contribute9,9,9,9,10 branches for the five labels.
That is96 per parent and192 in total. A completion with S outside the
frozen low bank is not covered. Higher S–52/S–L labels, different N
labels and other parent geometry remain outside this diagnostic.

All58 original facets remain fixed:52 projected originals plus six
fixed original whites. J,L,X,M,N,S add six new fixed facets, so R is
checked against58 explicit projected/product gates and six white
contacts,64 fixed facets in total. Original-side profile offsets and
all new profile offsets remain distinct.

## The whole polar line, not a label sweep

Write V for the three projected frame normals X,N,S and G for their
projected Lorentz Gram matrix. Every one of the192 saved frames has
signature(2,1) and rank3. Its orthogonal complement is therefore a
one-dimensional positive line. Since R has profile0, its physical
norm is its projected norm2, and its three required physical zero
products are also projected zero products. There are exactly two
opposite possible normals on that line.

Let c be the covector of projected products of M with X,N,S. Since
M also has profile0, its residual square is

    q = 2−c^t G^-1 c.

The positive one-dimensional complement gives the necessary identity

    (R·M)² = 2q.

The producer saves G, its inverse, inertia, c, q and both norm branches.
The checker independently constructs the null line of the three linear
equations using signed3×3 minors of V times the Lorentz metric. It
normalizes that cofactor vector to square2, obtains both signs and
compares the full vectors to the producer branches. It also recomputes
the residual identity and every retained native/center product.

Every saved frame is separated from singularity; no singular, null or
higher-dimensional orthogonal-complement case is silently called an
ideal corner. Such cases would be recorded unresolved by the producer.
There are zero unresolved frames here. The one-dimensional argument
covers any R–M label for each fixed frame without enumerating labels.
The numerical inequality2q<3 then rules out the proposed high finite
R–M label on these saved configurations.

All384 branches fail the finite-high R–M condition. Exactly one branch
per frame has a negative center product. Its R·M product is positive
and its only native blockers are L68 and M71. The opposite branch has
negative R·M of the same magnitude, a positive center product and many
native blockers. Neither sheet is dropped before checking. The full
branch records, blocker lists and per-N polar digests are saved.

## Accounting and unexecuted stages

The old deletion set is still30/44/63/64/65. The exact199-support ledger
loses19 finite-support certificates,12 compact two-face certificates
and one compact three-face certificate; three lost finite supports lie
on8. The source and checker reconstruct this debt independently of the
new normal calculations. This does not assert that all remaining old
ambient vertices survive a new numerical cut.

No R passes the necessary/native/center gates. The source retains the
full12-corner, target-cone/positive-definite-incidence and whole-interior-
face gates for an admitted completion, but none is executed here. No
new boundary ledger or open-end comparison is generated by this packet.
The preceding20 closest-native diagnostics remain24-finite/14-end
nonclosing cones; they are not improved by a failed polar candidate.

Together with the preceding complete low R–M bank, this covers the
missing high R–M possibility only for the declared20 N geometries and
low S domain. It does not prove an all-label impossibility for the
broader construction problem. A subsequent useful change must alter
one of those geometric premises or inspect a genuinely different
remaining S family; increasing an R–M label cutoff cannot change this
one-dimensional calculation.

## Verification

Native-first preflight90339/chunkb5c26d passed. Authoritative generation
99278/chunkd8afc1 and separate fresh replay12989/chunk96c7b6 both
terminate exit0/PASS. The final standalone saved checker chunk468132
passes192 cofactor lines,384 full branch records, all selected-domain
indices/digests and the exact19/12/1 old debit. The root agent completed
a full mathematical producer/report/helper read (919ae3), finding no
correction. Its scope was a source/report/helper read, with no execution
or independent bank or physical-interval replay. The packet is final
and frozen after direct-binding, AST and local-link checks.

The saved checker uses NumPy and the standard library, with no project
imports. It independently reconstructs the R polar lines from the saved
N/S frames. It checks the saved N/S norm equations and complete bound
provenance, but does not independently solve the preceding N/S norm
quadratics, regenerate their old low banks, reconstruct cones, repeat
the spherical classifier or provide physical interval certificates.
No exact new normal, ambient H9 volume or new example is claimed.

    python3 scripts/fresh_oblique_high_rm_right_polar_20261018.py --verify
    python3 scripts/check_oblique_high_rm_right_polar_saved_20261018.py
