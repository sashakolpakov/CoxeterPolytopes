# Root acceptance: compact dimension at most22

20 October 2026. **Full written composition PASS.** Every compact
hyperbolic Coxeter polytope has dimension at most22, at the internally
reviewed written-proof level and the explicitly stated external-input
boundary. This accepts the previously held assembly; its candidate
wording is historical. It does not claim external refereeing, a new
example, complete geometric Lean, or a new executable enumeration.

Accepted source:
`compact_dimension23_exclusion_composition_candidate_20261020.md`,
SHA256 `8a691e4fb5e51f3256a02c68abd70a6bafad3e9dea99ddc45d345585b4d86edd`.
The separate complete accounting review is
`review_compact_dimension23_exclusion_composition_20261020.md`.
I read both entire documents, checked their rank reductions and input
interfaces, and independently recomputed the following exhaustion.

The accepted old ceiling23 excludes all n>=24 without requiring
arbitrary metric faces to be Coxeter. Its component caps and published
ambient29 input are independent of the new dimension23 argument. It
therefore suffices to exclude n=23 at one actual farthest vertex.

Call A_r,r>=7 and B_r/D_r,r>=9 dominant. If none occurs, finite
spherical classification and the previously proved ordinary rates give
energy greater than22/14 after discarding the possible single short A1,
or greater than9/100+21/14 with the sole rank-two defect. Thus a
dominant exists. Every dominant has energy>1/2, so it is unique before
the spectator rates are used. The exceptional possibilities are a single
short A1, or one of the two isolated rank-two defects, never both.
J discards only its short rank; its light A1 remains nonshort. All costs
are orthogonal projections of the same actual vector of norm squared<1.

For A the half bound forces r>=16, and the component cap forces r<=20.
The respective available bounds at ranks16,17,18,19,20 are
4/7,5/7,5/7,6/7,6/7, each strictly exceeded. After the worst short-rank
discard their totals are respectively1,15/14,1,15/14,1, again
strictly exceeded. A rank-two defect gives a larger budget by13/700;
fully ordinary spectators also give a stronger budget. The uniform
A>=8 theorem supplies the rank16 estimate without an extra A16 input.

For B/D the9/14 estimate forces r>=18. Generic prefix M=r-6
satisfies the D rank requirement exactly and the B requirement with
room. At r18,19,20 the lower floors .8689525,.88864,.90715
plus the weakest spectator budgets4/14,3/14,2/14 all exceed1.
At r21 the floor .9240775 plus the rank-two defect's>.09 exceeds1.
For two A1 spectators, at least one is nonshort and pays>1/10;
for an ordinary rank-two component its cost>1/7 is stronger. This
handles J as well. The genuine1/10 estimate is needed here; replacing
it with1/14 would leave a gap.

At r22 the prefix floor .9401975 excludes a nonshort A1 remainder.
The only short remainder is isolated I: J and a rank-two defect cannot
fit. The reviewed B22/D22-plus-I theorem excludes exactly this case.
At r23 there is no spectator. B23 is capped; the whole-D23 theorem
has precisely the required scope. No arbitrary-spectator extension of
whole D23 was used.

Thus all possible finite farthest forests in dimension23 are excluded.
The original compact23 tail then gives the new unconditional compact
ceiling22; arithmeticity is nowhere assumed. Here unconditional means
the internal geometric lemmas have written proofs and all case families
are closed, rather than a conditional subclass result. The named
published classification and compact29 theorem remain external inputs.

The three new component proofs have complete root and separate-agent
written PASS. D23's corrected source a77f88f5 changes only the premise
of an exact rational estimate; the original source and review are
preserved. The outside composition review checked the delta and all
current source pins. The new uniform A theorem additionally has a third
adversarial review. There is no unresolved case at these boundaries.

The finished compact TeX/PDF and its compiled accounting kernel still
state23. Preserve that exact publication/formalization snapshot; the
new22 result is supporting work for the same paper and has not yet
been integrated or compiled in Lean. The other five papers, all existing
PDFs, arithmetic20 and general finite-volume57 are unchanged. Do not
create a manuscript for this increment. The gap to compact examples8
and the broader dimension/example goal remain open.
