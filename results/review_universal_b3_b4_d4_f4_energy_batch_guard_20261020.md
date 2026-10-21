# Review of four universal finite-component energy bounds

20 October 2026. Full source and paper review only; no checker, census,
optimizer, Lean or TeX process was run.

Reviewed in full:
`results/compact_universal_b3_b4_d4_f4_component_energy_20261020.md`,
SHA256 `2eb274288127f1bd2b7df36499bec1da0fedc607d3c7ac8bd7734ef8f80c44db`.
Supporting universal small-A source:
`results/compact_universal_small_a_component_energy_20261020.md`,
SHA256 `801581c4b855289aca7982b2d0f4bb6f5849f151a78f40b86c9d403e84f415e6`,
also previously read in full and reviewed separately.

**PASS within the exact stated scope.** At an actual vertex farthest
from a strict interior point of any compact hyperbolic Coxeter polytope,
B3 component energy exceeds 1/4, and B4, D4 and F4 component energies
exceed 1/3. All other finite components may be arbitrary. Combined with
the pinned small-A proof, this excludes whole farthest-vertex forests
with factors A1,A2,A3,A4,B3,B4,D4,F4 in dimension at least 14. It does
not cover H3,H4, the general mixed dihedral interfaces, or large A/B/D
components, and is not an unrestricted dimension ceiling.

I independently checked all four inverse Gram matrices and all displayed
projection calculations. The adjacent dual correlations meet the
respective pair/logarithmic-concavity thresholds. The equation at energy
1/3 reduces at correlation 1/sqrt(2) to 51k>34, so connected adjacency
indeed forces one ACTUAL opposite wall for the entire component.

The multiple-contact arguments use finite retained white diagrams at
actual endpoints. In each path type, any pair of contacts other than
the two old endpoints leaves a cycle after deleting an omitted old
node. In D4 every connecting path omits a node. The three permitted
endpoint-pair quadratic expressions are correct and have the asserted
strict lower bounds. None uses the false assertion that every
nonpositive principal subset is forbidden in a compact chamber.

Every single-contact diagonal, minimum normalized-derivative ratio,
and resulting native-label cutoff checks. The q=1 cases give a
nonzero positive null combination of actual normals; all-wall
acuteness places the resulting ideal point in every halfspace, so
compactness excludes it even with additional spectators. In the q<1
cases the stated height bounds are below 1/4. Small-height isolation
then removes all spectator contacts, giving h=q-1<=0 as required.

The B3 vector K^-1*1 has positive entries, and its covering coefficient
9-2sqrt(6)-2sqrt(2) exceeds 5/4. The argument correctly uses intrinsic
energy for B4,D4,F4 rather than asserting an unavailable uniform
covering coefficient. Consequently both exceptional height-propagation
arguments still charge every other component more than its rank/12.
The baseline charges dominate rank/14; in the isolated-A1 debit case
at dimension 14, an odd-rank component must be A3 or B3, with surplus
1/28, which exceeds the remaining debit 1/35. Thus the expanded forest
conclusion follows without a new spectator inventory.

The underlying farthest-point, shadow and cross-component mechanisms
were previously co-developed with the author and root. This is an
outside-author audit of the new component calculations and their full
composition, not an independent-origin claim or formal verification.
No new flaw was found.
