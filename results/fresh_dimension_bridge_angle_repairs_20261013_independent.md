# Independent H20 cap-bank audit and common escape

**PASS.** Separate exact arithmetic replays all 107 cap entries and proves that one strict future spacelike vector survives **all 82 base walls and all 107 caps simultaneously**. Therefore no subset of this cap bank can produce finite volume while retaining the 82 base walls.

The source is [the exact cap bank](../data/enumerations/fresh_dimension_bridge_angle_repairs_20261013_exact.json), tied by byte hashes to its scout and the frozen 82-wall construction. The audit imports only the earlier independent $\mathbb Q(\sqrt6)$ verifier and implements a separate recursive quadratic-tower arithmetic. It uses no numerical geometry, symbolic algebra package, or construction producer.

The replay verifies:

- the original 70-wall lift and all 12 subsequent polar additions;
- every new cap's positive norm, strict common-center pairing and advertised target cut;
- all 8,774 cap–base pair relations, with exact finite-angle polynomial identities and isolating intervals;
- the 19 orthogonalities and prescribed anchor angle for each cap, with rank 20 checked for all 40 distinct linear systems;
- that all 14 distinct target rays are still feasible in the frozen 82-wall cone;
- strict negativity of all 189 pairings with the common escape.

The cap-entry counts by anchor label are 14 for label 2, 7 for 3, 11 for 4, 11 for 5, 20 for 6, 12 for 8, 16 for 10 and 16 for 12. Entries may represent repeated mirrors; no mirror-deduplication claim is made.

The common escape is the saved strict escape of the 82-wall source. Its exact norm and common-center pairing are

\[
v^2=\frac{3849}{32768}>0,\qquad
(v,h)=-\frac{14839}{128}<0.
\]

The center has norm $-1518$. Every cap pairing is checked by exact nested-radical comparisons and is strictly negative. An interpolation with the strict timelike center reaches an interior point of the ideal boundary, proving infinite volume of the full halfspace intersection and every larger intersection obtained by taking a subset of the caps. No compatibility assumption among the new caps is needed for this obstruction.

The [audit artifact](fresh_dimension_bridge_angle_repairs_20261013_independent.json) stores all 107 exact common-escape pairings, the source hashes and all counts. Finite-angle checks avoid division in radical towers, so coincident or dependent radicals cause no field-independence assumption. Explicit dependent-radical sign cases pass, and a changed-coordinate fake normal is rejected.

```sh
python scripts/fresh_dimension_bridge_angle_repairs_20261013_independent.py --verify
```

Checker: [independent cap audit](../scripts/fresh_dimension_bridge_angle_repairs_20261013_independent.py). This closes this finite bank only. It does not exhaust possible caps, classify new–new angles, or exclude a completion using different normals or replacing base walls.
