# Dimension-21 arithmetic weight diagnostics

**No improved dimension bound was obtained.** Nine enlarged weight models have exact dual lower bounds above the required threshold 40. These are obstructions to the stated weight families, not nonexistence results for dimension-21 polytopes. The [synthesis checker](../scripts/fresh_sqrt5_dimension21_diagnostic_synthesis_20261016.py) and [exact output](fresh_sqrt5_dimension21_diagnostic_synthesis_20261016.json) independently recheck literal rational weak duality, inspect recurring support, and enumerate finite endpoint extensions. They rely on the frozen producer replays for reconstruction of face and forest coefficients.

All nine producers retain the nonnegative baseline, complete finite-queue singular obligations, inherited parameter-marked and dotted obligations, complete common-spectator enumeration, and exact finite-forest separation. All sources, certificates and reports are frozen after generation and saved replay.

| Fresh model prefix, with common prefix `fresh_sqrt5_dimension21_` and suffix `_20261016` | Certified lower bound exceeds |
|---|---:|
| [h3_cross](fresh_sqrt5_dimension21_h3_cross_20261016.md) | 40.25 |
| [e7_context](fresh_sqrt5_dimension21_e7_context_20261016.md) | 40.36 |
| [e8_context](fresh_sqrt5_dimension21_e8_context_20261016.md) | 40.07 |
| [e7e8_context](fresh_sqrt5_dimension21_e7e8_context_20261016.md) | 40.26 |
| [e7e8_exact](fresh_sqrt5_dimension21_e7e8_exact_20261016.md) | 40.23 |
| [a4_cross](fresh_sqrt5_dimension21_a4_cross_20261016.md) | 40.08 |
| [a4_exceptional](fresh_sqrt5_dimension21_a4_exceptional_20261016.md) | 40.14 |
| [a4_pairs](fresh_sqrt5_dimension21_a4_pairs_20261016.md) | 40.10 |
| [a4_pair_exceptional](fresh_sqrt5_dimension21_a4_pair_exceptional_20261016.md) | 40.07 |

These are nonoptimal dual certificates obtained at different separation stages; their numerical values should not be compared as model optima. The two largest models each have 472,694 variables and 120,996 common-spectator cases. No numerical optimization remains active.

The repeated obstruction involves exceptional components. `D5 + 2E8` and `A1 + D5 + E7 + E8` each appear among the eight largest forest multipliers in six models; `A1 + D5 + D7 + E8` appears in five. The face rows `old112:967`, `982`, `1107`, `958`, `984`, and `dotted:107` appear among the eight largest face supports in every model. They repeatedly couple E7/E8 vertices with D-type and I2(5) corners. Nonzero dual support does not assert global geometric realizability or optimality.

Merely adding further disconnected cross families is insufficient with the old contexts fixed. For each A4 model, the checker removes every dual face term containing a previously unpaid disconnected marked pair. The resulting dual remains above 40. In the final model, removing all nine D4–D4 terms leaves

\[
\frac{40672934387741}{1015000000000}>40.07.
\]

Any added nonnegative variable that vanishes on the retained faces and contributes nonnegatively on forests preserves this dual: its new dual column has nonpositive total. This rules out that specific extension strategy. It does not cover refined contexts, richer dependence within already paid pair families, baseline changes, or edge-dependent corrections.

The smallest concrete weight enlargement suggested by the census is to distinguish the **source root orbit** in existing A4 crosses. A4 has two root orbits, each of size two. Of the 75 destination-root orbits in the complete input, 44 occur with only one source orbit; 31 occur with both. Replacing `v_destination/(4d)` by `v_source,destination/(4d)` contributes `v/2` per destination root for each size-two source orbit. Setting an absent source bonus to zero therefore halves those 44 portions of vertex cost while preserving these input face payments. This uses the existing counts and remains additive; it has not been optimized or checked against a new global coverage binder. H3 similarly has three size-one source orbits, but all three occur for each of its six destination orbits.

The missing structural information is compatibility between actual neighboring vertices. The finite-diagram census has 68 component types and 311,574 rank-21 forests. E8 cannot be absorbed into a connected finite extension by adding one root; E7 can only extend to E8. Two specific rank-20 edge diagrams illustrate the resulting finite endpoint restrictions:

| Edge diagram | Possible finite endpoint types |
|---|---|
| `D4 + 2E8` | `A1 + D4 + 2E8`; `D5 + 2E8` |
| `D5 + E7 + E8` | `A1 + D5 + E7 + E8`; `D5 + 2E8`; `D6 + E7 + E8`; `E6 + E7 + E8` |

Every examined deletion still permits the original vertex type as a finite reextension. Thus finite-type lists alone exclude none of these vertices. A useful next structural computation is to classify which paired endpoint completions admit actual Lorentzian arithmetic edge Grams over Q(sqrt(5)), then impose their incidence constraints on vertex-type frequencies. The present per-vertex maximum contains no such global incidence accounting. Finite completion, Lorentzian signature, and finite-volume realization must remain separate checks.

Synthesis generation and saved replay passed, including 60 recursive manifests and 1,391 source references. No separate global arithmetic theorem is asserted.

```sh
python3 scripts/fresh_sqrt5_dimension21_diagnostic_synthesis_20261016.py --verify
```
