# Review of the changed H20 assembly's cusp diagnostic

**PASS:** source-level mathematical review and a fresh saved-result replay of [the cusp checker](../scripts/fresh_h20_changed_assembly_20261013_cusps.py). This review reuses its arithmetic implementation; it is not a second field-arithmetic implementation or a new cusp enumeration.

The script reconstructs the original and reflected source rays, primitive-normalizes and deduplicates them to 204 directions, and checks their exact null norm, future orientation and all 179 wall signs. Thus its feasibility statements do not rely only on the original H19 source incidences.

For each of the 17 all-short cases, it checks actual normal-span rank 19 and decomposes the incident Gram into components. Every component has a strictly positive kernel, rank one below its size, and an explicit positive weighted normal sum equal to a positive multiple of the cusp ray. The component ranks sum to 18. Positive semidefiniteness follows independently from incidence at a future null ray in the certified Lorentzian form; these vectors lie in its positive-semidefinite orthogonal hyperplane.

The stated deficiency is substantive: an H20 horosphere has tangent dimension 19, but the incident normals have quotient rank 18. The explicit vertical vector has norm 3, is orthogonal to the cusp ray and to every incident wall, and is therefore a nonzero recession direction on the horosphere. The saved strict timelike interior supplies a strict point in the cusp link. Nonincident walls have strictly negative pairing with the cusp and do not meet a sufficiently deep horoball, so they cannot remove this recession direction from the asymptotic link.

Consequently these 17 links are unbounded. The artifact carefully distinguishes actual normal-span rank from quotient rank and makes no boundedness claim for the other 187 sampled links. Its Lorentz signature, common strict interior and root definitions depend explicitly on the frozen 179-wall assembly certificate.

```sh
python scripts/fresh_h20_changed_assembly_20261013_cusps.py --verify
```

Result: `PASS_204_INHERITED_CUSP_SIGNS_AND_17_DEFICIENT_LINKS`. The [cusp artifact](../data/enumerations/fresh_h20_changed_assembly_20261013_cusps.json) and [main report](fresh_h20_changed_assembly_20261013.md) match this scope.
