# Review of the symmetric 862-facet branch

Source/proof review and a fresh run of

```sh
python scripts/fresh_h22_one_side_20261013.py --verify
```

both passed. This review reuses the existing interval implementation and the hash-bound 835/839 norm and contact certificates; it does not independently reconstruct every predecessor cap.

The [checker](../scripts/fresh_h22_one_side_20261013.py) correctly takes the 25 existing cap types, chooses their positive added-coordinate representatives, and adds both signs. Reflection preserves `diag(Q,2)` and explicitly permutes the 812 rational parent walls. The literal index map includes every old 839 facet. All 1,225 signed-cap pairings are covered: 1,217 rigorous strict dotted gates, six exact inherited right contacts and two exact parallel contacts from null old projections. Positive cap norms and a strict common center, together with the complete acute pair gate, establish 862 actual facets.

The norm-1,586 strict spacelike escape is checked against all 812 parents and all 50 signed caps; the chamber remains infinite volume. The same literal old cusp is strict on every cap. The two previously bounded links retain exactly their paired incident caps, and all added nonincident cap signs are strict.

The old compact quadrilateral also survives: every newly added mirror is dotted to its supporting old cap 17. Hence none intersects that supporting hyperplane, and the compatible orientation places the entire facet on its allowed side. Reflection supplies the second copy. This preservation is a geometric consequence of the complete pair gate, not a new whole-chamber vertex census.

The general one-side reduction, including the pure-vertical exception, is justified in the [extension feasibility report](fresh_h22_extension_feasibility_20261013.md). Neither the reduction nor this application proves finite volume or a new-dimension example.

Reviewed source SHA256: `e4b6ca471a9f311cf87186943314e50b3e2fd7d61ef80a17f048d0628dc025d5`.

Replayed [certificate](fresh_h22_one_side_20261013.json) SHA256: `908f40b81dc0259d79058f2139b06d2cfbfbaae3ba3e1b91c9a8229cd0805823`.
