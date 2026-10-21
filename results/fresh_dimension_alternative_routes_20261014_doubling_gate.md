# Why doubling the H6 source cannot open the H3-to-H4 lift

**Exact scoped obstruction.** No finite sequence of legitimate single-wall Coxeter doublings of the known ten-facet compact H6 source can produce a facet right-angled to all its intrinsic neighbors. Combined with the [H3 face-lift gate](fresh_dimension_alternative_routes_20261014.md), this rules out using those doubles to introduce a nonorthogonal H3-to-H4 source-face facet. It does not rule out all compact H9 extensions or more general gluing, reflection subgroups, rolling, or angle changes.

The [exact checker](../scripts/fresh_dimension_alternative_routes_20261014_doubling_gate.py) and [certificate](fresh_dimension_alternative_routes_20261014_doubling_gate.json) verify the finite diagram premise and independently form both first doubled Gram matrices. Each first double has 12 facets, all pair relations remain Coxeter, and each facet still has the required protected neighbor. Generation and saved replay pass. The assertion for arbitrarily many finite doublings follows from the invariant below, not from a bounded experiment.

## A general invariant under legitimate doubling

Call a facet protected if it has a finite odd-label neighbor. Suppose every facet of a compact Coxeter polytope is either protected, or has a finite nonright neighbor that is protected.

A protected facet cannot be a doubling wall: doubling a finite angle `pi/m` with odd m does not yield a Coxeter angle. A legitimate doubling wall has only even meeting labels. The double consists of the old chamber and its reflection; its facets are the surviving old facets and their reflected copies, with orthogonal copies merged. This is the usual doubling construction. [Allcock, Lemma 1.3](https://arxiv.org/pdf/0903.0138).

Every odd-label ridge between old protected facets survives doubling any other facet. Indeed the relative interior of the old codimension-two ridge has a point off the doubling wall; near this point the double is identical to the original chamber. Compact Coxeter polytopes are simple, so the third facet cannot contain the entire ridge. The same observation preserves any specified finite nonright ridge between a surviving unprotected facet and a protected facet. Reflection gives the identical assertion on the other half. Merging coplanar facets does not erase the old ridge's local geometry.

Consequently all surviving or reflected protected facets remain protected, and every other new facet retains a nonright protected neighbor. This proves the invariant by induction. In particular no all-right facet can appear at any finite stage.

## Literal source and both first doubles

The source's facets 1–6,9,10 all have an odd finite neighbor. The only remaining facets, 7 and 8, each meet facet 2 at `pi/4`; facet 2 is protected by its `pi/5` ridge with facet 1 and its `pi/3` ridge with facet 3. Thus the invariant holds at the start. The only eligible first doubling walls are exactly 7 and 8.

The checker uses the published finite diagram and the explicit unit-normal divergent entries

```text
G(6,7)=G(8,9)=-(2sqrt(2)+sqrt(10)),
G(7,8)=-(17+8sqrt(5)).
```

For doubling wall k, an old normal f reflects to `f-2(f,k)k`; right neighbors are retained once and every other surviving facet receives its reflected copy. This gives the two literal 12-by-12 exact Gram matrices saved in the certificate. No numerical angle fitting or optional facet deletion is used. [Burcroff, Figure 5](https://arxiv.org/pdf/2201.03437), [machine-readable source diagram](https://raw.githubusercontent.com/johnamcleod/10-facets-6-d-polytopes/main/verify_polytope.py).

The useful next change would therefore have to alter the protected odd-ridge structure or use a different intrinsic source. Repeating valid doublings of this source, even without a preset depth limit, cannot provide the all-right facet required by the coupled H3-face mechanism.
