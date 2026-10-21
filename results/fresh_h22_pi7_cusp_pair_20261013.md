# Two bounded cusp links in a coherent 837-wall H22 chamber

Two exact new facets repair the two unbounded ideal links identified by the [835-wall edge frontier](fresh_h22_pi7_edgewalk_20261013.md). The resulting chamber has **837 compatible actual facets and two certified bounded Euclidean cusp links**. Its global volume is still infinite: the previous strict spacelike escape survives both additions.

The [exact certificate](fresh_h22_pi7_cusp_pair_20261013.json) and [checker](../scripts/fresh_h22_pi7_cusp_pair_20261013.py) pass.

The new walls are the last-coordinate reflections of cap indices 4 and 5 (candidates 271 and 306). In the established coordinate convention the form is \(\operatorname{diag}(Q,2)\). Each original cap has last coordinate −1 and null old projection, so its reflected partner has last coordinate +1, the same norm two, and pairing −2 with the original. They form an exact parallel pair.

Reflection in the last coordinate preserves the form, fixes the strict center, and permutes all 812 rational parent walls. Their angle relations are therefore inherited exactly. Each new cap is strictly dotted to the other 22 old radical caps. The two new caps are also strictly dotted to one another. These checks establish all new Coxeter pair relations and preserve the prior chamber facets.

The two ideal null rays have last coordinate zero, so their pairings with each reflected cap agree with the original cap pairings. No other cap becomes incident. At each endpoint the 22 active rational parent walls split into affine components of ranks 8 and 12. The formerly isolated cap now has its parallel partner, adding an orthogonal affine A1 component. The complete active Gram therefore has inertia \((21,0,3)\).

This proves more than an affine-subdiagram match. The active quotient normals span the entire 21-dimensional Euclidean horosphere, and every component has a strictly positive kernel vector. If a vector lies in the recession cone of all active halfspaces, the positive relations force its pairing with every normal to vanish. Spanning then makes the vector zero. Thus each complete horospherical link is bounded. The previously certified null rays and all 837 wall inequalities ensure these are actual ideal endpoints of the chamber.

The old rational strict escape remains strictly feasible for every wall. Its norm is

\[
\frac{3956990142533}{12500000000000}>0.
\]

Consequently these are two local cusp completions, not a finite-volume H22 polytope. The coupled finite-edge construction can use the frozen 837-wall certificate as its compatibility base.

Reproduce from the repository root:

```sh
PYTHONPATH=src:. python scripts/fresh_h22_pi7_cusp_pair_20261013.py --verify
```

The checker uses exact reflection and Gram identities plus rigorous intervals for strict radical inequalities. It inherits the verified 835-wall roots and the separate exact ideal-ray certificates through source hashes. It does not assert that all ideal vertices or all outgoing directions of the enlarged chamber have been enumerated.
