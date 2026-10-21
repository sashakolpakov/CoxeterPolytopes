# H4-count compensation reduces but does not remove the model obstruction

Separating H4 compensation by the actual number of H4 components lowers the certified restricted-model obstruction to

**78271685370658/2232999997767 ≈ 35.05225502 > 35.**

All inherited bonus coefficients remain free. Each H4 pair-orbit compensation now has separate coefficients for actual H4 counts 1,2,3,4, in addition to the distinction between a facet mark inside or outside the paid H4. Compensation remains zero at count 5. This is a genuine enlargement of the [preceding free-bonus model](fresh_sqrt5_marked_free_bonus_proposal_20261016.md), with 30,260 variables.

The finite LP relaxation reaches numerical lower bound 35.05436701 after four solves and 4,461 marked-vertex cuts. The exact recovered dual has 1,051 terms and proves the rational lower bound above. It is a lower bound for the complete specified model, not a claim to its fully converged optimum.

The model retains all 1,979 old baseline-zero quadrilateral rows in 71,036 common contexts and the same 163 selected triangle families in 122,416 marked-context cases. Their H4-count distinctions produce 45,644 affine face profiles. Complete marked finite-forest separation remains in place, and every exact dual vertex is separately verified by literal root-pair summation.

The distance-two discount is 0 and the adjacent-3 discount is approximately 0.26359175 at the numerical obstructing LP. Both discount upper-bound dual marginals are 0, so loosening these inactive caps alone does not address this finite LP obstruction.

The dominant vertex remains 4H4+D5, with the facet root inside an H4, carrying approximately 59.36% of vertex-dual mass. The largest triangle terms now explicitly use core51_tail4 with external counts (A1,I2(5),H3,H4)=(1,0,0,3), and core46_tail5 with (0,0,0,3). Their facet marks often lie on a white A-type tail, whereas the obstructing vertex marks an H4. This identifies the next specific missing context: whether the facet root is in another H4 or in a different component, together with A1 absence/presence.

The [source](../scripts/fresh_sqrt5_marked_h4_count_free_bonus_20261016.py) generates and freshly replays the [certificate](fresh_sqrt5_marked_h4_count_free_bonus_20261016.json), including the exact dual, variable schema, literal witnesses, and numerical history.

**Scope:** an exact obstruction to this restricted payment model only. There is no new dimension bound, no compact arithmetic example, and no complete coverage of newly affected triangles or paid-corner quadrilaterals.
