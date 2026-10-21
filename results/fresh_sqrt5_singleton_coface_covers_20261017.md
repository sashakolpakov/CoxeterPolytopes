# Singleton mixtures already cover every individual coface slice

For all 34,520 exported individual-coface slices, the sum of the current
nonnegative base demands is at most the slice capacity. Place each demanded
base alone with its demanded mass and assign the remaining capacity to
the empty configuration. This is an explicit exact marginal cover.

The [checker](../scripts/fresh_sqrt5_singleton_coface_covers_20261017.py)
reconstructs the demands from actual provenance and rational face masses;
the [certificate](fresh_sqrt5_singleton_coface_covers_20261017.json)
stores every singleton mixture and its nonnegative unused capacity. Of the
870 positive slices, 833 have one positive block, 35 have two and two have
three. No optimizer or producer import is used. Generation and fresh
saved-data replay both pass.

This strengthens the interpretation of the
[complete disjoint-base cone calculation](fresh_sqrt5_complete_coface_base_cones_20261017.md).
Further incompatibilities among two or more bases inside a single exported
slice cannot exclude these particular marginals, provided each individual
base option remains allowed. Singleton configurations avoid every such
joint incompatibility. This assertion concerns the independently averaged
slice marginals; it does not construct a simultaneous global distribution.

The [cross-coface cap capacity](fresh_sqrt5_cross_coface_cap_capacity_20261017.md)
does exclude the old witness, by coupling different cofaces of the same
simplex. Its geometric obstruction is not represented by separate
singleton covers. No new dimension bound follows from either diagnostic.
