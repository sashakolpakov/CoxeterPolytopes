# Separate endpoint caps with a finite mutual corner

The [numerical scout](../scripts/fresh_nonarithmetic_s6_split_bridge_scout_20261017.py)
and [saved data](fresh_nonarithmetic_s6_split_bridge_scout_20261017.json) explore
two walls separately terminating S6+[30,44] and S6+[30,52], with a mutual
finite corner on S6+[30]. Neither wall must terminate both old edges, so the
earlier common-cap obstruction does not apply.

The retained chamber has63 walls. The tested endpoint profiles give50 surviving
branches for30/44 and26 for30/52, each with one continuous parameter interval.
All unprescribed old contacts are restricted to be strictly dotted. Sampling
the first parameter and solving the second from a prescribed mutual Coxeter
angle gives6228 numerical finite-corner hits;100 are saved. The label grids,
parameter samples and floating-point tolerances are explicit in the source.
This is not an exhaustive continuous or all-label classification.

The simple zero-contact pair at first parameter7 and mutual angle pi/4 was
subsequently [reconstructed exactly](fresh_nonarithmetic_s6_split_pair_20261017.md).
That exact package is the authority for the new walls; it does not depend on
numerical equality or this scout's parameter bounds.

Numerical generation and fresh saved replay passed. Neither a scout hit nor
the exact local pair establishes finite volume.
