# Exact deletion costs before a larger B-region repair

The control remains the [frozen64-wall chamber](fresh_nonarithmetic_replace47_close12_20261017.md), with198 finite vertices and15 open ends in its complete S6 boundary. No new chamber is promoted here.

The [source](../scripts/fresh_nonarithmetic_b_release_costs_20261017.py) reconstructs every compact face of dimension at least2 from the complete finite-vertex incidence list. It examines every subset of a finite vertex's nine active walls. Such a face is compact precisely when every incident edge at every one of its finite vertices has a second finite endpoint. For sufficiency, take a bounded projective truncation containing all those vertices. They form a nonempty closed-degree component of its graph; connectivity excludes any additional truncation vertices or unbounded direction. This uses simplicity at finite vertices and does not require an ambient exterior-ray census.

The resulting compact-face counts in dimensions2 through8 are558,531,364,174,56,11,1. The [exact result](fresh_nonarithmetic_b_release_costs_20261017.json) saves their maximal supports and the affected finite supports for each proposed deletion.

| Wall released | Finite vertices lost | Compact faces touched by dimension | Named S6 faces touched |
|---|---:|---|---|
|7|131|520 two-faces,519 three-faces,362 four-faces,174 five-faces,56 six-faces,11 seven-faces,1 eight-face|8,10,11,12,55|
|30|13|10 two-faces and1 three-face|10,44,55,63,64|
|52|5|4 two-faces|8,12,64,65|

In particular, releasing7 sacrifices the compact H8 face. Releasing30 sacrifices the compact tetrahedron with support `[0,1,2,3,51,55]`. Releasing52 loses the8-heptagon,12 and64 quadrilaterals, and65 pentagon, but no higher-dimensional compact face. A future replacement may rebuild them; these counts describe the old certificates that cease to apply, not a theorem against replacement geometry.

Generation and fresh saved replay passed. This is an exact combinatorial audit over the previously certified198 finite vertices, with its own source and input hashes. It is a primary computation; no new independent arithmetic review is claimed.
