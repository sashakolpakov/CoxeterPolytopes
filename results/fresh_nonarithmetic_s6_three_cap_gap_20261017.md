# Exact angle gap in the three-cap near-hit

The only near-hit from the [bounded three-cap search](fresh_nonarithmetic_s6_three_cap_bridge_scout_20261017.md)
is not a Coxeter cluster. The [exact checker](../scripts/fresh_nonarithmetic_s6_three_cap_gap_20261017.py)
and [certificate](fresh_nonarithmetic_s6_three_cap_gap_20261017.json) prove

    -2cos(pi/28) < A.B < -2cos(pi/27).

Here C is right to S6 and8, has labels15 to44 and4 to52, and uses the branch
directed negatively along the projection of old30. The caps A,B are the
negative normalized duals of52 and44 in the four-root frame[30,44,52,C].
They are respectively right to[30,44,C] and[30,52,C], as well as S6.

All four old side roots lie in S6-perp. The calculation therefore uses literal
three- and four-root Grams, with exact norm2, right0 and prescribed cosine
entries. Only the two old dotted products8/30 and44/52 are imported as outward
rational intervals from the frozen root recipes. Projection gives C.30, and
the inverse four-root Gram gives A.B directly. This avoids the interval
cancellation in an expanded physical-coordinate calculation.

Generation and a fresh saved replay passed. This excludes the displayed
configuration only; it does not exclude other labels, more walls, or changed
old walls. No finite-volume example is claimed.
