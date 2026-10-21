# Complete necessary contact profiles for a wall crossing the S75 plane

The [exact checker](../scripts/fresh_h20_s75_projection_gate_20261014.py) and [certificate](fresh_h20_s75_projection_gate_20261014.json) pass generation and saved replay. They reconstruct the rational inverse Gram and classify all possible support contacts of any acute Coxeter wall whose projection onto the retained S75 Lorentz-three complement is spacelike. They do not test the remaining chamber walls or assert realizability of the surviving profiles.

Let the support roots have square two, Gram `G`, and entrywise nonnegative inverse `J`. For a proposed normal `r` of square `q>0`, set

`c_i = -(r·s_i)/sqrt(2q) >= 0`.

Its projected square is `q(1-2c^T J c)`. A wall inducing a geodesic on the plane therefore requires the strict inequality `c^T J c < 1/2`. Since every nonright Coxeter contact has `c_i>=1/2`, including parallel and dotted cases, each such contact costs at least `J_ii/4` in this normalized inequality. Equivalently, it costs at least `q J_ii/2` in the physical projected norm.

Thus every `J_ii>=2` forces a right angle, and any other nonright finite contact must satisfy

`cos²(pi/m) < 1/(2J_ii)`.

Equality gives a null projection and does not define a geodesic in the hyperbolic plane. Nonnegative cross terms can strengthen these bounds; omitting them is only a necessary relaxation.

For the exact S75 support, the domains are:

| Walls | Inverse diagonal | Allowed individual finite labels |
| --- | --- | --- |
| 1,10 | 10/11 | 2,3,4 |
| 2,9 | 18/11 | 2,3 |
| 12 | 1 | 2,3 |
| 31,87 | 7/4 | 2,3 |
| 17 | 1/2 | Every finite label |
| 3,4,5,6,7,8,13,14,15,75 | At least2 | Right only |

In particular, **wall75 has inverse diagonal5 and must be right**. No support contact can be parallel or dotted because every diagonal is at least1/2.

The sum of the three smallest diagonal costs at label3 is `51/88>1/2`. Hence at most two support contacts can be nonright. If wall17 has another nonright neighbor and its own label is at least5, even the cheapest other contact gives cost greater than `95/176>1/2`, using `cos²(pi/5)>5/8`. Therefore only labels3 and4 need testing for paired wall17; this is an all-label argument, not a cutoff.

The complete necessary list consists of the all-right profile, nine single-contact profiles read from the table excluding wall17, the isolated wall17 tail with arbitrary finite label `m>=3`, and these seven two-contact profiles:

- `(1,12)` and `(10,12)`, both contacts label3.
- `(17,1)` and `(17,10)`, with label3 or4 at17 and label3 at the other wall.
- `(17,12)`, both label3.

The exact two-contact census checks 52 candidates in `Q(sqrt2)`; seven satisfy the strict projection condition and45 fail. Thus there are17 discrete profiles including the all-right profile, plus one arbitrary finite-label tail. Other old-wall angles, common-interior conditions and pair compatibility between new caps remain separate obligations.

```
python3 scripts/fresh_h20_s75_projection_gate_20261014.py --verify
```

The input is the frozen [independent whole-face certificate](fresh_h20_bridge_faces_20261014_review.md). No frozen predecessor files were changed.
