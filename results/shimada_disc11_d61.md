# The common physical \(D_{61}\) K3 eigensection

The two routes

\[
D_{30}\xrightarrow{(22,0)}D_{61},\qquad
D_{27}\xrightarrow{(23,0)}D_{61}
\]

produce exactly the same marked chamber: their cumulative row isometries and
their complete transported 32-wall sets agree. This is stronger than equality
of the representative label.

Exact restriction over \(\mathbb Q(\phi)\) gives

\[
\boxed{32\text{ walls}\longrightarrow25\text{ rays}\longrightarrow11\text{ facets}}
\]

with seven root facets and four internal facets. Eleven exact timelike face
witnesses and fourteen exact positive redundancy identities pass the generic
half-space kernel.

The 18 ambient `-2` walls contain 17 integral field roots, with shell counts
eight of norm \(2\) and nine of norm \(4-2\phi\). They span the full field
rank 10. Wall `(16,0)` is the unique non-field root and is redundant.

The internal first star is

| wall | physical target | status |
|---:|---:|---|
| `(18,0)` | \(D_{84}\) | new |
| `(19,0)` | \(D_{30}\) | exact return |
| `(20,0)` | \(D_{27}\) | exact return |
| `(23,0)` | \(D_{37}\) | new |

Both returns are checked by cumulative isometry and full wall-set equality;
the return walls are the negatives of the corresponding entry walls.

Reproduce with:

```sh
PYTHONPATH=src:. python3 scripts/analyze_shimada_disc11_d61.py /path/to/disc11.txt
```
