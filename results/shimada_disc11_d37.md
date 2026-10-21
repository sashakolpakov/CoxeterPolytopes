# Exact high-root continuation through physical \(D_{37}\)

The two genuinely new neighbors of \(D_{61}\) were ranked exactly by their
restricted root facets. Physical \(D_{84}\) has six: the verifier supplies
six exact face witnesses and positive exact dependences for all six remaining
restricted root classes. Physical \(D_{37}\) has nine. Thus the prescribed
high-root continuation selects \(D_{37}\).

Its complete exact reduction is

\[
\boxed{31\text{ walls}\longrightarrow22\text{ rays}\longrightarrow11\text{ facets}}
\]

with nine root facets and only two internal facets. All eleven faces have
exact timelike witnesses, all eleven remaining rays have exact positive
dependences, and the generic \(\mathbb Q(\phi)\) kernel passes.

The nine root classes are

\[
\begin{gathered}
\{(0,0),(13,0)\},\ \{(1,0),(5,0)\},\ \{(2,0),(19,0)\},\\
\{(3,0),(12,0)\},\ \{(4,0),(16,0)\},\ \{(6,0),(17,0)\},\\
\{(7,0),(18,0)\},\ \{(8,0),(14,0)\},\ \{(9,0),(15,0)\}.
\end{gathered}
\]

Among 20 ambient `-2` walls, 19 are integral field roots, with shell counts
nine of norm \(2\) and ten of norm \(4-2\phi\), spanning field rank 10.
Wall `(11,0)` is the unique non-field root; field root `(10,0)` is also
redundant.

The two exact crossings are

\[
(20,0)\longrightarrow \text{the same physical }D_{61},\qquad
(21,0)\longrightarrow \text{new physical }D_{22}.
\]

Thus the high-root branch has a single new forward endpoint, \(D_{22}\).

Reproduce with:

```sh
PYTHONPATH=src:. python3 scripts/analyze_shimada_disc11_d37.py /path/to/disc11.txt
```
