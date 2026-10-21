# Katsura--Kondo--Shimada characteristic-5 computation data

These are verbatim computation-data files for:

T. Katsura, S. Kondo, and I. Shimada, *On the supersingular K3 surface in
characteristic 5 with Artin invariant 1*, Michigan Math. J. 63 (2014),
803--844, DOI `10.1307/mmj/1417799227`, arXiv `1312.0687`.

Official sources:

- paper: <https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/KKSssK3char5/KKSssK3char5.pdf>
- data landing page: <https://home.hiroshima-u.ac.jp/ichiro-shimada/ComputationData.html>
- `compdataB.txt`: <https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/KKSssK3char5/compdatassK3char5/compdataB.txt>
- `compdataChams.txt`: <https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/KKSssK3char5/compdatassK3char5/compdataChams.txt>

Retrieved from the official site on 2026-08-06.  SHA-256 checksums:

```text
956391d546707c958947d55c9bfc43ae1ae7952f73a01d4a4d1a0505169e5a1e  compdataB.txt
5173c1c2974c9db1a78c3bf19f6cdc2ad4f08c10be3df9e99cfdd85c6dce7805  compdataChams.txt
```

`GramSX` is the rank-22 Neron--Severi Gram matrix.  In
`compdataChams.txt`, `walls[i]` is the complete published wall list for the
highlighted induced chamber `D_i`; `SAutXD[i]` and `TAutXD[i]` factor the
published period-preserving group `Aut_X(D_i)`.

Scope is important.  The paper distinguishes the full finite lattice chamber
stabilizer

```text
Aut(D) = {g in O^+(S_X) : Dg = D}
```

from the geometric/period-preserving subgroup

```text
Aut_X(D) = Aut(D) intersect G_X.
```

The archived matrices are named for and documented as `Aut_X(D_i)`; the full
`Aut(D_i)` matrices are not supplied.  The paper also reports 12 and 328
further chamber classes at levels 2 and 3, and more than 6000 found at level
4, but it does not give complete wall matrices for those chambers.  Therefore
the repository's certificate is deliberately limited to the three published
groups and is not a census of the whole surface.

