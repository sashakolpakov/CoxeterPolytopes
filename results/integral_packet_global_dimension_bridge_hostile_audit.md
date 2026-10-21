# Hostile audit: integral-packet global-dimension bridge

**Verdict (2026-08-10): PASS, with the stated conditional scope.**  The
faithful-rank inequality and bounded-packet finiteness theorem are valid under
their explicit global hypotheses.  The star and blowup families correctly
show why bounded local four-node data or bounded individual curve degrees do
not supply those hypotheses.  No universal faithful compiler for all
reflection-maximal arithmetic finite-volume chambers is proved.

## Pinned artifacts

| artifact | SHA-256 |
|---|---|
| `results/integral_packet_global_dimension_bridge.md` | `1323b7243579d24f05a10c0e23d916808c17594fab10605c2470280e7b780920` |
| `scripts/certify_integral_packet_dimension_bridge.py` | `80ed3814c8f00d7d75f0d06b4fb2cdd1049d540ac708fd504672a5a5f6cc6852` |
| `data/enumerations/integral_packet_dimension_bridge.json` | `61f5721f00e348b204af49ec4ac88e1e001909a44f01854022209184230451e9` |
| `tests/test_integral_packet_dimension_bridge.py` | `49d51380c9a99bee8e1c0ca9aeaf789624da83a894c2102730898adf52e9d245` |

## Rank, degree, and finite enumeration

For a certified rank-`n+1` `k`-form `L_k`, restriction of scalars has
rational dimension `d(n+1)`.  The stipulated injection

```text
Res_{k/Q}(L_k) -> H_Q
```

therefore gives `d(n+1)<=rank(H_Z)`, and an independently obtained ceiling
`rank(H_Z)<=B` gives `n<=floor(B/d)-1`.  This is linear algebra, not a
classification or volume argument.  The rank-`n+1` `k`-form and the global
injection are part of the faithful-packet certificate: neither follows from
a collection of local cyclic-product fields.  In particular, a naive
unit-normal span should not be substituted without certifying its `k`-form.

The separate facet bound is also exact.  If distinct irreducible facet curves
satisfy `sum A.D_i<=D` for one integral ample `A`, every summand is a positive
integer.  Hence the facet number is at most `D`, and a full-dimensional
hyperbolic chamber has `n+1` or more facets, so `n<=D-1`.  A per-curve bound
does not imply this total-degree hypothesis.

For fixed `B,D,C,t`, the displayed exponent

```text
R(R+1)/2 + t R^2 + D R
```

correctly counts the entries of a symmetric pairing, `t` endomorphisms, and
an ordered `D`-slot facet list at rank `R`; unused slots can be zero padded.
Summing over `R<=B` and applying the stipulated finite grammar leaves a finite
set.  Exact eigenvalue isolation, Gram/inertia checks, irredundancy, and the
full finite/ideal link test are terminating at bounded rank and facet count.
The returned cyclic products therefore generate only finitely many fields.
This conclusion would fail if a repeated eigenspace, moving scalar,
cross-ratio, or wall suffix could be stored outside that fixed finite
grammar; the theorem explicitly excludes all four.

The interval probe only makes the coefficient part effective:
`[Q(alpha_j):Q]<=23`, hence `[k:Q]<=23^t` for `t` generators (and `d<=23`
when one is primitive).  The report correctly does not infer a dimension
bound from this field bound alone.

## Exact star obstruction and deletion gate

For `G_m=-2I+Adj(K_{1,m})`, a Schur-complement calculation gives

```text
chi_m(lambda)=(lambda+2)^(m-1)*((lambda+2)^2-m).
```

Thus the eigenvalues are `-2` with multiplicity `m-1` and
`-2 +/- sqrt(m)`.  For every `m>=5`, `G_m` is nonsingular of inertia
`(1,m,0)` and rank `m+1`.

I independently enumerated its induced subgraphs on at most four vertices.
Without the centre they are the four edgeless graphs `s K_1`, `1<=s<=4`;
with the centre they are `K_{1,j}`, `1<=j<=3`, with the centre-only case
already equal to `K_1`.  These are exactly the seven archived types, and all
are negative definite.  The predecessor artifact omitted the disconnected
edgeless types; the four hashes above pin the repaired exhaustive ledger.

Deleting the centre gives `-2I_m`; deleting a leaf gives `G_{m-1}`.  Hence
`m=5` has an elliptic centre deletion and parabolic leaf deletions, whereas
every `m>=6` has an indefinite leaf deletion.  The report uses this only as
a necessary simplex finite/ideal-vertex gate.  It does not claim that the
unbounded tail consists of finite-volume chambers, and therefore establishes
exactly the advertised local-Gram-to-global-rank no-go without a false
polytope promotion.

## Surface analogy and the open descent lemma

For distinct points `p_i`, the exceptional curves on
`Bl_{p_1,...,p_r}(P^2)` have matrix `-I_r` and the Picard rank is `r+1`.
Writing `pi` for the blowdown, `O(-sum E_i)` is `pi`-ample; after tensoring by
`pi^*O(a)` for sufficiently large `a`, the class
`A=aH-sum E_i` is ample and satisfies `A.E_i=1`.  This validates the claimed
bounded individual degree.  The polarization and its global Hilbert data
vary with `r`, exactly as the report warns, so this is not a counterexample
to bounded-Hilbert or bounded-total-degree finiteness.

The fixed-Hilbert faithful descent lemma is noncircular as a **conditional
completion target** only if its constants and finite grammar are derived
uniformly, independently of `n`, any chamber list, and the desired finiteness
conclusion.  The report states that this construction is open.  Its six
clauses form a sufficient package for effective chamber enumeration.  They
are not all logically minimal for a dimension bound alone: a faithful global
injection plus an independent uniform rank ceiling `B` already suffices, and
the total-degree ceiling `D` is an alternative direct route.  Coordinate,
generator, and complete-return clauses are what upgrade that bare bound to
an effective enumeration of fields and chambers.

Reflection-maximality is not used in Theorem 2 and is not recovered by the
star or blowup examples.  Any universal compiler must separately prove
exact chamber return, cusp coverage, and preservation of maximality.  The
frozen report keeps this boundary explicit.

## Replay

```text
python3 scripts/certify_integral_packet_dimension_bridge.py --check
PASS: faithful-rank inequality frozen; star local-to-global obstruction exact

python3 -m pytest -q tests/test_integral_packet_dimension_bridge.py
6 passed in 1.09s
```

Both Python files compile.  A control-byte scan of all four frozen artifacts
was empty.  `/usr/bin/time -l` reported `0.80 real`, then returned wrapper
status 1 solely because the sandbox denied `sysctl kern.clockrate`; a direct
replay returned status 0.  No candidate producer, archive, theorem, or test
was modified by this audit.
