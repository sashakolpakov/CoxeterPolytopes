# Every finite endpoint on the original 62 germs misses the escape

The [exact interval certificate](fresh_h9_h8_disjoint_walk_20261014_universal_gate.json) proves that every compatible finite endpoint cap on any of the changed fifteen-wall seed's 62 open finite-corner germs misses its literal strict escape. It covers all 1,076 local spherical profiles, all 64 unrestricted dihedral tails, all real endpoint parameters and all integer Coxeter labels on the remaining walls. This strengthens the bounded first-stage scout; it does not cover the later 71 germs created by the quartet.

For each support the checker reconstructs the normalized Gram and its inverse exactly. The saved finite point `X` and outgoing tangent `Y` give a future orthonormal Lorentz frame `U,V`. With support projection `p` and `ρ=2−p²`, every outgoing finite endpoint cap normalized to norm two has

`r=p+hU+√(h²+ρ)V`, with `h>0`.

The proof includes the larger closed domain `h≥0`. An isolated dihedral tail is covered by the continuous interval `√3≤k≤2`, with `p=k p₀` and `ρ=2−k²/2`. This includes every integer tail label, without a cutoff.

Every wall or escape pairing is `b+a h+c√(h²+ρ)`. Exact outward intervals establish `a≤0` and `a+c<0`, so these expressions decrease with `h`. The escape also has `a<0<c`; its linear upper bound gives a finite certified cutoff beyond which every cap misses it.

On the remaining closed parameter boxes, a terminal box either pairs strictly negatively with the escape or violates a necessary old-wall angle condition. The angle test keeps labels 2 through 15 individually and retains the whole continuous interval containing labels 16 and above, parallel and dotted contacts. It therefore rejects no valid high-label relation. All arithmetic is outward-rounded at 128 bits, with rational Machin/Taylor enclosures for the angle constants. The complete proof uses 4,964 boxes for 1,140 finite/tail profile families.

The [checker](../scripts/fresh_h9_h8_disjoint_walk_20261014_universal_gate.py) rebuilds its coefficients from the exact [frontier](fresh_h9_h8_disjoint_walk_20261014.json). It reuses the frozen interval engine from the earlier mixed-D4 search; this is an exact extension of that method, not an independent reimplementation. Source hashes bind both the inputs and the engine.

This obstruction does not prohibit adding these caps as intermediate steps. Ideal endpoints, additional walls without these finite endpoints, and later supports remain separate questions. In particular, it does not establish nonexistence of a finite-volume completion.

```sh
python scripts/fresh_h9_h8_disjoint_walk_20261014_universal_gate.py --verify
```
