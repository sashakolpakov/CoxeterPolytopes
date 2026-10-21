# Exact native obstruction to the flipped disk's Z duty

**Status: FINAL/FROZEN; exact generation and fresh saved replay PASS.** This certificate excludes Z on the 20 explicitly reconstructed algebraic A/C/X/Y/R prefixes used by the [old-30 connector scout](fresh_old30_connector_flip_20261017.md). It does not exclude continuous changes to those prefixes or other disk incidences. The 63-wall, 199-finite-vertex, 10-open-end infinite-volume control is unchanged.

The flipped disk requires both finite corners `S6+X+Y+Z` and `S6+Y+R+Z`. Since `S6=H4+A2`, every added root in such a spherical extension is right to H4 and has one of the seven allowed A2 profiles. The inherited finite-diagram classifier enumerates the two A2-plus-three-root corners simultaneously. A finite Coxeter bond with label at least 6 belongs to an isolated I2 component, so replacing it by the symbol 6 preserves this finite-domain classification. This gives exactly 295 low-label recipes and two single-high families; there is no two-high case.

Each prefix is defined by exact norm-two/contact equations in this order:

- A: frame `[7,8,53]`, zero A2 profile, labels `[2,2,3]`.
- C: frame `[7,30,A]`, zero A2 profile, labels `[2,2,AC]` from its saved recipe.
- R: frame `[8,44,63]`, zero A2 profile, labels `[2,R44,2]` from its saved recipe.
- X: frame `[8,A,C]`, its saved profile and three assigned labels.
- Y: frame `[8,R,X]`, its saved profile and three assigned labels.

The JSON records the selected algebraic sign for each step. Comparison with the saved numerical product vector chooses this literal sign and is only a discovery sanity check. The obstruction does not use a numerical vector as an exact root or certify all native contacts of the prefix itself.

All auxiliary old targets are projected off the entire six-root S6, using its literal H4+A2 Gram matrix. This matters because old4,9,10,11 need not be right to H4. For new roots right to H4, the native product is recovered from the residual product by adding the A2 projection term. Thus every rejection is a native norm-two Coxeter gate, not a projected-angle test.

For Z, invert the projected Gram of `[X,Y,R]`, solve its prescribed contact projection, and use an explicitly positive complementary kernel line. The latter certifies that the nondegenerate three-frame is Lorentzian before any remaining-norm sign is used. The norm equation then has the two recorded square-root branches; both are tested. The code also covers a zero remaining norm conservatively and stops if any required branch cannot be excluded.

Each high magnitude `2 cos(π/m)`, `m≥6`, is relaxed to the complete closed interval `[√3,2]`, including the limit. All 297 cases certify without subdivision. The 594 branch certificates comprise 390 strict obtuse products and 204 strict consecutive-angle gaps. Their first witnessing walls are old4 (220), old7 (151), old8 (126), old10 (57), and old30 (40). A consecutive-angle gap excludes every finite label because its two thresholds are consecutive; searching thresholds through label64 is a search for rejection witnesses, not a cutoff on candidate labels.

Both high families occur at zero-based prefix7, profile0, with labels `[H,2,2]` and `[H,2,3]` to `[X,Y,R]`. Over their entire relaxed intervals, the negative kernel branch has positive product with old30 (lower bounds greater than 1.095 and 0.990 respectively); the positive branch has positive product with old8 (lower bounds greater than 1.999 and 2.065). The exact rational enclosures, rather than these rounded displays, are the certificate. Neither a common-center assumption nor a choice of the future branch is needed for these exclusions.

The calculation reuses the frozen outward rational physical-interval engine and finite-diagram classifier. Generation69439 and fresh saved replay91267 both terminated with exit0/PASS. Syntax, all13 direct source bindings and five local report links pass. A separate saved-output scan also checks both signs and the strict rational inequality in every one of the594 rejection certificates; it does not reconstruct their physical interval enclosures. Root's independent full source/proof read passed, explicitly checking the full S6 projection on both old arguments, positive-complement Lorentzian-frame argument, both norm branches, finite-diagram high-label relaxation and restriction to20 literal prefixes. Root did not independently rebuild the physical intervals or numerical banks. Generation and saved replay check deterministic reconstruction and every interval rejection; they are not independent arithmetic implementations. The statement is exact for these20 literal prefixes and this mandatory Z duty. It does not establish completeness of the earlier numerical prefix search, exclude continuously varied prefixes, or yield an ambient volume result.

[Source](../scripts/fresh_old30_connector_native_gate_20261017.py) · [Saved certificates](fresh_old30_connector_native_gate_20261017.json) · [Inherited exact frame engine](fresh_new_cavity_construction_gate_20261017.md) · [Bounded scout and predecessor correction](fresh_old30_connector_flip_20261017.md).
