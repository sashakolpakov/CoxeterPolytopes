# Golden odd-Leech section H23 independent hostile audit

**Verdict: PASS, with one repaired schema label.**  A clean-room exact
replay confirms both `Q(sqrt(5))` first-prefix screens on

`S23=(1^24)^perp` inside the catalogue odd-Leech lattice:

- `<-phi> + S23`, with first packet size 553,104;
- `<-2phi> + S23`, with first packet size 53,682.

The audit imports no project module.  It reconstructs the extended Golay
code from twelve standalone catalogue masks, reconstructs the odd-Leech
Construction-B lattice by HNF, and then forms the primitive zero-sum
section.  The resulting section independently has rank 23, determinant 3,
Smith form `1^22,3`, minimum 3, a 2,576-vector minimal shell spanning rank
23, and the same frozen minimal-shell digest.

The section Gram has rank 22 modulo 3 and therefore a one-dimensional
pairing radical.  This revealed a semantic defect in the production JSON:
the old key `radical_mod_three_rank: 22` was really storing the **Gram**
rank, not the radical dimension.  The generator and prose were already
mathematically correct.  The schema is now repaired to record both

- `gram_rank_mod_three: 22`, and
- `pairing_radical_dimension_mod_three: 1`.

The hostile replay goes beyond the prose exclusion for the radical.  A full
Construction-B dynamic program finds zero nonzero radical vectors of raw
square at most 48, independently proving radical minimum at least 7.  It
therefore validates the 3-gated height-zero exclusions.

The same full dynamic program, over all 4,096 reconstructed Golay words,
enforces coordinate sum zero, exact square, code congruence, and the even
residual-sum condition.  It obtains the complete coset-resolved shell counts

- norm 4: `(53,682, 0)` in the even and odd cosets;
- norm 5: `(0, 553,104)` in the even and odd cosets.

Standalone `Q(phi)` arithmetic confirms the canonical local mirror-norm
alphabet `{1,2,3,6}`, the ungated/gated height-zero type counts
`(0,0,1,1)/(0,0,0,0)`, all seven controls through the `a=phi` packet, and
both controls through the `a=2phi` packet.  Each prefix has exactly one
nonempty type: `(u^2,v^2,u.v)=(0,5,0)` and `(0,4,0)`, respectively.

Every displayed packet root is primitive because its time coordinate is a
unit.  Direct calculation confirms root norm two and reflection integrality
in both the time and full spatial lattice coordinates.  The section minimum
gives distinct-vector dot bounds 3 and 2, so both complete simultaneous
packets have root-pairing upper bound `-2phi<0` and are mutually acute.

Finally, the replay confirms both strict hyperideal witnesses:

- for `a=phi`, `z=(1,(2/3)e_4)` has square `16/9-phi>0` and closest wall
  pairing `-1`;
- for `a=2phi`, `z=(1,(4/3)e_3)` has square `16/3-2phi>0` and closest wall
  pairing `-2+(2/3)phi<0`.

Thus both exact finite prefixes are neither finite-volume nor compact.  This
does not decide the eventual chambers: later mirrors and other timelike
coefficient ideals remain open.

Reproduce with:

```sh
python3 scripts/independent_hostile_audit_golden_odd_leech_section_h23.py
PYTHONPATH=src python3 scripts/audit_golden_odd_leech_section_h23.py --verify
PYTHONPATH=src pytest -q tests/test_golden_odd_leech_section.py tests/test_golden_odd_leech_section_h23_independent_hostile_audit.py
```
