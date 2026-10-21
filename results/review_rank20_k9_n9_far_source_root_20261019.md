# Root source review: one cubic arithmetic batch

19 October 2026. **Full source review PASS for a guarded research run.**
I read both complete sources and the complete mathematical design and
handoff. This review does not report a completed enumeration or absence
result.

* Pure helper `scripts/exact_totally_real_field_20261019.py`, SHA256
  `fbd896a6160f232e0f6b4e8a715c990b05498d74be5c846e91aa699c95d598d2`.
* Fixed producer `scripts/fresh_rank20_k9_n9_endpoint3_far_20261019.py`, SHA256
  `2310f437a8210d34891a3e919965ffa35b5cde1a5102a405cd6a64826706ba5d`.
* Design pin `b1dcfe5e74d1372d2250cdbfc3f5dc256bacc7b704aada73760cdad5e3e2dc26`;
  height theorem pin `9fc0583124cba17afed393bb569f83835e8f2c4be8ce5cfa64839f0369254516`.

The helper uses rational polynomial Euclid, sign-preserving Sturm
normalization, exact root counts/refinement, interval Horner signs,
Newton traces and rational trace duals. The fixed cubic is irreducible
by its rational-root test, and the three intervals certify all real
embeddings in the right order. Its asserted integral basis remains the
explicit external real-cyclotomic ring input.

The producer enumerates every integer trace slice and every integer
tail of each complete outer coefficient box. Recovering c0 from the
trace equation introduces no cutoff. The rational nonphysical enclosure
contains the full `|9A_i+1|<sqrt(10)` interval. I checked the expressions

```
qc = 4/9 + Z^2/18,  C = 18A^2+4A-12,
delta = 40/C,  (Z^2/180-11/36)*delta = 1.
```

All branches use exact field signs and power-basis integrality. The
three-embedding z test enforces the declared actual-ground-field
restriction. Every surviving necessary-condition vector is retained;
complete nonempty and complete empty statuses are distinct. An exit-zero
bounded wrapper run would mean that this traversal finished, and would
not itself imply no hits or global geometric existence.

The producer writes only a new bounded certificate path, hashes source
dependencies before and after work, and records incomplete status on
errors or resource/visit caps. It spawns no process. The root invocation
must still use the shared-lock 59-second/1-GiB wrapper. A separate full-box
replay has not been implemented and remains required for an absence
claim.

## Actual launch status

The first scheduled generation attempt was deferred BEFORE scientific
child launch: disk 184.8 GiB, available memory 5.0 GiB, below the 6-GiB
launch gate. No result certificate or generation ledger was created.
No cap was relaxed. The source is ready; execution remains queued.
