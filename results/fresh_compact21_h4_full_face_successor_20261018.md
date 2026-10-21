# Exact full native-face primal below 40; full forests still pending

The single authorized full-face solve found a rational rule satisfying every
old native face inequality and the existing forest support at

    z = 49986480339/1250000000 = 39.9891842712,
    H4 transfer t = 3/4.

It uses the separately proved strengthened A5 source capacity one. The former
capacity-two packet remains unchanged and valid. **No dimension improvement is
claimed here:** the exact vector still needs complete finite-forest separation
and final universal-coverage composition.

The streaming reconstruction processed all 1,979 native faces and all 102,917
common spectator cases. It deduplicated only exactly identical rational rows
and right-hand sides, obtaining exactly the original 51,443 face inequalities.
The model then added 510 existing forest rows and t<=1, for 51,954 total rows,
30,750 occurring columns and 187,756 nonzero coefficients. Missing bonus
coordinates are zero. All 93,045 old bonus coordinates retain their original
global indices; removing unused columns in the numerical solve changes no
selected inequality.

The worker retained 13,104 positive rational bonus values, 2,791 exact dual
multipliers and all 2,791 nonzero numerical dual multipliers. Every literal
deduplicated row was verified rationally after conservative primal rounding.
The minimum exact row slack is zero. Exact weak duality was checked separately,
with lower value 39989184274324/1000000000253. The numerical objective was
39.98918426937359; it is a proposal, not proof. The old native face coverage
and the reviewed D5/H4 face-preservation lemmas are theoretical dependencies.

Generation passed in terminal chunk `f236db`. The guard recorded 8.679 seconds,
one thread, 258,949,120 bytes combined peak RSS and 282 bytes stdout. Matrix
construction took 3.73 seconds. The solver limit was 42.27 seconds, below the
authorized 45-second cap. Total wall and RSS caps were 59 seconds and 1 GiB.
Initial free disk was 197.96 GiB and available RAM was 7.22 GiB. The process
exited, released the global lock and removed its temporary output. The final
necessary sparse artifact is 443,626 bytes, below its 8 MiB hard file cap.

```sh
python3 -B scripts/fresh_compact21_h4_full_face_successor_20261018.py --verify
```

This solver-free replay streams the same native inventory, checks its coverage
ledger and source hashes, rebuilds only the small sparse rational rows, and
checks both saved certificates. It imports no numerical package. A fresh
separate-process replay remains a distinct verification step from the worker's
generation check.

Verification status: exact internal full native-face verification and exact
verification of the existing forest support. No complete finite-forest bound,
Lean proof, external review or dimension-20 theorem is asserted in this packet.
