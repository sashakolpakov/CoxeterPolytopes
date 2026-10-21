# Hostile audit: six-carrier integral-packet calibration

**Verdict (2026-08-10): PASS after the direct-input provenance repair.**
All six rows have the claimed instancewise full-trace injection and packet
rank.  The coordinate ceilings and interval-generator counts recompute, and
the only two finite `D` entries use actual integral ample classes.  These are
carrier- and basis-dependent certificates for six frozen examples, not
uniform constants and not a theorem about general reflection-maximal
chambers.

## Pinned artifacts

| artifact | SHA-256 |
|---|---|
| `results/six_carrier_integral_packet_rank_calibration.md` | `9b888ad0a9fd92fb2a1cd792f0299177bb5ca82bc2b9b212fe81f3030b076fc4` |
| `data/enumerations/six_carrier_integral_packet_rank_calibration.json` | `e36f558b0d40e817c718ccb4760687240dab89a3e2f8a12980f5cd28f816eabd` |
| `tests/test_six_carrier_integral_packet_rank_calibration.py` | `4c4762f427cd9756137f237d15489b6e5c76b528a741b602963338af176f44d7` |

The ledger now also pins the three direct Bugaenko inputs which the test uses:

| direct input | SHA-256 |
|---|---|
| `src/coxeter_polytopes/known_examples/bugaenko.py` | `bd08aa582cce4ea28cf9212ae3ead602d1f81f87ab9003eaca3316865ef03f5c` |
| `scripts/bugaenko_k3_eigen_nef.py` | `e429b97676d63ad00536110569ebc220ac8eaeead505cb385785f3bb52bf3554` |
| `scripts/bugaenko_trace_lattice.py` | `1d7900fd72f41a7c1ff346ba17aa8c6170d24db4b7b3504acb79274c4b4c85b3` |

The predecessor ledger omitted these three direct pins while claiming pinned
inputs.  The hashes above identify the repaired provenance state.

## Row audit

| row | hostile check |
|---|---|
| Bugaenko `H8/11` | The `20 x 18` divisor embedding has rank `18` and satisfies `E^t I_(1,19) E=4G_N`; the common positive factor is recorded and does not change the chamber.  The rational Chow projector has rank `18`, so `B=18=d(n+1)` is a genuine global subpacket, not the ambient Picard rank `20`.  In its integral trace basis the embedded pairing has height `8`, multiplication by `2phi` has height `4`, and the eleven facet vectors have maximum height `48`, giving `C=48`.  The action satisfies `a^2-2a-4=0`.  Its two real eigenvalues are `1+sqrt(5)` and `1-sqrt(5)`, in `[2,14]` and `[-2,1]`, so `t=1`.  The `2phi` action is integral on the rank-18 packet and comes from the rational divisor-Tate action; the separate 40-term integral ambient correspondence acts by `4phi`.  Neither is called monodromy.  `D` is correctly absent because the weak transported divisors are not certified effective irreducible curves. |
| toric `H6/9` | The nine walls span `NS` with Smith factors all one, so `B=7=d(n+1)`; the odd-unimodular pairing and wall coordinates give `C=1`, and `t=0` for `Q`.  The archived chamber is the complete nef cone and `A=7H-2(E_1+E_2+E_3)-(E_4+E_5+E_6)` is in its strict interior, hence is ample.  Its facet degrees are six ones and three twos, so `D=12`. |
| Coble II `H10/30` | The walls span the full rank-11 lattice integrally, giving `B=11`; the pinned basis gives `C=17`, and `t=0`.  No ample total-degree packet is archived, so `D` is correctly null. |
| Coble I `H11/26` | The physical walls have rational rank `12`, and the pinned arithmetic walls span their full lattice; this supports exactly the stated identity with `NS_Q` and `B=12`.  The full-NS basis gives `C=21`, and `t=0`.  Raw plane degree is not promoted to an ample degree after blowup, so `D` is correctly null. |
| Coble VI `H14/60` | The physical wall matrix has rank `15` in `NS_Q=I_(1,14)`, giving `B=15`; the pinned basis gives `C=8`, and `t=0`.  No total ample degree is claimed. |
| Coble MII `H17/48` | The forty-eight physical curves span `Pic` with index one, so `B=18`; pairing and curve coordinates have `C=1`, and `t=0`.  Negative-curve exhaustion identifies their chamber with the full nef cone.  The strict interior class `H=5F_1+5F_2-sum E_ij` is therefore ample.  It has degree one on the eight boundary and sixteen exceptional curves and degree six on the twenty-four graph curves, giving `D=24+144=168`. |

For the five rational-surface rows, `rho=b2` has its literal blowup meaning.
The Bugaenko `rho=20` is instead the chosen weak ambient carrier and is not
claimed minimal or intrinsic to the abstract chamber.

## Scope and noncircularity

None of the displayed values is a uniform consequence of reflection
maximality.  `B` depends on the already supplied global carrier, `C` depends
on its pinned integral basis, and `D` depends on both the carrier and a chosen
ample polarization.  The Bugaenko equality `B=18` is the minimal possible
rank once that full degree-two trace module is required, but it does not
bound ranks for other chambers.  Likewise `12` and `168` are valid upper
certificates for their two fixed rows, not optimal or class-uniform degree
bounds.

No monodromy rank is manufactured.  The Bugaenko operator is a Chow
correspondence action; the toric group of order `216` is a marking-cover deck
group; the other four rows supply no monodromy packet.  The report therefore
correctly leaves all four general gates `B,D,C,t` open and derives no global
dimension bound from the sample maxima.

## Replay

```text
python3 -m pytest -q tests/test_six_carrier_integral_packet_rank_calibration.py
6 passed in 0.47s
```

The focused test rechecks all declared source hashes, six rank identities,
wall ranks/spans, both total-degree vectors, all six coordinate ceilings, the
golden interval generator, and the absence of monodromy ranks.  The test file
compiles, and a control-byte scan of the report, ledger, and test was empty.
No candidate artifact or pinned source was modified by this audit.
