# A twenty-wall arithmetic H9 prefix over Q(sqrt(2))

The [exact search](../scripts/fresh_silver_e8_h9_20261015.py) on
L=⟨−(1+sqrt(2))⟩⊕E8⊕⟨1⟩ gives a complete Vinberg prefix through height30
with twenty walls. Generation and saved replay pass. The
[certificate](fresh_silver_e8_h9_20261015_h30.json) includes all seventeen
time controls, full root vectors and a strict timelike interior point.

The lattice is unimodular over Z[sqrt(2)]. The primitive reflection norm
classes, modulo unit squares, are1,2+sqrt(2),2; all three are included.
The zero-height chamber is E8+A1. The enumeration uses exact simultaneous
bounds in both embeddings and dominant E8 covectors.

The current prefix has infinite volume. An exact feasible spacelike ray is

    (1, 2s, 4s, 6s, 5s, −1+9s/2, −2+4s, −1+2s, 3s, 0),  s=sqrt(2),

with squared norm10−7s>0 and nonpositive pairing with all twenty walls.
Combining it with the strict interior point gives a strict spacelike escape.
The geometry routine found this witness after three distinct feasible rays;
it did not enumerate every vertex or cusp.

This does not prove nonreflectivity of the full integral reflection group:
later Vinberg walls may cut the escape. Nor is an integral isometry or
inequivalence to earlier hosts asserted. No finite-volume example or new
dimension has been obtained. Independent review of this prefix is pending;
the immediate research priority is the arithmetic dimension-bound work.

```sh
env PYTHONPATH=src python3 scripts/fresh_silver_e8_h9_20261015.py --verify
```
