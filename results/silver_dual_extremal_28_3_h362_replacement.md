# Silver DualExtremal(28,3): first cutter at H362 and replacement

The complete eight-class chronology through `H=362 + 256*sqrt(2)` has
103 controls on 70 heights. After the
audited first packet, the only nonempty necessary levels before the endpoint are
`124+88*sqrt(2)`, `181+128*sqrt(2)`, and `186+132*sqrt(2)`.
The continuation also repairs the former fixed `[-3,3]` binary minimum gate:
exact determinant coefficient bounds remove six spurious types, including the
otherwise nonempty-looking level `248+176*sqrt(2)`.

At the endpoint, the q=2 type `(362,182,256)` has exact augmented maximum
`d.(u+sqrt(2)v)=53+37*sqrt(2)` with determinant zero. The explicit vectors in
the JSON satisfy `u=2d+v`, so this is a realized primitive reflective root. Its
pairing with `(1,d/2)` is `1/2 + 1/2*sqrt(2)>0`.
Exact four-vector PSD cross-Gram bounds, with the two sole positive relaxed
extensions rejected by explicit norm-6 and norm-2 combinations, make this root
acute to every earlier necessary type. The elementary equal-height norm-two
difference lemma handles its simultaneous packet.

The cutter does not close the cone. The replacement `z'=(1,(12/25)d)` has exact norm
`1777/625 - 2*sqrt(2)>0`; its tight necessary-type margin is
`14/25 + 6/25*sqrt(2)>0`. Thus the accepted prefix through H362 still
has an open hyperideal patch and infinite volume.

This is a bounded result, not Vinberg termination or global nonreflectivity.
