# Narrow independent review of the release8 numerical deformation

**PASS: conceptual and source review, numerical scope only.** This note reviews the [source](../scripts/fresh_nonarithmetic_release8_ring_deformation_20261017.py) and [saved numerical output](fresh_nonarithmetic_release8_ring_deformation_20261017.json). It does not independently repeat the numerical root solve, contact checks or boundary enumeration, and does not certify a new exact chamber.

The exact bytes reviewed are:

- source SHA256: `a9f254b10ee275b8bc62981e87c86803658091191e63c66195e5c2e4b45c2cc1`;
- data SHA256: `e6d54626bde5be2579df6585d443e7227d5249dafc36b2b137c87fa892eb3e54`.

The free parameter is the native product of the replacement scaffold12 with old8; its contacts with7 and52 remain prescribed. The subsequent norm-constrained constructions subtract the A2 projection offsets from their native contact equations. Each has at most two real branches, with numerical branch choice by continuity from the saved stage25 configuration. The scalar solve makes S and T right, within floating tolerance. These operations are valid discovery mechanisms; the source correctly keeps all their geometric conclusions numerical.

The positive scaffold12·8 product requires dropping one of those two walls. Both actual selections are evaluated separately. Retaining scaffold12 and deleting8/old12/47/53 gives18 open ends. Retaining8 and deleting12/47/53 gives15, with no scaffold in that chamber. The invalid selection containing both is not promoted as an11-end chamber. All retained old native contacts, all new-new contacts and the common-center sign are checked for each selected system. Deleted finite supports are debited from the exact old198-support census; the source does not present that debit as a complete new ambient census.

The restoration search for the second selection includes native tests against all six constructed walls, not merely the old bank. Its finite alphabet and both the one-cap and split-cap duties are explicit. The53 finite-angle gaps and23 dotted pairs among the bounded2×38 split pools justify failure of that finite mutual-contact bank, not an all-parameter theorem.

No mathematical scope issue was found. The result is a reproducible numerical diagnostic: one continuous parameter closes the selected central contact, but the valid deletion costs leave18 or15 open ends. There is no interval-certified deformation root, exact new native Coxeter system, complete ambient cone, or finite-volume construction in this package. This review is final/frozen for the hashes above.
