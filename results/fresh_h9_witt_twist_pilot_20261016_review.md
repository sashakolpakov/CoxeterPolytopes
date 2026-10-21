# Independent complete height20 pilot on the local-Witt host

**PASS, generation and fresh saved replay.** The
[independent checker](../scripts/fresh_h9_witt_twist_pilot_20261016_review.py)
and [certificate](fresh_h9_witt_twist_pilot_20261016_review.json) reconstruct
all13 controls and16 walls of the height20 Vinberg prefix, using all four
primitive norm classes and the complete B7+B2 zero sector from the
[independent setup review](fresh_h9_witt_twist_20261016_review.md).

This reviewer uses project Golden arithmetic, with no producer or producer
helper imports. Its coefficient boxes are deliberately larger. It rebuilds
the control heights from the norm equation at both embeddings, and checks
the resulting coordinate alphabets against the saved data only after their
independent construction.

The search algorithm differs substantially from the producer. It enumerates
ordered spatial blocks of sizes3,4,2, hashes their exact energies, and joins
them into complete norm solutions before applying chronological wall tests.
It does not reuse the producer's reverse-coordinate recursion, remaining
coordinate lower bound, or incremental wall pruning. There are282 three-
coordinate candidates,632 four-coordinate candidates and58 two-coordinate
candidates:972 total block candidates. The joins give32 complete ordered
norm vectors across all controls; wall filtering recovers exactly the seven
positive-height walls. Every same-height pair is also checked.

The global boxes follow from q<5, q′<10, φ>8/5, −φ′>3/5, β>2 and β′>4.
At height20 they give x<8, |x′|<5, y<11, |y′|<4, z<8 and |z′|<2.
The independent coefficient conversion uses strict rational bounds for φ
and √5 and adds two units of padding. Physical dominance is enforced within
the seven- and two-coordinate blocks; no positivity restriction is imposed
at the conjugate embedding.

For every reconstructed vector the checker verifies the exact norm and
reflection integrality. It then checks the complete16-by16 Gram matrix,
all120 pair types, a strict future interior point, and all16 actual facet
witnesses. These tests establish a complete bounded prefix with genuine
facets. They do not on their own establish finite volume or nonreflectivity.

The producer's complete cone census and corner-map exploration are separate
claims. The two-map nonreflectivity proof, when independently checked, needs
only its literal finite corner witnesses and this complete bounded prefix;
it need not assume that every finite corner was found.

```sh
python3 scripts/fresh_h9_witt_twist_pilot_20261016_review.py
python3 scripts/fresh_h9_witt_twist_pilot_20261016_review.py --verify
```
