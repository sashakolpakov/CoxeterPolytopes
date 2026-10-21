# Silver rank-24 modular sequence: the complete `n=4` prefix

The literal Nebe--Sloane endpoint Grams define `G_n=G_0+n(G_1-G_0)`.
Exact rational matrix arithmetic gives
`D G_0^(-1) D = D-G_0/5`, hence
`det(G_n)=(n^2+5n+5)^12`.  For `n=4`, `G_4` has Smith invariants
`1^12,41^12`, minimum 24, and 37,800 signed minimal vectors.

For the admissible matched coefficient `a=-3+5*sqrt(2)`, the complete local
alphabet is the dyadic choice times the two independent split-41 choices:
`{1,2}*{1,7+2*sqrt(2)}*{1,7-2*sqrt(2)}`.  The mod-41 pairing-radical
preimage has index `41^12`; after dividing its Gram by 41 it again has
determinant `41^12`, minimum 24, and 37,800 signed minimal vectors.  This
forces divisor-square normalization on every split-prime gated target.

The complete eight-class chronology through `H=5000`
contains 196 controls on 132
heights.  Its first nonempty necessary packet is the ungated norm-two level
`H=806 + 570*sqrt(2)`, with the two binary Gram
types `(798,408,570)` and `(794,410,570)`.

For `h=G_4^(-1)e_1`, `h^2=24/41`.  Exact augmented PSD enumeration plus the
complete coordinate bounds on the norm-24, 26, and 28 shells gives directional
maxima `17+12*sqrt(2)` and `14+10*sqrt(2)`.  Therefore
`z=(1,3h)` has Lorentz norm `339/41 - 5*sqrt(2)>0`
and margins `7 + 4*sqrt(2)` and
`16 + 10*sqrt(2)`, both positive.  The whole bounded
prefix retains an open hyperideal patch and has infinite volume.

This is a bounded negative design result, not Vinberg termination or a proof
of global nonreflectivity.  The broader sequence ranking used only the `q=1,2`
necessary relaxation and is intentionally not promoted to a theorem here.

Sources: [P24.9](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/P24.9.html) and [L_24.11](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/L_24.11.html).
