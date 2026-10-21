# Exact failed faces and maximizing marked forests

The [source](../scripts/fresh_sqrt5_h3_center_selected_diagnostic_20261016.py)
and [data](fresh_sqrt5_h3_center_selected_diagnostic_20261016.json) reconstruct
306 maximizing witnesses of the
[selected exact budget-35 rule](fresh_sqrt5_marked_h3_center_sparse_recovery_20261016.md).
Each witness includes its rank-21 forest, marked component and root, marked
tree code, full marked-component edges, and literal payment profile. Every
profile is recounted directly from root pairs and pays exactly 35. These are
one chosen witness per maximizing DP context/mode front, not an enumeration of
all tied maximizers.

Representative forests include E7+H4+2D5, E8+2H4+D5, and B4+3H4+D5. Of the 306
witnesses, 103 mark an A1 component, 81 mark H3, and 51 mark I2(5). The data also
record all remaining marked-component types and root positions.

Two exact newly affected face failures are preserved:

- Triangle row 136 has seven roots and corners H4+A2, D6, A3+I2(5)+A1. Common
  external counts (A1,I2(5),H3,H4,A4,D4,D5) are (2,1,1,0,0,0,0); facet root is
  white root 3. Its payment is exactly 2 minus the positive adjacent discount.
- Quadrilateral class 4202 uses common counts (0,2,0,0,0,0,0) and white facet
  root 4. Its payment is exactly 1 minus the same discount.

The data contain literal Gram-label rows, sparse profiles, full rational
payments and shortfalls, and the names and values of every involved variable.
All compensating bonus coefficients in these two profiles are exactly zero.
The approximate payments are 1.364071469188652 and 0.364071469188652.

Generation and fresh saved replay pass. This is an exact payment and witness
extraction diagnostic. Arithmetic geometric admissibility is checked in the
separate face producers or subsequent audits, not established by this script.
These failures concern the attempted compact arithmetic bound of 20; they do
not affect the frozen compact bound of 21.
