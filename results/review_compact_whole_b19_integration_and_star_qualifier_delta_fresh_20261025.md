# Whole B19: corrected append and final affine-star qualifier

25 October2026. SOURCE ONLY. No import/test/proof/build job or manuscript
edit. This supplements the frozen full review `15ded4f6` without
rewriting that historical report.

I read the full integration receipt and the complete extracted-append
diff against original 7e lines904--3001. The only mathematical
differences in that append are the previously reviewed c2bb Bwindow,
25aa Bhead and 78fca common-list sections; there is one final blank
line. Every other B19 proof body is unchanged. The patch has only
this append and does not remove or rewrite the existing A19 or
higher-tail proofs.

Integration patch
`results/compact17_whole_b19_corrected_working_integration_20261025.patch`,
SHA256 `a504925658613e28c775c20dcfea0dd82fde4011f3478e4bed3dae19270f696d`.
Receipt same basename `.json`, SHA256
`32e0856aac8a63dc54230778ee7cb498bff27119265d685ab44780e01736dd6f`.
It binds main483453fb to main3f71ad54 and the three exact section
pins of the full review.

Root then correctly observed one further explicit premise needed in
the F3 star sentence: leaves c,n16,n18 are mutually orthogonal only
when `c.n16=0`. I confirmed this directly. A full bridge contacting
the near-left16 root is already directly strong by the bridge
ratio>2, so it does not use the star. The one-sided right-short row
and the far-left bridge have zero contact at16 and have the stated
positive affine star. This refinement avoids giving the star a
broader domain than its actual Gram matrix. The frozen full report's
F3 discussion must be read WITH this final qualifier, not as an
assertion about arbitrary bridge contacts.

I read the entire one-hunk final patch:
`results/compact17_b19_star_domain_qualifier_working_integration_20261025.patch`,
SHA256 `ded5581c6cc1394253ae2a3b78603d3cd5e6593d1504f6d7aba899640d1a8e78`.
It only adds “and zero contact at sixteen” and reflows those three
lines. All prices, numbers, lists and other proof bytes are unchanged.
Receipt same basename `.json`, SHA256
`6ef8968f879bc6ec124c8e70a7890fc17760cfc311100ef201acfbcac33391f3`,
binds the preserved main3f71 preimage to main43198f86, 744887bytes:
`43198f86c20d19c27e62a77fb8e966810bbf9f5cd2393c975da4e12d34e0d87a`.

Verdict: these exact append and final-qualifier deltas are supported
at the full review's inherited geometric/prefix interfaces. No lost
case or changed bound. This does not claim exact20/D19 coverage,
a complete compact manuscript, a compiled PDF or Lean verification.
