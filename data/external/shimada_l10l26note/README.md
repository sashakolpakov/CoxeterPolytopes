# Shimada's exact Enriques configuration data

These files come from Ichiro Shimada's official computational-data page:
<https://www.math.sci.hiroshima-u.ac.jp/~shimada/K3andEnriques.html>.

- `L10L26Note.pdf`, *A note on configurations of (-2)-vectors on Enriques
  surfaces*: SHA-256
  `b88acf9d1f4602724117e1270431a3cf0afe31e4214ab820ad2796125be2fc3e`.
- `L10L26NoteCompData.pdf`, data specification: SHA-256
  `a2744b6f5826ee23e8170aa7b90beea9fb836e46078a465ccee0ba2349840361`.
- `L10L26compdata2.txt`, exact GAP records: SHA-256
  `218c8188e467bde29369cdeda3771d767c3f70105d8fb56def8d698b90e75a6c`.

The note gives the dictionary

```text
I=12A, II=12B, III=20B, IV=20F,
V=20A, VI=20E, VII=20C=20D.
```

The raw GAP records contain `GramL10`, the ten-coordinate `walls`, and the
complete wall `configuration` matrix.  The project parser reads only balanced
integer arrays and quoted record names; it does not execute the GAP file.
