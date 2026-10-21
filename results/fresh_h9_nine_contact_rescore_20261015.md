# Three saved banks do not cut the latest escape

The [standalone rescorer](../scripts/fresh_h9_nine_contact_rescore_20261015.py) and [certificate](fresh_h9_nine_contact_rescore_20261015.json) test three frozen numerical banks against the exact 30-wall chamber and its latest strict escape.

| Bank | Entries | Compatible entries | Escape cutters |
|---|---:|---:|---:|
| Earlier 28-wall complete-line proposals | 68,111 | 66,912 | 0 |
| 29-wall complete-line proposals | 68,526 | 68,004 | 0 |
| Nine new 29-wall finite-corner germs | 2,173 | 1,910 | 0 |

The reconstruction uses a 100-digit physical frame followed by NumPy pair gates. It inherits old numerical gates and checks new finite labels only through 256, with the explicit tolerances in the certificate. Candidate overlap between banks is not removed, so the rows must not be interpreted as a count of distinct walls across all banks. Generation and saved replay pass.

This is a bounded numerical diagnostic. It does not certify an all-label exclusion, include every continuous parameter, or search the new frontiers created by wall 30.

```sh
python scripts/fresh_h9_nine_contact_rescore_20261015.py --verify
```
