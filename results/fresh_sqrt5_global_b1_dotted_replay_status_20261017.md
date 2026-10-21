# Full saved dotted attachment replay completed

The full saved replay of
[the dotted attachment closure](fresh_sqrt5_global_b1_dotted_fast_20261017.json)
finished successfully on2026-10-17. The existing process handle27242 returned
terminal exit status0. Its final output was:

```
PASS_SAVED_DOTTED_B1_ATTACHMENT_CLOSURE 25535 25535 True
```

This completes the previously pending replay of all25,535 saved states.
It was not restarted from a partial log. The command was:

```
python3 scripts/fresh_sqrt5_global_b1_dotted_fast_20261017.py --verify
```

The preserved [source](../scripts/fresh_sqrt5_global_b1_dotted_fast_20261017.py)
has SHA256 `8afd90e0cbf2d62cb3ee6af1832d13ccb1a38504d180fca5c1de677019e8f675`.
The JSON has SHA256 `4b4168d24a3c998262dd0ff729a5332d97d1eddc570dbc21e0f310f3001944aa`.
The [saved database](../data/enumerations/fresh_sqrt5_global_b1_dotted_fast_20261017.sqlite)
has SHA256 `be3a60b3665bdf114ebcb859810f2abcca160aeaebdd423df4bc427cf8fa1ad0`,
matching the producer manifest. No producer or database was changed.

The completed closure remains conditional on its specified pruning assumptions,
including adjacent discounts at most one half, zero distance-two discount and
nonnegative bonuses. A complete replay does not make this an unrestricted face
catalogue or a dimension theorem. The separate finite attachment replay is
still active; this note makes no claim about its completion. Earlier frozen
checkpoints correctly record the replay as pending at their creation time.
