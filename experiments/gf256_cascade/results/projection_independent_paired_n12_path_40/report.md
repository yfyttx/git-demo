# GF(256) cascade: independent versus coupled outer projections

## Paired setup

- Fixed 12-symbol `path` inner generator and fixed TIT outer matrices (hashes in `config.json`).
- Depolarizing channel `p_D=0.07125`, hence marginal `f_m=0.0475`.
- Forty identical channel frames for all modes: seeds `20272000`–`20272007`, five frames per seed.
- The X decoder and its conditional Z channel prior are the same in all modes. The independent mode runs each of the eleven outer soft-BP decoders once for 20 iterations, with no inter-projection Turbo exchange. It then uses the same physical-symbol reconstruction as the coupled decoder. It has no robust or adaptive fallback.
- The coupled `conditioned_robust` and `conditioned_robust_adaptive` modes use their existing fallback schedules. Thus this comparison tests the practical independent baseline against the existing decoders; it does not isolate message coupling from all differences in compute budget.

| Z decoder | Logical failures / 40 | Z-syndrome failures | Undetected logical failures | CPU seconds / frame |
| --- | ---: | ---: | ---: | ---: |
| Coupled robust | 9/40 (22.5%) | 8 | 1 | 24.92 |
| Coupled robust adaptive | 4/40 (10.0%) | 3 | 1 | 30.67 |
| Independent outer projections | 40/40 (100%) | 40 | 0 | 4.29 |

The independent decoder failed on all 31 frames where coupled robust succeeded and all 36 frames where coupled robust adaptive succeeded. No independent-only success occurred. The exact two-sided paired McNemar p-value versus adaptive is `2.91e-11`. This is a 40-frame pilot at one matrix and one channel point, not a general FER curve.

## One-frame projection diagnosis

For seed `20272000`, trial `0`, the independent decoder had `bp_valid=0/11`: all eleven decoded outer projections missed their syndromes and the complete `H_X` estimate left 1,297 checks unsatisfied. With 100 independent outer BP iterations, the same frame still had `bp_valid=0/11` and 1,293 unsatisfied checks. The coupled adaptive decoder correctly recovered all eleven projections on this frame after its retry and passed the full syndrome. The sampled physical X error was decoded exactly in all three runs. These diagnostics use the true Z error only **after** candidate selection to count differences; it does not enter decoding.

## Interpretation

This simple independent-projection decoder is substantially worse at this operating point. The evidence supports keeping information exchange through shared inner symbols. It does not prove that every independent decoder must fail: stronger per-projection postprocessing, a different code, or a lower `p_D` could alter the result. Candidate selection must still be evaluated by logical failure, not syndrome satisfaction alone.

Raw counts and pairwise transitions are in `summary.csv`, `pairwise.csv`, and `failure_transitions.csv`; per-batch logs and executable/matrix hashes are retained in this directory.

To repeat with a newly built simulator, choose a **new** run directory because
the saved configuration pins the exact executable hash:

```bash
python3 experiments/gf256_cascade/run_paired_decoder_study.py \
  --p-depolar 0.07125 --frames-per-point 40 --batch-size 5 \
  --seed-start 20272000 --workers 4 \
  --modes conditioned_robust conditioned_robust_adaptive conditioned_projection_independent \
  --matrix-dir experiments/gf256_cascade/results/g2_basis_path_n12_300frames/path \
  --inner-length 12 --allow-custom-matrices \
  --run-dir experiments/gf256_cascade/results/projection_independent_repeat
```
