# Z-side projection ablation, GF(256), inner length 12

## Fixed conditions

TIT outer matrices and `H_Z` are held fixed within each basis. The star and path `G_2` representations have the same row space. Matrix seed: `20261018`; depolarizing probability: `p_D=0.07125`, so `f_m=0.0475`. The original 300-frame baseline uses channel seeds `20268300`–`20268359`, five frames per seed, and decoder mode `conditioned_robust`.

Baseline logical failures: star 39/300, including 38 Z-syndrome failures; path 45/300, including 43 Z-syndrome failures. X decoding was exact on all 300 frames for both bases.

## Removing staged Z projection decoding

The experimental mode `conditioned_fullhx_bp` retains the baseline HZ/X decoder and the same X-conditioned physical Z channel prior, but runs GF(256) BP directly on the complete `H_X` Tanner graph. The default is 20 iterations; `GF256_FULLHX_BP_ITERS` can override this. It does not use the structured Z decoder or its fallback.

For channel seed `20268312`, star frame 0 remained a logical failure: the number of unsatisfied `H_X` checks increased from 13 under `conditioned_robust` to 1429 under direct BP. Path frames 0, 1, and 2 all failed under direct BP, with 671, 799, and 754 unsatisfied checks; baseline path succeeded on frames 0 and 1 and failed on frame 2. This is a four-frame diagnostic, not a population FER estimate. It rules out this 20-iteration direct-BP implementation as a drop-in replacement on these frames.

## Repairing syndrome versus logical class

The existing `conditioned_robust_joint_projection` mode was replayed on the first 10 path frames selected **because the baseline Z syndrome failed**. It made all 10 final Z syndromes valid, but rescued only 1 logical failure. The other 9 became undetected logical failures. This conditional 10-frame sample does not estimate FER; it demonstrates that satisfying the syndrome is insufficient.

## Failure-triggered adaptive Turbo schedule

`conditioned_robust_adaptive` uses the baseline result unless its Z syndrome fails; only then does it try additional outer/Turbo BP iterations with lower damping. The replay utility `experiments/gf256_cascade/replay_z_failures.py` identifies the original Z-syndrome failures and reruns the same seed and trial index. `GF256_CASCADE_SKIP_TARGET_DIAG=1` disables diagnostic-only oracle computations during replay; it does not enter the decoder's decision path.

The first 10 baseline Z failures yielded 6 rescues for star and 9 for path. These are selected failures, not random 10-frame FER estimates.

The full replay results are:

| Basis | Original logical failures | Original Z-syndrome failures rescued | Adaptive logical failures | Remaining Z-syndrome failures | Undetected logical failures |
|---|---:|---:|---:|---:|---:|
| star | 39/300 | 21/38 | **18/300** | 17 | 1 |
| path | 45/300 | 36/43 | **9/300** | 7 | 2 |

The original syndrome-valid logical failures (1 for star, 2 for path) were unchanged. The adaptive schedule produced no new syndrome-valid logical failures among the replayed baseline Z failures. The complete summaries are in `replay_adaptive_star_all/summary.json` and `replay_adaptive_path_all/summary.json`.

These 300-frame totals are exactly derived from targeted replay because `conditioned_robust_adaptive` follows the identical code path as `conditioned_robust` until the Z-syndrome failure branch. No initially successful frame is replaced. This is a **development-set comparison**: the schedule was chosen after inspecting these same original 300 frames.

## Independent channel seeds: path basis

With the same path matrix and `p_D=0.07125`, a fresh paired experiment used seeds `20269000`–`20269009`, five frames per seed (50 frames total). Both modes decoded each identical channel frame. Full source rows and transitions are in `fresh_path_50/summary.csv`, `pairwise.csv`, and `failure_transitions.csv`.

| Decoder | Logical failures | Z-syndrome failures | Undetected logical failures | CPU seconds/frame |
|---|---:|---:|---:|---:|
| `conditioned_robust` | 8/50 (0.16) | 8 | 0 | 25.52 |
| `conditioned_robust_adaptive` | 2/50 (0.04) | 2 | 0 | 29.12 |

The adaptive mode rescued 6 paired frames and regressed on 0. Exact two-sided McNemar p-value: 0.03125. CPU time increased by about 14.1%. The 95% Wilson intervals are 0.083–0.285 for 8/50 and 0.011–0.135 for 2/50. This fresh-seed result supports the improvement at this single matrix and channel point, but 50 frames yield wide FER intervals. More seeds and multiple `p_D` points are needed before making a general performance claim.
