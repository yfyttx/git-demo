# PRA 103, 022617 figure and metric benchmark for this project

Reference: P. Fuentes et al., “Design of low-density-generator-matrix-based quantum codes for asymmetric quantum channels,” *Physical Review A* **103**, 022617 (2021), DOI: 10.1103/PhysRevA.103.022617. The supplied file is the accepted manuscript.

## What the paper measures and plots

The paper reports two post-decoding quantities:

- **QBER**: fraction of physical qubits that retain an error after decoding.
- **WER**: fraction of blocks with at least one residual physical-qubit error.

Its principal performance plots use a logarithmic error-rate axis and gross Pauli flip probability (p=p_X+p_Y+p_Z) on the horizontal axis. For asymmetric noise it fixes the asymmetry ratio α and varies (p). Figures 4–8 compare CSS schemes and constituent X/Z decoders; Fig. 6 explicitly plots X component, Z component, and the full CSS scheme together. Figs. 9–10 vary the underlying LDGM degree and block length. Fig. 11 compares several channel asymmetries. Fig. 12 summarizes achievable rate at a target WER of (10^{-3}) against the Pauli-channel hashing bound.

For the paper's depolarizing case, α=1 and (p_X=p_Y=p_Z=p/3), so the marginal X and Z flip probability is (f_m=2p/3). In this project the equivalent gross channel parameter is (p_D); thus a depolarizing-channel curve should use (p_D) on the x-axis, or explicitly convert from (f_m). The paper's asymmetric curves require a different channel generator and should be presented as a separate experiment.

## Metric mapping for a degenerate quantum code

The paper's physical QBER/WER are useful secondary metrics, but they are not the same as logical FER. A quantum decoder may return a different physical error vector that differs from the true error by a stabilizer; this is physically nonzero residual support but still logically successful. Therefore this project should report all levels explicitly:

1. **Primary quantum metric — logical FER**: fraction of frames for which either residual Pauli component is outside the corresponding stabilizer subgroup. Report X-logical, Z-logical, and union (overall) rates.
2. **Syndrome failure**: fraction whose decoded estimate fails its complete measured syndrome. Report X and Z separately.
3. **PRA-style physical QBER**: residual physical-qubit fraction after applying the estimated correction. Report X, Z, and Pauli-union residual fractions.
4. **PRA-style physical WER**: fraction of frames with any residual physical support, again separated into X, Z, and union. Label this “physical residual WER”; do not call it logical FER.
5. **Undetected logical failure**: syndrome passes but the residual is outside the stabilizer subgroup. Keep this visible because syndrome convergence alone is insufficient.

The paper's figures do not show confidence intervals or state a per-point Monte Carlo frame budget in the inspected simulation section. For a submission-quality reproduction, this project should add uncertainty intervals and a prespecified stopping rule rather than copying that omission.

## Current matched reference point from this project

Fresh reruns use the same matrix, same (p_D=0.07125), same seeds, 100 iterations, conditioned Z prior, and no X fallback. Two 500-frame runs were combined only after checking that matrix, channel, decoder, and protocol matched.

| Quantity | Result, 1000 frames |
|---|---:|
| Code parameters | (N=5238, k=1552, R_Q=8/27approx0.2963), GF(256) outer symbols with binary companion-block expansion |
| Channel | Depolarizing, (p_D=0.07125, f_m=0.0475) |
| Overall logical FER | (69/1000=6.9\%), Wilson 95% CI ([5.49\%,8.64\%]) |
| X / Z syndrome failures | (68/1000) / (68/1000) |
| Undetected logical failures | (0/1000) in this run |
| Residual QBER X / Z / Pauli union | 0.325% / 4.849% / 5.108% |
| Physical residual WER X / Z / union | 6.8% / 100% / 100% |

The 100% Z physical residual WER does **not** mean all frames are logical failures: stabilizer-equivalent physical residuals can be logically harmless. That is exactly why the manuscript should keep physical QBER/WER and logical FER in separate columns/curves. This one point is an audit reference, not a performance curve or evidence of paper-ready performance.

## First fixed-matrix channel scan

To follow the paper's gross-flip-probability sweep while keeping the method fixed, a six-point scan was run with the same `20261019` matrix, 100 BP iterations, conditioned Z prior, and no X fallback. The initial 0.0600 and 0.0650 points used 200 frames; the 0.0675 point used 1000 frames; fresh-seed extensions raised the 0.0685 and 0.0705 points to 4200 and 2700 frames; the 0.07125 point combines the two compatible 500-frame runs above.

| (p_D) | Frames | Overall logical failures | Logical FER | X / Z syndrome failures | Residual QBER X / Z / union |
|---:|---:|---:|---:|---:|---:|
| 0.06000 | 200 | 0 | 0% | 0 / 0 | 0% / 4.018% / 4.018% |
| 0.06500 | 200 | 1 | 0.5% | 1 / 1 | 0.034% / 4.291% / 4.319% |
| 0.06750 | 1000 | 9 | 0.9% | 9 / 6 | 0.043% / 4.532% / 4.569% |
| 0.06850 | 4200 | 85 | 2.024% | 84 / 70 | 0.093% / 4.608% / 4.685% |
| 0.07050 | 2700 | 116 | 4.296% | 115 / 106 | 0.206% / 4.785% / 4.951% |
| 0.07125 | 1000 | 69 | 6.9% | 68 / 68 | 0.325% / 4.849% / 5.108% |

The sweep shows a rising logical-failure trend over the tested points. The two low-noise points remain small pilots: one event at 0.065 has a Wilson 95% interval of about 0.09–2.78%, and zero events at 0.060 has a one-sided exact 95% upper bound of 1.49%. The transition-region points have 9, 85, and 116 events, with respective Wilson 95% intervals 0.47–1.70%, 1.64–2.50%, and 3.59–5.13%. The Z-side physical residual QBER remains around 4–5% while logical FER is much lower, consistent with many residual vectors being stabilizer-equivalent. The scan locates the rising region but is not yet a final low-FER threshold estimate; low-noise points and independent matrix-seed replication remain missing.

## Distinct-matrix check at the transition point

A second outer GF(256) matrix pair of the same dimensions was expanded with the same `[9,8,2]` inner code and decoder. The constructor verified CSS orthogonality, full expected ranks, and the same `(N,k,R_Q)=(5238,1552,8/27)`. At `p_D=0.0685`, 11/1000 frames failed logically (1.1%; Wilson 95% CI 0.62–1.96%), with 11 X-syndrome failures, 9 Z-syndrome failures, and no undetected logical failures. Residual QBER X/Z/union was 0.055%/4.561%/4.608%.

The second matrix was then swept at three additional points, with 1000 frames per point:

| `p_D` | Frames | Logical failures | FER (Wilson 95% CI) | X / Z syndrome failures | Residual QBER X / Z / union |
|---:|---:|---:|---:|---:|---:|
| 0.06500 | 1000 | 3 | 0.30% (0.10–0.88%) | 3 / 1 | 0.013% / 4.375% / 4.387% |
| 0.06850 | 1000 | 11 | 1.10% (0.62–1.96%) | 11 / 9 | 0.055% / 4.561% / 4.608% |
| 0.07050 | 1000 | 41 | 4.10% (3.04–5.51%) | 41 / 37 | 0.200% / 4.800% / 4.960% |
| 0.07125 | 1000 | 60 | 6.00% (4.69–7.65%) | 60 / 57 | 0.312% / 4.867% / 5.115% |

The two matrices show a similar increasing trend; their intervals overlap at the common points. At `p_D=0.0685`, matrix 1 gives 85/4200 (2.024%; Wilson 95% CI 1.640–2.496%) versus 11/1000 for matrix 2. The point estimates differ, but the intervals overlap; these data do not establish that one matrix is better. Each second-matrix point has only 3–60 events, so this is still a replication pilot rather than a precision performance comparison. The matrix-level comparison also uses a single realization per design.

Internal plots (not manuscript-final) and data: `results/pra_metric_reference_sweep_20261118/`. `pra_metric_reference_pilot.pdf` shows the first matrix's X/Z/overall logical frame failures with intervals and residual physical QBER with frame-cluster bootstrap intervals. `pra_metric_two_matrix_fer.pdf` compares overall logical FER on the two outer-matrix realizations. Figure source and command are `plot_pra_metric_reference_sweep.py`.

## Figure plan modeled on the reference

### Main performance figure (analogous to Figs. 4 and 6)

- X-axis: gross depolarizing probability (p_D); log y-axis.
- Panel A: overall logical FER with Wilson 95% intervals; include X-logical and Z-logical curves as secondary series only if legible.
- Panel B: X/Z syndrome-failure rates and undetected logical failure, to distinguish decoder nonconvergence from wrong logical class.
- A separate supplemental figure can show PRA-style QBER and physical residual WER; caption must say these are physical residual metrics and not logical failure rates.
- Use fixed matrix, fixed decoder, and prespecified iteration/stopping budget across points. Record number of frames and failures at every point; for low FER, stop after a predeclared maximum frame count or target number of failures and report exact binomial intervals.

### Ablation figures (analogous to Figs. 5, 7–10)

- Decoder ablation: baseline, conditional prior, full-matrix fallback, candidate-list/OSD refinement; paired same-noise comparisons and CPU per frame.
- Inner-code ablation: compare (H_2,G_2) candidates only after matching physical length, quantum rate, outer matrix, and decoder budget. Include girth/short-cycle counts as explanatory properties, not as a substitute for FER.
- Length scaling: at fixed (R_Q), inner-code family, channel, and decoder, compare at least two physical block lengths.

### Asymmetric-channel figure (analogous to Fig. 11)

Only after an asymmetric Pauli channel is implemented and independently validated, plot logical FER against gross (p) in separate panels for fixed α values (for example α=1, (10^2), and a larger value). Use (p_X,p_Y,p_Z) from one explicitly stated channel convention; do not label the current depolarizing results as asymmetric.

### Rate benchmark (analogous to Fig. 12)

For each code rate, estimate the largest channel (p_D) (or fixed-α gross (p)) meeting a predeclared logical-FER target. This is the quantum-appropriate analogue of rate at WER (10^{-3}); if also reproducing the paper's physical WER target, show it as a distinct benchmark.

## Evidence gap and publication-readiness assessment

This update has produced a first fixed-matrix multi-point logical-FER curve, rather than only a single reference point. It now supports a descriptive claim that FER rises from 0.9% at `p_D=0.0675` (9/1000) to 2.024% at `0.0685` (85/4200), 4.296% at `0.0705` (116/2700), and 6.9% at `0.07125` (69/1000). The 0.0685 and 0.0705 points have useful event counts, but the .0600/.0650 points remain under-sampled and all points use one outer-matrix seed. The new 9-bit construction also has significant X-side syndrome nonconvergence (e.g. 84/4200 at .0685, 115/2700 at .0705), so the curve characterizes this implementation and is not evidence that the construction has reached its achievable decoding performance.

This is a credible prototype and debugging/evaluation package, and now includes a first multi-point curve on two distinct outer-matrix realizations. It is not yet a complete publication-grade evidence set for a new method: the second curve is a pilot with only 1000 frames per point, low-FER sampling remains sparse, and there is only one matrix realization per design. Still missing are adequately powered, prespecified low-FER comparisons; a matched baseline or ablation on exactly the same code and noise frames; rate/length scaling; and an asymmetric-channel study if asymmetry is central to the contribution. The old direct Kasai runs are a different construction and decoder, so they cannot fill these gaps for the current cascade.

Next: freeze the present decoder as a baseline; increase frame counts on both matrices around the transition region; improve and evaluate X-side recovery without selecting methods on the final holdout; then run a paired, untouched comparison against the frozen baseline. After that, add the rate/length and channel-asymmetry figures required by the intended contribution. The present figures should be labeled exploratory in talks and drafts until these checks are complete; no experimental workflow can guarantee acceptance.

## Traceability

- Metric implementation: `experiments/gf256_cascade/run_companion_staged_xz.py`
- Two 500-frame source runs: `results/companion_pra_metrics_seed20261125.json`, `results/companion_pra_metrics_seed20261126.json`
- Expanded transition-point runs: `results/companion_pra_sweep_pD0p0685_seed20261200.json`, `results/companion_pra_sweep_pD0p0705_seed20261201.json`
- Distinct outer-matrix check: `results/companion_binary_n9_k8_outer_path2/manifest.json`, `results/companion_pra_matrix2_pD0p0685_seed20261300.json`
- Additional second-matrix points: `results/companion_pra_matrix2_pD0650_seed20261301.json`, `results/companion_pra_matrix2_pD0705_seed20261302.json`, `results/companion_pra_matrix2_pD07125_seed20261303.json`
- Two-matrix comparison figure/data: `results/pra_metric_reference_sweep_20261118/pra_metric_two_matrix_fer.pdf`, `two_matrix_fer_data.csv`
- Aggregation: `results/companion_pra_metrics_1000_20261118.json`
- Decoder failure diagnosis: `X_FAILURE_MAP_DIAGNOSIS_20261118.md`
