# Provenance of published values

Maps each figure, table and headline number in the manuscript to the
script that produces it. Values are current; superseded figures from
earlier drafts are not reproduced here.

Nominal parameters throughout: θ = [σ_u, σ_s, d, r] = [0.4, 0.4, 4,
0.002]. Total sensory lag (1+d)·Δt = 50 ms at nominal d. 30 synthetic
subjects, N = 200 trials each, N_pred = 20000, nRep = 400. All results
produced under MATLAB R2025b.

**Subject 1 unless stated.** minD1, minD2, and the arrays in
`thesis_fig_data.mat` (`prof1`, `prof2`, `D_su_r`, `D_ss_d`) are
subject 1, not aggregates over the 30.

**Provenance markers.** (O) marks Octave 8.4.0 output; unmarked values
are MATLAB. Octave-derived values should be regenerated before
submission.

## Figures

| manuscript | script | notes |
|---|---|---|
| Fig. 1 | `figures/fig3_tier1_escape.m` | reads `thesis_fig_data.mat` |
| Fig. 2 | `figures/fig4_compensation.m` | recomputes clouds; null vectors hardcoded, see below |
| Fig. 3 | `figures/fig5_tier2_profiles.m` | reads `thesis_fig_data.mat` |
| Fig. 4 | `figures/fig6_widths_vs_N.m` | values inline in script |

`figures/fig1_endpoint_clouds.m` is not a manuscript figure. It produces
the aspect ratio (≈3.3) cited in methods and is retained for that value.

Figure scripts must be regenerated for submission at 174 mm width,
`-depsc` rather than `-dpdf`, lowercase panel letters, and descriptive
titles moved from the axes into the captions.

## Calibration constants

| quantity | value | source |
|---|---|---|
| tier 1 threshold, 3 statistics | 7.9219 | `stats/hotelling.m` |
| tier 2 threshold, 18 statistics | 30.6630 | `stats/hotelling.m` |
| single-parameter profile threshold | 3.8415 | χ², 1 dof |
| min D, tier 1, subject 1 | 0.1109 | |
| min D, tier 2, subject 1 | 12.1407 | |

30.6630 is the Hotelling value for 18 statistics at nRep = 400, which
independently confirms the tier 2 summary vector's dimension.

## Structural result

| quantity | value | script |
|---|---|---|
| cond(Σ), 3 genuine statistics | 2.60e5 | |
| cond(Σ), 5 reported statistics | 4.83e12 | |
| matrix equivariance, Sim 2 | all exactly 0 | `robustness/sym_check.m` |
| matrix equivariance, Sim 1 | fails on Q and x₁ only | `robustness/sym_check.m` |
| gain equivariance, ΔL / ΔK | 2.4e-14 / 1.4e-15 | `robustness/sym_check.m` |
| pathwise, d = 0 / d = 4 | 5.44e-15 / 8.66e-15 | `robustness/path_sym2.m` |
| largest cross-block correlation | 0.070 (SE 0.050) | |

CONFLICT: a second source gives the pathwise residual as 2.988e-16. Both
are zero to machine precision, but one value must be chosen.

This is the justification for scoring 3 statistics rather than 5. Tier 2
G2 gate: `tier2/gate_tier2.m`.

## Tier 1

| quantity | value | script |
|---|---|---|
| rank of F | 3 of 4 | `tier1/sizing.m` |
| accepted region volume | 1.22% of box | `tier1/escape.m` |
| subjects below threshold on outer shell | 28 of 30 | `tier1/escape.m` |
| min D on outer shell | 20% of threshold | `tier1/escape.m` |
| median walls reached | see note | `tier1/escape.m`, `tier1/walls.m` |
| correlation σ_u / r | +0.838 | |
| correlation σ_s / d | −0.65 | see caveat |

**cond(F) is NOT reportable at tier 1.** The fourth eigenvalue is
exactly zero in exact arithmetic, so the ratio is set by floating-point
noise: across three seeds it returns 5.98e18, Inf, and 2.39e18, with
fourth eigenvalues 1.53e-14, 1.28e-13, 3.18e-13. Report rank and the
eigenvalue spectrum (largest ≈2.6e3, smallest ≈1e-14), never the ratio.
Tier 2's cond(F) is stable across the same seeds (1.59–1.62e3) and may
be quoted.

**Wall count convention unresolved.** `escape.m` counts 6 walls, median
3; `walls.m` counts 8, median 5. Four parameters give 8 box faces, so
the 6-wall count excludes two. Choose one, state it in the text, and
make the figure agree.

**σ_s/d correlation is qualitative only.** No subject's region closes in
both parameters, so it is measured on a set the sweep box truncates.
Dropping wall-touching points shifts it by roughly 0.1 (O, illustrative
— regenerate before quoting digits). Report to two decimals and never
present the tier 1 / tier 2 near-equality as evidence.

## Tier 2

| quantity | value | script |
|---|---|---|
| rank of F | 4 of 4 | `tier1/sizing.m`, `robustness/robustness.m` |
| cond(F) | 2.07e3 | `tier1/sizing.m` |
| closure margin | 237% of threshold | `tier2/escape_t2.m` |
| subjects below threshold | 0 of 30 | `tier2/escape_t2.m` |
| walls reached | 0 | `tier2/escape_t2.m` |
| correlation σ_u / r | +0.397 | |
| correlation σ_s / d | −0.66 | |

`gate_tier2.m` does not compute F. Correlations are of grid indices,
equivalent to correlating log-parameters on the geometric axes and d on
the linear delay axis.

## Interval widths at N = 200

Full factors, hi/lo. Delay in 10 ms steps. Every width in the manuscript
is labelled GRID or CURVATURE at point of use.

| parameter | GRID | CURVATURE | grid/curvature |
|---|---|---|---|
| σ_u | 1.11 | 1.117 | 94% |
| σ_s | 1.33 | 1.409 | 83% |
| r | 1.23 | 1.330 | 73% |
| d | see below | 1.71 steps | |

Ratios are of log-widths, 73–94%; curvature is the conservative bound.
Source: `figures/fig6_widths_vs_N.m`; curvature values from
`tier2/nsweep_fisher.m`.

**Widths scale as 1/√N in log-width for the three multiplicative
parameters, and directly in steps for the delay.** The factor itself
does not scale as 1/√N — plotting it on log-log gives a slope near −0.1.

Delay, per subject: no interval spans more than one 10 ms step. All 30
contain the true value — 12 extend below, 3 exact, 15 above. Realised
coverage 30/30 exceeds nominal 95% because quantisation makes an integer
parameter's coverage jump rather than track smoothly. Script:
`tier1/delay_report.m` (currently reports in 1+d with truth 5; the
manuscript uses d with truth 4).

## Coverage

Two independent draws, 93.2% and 94.5%, pooling to 93.9% over 1200
subjects. Both draws are `tier2/coverage_gate2.m` under different seeds.
`tier1/coverage_gate.m` is the tier 1 three-statistic gate and returns
95.3% — a separate quantity.

## Compensation direction

Null space of the 3×4 Jacobian at nominal, columns scaled by nominal
values. Script: `robustness/compensation_direction.m`, line 38,
`v = null(Jrel)`.

    v_rel = [−0.084, −0.401, +0.906, −0.105]

Independently in log coordinates, delay column rescaled by (1+d)·ln10:

    v_log = [−0.0969, −0.4846, +0.8618, −0.1142]     (O)

The two directions agree to **5.5°**. Report the angle, not
componentwise agreement: componentwise they differ by up to 21% on σ_s,
but they are differently normalised so componentwise comparison is not
meaningful.

`fig4_compensation.m` hardcodes both vectors and computes neither.
`v_log` is Octave-derived; the rescaling that produces it appears in
`fig_rank.m` and `check_fig3.m`. Confirm which generated the reported
digits and regenerate in MATLAB.

Displaced point: σ_s halved, delay 4 → 9, σ_u −8%, r −13%. Endpoint
scatter: SD_irrel 1.4852 vs 1.4834 cm, a 0.1% difference. Aspect ratio
3.32 vs 3.42 is a 3% difference against a ±2% Monte Carlo band — quote
SD_irrel, and if the AR is quoted say "within Monte Carlo error" rather
than "near-identical".

Two-stage check, reported as a negative result: mean control magnitude
differs by 4.7% at peak (0.3669 vs 0.3843), max|L| by 4.4%; state
estimation error runs ~15% lower at the displaced point throughout, in
one direction with no crossing.

CAUTION: 4.7% sits close to the ±2% MC band. Confirm it is outside noise
and lead the negative result with the estimation-error finding, which is
larger and directionally wrong for a trade.

## Go/no-go test

Subject 1: 494 tier-1-accepted vectors, 54 on a low wall, 45 low-σ_s and
9 high-σ_s. Eight sampled points, all accepted at tier 1, all rejected
at tier 2. Weakest 4.3× threshold (low-σ_s arm), strongest 305×
(high-σ_s arm) — a seventy-fold spread.

Wall points must be split by σ_s relative to truth before sampling;
sorting on discrepancy alone lands on whichever arm dominates.

## Robustness

Tier 1 condition numbers are omitted deliberately (see Tier 1 above).

| variant | tier 1 | tier 2 cond(F) |
|---|---|---|
| Sim 2, P = 1 (adopted, Q scale 1/3) | rank 3 of 4 | 2.07e3 |
| Sim 2, P = 2 (Q scale 1/4) | rank 3 of 4 | 2.31e3 |
| Sim 1, P = 1 | rank 3 of 4 | 2.82e3 |

Script: `robustness/robustness.m`. Under Q scale 1/4, tier 2 widths move
by at most 10%: σ_u ×1.117 → ×1.114, σ_s ×1.409 → ×1.441, delay 1.71 →
1.88 steps, r ×1.330 → ×1.291.

Task geometry: Sim 1 recovers r better than Sim 2 (×1.218 vs ×1.330) but
σ_s substantially worse (×2.130 vs ×1.409).

Seed stability of the ∂mean_t/∂θ row, relative to ∂mean_n/∂θ as control:
Sim 2 (proved zero, calibration) 160×, Sim 1 152×. Script:
`robustness/sim1_5stat.m`.

## Registered predictions

P1 pilot: 7 of 8 regions closed with σ_s and d held fixed, against 9 of
30 with them free. Scripts: `bias-investigation/stageB_pilot.m`,
`stageB_pilot2.m`. This is the evidence that fixing noise parameters
manufactures identifiability, and it is load-bearing for the discussion.
