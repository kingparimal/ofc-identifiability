# Two published discretisations of the Todorov & Jordan (2002) muscle filter are not the same plant

MATLAB code accompanying "A specification discrepancy in the muscle filter of the Todorov–Jordan optimal feedback control model, and the convention that reproduces the published behavior<img width="468" height="73" alt="image" src="https://github.com/user-attachments/assets/59df3e8f-0ca2-4d4c-80b1-8a276e416a65" />
" (submitted to *eNeuro*).

Todorov & Jordan (2002) and Todorov (2005) describe the same
second-order muscle filter but discretise it differently. The two
conventions produce plants with DC gains differing by a factor of 20.4,
and because the effort term of the cost function is quadratic in the
control signal, this squares to a factor of ~418 in the effective effort
weight. Under the 2005 convention with the 2002 parameter values, the
modelled movement never decelerates and does not reach the target.

This repository reproduces that comparison and the manuscript's single
figure.

## Requirements

MATLAB R2025b. No toolboxes.

## Reproducing the figure

    addpath(genpath('.'))
    fig2_filter_conventions

Writes `Figure-1.tif` (400 dpi), `Figure-1.eps` (vector) and
`Figure-1.pdf` to the working directory, and prints the validated
quantities to the console. Runs in a few seconds.

The script asserts all six published values at 5% tolerance before it
plots, so it cannot silently produce a figure from a wrong model build.

## The two conventions

Per axis, with tau1 = tau2 = 40 ms and dt = 10 ms:

| | force update | input coefficient | DC gain per stage |
|---|---|---|---|
| 2002, exponential | `f(t+dt) = exp(-dt/tau)*f(t) + g(t)` | 1 | 4.5208 |
| 2005, Euler | `f(t+dt) = (1-dt/tau)*f(t) + (dt/tau)*g(t)` | dt/tau | 1 |

Second-order cascade: 20.438 against 1.

Three consequences follow from the same factor in `B`:

- **Effort weight, 418x.** The effort term is `r*|u|^2`, so a 20.4x
  change in the control needed for a given force squares.
- **Control-dependent noise scale, 16x.** `C = sigma_u*B`, and the noise
  term in the gain recursion goes as `sigma_u^2 * B'*S*B`, so the same
  nominal sigma_u means different things in the two plants.
- **Peak control magnitude, 6.3x.** Measured, see below.

## Validated values

Deterministic Sim 2 (interception), no noise realisations, no added
sensory delay, 50 steps at dt = 10 ms.

| quantity | 2002 (exponential) | 2005 (Euler) |
|---|---|---|
| terminal speed (m/s) | 0.0252 | 0.3638 |
| endpoint bias (cm) | 0.0776 | 4.7551 |
| peak \|u\| | 0.3178 | 1.9944 |

The movement starts with the two masses 20 cm apart. Under the 2002
convention it closes to within 0.08 cm and stops; under the 2005
convention it is still 4.76 cm short and still moving at 0.36 m/s when
the trial ends.

## Gain convention

Deterministic comparisons use **noise-free gains**: `C` is zeroed before
`solveLQGND`, so the controller is solved without the control-dependent
noise term, and the closed loop is then rolled forward with no noise
realisations.

This is not merely a convenience. `C = sigma_u*B`, and the Euler `B`
carries a factor dt/tau1 = 0.25 where the exponential `B` carries 1, so
holding sigma_u fixed across the two plants would vary the noise
contribution to the gain recursion by a factor of 16 in lockstep with
the discretisation under test. Zeroing `C` before solving removes that
confound and isolates the filter convention.

`fig2_filter_conventions.m` exposes a `GAINS` switch so the comparison
can be reproduced both ways. `'noisefree'` is the published setting;
`'nominal'` retains sigma_u = 0.4 in the control problem.

## Layout

    buildModel2D.m              plant construction, 2002 convention
    solveLQGND.m                coordinate-descent LQG solver
    fig2_filter_conventions.m   builds both plants, validates, plots
    simulink/                   independent re-implementation

The Euler plant is constructed inside the figure script by overwriting
`A` and `B`, so `buildModel2D.m` is unmodified from the identifiability
study and the two repositories can be compared directly.

## Relation to the other study

The same model underpins a companion analysis of parameter
identifiability <!-- FILL: citation or DOI once available -->. That
repository is separate; the shared model files here are a frozen
snapshot so this paper's figure is reproducible independently.



## Licence

MIT. See `LICENSE`.
