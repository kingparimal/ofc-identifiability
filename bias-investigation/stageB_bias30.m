%% stageB_bias30.m -- 30-subject re-check of the possible upward bias in r.
%
% Same grid geometry as stageB_pilot2 (which the grid cost is dominated by, and
% which does not depend on the subject -- see IDENTIFIABILITY_RESULTS sec 5.2,
% so 30 subjects cost essentially no more than 8).
%
% Two changes from stageB_pilot2 that matter:
%   1. nSubj = 30 rather than 8.
%   2. best-fit stored as a grid INDEX as well as a value, so the comparison
%      against truth is done on indices and cannot be corrupted by the fact
%      that logspace places the r = 0.002 node 4.34e-19 above 0.002.
%
% Also reports the 3-statistic scoring alongside the 5-statistic one. The two
% zero-signal statistics (mean_t, cov_nt) are exactly zero in expectation for
% Sim 2, so the 3-stat version is the one with power; if the two disagree on
% the bias that is itself worth knowing.

clear; clc;
task   = 'sim2';
theta0 = [0.4; 0.4; 4; 0.002];
Nobs   = 200;
Npred  = 20000;
nSubj  = 30;
nRep   = 400;

nSU = 15;  su_range = linspace(0.30, 0.50, nSU);
nR  = 15;  r_range  = logspace(log10(2e-4), log10(2e-2), nR);
[~, iTrueSU] = min(abs(su_range - theta0(1)));
[~, iTrueR ] = min(abs(r_range  - theta0(4)));
fprintf('truth nodes: sigma_u %d/%d, r %d/%d\n', iTrueSU, nSU, iTrueR, nR);
fprintf('r node - truth = %.3e (nonzero: compare INDICES, not values)\n\n', ...
        r_range(iTrueR) - theta0(4));

Sigma5 = statsCovariance(theta0, Nobs, nRep, 60000, task);
idx3   = [1 3 4];
Sigma3 = Sigma5(idx3, idx3);
thr5 = hotelling(5, nRep);
thr3 = hotelling(3, nRep);
fprintf('thresholds: 5-stat %.3f   3-stat %.3f\n\n', thr5, thr3);

%% ---- predicted statistics on the grid (subject independent) -------------
Spred = nan(5, nSU, nR);
t0 = tic;
for i = 1:nSU
    for j = 1:nR
        th = [su_range(i); theta0(2); theta0(3); r_range(j)];
        Spred(:,i,j) = ofcSummary(th, Npred, 12345, task);
    end
    fprintf('  grid row %d/%d (%.0f s)\n', i, nSU, toc(t0));
end
save('bias30_grid.mat', 'Spred', 'su_range', 'r_range', 'theta0', 'Npred', 'task', '-v7');
fprintf('\n');

%% ---- score every subject both ways --------------------------------------
bR5 = nan(nSubj,1); bR3 = nan(nSubj,1);
bS5 = nan(nSubj,1); bS3 = nan(nSubj,1);
for k = 1:nSubj
    s_obs = ofcSummary(theta0, Nobs, 70000 + k, task);
    D5 = nan(nSU, nR); D3 = nan(nSU, nR);
    for i = 1:nSU
        for j = 1:nR
            D5(i,j) = discrepancy(s_obs,       Spred(:,i,j),       Sigma5);
            D3(i,j) = discrepancy(s_obs(idx3), Spred(idx3,i,j),    Sigma3);
        end
    end
    [~, m] = min(D5(:)); [bi, bj] = ind2sub(size(D5), m); bS5(k) = bi; bR5(k) = bj;
    [~, m] = min(D3(:)); [bi, bj] = ind2sub(size(D3), m); bS3(k) = bi; bR3(k) = bj;
end
save('bias30_fits.mat', 'bR5', 'bR3', 'bS5', 'bS3', 'iTrueR', 'iTrueSU', '-v7');

%% ---- report -------------------------------------------------------------
report('r, 5-stat',      bR5, iTrueR,  nSubj);
report('r, 3-stat',      bR3, iTrueR,  nSubj);
report('sigma_u, 5-stat [CONTROL]', bS5, iTrueSU, nSubj);
fprintf(['\nFirst 8 subjects only, for continuity with the recorded pilot:\n']);
report('r, 5-stat, n=8', bR5(1:8), iTrueR, 8);
