%% stageB_pilot.m  --  two-parameter pilot: does the r / sigma_u ridge exist?
%
% The cheapest test of prediction P1, run before committing to the full
% four-parameter grid. sigma_s and d are held at their nominal values and
% treated as known; only the effort weight r and the motor noise scale
% sigma_u are varied.
%
% P1 (registered in the proposal): r and sigma_u will form a RIDGE -- a
% diagonal confidence region rather than a compact one -- because their bias
% curves are nearly superimposed in the Stage 4 sweep.
%
% P2: the ridge will be PARTIALLY broken by the aspect ratio, since the two
% parameters move it in opposite directions, so the region should be finite
% but elongated rather than unbounded.
%
% WHAT WOULD FALSIFY P1: a compact, roughly circular region. That would mean
% the two parameters are separable from endpoint statistics alone at N = 200,
% and the confound Todorov describes in prose is not a practical problem.
%
% Cost: nGrid^2 solves at high N, plus one subject and one covariance at the
% realistic N. A few minutes.

clear;  clc;
task   = 'sim2';
theta0 = [0.4; 0.4; 4; 0.002];      % sigma_u, sigma_s, d, r
Nobs   = 200;                        % realistic trial count
Npred  = 20000;                      % for the predicted statistics
nGrid  = 9;

% Ranges from the Stage 4 sweep (REFERENCE_NUMBERS sec 2.6).
su_range = linspace(0.1,  0.7,   nGrid);
r_range  = linspace(5e-4, 4e-3,  nGrid);

fprintf('=== Stage B pilot: r vs sigma_u ===\n');
fprintf('    %s, N_obs = %d, %dx%d grid, sigma_s and d held at nominal\n\n', ...
        task, Nobs, nGrid, nGrid);

%% ---- the synthetic subject -------------------------------------------
% One experiment, at the true parameters, with a realistic number of trials.
s_obs = ofcSummary(theta0, Nobs, 777, task);

% Sampling covariance at that trial count.
Sigma = statsCovariance(theta0, Nobs, 400, 50000, task);

% Threshold. chi2inv(0.95,5) = 11.07 for a known Sigma. Because Sigma is
% ESTIMATED from nRep replicates the statistic is Hotelling-like rather than
% exactly chi-square, and the correct threshold is slightly larger. With
% nRep = 400 and p = 5 the correction is about 2%; both are printed so the
% choice is visible rather than buried.
p    = 5;  nRep = 400;
thr  = chi2inv_local(0.95, p);
thrH = (p*(nRep-1)/(nRep-p)) * finv_local(0.95, p, nRep-p);
fprintf('  threshold  chi2  %.2f      Hotelling-corrected  %.2f\n\n', thr, thrH);

%% ---- the grid --------------------------------------------------------
D = nan(nGrid, nGrid);              % rows: sigma_u, cols: r
AR = nan(nGrid, nGrid);
t0 = tic;
for i = 1:nGrid
    for j = 1:nGrid
        th = [su_range(i); theta0(2); theta0(3); r_range(j)];
        [s_pred, info] = ofcSummary(th, Npred, 12345, task);
        D(i,j)  = discrepancy(s_obs, s_pred, Sigma);
        AR(i,j) = info.AR;
    end
    fprintf('  row %d/%d  (%.0f s elapsed)\n', i, nGrid, toc(t0));
end
fprintf('\n');

%% ---- read off the result ---------------------------------------------
inside = D <= thrH;
nIn    = sum(inside(:));

[~, k]    = min(D(:));
[bi, bj]  = ind2sub(size(D), k);

fprintf('--- result ---\n');
fprintf('  best-fit grid point : sigma_u %.3f, r %.2e   (D = %.2f)\n', ...
        su_range(bi), r_range(bj), D(bi,bj));
fprintf('  truth               : sigma_u %.3f, r %.2e\n', theta0(1), theta0(4));
fprintf('  grid points inside 95%% region : %d of %d (%.0f%%)\n', ...
        nIn, numel(D), 100*nIn/numel(D));

% Is the truth inside its own confidence region? It should be, 95% of the
% time. A single subject can miss; this is a sanity check, not a test.
[~, ti] = min(abs(su_range - theta0(1)));
[~, tj] = min(abs(r_range  - theta0(4)));
fprintf('  truth inside region : %s  (D = %.2f)\n', ...
        string(inside(ti,tj)), D(ti,tj));

%% ---- quantify the shape ----------------------------------------------
% A ridge is an ELONGATED region. Measure that rather than eyeballing it:
% take the grid points inside the region, normalise each axis to its range
% so the two parameters are comparable, and compare the spread along the
% principal axes.
if nIn >= 3
    [I, J] = find(inside);
    x = (su_range(I)' - min(su_range)) / range(su_range);
    y = (r_range(J)'  - min(r_range))  / range(r_range);
    Creg = cov([x y]);
    ev   = sort(eig(Creg), 'descend');
    elong = sqrt(ev(1)/max(ev(2), eps));
    [V,~] = eig(Creg);
    ang = atan2d(V(2,end), V(1,end));
    fprintf('  region elongation (normalised) : %.2f\n', elong);
    fprintf('  principal axis                 : %.0f deg\n', ang);
    fprintf('\n');
    if elong > 2.5
        fprintf('  >>> RIDGE. P1 supported: the region is elongated, so r and\n');
        fprintf('      sigma_u trade off against each other at N = %d.\n', Nobs);
    else
        fprintf('  >>> COMPACT. P1 not supported at this trial count: the two\n');
        fprintf('      parameters are separable from endpoint statistics.\n');
    end
else
    fprintf('  region too small to characterise (%d points). Refine the grid.\n', nIn);
end
fprintf('\n');

%% ---- figure ----------------------------------------------------------
figure('Color','w','Position',[100 100 980 400]);

subplot(1,2,1); hold on; grid on;
contourf(r_range, su_range, D, 20, 'LineColor','none');
contour(r_range, su_range, D, [thrH thrH], 'w-', 'LineWidth', 2.5);
plot(theta0(4), theta0(1), 'w+', 'MarkerSize', 14, 'LineWidth', 2.5);
plot(r_range(bj), su_range(bi), 'ro', 'MarkerSize', 9, 'LineWidth', 2);
colorbar; xlabel('effort weight  r'); ylabel('motor noise  \sigma_u');
title(sprintf('Discrepancy, N = %d  (white = 95%% region)', Nobs));
legend('', '95% boundary', 'truth', 'best fit', 'Location','northeast');

subplot(1,2,2); hold on; grid on;
contourf(r_range, su_range, AR, 20, 'LineColor','none');
contour(r_range, su_range, AR, [3.35 3.35], 'w-', 'LineWidth', 2);
plot(theta0(4), theta0(1), 'w+', 'MarkerSize', 14, 'LineWidth', 2.5);
colorbar; xlabel('effort weight  r'); ylabel('motor noise  \sigma_u');
title('Aspect ratio (white = nominal 3.35)');

% The second panel is the mechanism behind P2. If the aspect-ratio contours
% run ACROSS the ridge rather than along it, that statistic is what keeps
% the region finite. If they run parallel to it, the aspect ratio is not
% breaking the confound and the region is bounded by something else.

%% ---------------------------------------------------------------- helpers
function x = chi2inv_local(p, k)
    if exist('chi2inv','file') == 2, x = chi2inv(p, k); return, end
    lo = 0; hi = 1;
    while gammainc(hi/2, k/2) < p, hi = hi*2; end
    for i = 1:200
        m = (lo+hi)/2;
        if gammainc(m/2, k/2) < p, lo = m; else, hi = m; end
    end
    x = (lo+hi)/2;
end

function x = finv_local(p, d1, d2)
    if exist('finv','file') == 2, x = finv(p, d1, d2); return, end
    lo = 0; hi = 1;
    while betainc(d1*hi/(d1*hi+d2), d1/2, d2/2) < p, hi = hi*2; end
    for i = 1:200
        m = (lo+hi)/2;
        if betainc(d1*m/(d1*m+d2), d1/2, d2/2) < p, lo = m; else, hi = m; end
    end
    x = (lo+hi)/2;
end
