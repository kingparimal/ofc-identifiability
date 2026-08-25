%% stageB_pilot2.m  --  refined r / sigma_u pilot
%
% The first pilot found a ridge but could not measure it: the accepted region
% was thinner than one grid cell in sigma_u, so every accepted point fell in
% a single row, the covariance of those points was rank 1, and the elongation
% metric returned 1.8e7. The region was also OPEN at the upper r boundary, so
% the extent of the ridge was unknown.
%
% Three changes:
%
%   1. sigma_u resolved finely over a narrow range around the region;
%      r swept over a much wider range, logarithmically, to find its
%      upper bound if one exists.
%
%   2. The shape metric no longer uses eigenvalues of the accepted grid
%      points. It reports the extent of the region along each parameter as a
%      fraction of that parameter's swept range, and flags whether the region
%      touches a boundary. That is well defined however thin the region is,
%      and it distinguishes "narrow" from "unbounded", which the eigenvalue
%      ratio cannot.
%
%   3. Several synthetic subjects rather than one. The first pilot's best fit
%      landed at r = 3.56e-3 against a truth of 2.0e-3; whether that is the
%      ridge or that one subject's noise cannot be told from a single draw.
%      This is nearly free: the grid of PREDICTED statistics does not depend
%      on the subject, so it is computed once and scored repeatedly.

clear;  clc;
task   = 'sim2';
theta0 = [0.4; 0.4; 4; 0.002];
Nobs   = 200;
Npred  = 20000;
nSubj  = 8;

nSU = 15;   su_range = linspace(0.30, 0.50, nSU);          % fine, narrow
nR  = 15;   r_range  = logspace(log10(2e-4), log10(2e-2), nR);  % wide, log

fprintf('=== Refined pilot: r vs sigma_u ===\n');
fprintf('    sigma_u %.2f..%.2f (%d pts, step %.4f)\n', ...
        su_range(1), su_range(end), nSU, su_range(2)-su_range(1));
fprintf('    r %.1e..%.1e (%d pts, log spaced)\n', ...
        r_range(1), r_range(end), nR);
fprintf('    %d synthetic subjects at N = %d\n\n', nSubj, Nobs);

Sigma = statsCovariance(theta0, Nobs, 400, 60000, task);
p = 5;  nRep = 400;
thr = (p*(nRep-1)/(nRep-p)) * finv_local(0.95, p, nRep-p);
fprintf('  threshold %.2f\n\n', thr);

%% ---- predicted statistics on the grid (subject independent) ----------
Spred = nan(5, nSU, nR);
AR    = nan(nSU, nR);
t0 = tic;
for i = 1:nSU
    for j = 1:nR
        th = [su_range(i); theta0(2); theta0(3); r_range(j)];
        [s, info] = ofcSummary(th, Npred, 12345, task);
        Spred(:,i,j) = s;
        AR(i,j)      = info.AR;
    end
    fprintf('  grid row %d/%d (%.0f s)\n', i, nSU, toc(t0));
end
fprintf('\n');

%% ---- score each subject ----------------------------------------------
suLo = nan(nSubj,1); suHi = nan(nSubj,1);
rLo  = nan(nSubj,1); rHi  = nan(nSubj,1);
bestSU = nan(nSubj,1); bestR = nan(nSubj,1);
openTop = false(nSubj,1);  openBot = false(nSubj,1);
nIn = nan(nSubj,1);

for k = 1:nSubj
    s_obs = ofcSummary(theta0, Nobs, 70000 + k, task);

    D = nan(nSU, nR);
    for i = 1:nSU
        for j = 1:nR
            D(i,j) = discrepancy(s_obs, Spred(:,i,j), Sigma);
        end
    end
    if k == 1, D1 = D; end

    inside  = D <= thr;
    nIn(k)  = sum(inside(:));
    if nIn(k) == 0, continue, end

    [~, m]  = min(D(:));  [bi, bj] = ind2sub(size(D), m);
    bestSU(k) = su_range(bi);  bestR(k) = r_range(bj);

    [I, J]  = find(inside);
    suLo(k) = su_range(min(I));  suHi(k) = su_range(max(I));
    rLo(k)  = r_range(min(J));   rHi(k)  = r_range(max(J));
    openTop(k) = any(J == nR);   openBot(k) = any(J == 1);
end

%% ---- report ----------------------------------------------------------
fprintf('--- per-subject 95%% regions ---\n');
fprintf('  %3s %6s %-18s %-22s %s\n', 'k', 'nIn', 'sigma_u range', 'r range', 'open?');
for k = 1:nSubj
    if nIn(k) == 0
        fprintf('  %3d %6d  (empty)\n', k, 0);  continue
    end
    tag = '';
    if openTop(k), tag = [tag 'TOP ']; end
    if openBot(k), tag = [tag 'BOTTOM']; end
    fprintf('  %3d %6d  %.3f - %.3f      %.2e - %.2e   %s\n', ...
            k, nIn(k), suLo(k), suHi(k), rLo(k), rHi(k), tag);
end

fprintf('\n--- summary ---\n');
suW = (suHi - suLo) / range(su_range);
rW  = (log10(rHi) - log10(rLo)) / (log10(r_range(end)) - log10(r_range(1)));
fprintf('  sigma_u: region spans %.0f%% of swept range (median over subjects)\n', ...
        100*median(suW, 'omitnan'));
fprintf('  r      : region spans %.0f%% of swept range (log scale, median)\n', ...
        100*median(rW, 'omitnan'));
fprintf('  best-fit sigma_u: %.3f +/- %.3f   (truth %.3f)\n', ...
        mean(bestSU,'omitnan'), std(bestSU,'omitnan'), theta0(1));
fprintf('  best-fit r      : %.2e (median), spread %.2e to %.2e   (truth %.2e)\n', ...
        median(bestR,'omitnan'), min(bestR), max(bestR), theta0(4));
fprintf('  subjects whose region is OPEN at high r : %d of %d\n', ...
        sum(openTop), nSubj);

fprintf('\n');
if median(rW,'omitnan') > 3*median(suW,'omitnan')
    fprintf(['  >>> RIDGE along r. sigma_u is pinned; r is not. The region is\n' ...
             '      %.0fx wider in r than in sigma_u as a fraction of the\n' ...
             '      swept range.\n'], median(rW,'omitnan')/median(suW,'omitnan'));
    if sum(openTop) > 0
        fprintf(['      %d subject(s) have a region reaching the top of the r\n' ...
                 '      sweep -- r is bounded below but not above within it.\n'], ...
                 sum(openTop));
    end
else
    fprintf('  >>> Region comparably tight in both. Refine or widen further.\n');
end
fprintf('\n');

%% ---- figure ----------------------------------------------------------
figure('Color','w','Position',[100 100 980 420]);

subplot(1,2,1); hold on; grid on;
contourf(log10(r_range), su_range, min(D1, 200), 25, 'LineColor','none');
contour(log10(r_range), su_range, D1, [thr thr], 'w-', 'LineWidth', 2.5);
plot(log10(theta0(4)), theta0(1), 'w+', 'MarkerSize', 14, 'LineWidth', 2.5);
plot(log10(bestR(1)), bestSU(1), 'ro', 'MarkerSize', 9, 'LineWidth', 2);
colorbar; xlabel('log_{10} effort weight  r'); ylabel('motor noise  \sigma_u');
title(sprintf('Subject 1, N = %d  (white = 95%% region)', Nobs));

subplot(1,2,2); hold on; grid on;
for k = 1:nSubj
    if nIn(k) == 0, continue, end
    plot(log10([rLo(k) rHi(k)]), [k k], '-', 'LineWidth', 3, 'Color', [.4 .5 .8]);
    plot(log10(bestR(k)), k, 'r.', 'MarkerSize', 18);
end
xline(log10(theta0(4)), 'k--', 'LineWidth', 2);
xlabel('log_{10} r'); ylabel('subject'); ylim([0 nSubj+1]);
title('Range of r consistent with each subject''s data');
legend('95% interval', 'best fit', 'truth', 'Location','best');

%% ---------------------------------------------------------------- helper
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
