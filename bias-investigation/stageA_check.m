%% stageA_check.m  --  Stage A gate
%
% Four tests, each able to fail. Nothing in Stage B is meaningful until all
% four pass, because a miscalibrated discrepancy produces confidence regions
% of the wrong size and a false identifiability map.
%
%   A1  reproducibility and the sensitivity baseline
%   A2  the summary vector is not rank-deficient
%   A3  sampling covariance scales as 1/N
%   A4  COVERAGE -- the decisive test
%
% Runtime: a few minutes, dominated by A4.

clear;  clc;
task  = 'sim2';
theta0 = [0.4; 0.4; 4; 0.002];      % nominal: sigma_u, sigma_s, d, r

fprintf('=== Stage A gate ===\n');
fprintf('    task %s, nominal theta = [%.2f %.2f %d %.4f]\n\n', ...
        task, theta0(1), theta0(2), theta0(3), theta0(4));

%% ---------------------------------------------------------------- A1
% Reproducibility, and agreement with the recorded Stage 4 / sec 2.4 values.
% The d = 4 Sim 2 row of REFERENCE_NUMBERS sec 2.4 gives SD_rel 0.4471 cm,
% SD_irrel 1.4965 cm, AR 3.347. Those came from a different random stream, so
% agreement is expected within Monte Carlo error (about +/-2% at N = 5000),
% not exactly.
fprintf('--- A1  reproducibility and baseline ---\n');

[sA, iA] = ofcSummary(theta0, 5000, 1, task);
[sB, iB] = ofcSummary(theta0, 5000, 1, task);
same = isequal(sA, sB);
fprintf('  same seed reproduces exactly       : %s\n', vd(same));

[~, iC] = ofcSummary(theta0, 5000, 2, task);
fprintf('  seed 1  SD_rel %.4f cm  SD_irrel %.4f cm  AR %.3f\n', ...
        iA.sd_rel*100, iA.sd_irrel*100, iA.AR);
fprintf('  seed 2  SD_rel %.4f cm  SD_irrel %.4f cm  AR %.3f\n', ...
        iC.sd_rel*100, iC.sd_irrel*100, iC.AR);
fprintf('  recorded (sec 2.4)  0.4471          1.4965            3.347\n');

relAR = abs(iA.AR - 3.347) / 3.347;
ok1 = same && relAR < 0.05;
fprintf('  AR within 5%% of recorded            : %.1f%%   %s\n\n', ...
        100*relAR, vd(ok1));

%% ---------------------------------------------------------------- A2
% The summary vector must not contain derived quantities. Five independent
% statistics -> Sigma full rank. If aspect ratio or average SD were included
% alongside the SDs they are computed from, Sigma would be singular.
fprintf('--- A2  summary vector is full rank ---\n');

[Sig, S] = statsCovariance(theta0, 500, 200, 10000, task);
c = cond(Sig);
ok2 = c < 1e10;
fprintf('  cond(Sigma)                        : %.3e   %s\n', c, vd(ok2));

% Demonstrate the failure mode deliberately: append AR and avgSD.
sd_n = sqrt(S(3,:));  sd_t = sqrt(S(4,:));
Sbad = [S; sd_t./sd_n; (sd_n+sd_t)/2];
cbad = cond(cov(Sbad'));
fprintf('  cond with AR and avgSD appended    : %.3e  <- singular, as expected\n', cbad);
fprintf('  (this is why the reported five statistics are only three)\n\n');

%% ---------------------------------------------------------------- A3
% Sampling covariance should scale as 1/N. Doubling N should halve the
% variance of each summary statistic. A departure means the replicates are
% not independent, or N is not doing what it should.
fprintf('--- A3  sampling covariance scales as 1/N ---\n');

Ns = [250 500 1000];
v  = zeros(numel(Ns), 5);
for k = 1:numel(Ns)
    Sk = statsCovariance(theta0, Ns(k), 200, 20000 + 1000*k, task);
    v(k,:) = diag(Sk)';
end
ratio = v(1,:) ./ v(3,:);          % expect ~4 (N goes 250 -> 1000)
fprintf('  var(N=250)/var(N=1000), expect ~4.0:\n');
lbl = {'mean_n','mean_t','var_n','var_t','cov_nt'};
for j = 1:5
    fprintf('        %-8s %.2f\n', lbl{j}, ratio(j));
end
ok3 = all(ratio > 2.5 & ratio < 6);
fprintf('  all in [2.5, 6]                    : %s\n\n', vd(ok3));

%% ---------------------------------------------------------------- A4
% COVERAGE. The decisive test, and the one that can genuinely fail.
%
% Generate many synthetic subjects at the SAME known theta. For each, compute
% the discrepancy of its summary vector against the expected value. If the
% Gaussian approximation holds and Sigma is right, those discrepancies follow
% chi2(5), so about 95% fall below chi2inv(0.95,5) = 11.07.
%
% If coverage comes out at 80%, every confidence region in Stage B is too
% small and the identifiability map is wrong in the optimistic direction --
% it would claim parameters are recoverable when they are not.
fprintf('--- A4  COVERAGE (the decisive test) ---\n');

Ncov  = 200;                       % realistic trial count
nSub  = 300;                       % synthetic subjects
fprintf('  N = %d, %d synthetic subjects\n', Ncov, nSub);

% Expected summary vector: high-N run at the true theta.
s_true = ofcSummary(theta0, 200000, 999, task);

% Sampling covariance at the realistic N.
Sig_cov = statsCovariance(theta0, Ncov, 400, 30000, task);

D = zeros(nSub,1);
for k = 1:nSub
    s_k  = ofcSummary(theta0, Ncov, 40000 + k, task);
    D(k) = discrepancy(s_k, s_true, Sig_cov);
end

thr95 = chi2inv_local(0.95, 5);
thr68 = chi2inv_local(0.68, 5);
cov95 = mean(D <= thr95);
cov68 = mean(D <= thr68);

fprintf('  nominal 95%% -> actual %.1f%%   (threshold %.2f)\n', 100*cov95, thr95);
fprintf('  nominal 68%% -> actual %.1f%%   (threshold %.2f)\n', 100*cov68, thr68);
fprintf('  median D %.2f (expect ~%.2f),  mean D %.2f (expect 5)\n', ...
        median(D), chi2inv_local(0.5,5), mean(D));

ok4 = abs(cov95 - 0.95) < 0.05 && abs(cov68 - 0.68) < 0.08;
fprintf('  coverage within tolerance          : %s\n\n', vd(ok4));

figure('Color','w','Position',[100 100 520 380]);
histogram(D, 30, 'Normalization','pdf', 'FaceColor',[.6 .6 .65]); hold on; grid on;
xg = linspace(0, max(D), 300);
plot(xg, chi2pdf_local(xg,5), 'r-', 'LineWidth', 2);
xline(thr95, 'k--', 'LineWidth', 1.5);
xlabel('discrepancy D'); ylabel('density');
legend('observed', '\chi^2_5', '95% threshold', 'Location','best');
title(sprintf('Coverage check: %.1f%% below the 95%% threshold', 100*cov95));

%% ----------------------------------------------------------------
allok = ok1 && ok2 && ok3 && ok4;
fprintf('=== STAGE A: %s ===\n\n', vd(allok));
if allok
    fprintf(['Machinery is calibrated. Proceed to the two-parameter pilot\n' ...
             '(r and sigma_u, 9x9 grid, others fixed) before the full grid.\n\n']);
else
    fprintf(2, ['Do NOT proceed to Stage B. A miscalibrated discrepancy\n' ...
                'produces confidence regions of the wrong size, and the\n' ...
                'identifiability map would be wrong in the optimistic\n' ...
                'direction -- claiming recoverability that is not there.\n\n']);
end

%% ---------------------------------------------------------------- helpers
function s = vd(ok)
    if ok, s = 'PASS'; else, s = 'FAIL'; end
end

function x = chi2inv_local(p, k)
% Avoids a Statistics Toolbox dependency. chi2 with k dof is gamma(k/2, 2).
    if exist('chi2inv','file') == 2
        x = chi2inv(p, k);  return
    end
    x = gaminv_local(p, k/2, 2);
end

function x = gaminv_local(p, a, b)
    lo = 0;  hi = 1;
    while gamcdf_local(hi, a, b) < p, hi = hi * 2; end
    for i = 1:200
        mid = (lo + hi)/2;
        if gamcdf_local(mid, a, b) < p, lo = mid; else, hi = mid; end
    end
    x = (lo + hi)/2;
end

function y = gamcdf_local(x, a, b)
    y = gammainc(x/b, a);
end

function y = chi2pdf_local(x, k)
    if exist('chi2pdf','file') == 2
        y = chi2pdf(x, k);  return
    end
    y = (x.^(k/2-1) .* exp(-x/2)) / (2^(k/2) * gamma(k/2));
end
