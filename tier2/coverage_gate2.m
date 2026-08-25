%% coverage_gate2.m -- does the TIER 2 scoring cover at 95%?
%
% Must pass BEFORE the 3-hour tier 2 grid is built. Needs no grid: coverage is
% P(D(theta_true) <= thr), evaluated at the truth point only. ~2 minutes.
%
% PASS BAND 93-97% (3 SE at n = 600). A scoring that does not cover produces
% regions that are WRONG, not merely wide -- and every tier 2 width would
% inherit it.

clear; clc;
theta0=[0.4;0.4;4;0.002]; Nobs=200; task='sim2'; nRep=400; nSubj=600;
tS=[16 23 30 37 44 51];
p = 3*numel(tS);

[Sigma2,~] = statsCovariance2(theta0,Nobs,nRep, 61000,task,tS);   % Sigma
[~,Sp]     = statsCovariance2(theta0,Nobs,nRep, 91000,task,tS);   % s_pred
[~,So]     = statsCovariance2(theta0,Nobs,nSubj,310000,task,tS);  % subjects
s_pred = mean(Sp,2);
thr = hotelling(p,nRep);

D=nan(nSubj,1);
for k=1:nSubj, D(k)=discrepancy(So(:,k), s_pred, Sigma2); end
E = p*(nRep-1)/(nRep-p-1);
fprintf('tier 2: %d statistics at times %s\n', p, mat2str(tS));
fprintf('cond(Sigma) = %.3e   threshold = %.3f\n\n', cond(Sigma2), thr);
fprintf('%-14s %8s %8s %8s %8s\n','scoring','cov95','medD','meanD','E[D]');
fprintf('%-14s %7.1f%% %8.2f %8.2f %8.2f\n','tier 2', ...
        100*mean(D<=thr), median(D), mean(D), E);
se = 100*sqrt(.95*.05/nSubj);
fprintf('\nbinomial SE %.1f%% -> pass band %.1f%% to %.1f%%\n', se, 95-3*se, 95+3*se);
fprintf('tier 1 comparison (same design): 3-stat gave 95.3%%\n');
