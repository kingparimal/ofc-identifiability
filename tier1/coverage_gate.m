%% coverage_gate.m -- does the 3-statistic scoring actually cover at 95%?
%
% Must pass BEFORE any projection from the full grid is interpreted. Needs no
% grid: coverage is P(D(theta_true) <= thr), evaluated at the truth node only.
%
% Uses statsCovariance to generate subjects because it solves the control
% problem once and rolls out nRep times -- 600 subjects for the cost of one solve.
% All three seed blocks are disjoint, so Sigma, s_pred and the subjects are
% independent.
%
% PASS: 93-97% at nominal 95. Below that, regions under-cover and every width
% in the study is an underestimate. Above, they over-cover and widths are
% conservative.

clear; clc;
theta0=[0.4;0.4;4;0.002]; Nobs=200; task='sim2'; nRep=400; nSubj=600;

[Sigma5,~] = statsCovariance(theta0,Nobs,nRep, 60000,task);   % Sigma
[~,Sp]     = statsCovariance(theta0,Nobs,nRep, 90000,task);   % s_pred
[~,So]     = statsCovariance(theta0,Nobs,nSubj,300000,task);  % subjects
s_pred = mean(Sp,2);

idx3=[1 3 4]; idx2=[2 5];
SigB = Sigma5; SigB(idx3,idx2)=0; SigB(idx2,idx3)=0;   % impose the exact block
Sig3 = Sigma5(idx3,idx3);
thr5 = hotelling(5,nRep); thr3 = hotelling(3,nRep);

D5=nan(nSubj,1); D5b=nan(nSubj,1); D3=nan(nSubj,1);
for k=1:nSubj
    D5(k) =discrepancy(So(:,k),      s_pred,       Sigma5);
    D5b(k)=discrepancy(So(:,k),      s_pred,       SigB);
    D3(k) =discrepancy(So(idx3,k),   s_pred(idx3), Sig3);
end
E5=5*(nRep-1)/(nRep-6); E3=3*(nRep-1)/(nRep-4);
fprintf('%-24s %8s %8s %8s %8s\n','scoring','cov95','medD','meanD','E[D]');
fprintf('%-24s %7.1f%% %8.2f %8.2f %8.2f  [thr %.3f]\n','5-stat, Sigma est', ...
        100*mean(D5<=thr5),  median(D5), mean(D5), E5, thr5);
fprintf('%-24s %7.1f%% %8.2f %8.2f %8.2f  [thr %.3f]\n','5-stat, Sigma block', ...
        100*mean(D5b<=thr5), median(D5b),mean(D5b),E5, thr5);
fprintf('%-24s %7.1f%% %8.2f %8.2f %8.2f  [thr %.3f]\n','3-stat  <-- PRIMARY', ...
        100*mean(D3<=thr3),  median(D3), mean(D3), E3, thr3);
fprintf('\nbinomial SE at p=0.95, n=%d: %.1f%%   -> pass band %.1f%% to %.1f%%\n', ...
        nSubj, 100*sqrt(.95*.05/nSubj), 95-3*100*sqrt(.95*.05/nSubj), ...
        95+3*100*sqrt(.95*.05/nSubj));
fprintf('recorded A4 (5-stat, n=300): 94.0%%\n');
