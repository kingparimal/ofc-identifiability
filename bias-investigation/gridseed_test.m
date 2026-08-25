%% gridseed_test.m -- is the apparent bias in the estimator, or in the grid?
%
% stageB_bias30 gave a CONTROL (sigma_u) more significant than the effect (r),
% in the opposite direction. The suspect is the prediction grid: every point is
% computed with one fixed seed (12345) at N_pred = 20000, so the whole surface
% carries a single Monte Carlo displacement shared by all subjects. It does not
% shrink as subjects are added.
%
% This rebuilds the grid under several seeds and re-scores the SAME 30 subjects
% each time. Subject noise is therefore identical across conditions and the only
% thing varying is the grid realisation -- a paired design.
%
% READING IT:
%   counts move around / control changes sign  -> the grid, not the estimator.
%                                                 Bias untestable at this N_pred.
%   counts stable across all seeds             -> a real effect. Report it.
%
% The pooled row scores against the mean of the grids, which has sqrt(nSeed)
% less displacement, and is the best available estimate of the truth.

clear; clc;
task='sim2'; theta0=[0.4;0.4;4;0.002];
Nobs=200; Npred=20000; nSubj=30; nRep=400;
gseeds = [12345 23456 34567 45678];

nSU=15; su_range=linspace(0.30,0.50,nSU);
nR =15; r_range =logspace(log10(2e-4),log10(2e-2),nR);
[~,iTS]=min(abs(su_range-theta0(1)));
[~,iTR]=min(abs(r_range -theta0(4)));

Sigma5=statsCovariance(theta0,Nobs,nRep,60000,task);
idx3=[1 3 4]; Sigma3=Sigma5(idx3,idx3); thr3=hotelling(3,nRep);

% subjects computed ONCE and reused, so the pairing is exact
Sobs=zeros(5,nSubj);
for k=1:nSubj, Sobs(:,k)=ofcSummary(theta0,Nobs,70000+k,task); end

nG=numel(gseeds);
Sg=nan(5,nSU,nR,nG);
t0=tic;
for g=1:nG
    for i=1:nSU
        for j=1:nR
            th=[su_range(i); theta0(2); theta0(3); r_range(j)];
            Sg(:,i,j,g)=ofcSummary(th,Npred,gseeds(g),task);
        end
    end
    fprintf('  grid seed %d built (%.0f s)\n', gseeds(g), toc(t0));
end
Sbar=mean(Sg,4);
save('gridseed.mat','Sg','Sbar','Sobs','su_range','r_range','gseeds','-v7');

fprintf('\n%-16s %-28s %-28s\n','grid seed','r  (above/at/below, p)','sigma_u CONTROL');
for g=1:nG+1
    if g<=nG, S=Sg(:,:,:,g); tag=sprintf('%d',gseeds(g));
    else,     S=Sbar;        tag='POOLED (mean)'; end
    bi=nan(nSubj,1); bj=nan(nSubj,1);
    for k=1:nSubj
        D=nan(nSU,nR);
        for i=1:nSU, for j=1:nR
            D(i,j)=discrepancy(Sobs(idx3,k),S(idx3,i,j),Sigma3);
        end, end
        [~,m]=min(D(:)); [bi(k),bj(k)]=ind2sub(size(D),m);
    end
    fprintf('%-16s %-28s %-28s\n', tag, fmt(bj,iTR), fmt(bi,iTS));
end

% how big is the grid displacement, in units of the 30-subject standard error?
fprintf('\n--- grid displacement vs subject noise ---\n');
sdSub = std(Sobs,0,2)/sqrt(nSubj);           % SEM across 30 subjects
sdGrid= std(reshape(Sg(:,iTS,iTR,:),5,nG),0,2); % spread of the truth node over seeds
nm={'mean_n','mean_t','var_n','var_t','cov_nt'};
for i=[1 3 4]
    fprintf('  %-7s grid spread %.3e   subject SEM %.3e   ratio %.2f\n', ...
            nm{i}, sdGrid(i), sdSub(i), sdGrid(i)/sdSub(i));
end
