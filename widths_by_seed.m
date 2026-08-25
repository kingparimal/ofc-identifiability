%% widths_by_seed.m -- are the INTERVAL WIDTHS robust to grid realisation?
%
% gridseed_test showed the argmin LOCATION is not: the sign-test p for r swung
% from 0.023 to 1.000 on grid seed alone, with subjects held fixed. That kills
% argmin-based bias testing at this N_pred.
%
% The widths are the actual results of the study (RQ1/RQ2), so they need their
% own check. Expectation: a displacement of the prediction surface MOVES the
% accepted region without much changing its SIZE, so log-widths should be far
% more stable across seeds than the counts were.
%
% Reads gridseed.mat, written by gridseed_test. No grid rebuild.
%
% READ IT AS:
%   log-widths stable to a few percent across seeds -> widths are trustworthy,
%       the headline results stand, and the bias question is simply closed as
%       not answerable at affordable N_pred.
%   log-widths swinging like the counts did          -> the pilot's factor 7.2
%       is a property of seed 12345 and the whole width result needs N_pred
%       raised before anything is written down.

clear; clc;
S=load('gridseed.mat');            % Sg, Sbar, Sobs, su_range, r_range, gseeds
theta0=[0.4;0.4;4;0.002]; Nobs=200; task='sim2'; nRep=400;
su_range=S.su_range; r_range=S.r_range;
nSU=numel(su_range); nR=numel(r_range); nSubj=size(S.Sobs,2); nG=numel(S.gseeds);

Sigma5=statsCovariance(theta0,Nobs,nRep,60000,task);
idx3=[1 3 4]; Sigma3=Sigma5(idx3,idx3); thr=hotelling(3,nRep);

fprintf('%-16s %14s %14s %10s %10s\n', ...
        'grid seed','su log-width','r log-width','ratio','open top');
W=nan(nG+1,2);
for g=1:nG+1
    if g<=nG, P=S.Sg(:,:,:,g); tag=sprintf('%d',S.gseeds(g));
    else,     P=S.Sbar;        tag='POOLED'; end
    wS=nan(nSubj,1); wR=nan(nSubj,1); op=0;
    for k=1:nSubj
        D=nan(nSU,nR);
        for i=1:nSU, for j=1:nR
            D(i,j)=discrepancy(S.Sobs(idx3,k),P(idx3,i,j),Sigma3);
        end, end
        in=D<=thr; if ~any(in(:)), continue, end
        [I,J]=find(in);
        wS(k)=log10(su_range(max(I))/su_range(min(I)));
        wR(k)=log10(r_range(max(J))/r_range(min(J)));
        if any(J==nR), op=op+1; end
    end
    W(g,:)=[median(wS(~isnan(wS))) median(wR(~isnan(wR)))];
    fprintf('%-16s %6.4f (x%.3f) %6.4f (x%.3f) %10.2f %10d\n', ...
            tag, W(g,1), 10^W(g,1), W(g,2), 10^W(g,2), W(g,2)/W(g,1), op);
end
fprintf('\n--- stability across seeds (the actual test) ---\n');
fprintf('  sigma_u log-width: spread %.4f, %.1f%% of its median\n', ...
        max(W(1:nG,1))-min(W(1:nG,1)), 100*(max(W(1:nG,1))-min(W(1:nG,1)))/median(W(1:nG,1)));
fprintf('  r       log-width: spread %.4f, %.1f%% of its median\n', ...
        max(W(1:nG,2))-min(W(1:nG,2)), 100*(max(W(1:nG,2))-min(W(1:nG,2)))/median(W(1:nG,2)));
fprintf('  ratio r/su       : %.2f to %.2f   (pilot recorded 8.5)\n', ...
        min(W(1:nG,2)./W(1:nG,1)), max(W(1:nG,2)./W(1:nG,1)));
fprintf('\nCompare: the argmin sign-test p ranged 0.023 to 1.000 on the same data.\n');
