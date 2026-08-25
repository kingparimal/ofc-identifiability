%% nsweep_t2.m -- TIER 2 profile widths vs trial count N.
%
% Produces the sample-size recommendation: how many trials are needed to
% recover each parameter to a stated precision, now that all four ARE
% recoverable (tier 1 could not bound any of them).
%
% No grid rebuild. Predicted statistics do not depend on the observation trial
% count; only s_obs and Sigma change with N.
%
% WHICH GRID. Widths scale as 1/sqrt(N), so no single grid spans the range:
%   N <= 100  -> grid4c_t2.mat  (coarse; regions are large and would run off
%                                the fine grid)
%   N >= 200  -> grid4e_t2.mat  (fine; sized for N = 200)
% The switch is a property of the grids, not the model, and is reported.
%
% FLAGS, printed per cell:
%   OPEN    region touches a boundary -> width is a LOWER BOUND
%   COARSE  fewer than 3 cells across -> width is at the grid resolution floor
%
% EXPECT at the ends of the range: N = 1000 and 5000 will read COARSE on the
% fine grid (at N = 5000 sigma_u's width is ~1 cell). Report those as bounds
% and extrapolate using the 1/sqrt(N) law rather than building another grid --
% tier 1's sigma_u followed that law to three figures (0.3876, 0.2907, 0.1938
% at N = 50, 100, 200; ratio 0.500 against a predicted 0.500).
%
% The delay is INTEGER, quantised at dt = 10 ms. On grid4c_t2 it is sampled
% every 2 steps, so its resolution there is 2 steps. It cannot narrow below
% 1 step at any N: that is physical, not a grid limit.

clear; clc;
theta0=[0.4;0.4;4;0.002]; task='sim2'; nSubj=30; nRep=400;
Nlist   = [50 100 200 500 1000 5000];
useFine = Nlist >= 200;
thrP = 3.841458820694124;                 % chi2inv(0.95,1)
axes_n = {'sigma_u','sigma_s','delay(1+d)','r'};

Gc = load('grid4c_t2.mat');  Gf = load('grid4e_t2.mat');
tS = Gf.tSample;
thrR = hotelling(size(Gf.Spred,1), nRep);
fprintf('region thr %.3f | profile thr %.3f | %d subjects | %d statistics\n\n', ...
        thrR, thrP, nSubj, size(Gf.Spred,1));
fprintf('%6s %-7s %-12s %14s %9s %6s %s\n', ...
        'N','grid','parameter','95%% profile','logwidth','cells','flag');

for a=1:numel(Nlist)
    N = Nlist(a);
    if useFine(a), G=Gf; gname='fine'; else, G=Gc; gname='coarse'; end
    sz=[numel(G.su_range) numel(G.ss_range) numel(G.d_range) numel(G.r_range)];
    av={G.su_range, G.ss_range, 1+G.d_range, G.r_range};

    Sigma = statsCovariance2(theta0,N,nRep,60000,task,tS);
    % Subjects are generated with ofcSummary2 directly, NOT statsCovariance2:
    % that function asserts nRep >= 10*p for a sound covariance estimate, which
    % is correct for Sigma but irrelevant for drawing 30 subjects. This also
    % guarantees the subjects are bitwise the ones score_grid4e_t2 uses
    % (seeds 70001..70030), so the N = 200 continuity check is exact.
    So = zeros(size(Gf.Spred,1), nSubj);
    for k=1:nSubj
        So(:,k) = ofcSummary2(theta0, N, 70000+k, task, tS);
    end
    P = reshape(G.Spred, size(G.Spred,1), []);

    plo=nan(nSubj,4); phi=nan(nSubj,4); pop=false(nSubj,4); ncell=nan(nSubj,4);
    nEmpty=0;
    for k=1:nSubj
        e = P - So(:,k);
        D = reshape(sum(e.*(Sigma\e),1), sz);
        if ~any(D(:)<=thrR), nEmpty=nEmpty+1; end
        pmin=min(D(:));
        for q=1:4
            Dp=D; for o=sort(setdiff(1:4,q),'descend'), Dp=min(Dp,[],o); end
            ok=find(squeeze(Dp)<=pmin+thrP);
            plo(k,q)=av{q}(min(ok)); phi(k,q)=av{q}(max(ok));
            pop(k,q)=(min(ok)==1)|(max(ok)==sz(q));
            ncell(k,q)=max(ok)-min(ok)+1;
        end
    end
    for q=1:4
        w = log10(phi(:,q)./plo(:,q));  w=median(w(~isnan(w)));
        cel = median(ncell(:,q));  op = sum(pop(:,q));
        flag='';
        if op > nSubj/2, flag=[flag 'OPEN ']; end
        if cel < 3,      flag=[flag 'COARSE']; end
        fprintf('%6d %-7s %-12s %6.4g..%-7.4g %9.4f %6.0f %s\n', ...
                N, gname, axes_n{q}, median(plo(~isnan(plo(:,q)),q)), ...
                median(phi(~isnan(phi(:,q)),q)), w, cel, flag);
    end
    if nEmpty>0, fprintf('%6s (%d empty joint regions)\n','',nEmpty); end
    fprintf('\n');
end
fprintf(['Expected: log-widths halve per 4x in N, except where flagged.\n\n' ...
         'CONTINUITY CHECK. The N = 200 row should sit within ONE GRID CELL of\n' ...
         'the pass-5 scoring (sigma_u 0.0456, sigma_s 0.1231, delay 0.0792,\n' ...
         'r 0.0910). It will not match exactly: Sigma is re-estimated here from\n' ...
         'a different seed block, and one cell is 0.0091 (sigma_u), 0.0274\n' ...
         '(sigma_s), 1 step (delay) and 0.0228 (r). Larger drift than that\n' ...
         'means something other than the Sigma draw has changed.\n\n' ...
         'The delay reads COARSE at every N: its region is ~2 cells and cannot\n' ...
         'narrow below 1 step, because d is integer at dt = 10 ms. That is\n' ...
         'physical quantisation, not a grid limit. Report it as +/- one step.\n']);
