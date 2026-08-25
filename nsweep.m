%% nsweep.m -- profile widths as a function of trial count N.
%
% Converts the identifiability result into a SAMPLE-SIZE RECOMMENDATION (sec 6
% item 2). No grid rebuild: the predicted statistics do not depend on the
% observation trial count. Only s_obs and Sigma change with N.
%
% WHICH GRID. Regions shrink as 1/sqrt(N), so no single grid serves the whole
% sweep:
%   N <= 200  -> grid4.mat  (pass 1, wide)  regions are large and would run off
%                                            the fine grid
%   N >= 500  -> grid4b.mat (pass 2, fine)  regions are small and would sit
%                                            inside one cell of the coarse grid
% The switch is a property of the grids, not of the model, and is reported.
%
% RESOLUTION FLAGS, printed per cell, because a width is only meaningful if the
% region is both closed and several cells across:
%   OPEN    region touches a sweep boundary -> width is a LOWER BOUND (sec 5.1)
%   COARSE  fewer than 3 cells across       -> width is at the resolution floor
%
% Subjects are generated through statsCovariance, which solves once and rolls
% out nSubj times. With seed0 = 70000 this reproduces the ofcSummary(70000+k)
% subjects used everywhere else, at a fraction of the cost.

clear; clc;
theta0=[0.4;0.4;4;0.002]; task='sim2'; nSubj=30; nRep=400;
Nlist   = [50 100 200 500 1000 5000];
useFine = Nlist >= 500;
thrR = hotelling(3,nRep);
thrP = 3.841458820694124;                 % chi2inv(0.95,1)
idx3 = [1 3 4];
axes_n = {'sigma_u','sigma_s','delay(1+d)','r'};

Gc = load('grid4.mat');  Gf = load('grid4b.mat');
fprintf('region thr %.3f | profile thr %.3f | %d subjects\n\n', thrR, thrP, nSubj);
fprintf('%6s %-6s %-12s %12s %8s %7s %s\n', ...
        'N','grid','parameter','95%% profile','logwidth','cells','flag');

for a=1:numel(Nlist)
    N = Nlist(a);
    if useFine(a), G=Gf; gname='fine'; else, G=Gc; gname='wide'; end
    sz=[numel(G.su_range) numel(G.ss_range) numel(G.d_range) numel(G.r_range)];
    av={G.su_range, G.ss_range, 1+G.d_range, G.r_range};
    % log step per axis (delay is integer -> count nodes instead)
    st=[log10(G.su_range(2)/G.su_range(1)), log10(G.ss_range(2)/G.ss_range(1)), ...
        NaN, log10(G.r_range(2)/G.r_range(1))];

    Sigma = statsCovariance(theta0,N,nRep,60000,task);
    Sig3  = Sigma(idx3,idx3);
    [~,So] = statsCovariance(theta0,N,nSubj,70000,task);

    P3 = reshape(G.Spred(idx3,:,:,:,:),3,[]);
    plo=nan(nSubj,4); phi=nan(nSubj,4); pop=false(nSubj,4); ncell=nan(nSubj,4);
    nEmpty=0;
    for k=1:nSubj
        e = P3 - So(idx3,k);
        D = reshape(sum(e.*(Sig3\e),1), sz);
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
        cel = median(ncell(:,q));
        op  = sum(pop(:,q));
        flag='';
        if op > nSubj/2,  flag=[flag 'OPEN ']; end
        if cel < 3,       flag=[flag 'COARSE']; end
        fprintf('%6d %-6s %-12s %6.4g..%-6.4g %8.4f %7.0f %s\n', ...
                N, gname, axes_n{q}, median(plo(~isnan(plo(:,q)),q)), ...
                median(phi(~isnan(phi(:,q)),q)), w, cel, flag);
    end
    if nEmpty>0, fprintf('%6s %s\n','', sprintf('(%d empty joint regions)',nEmpty)); end
    fprintf('\n');
end
fprintf(['Expected: log-widths fall as 1/sqrt(N), i.e. halving per 4x in N,\n' ...
         'EXCEPT where flagged OPEN (bounded by the sweep) or COARSE (bounded\n' ...
         'by the grid). Delay should stay saturated at every N -- if it ever\n' ...
         'closes, that N is the answer to "how many trials to measure delay".\n']);
