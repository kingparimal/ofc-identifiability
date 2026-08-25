%% escape.m -- does the 95%% region CLOSE inside the swept box, or escape it?
%
% Run after build_grid4c. Needs grid4c.mat, plus statsCovariance, ofcSummary,
% discrepancy and hotelling on the path. Takes a couple of minutes.
%
% WHY. r's width read x8.0 on pass 1, x22.6 on pass 2 and x129.8 on pass 3 --
% growing every time the grid grew. That is the signature of a region that is
% not contained by the box, in which case each "width" is a measurement of the
% analyst's sweep rather than of the model (sec 5.1, one level up: not the
% sweep window but the whole 4-D box).
%
% WHAT IT MEASURES. The minimum of D on the outer shell of the grid. If the
% region merely grazes a wall, that minimum sits AT the threshold. If the
% region is still descending as it exits, the minimum sits well BELOW it, and
% no finite widening will close it.
%
% READ IT AS:
%   min D on shell ~ threshold        -> region nearly contained; a wider box
%                                        would close it, and is worth building
%   min D on shell << threshold       -> region escapes; STOP BUILDING GRIDS.
%                                        Report every width as a lower bound
%                                        with the box stated.
%
% Also reports the region volume as a fraction of the box, and how many of the
% 6 walls it reaches. A thin sheet reaching several walls at once is a
% correlated escape along a ridge, not independent clipping of separate axes.

% Does the 95% region CLOSE inside the swept box, or escape through its walls?
% If the minimum of D on the outer shell is well BELOW threshold, the region is
% still descending as it exits -- no finite widening will close it.
G=load('grid4c.mat'); theta0=G.theta0; task=G.task; Nobs=200; nSubj=30;
sz=[numel(G.su_range) numel(G.ss_range) numel(G.d_range) numel(G.r_range)];
Sigma=statsCovariance(theta0,Nobs,400,60000,task); idx3=[1 3 4]; Sig3=Sigma(idx3,idx3);
thr=hotelling(3,400);
P3=reshape(G.Spred(idx3,:,:,:,:),3,[]);
[A,B,C,E]=ndgrid(1:sz(1),1:sz(2),1:sz(3),1:sz(4));
shell = A==1|A==sz(1)|B==1|B==sz(2)|E==1|E==sz(4);   % faces in su, ss, r
fprintf('threshold %.3f | shell is %.1f%% of the grid\n\n', thr, 100*mean(shell(:)));
mn=nan(nSubj,1); vol=nan(nSubj,1); nface=nan(nSubj,1);
for k=1:nSubj
  s=ofcSummary(theta0,Nobs,70000+k,task);
  e=P3-s(idx3); D=reshape(sum(e.*(Sig3\e),1),sz);
  mn(k)=min(D(shell)); vol(k)=100*mean(D(:)<=thr);
  w={D(1,:,:,:),D(end,:,:,:),D(:,1,:,:),D(:,end,:,:),D(:,:,:,1),D(:,:,:,end)};
  f=0; for q=1:6, t=w{q}; if any(t(:)<=thr), f=f+1; end, end
  nface(k)=f;
end
fprintf('min D on the outer shell, across 30 subjects:\n');
fprintf('  median %.3f | range %.3f .. %.3f | below threshold for %d/30\n', ...
        median(mn), min(mn), max(mn), sum(mn<=thr));
fprintf('  median as %% of threshold: %.0f%%\n', 100*median(mn)/thr);
fprintf('\nregion volume as %% of the swept box: median %.2f%% (range %.2f .. %.2f)\n', ...
        median(vol), min(vol), max(vol));
fprintf('number of the 6 walls the region reaches: median %.1f, max %d\n', ...
        median(nface), max(nface));
