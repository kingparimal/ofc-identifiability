%% walls.m -- WHICH walls does the region escape through?
%
% escape.m established that the region reaches 3 of 6 walls for every subject.
% Which 3 names the unidentified DIRECTION in the full 4-parameter space --
% the numerical form of the one-dimensional preimage family that sec 1.1
% predicts structurally (the map theta -> E[s] runs R^4 -> R^3).
%
% Run after build_grid4c. Needs grid4c.mat.

clear; clc;
G=load('grid4c.mat'); theta0=G.theta0; task=G.task; Nobs=200; nSubj=30;
sz=[numel(G.su_range) numel(G.ss_range) numel(G.d_range) numel(G.r_range)];
Sigma=statsCovariance(theta0,Nobs,400,60000,task); idx3=[1 3 4];
Sig3=Sigma(idx3,idx3); thr=hotelling(3,400);
P3=reshape(G.Spred(idx3,:,:,:,:),3,[]);
wname={'sigma_u LOW','sigma_u HIGH','sigma_s LOW','sigma_s HIGH', ...
       'delay LOW (d=0)','delay HIGH (d=10)','r LOW','r HIGH'};
hit=zeros(1,8); nface=nan(nSubj,1);
for k=1:nSubj
  s=ofcSummary(theta0,Nobs,70000+k,task);
  e=P3-s(idx3); D=reshape(sum(e.*(Sig3\e),1),sz);
  w={D(1,:,:,:),D(end,:,:,:),D(:,1,:,:),D(:,end,:,:), ...
     D(:,:,1,:),D(:,:,end,:),D(:,:,:,1),D(:,:,:,end)};
  f=0;
  for q=1:8
     t=w{q}; if any(t(:)<=thr), hit(q)=hit(q)+1; f=f+1; end
  end
  nface(k)=f;
end
fprintf('threshold %.3f | all 8 walls, including the two delay faces\n\n', thr);
fprintf('%-20s %10s\n','wall','subjects');
for q=1:8, fprintf('%-20s %7d/%d\n', wname{q}, hit(q), nSubj); end
fprintf('\nwalls reached per subject: median %.1f, range %d..%d\n', ...
        median(nface), min(nface), max(nface));
fprintf(['\nThe consistently-hit set names the escape direction. Read the sign\n' ...
         'against the pairwise correlations: sigma_s/delay -0.65 means low\n' ...
         'sigma_s pairs with long delay along the ridge.\n']);
