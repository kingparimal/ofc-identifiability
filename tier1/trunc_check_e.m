%% trunc_check_e.m -- is the PROFILE minimisation truncated by the box?
%
% A profile fixes one parameter and minimises over the other three. If that
% minimisation runs into a wall, the profile is a LOWER BOUND, not a width.
% This is the pass-2 failure mode: a too-narrow axis makes every OTHER profile
% too narrow as well, and the error is invisible because it looks like a number.
%
% Wall contact alone is not proof of distortion -- on pass 4 sigma_u touched for
% 30/30 subjects yet still measured 95% of its box-independent curvature
% prediction, so the minimum was flat there. Read this table TOGETHER with the
% predicted-vs-measured comparison below: truncation matters where wall contact
% and a shortfall against prediction COINCIDE. On pass 4 that was sigma_s (74%).
%
% Run after score_grid4e_t2. Needs grid4e_t2.mat.

G=load('grid4e_t2.mat'); theta0=G.theta0; task=G.task; Nobs=200;
tS=G.tSample; sz=[numel(G.su_range) numel(G.ss_range) numel(G.d_range) numel(G.r_range)]; nSubj=30;
Sig=statsCovariance2(theta0,Nobs,400,60000,task,tS);
thrP=3.841458820694124;
P=reshape(G.Spred,size(G.Spred,1),[]);
nm={'sigma_u','sigma_s','delay','r'};
hitwall=zeros(nSubj,4);
for k=1:nSubj
  s=ofcSummary2(theta0,Nobs,70000+k,task,tS);
  e=P-s; D=reshape(sum(e.*(Sig\e),1),sz); pmin=min(D(:));
  for q=1:4
    others=setdiff(1:4,q); w=0;
    for iq=1:sz(q)
      idx=repmat({':'},1,4); idx{q}=iq;
      sub=D(idx{:});
      if min(sub(:))>pmin+thrP, continue, end
      [~,m]=min(sub(:)); c=cell(1,4); [c{1},c{2},c{3},c{4}]=ind2sub(sz,m);
      for o=others
        if c{o}==1 || c{o}==sz(o), w=w+1; break, end
      end
    end
    hitwall(k,q)=w;
  end
end
fprintf('profile minimiser lands ON A WALL (of the other 3 params):\n');
fprintf('%-10s %28s\n','parameter','subjects with >=1 such slice');
for q=1:4
  fprintf('%-10s %20d/%d\n', nm{q}, sum(hitwall(:,q)>0), nSubj);
end
fprintf('\npredicted vs measured full log-widths:\n');
pr=[0.0480 0.1512 NaN 0.1224]; me=[0.0456 0.1118 NaN 0.1138];
for q=[1 2 4]
  fprintf('  %-9s predicted %.4f  measured %.4f  (%.0f%% of predicted)\n', ...
          nm{q}, pr(q), me(q), 100*me(q)/pr(q));
end
