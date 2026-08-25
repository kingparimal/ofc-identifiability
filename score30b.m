load('pilot_grid.mat');
Nobs=200; nSubj=30; nRep=400; [nSU,nR]=size(AR);
[~,iTrueSU]=min(abs(su_range-theta0(1)));
[~,iTrueR ]=min(abs(r_range -theta0(4)));
fprintf('truth grid indices: su %d/%d, r %d/%d\n',iTrueSU,nSU,iTrueR,nR);
fprintf('grid step in r = factor %.4f\n\n', r_range(2)/r_range(1));
Sigma5 = statsCovariance(theta0,Nobs,nRep,60000,task);
idx3=[1 3 4]; Sigma3=Sigma5(idx3,idx3);
thr5=hotthr(5,nRep); thr3=hotthr(3,nRep);
lab={'5-stat @ 11.30','3-stat @ 7.92'};
for v=1:2
  if v==1, ii=1:5; Sg=Sigma5; thr=thr5; else ii=idx3; Sg=Sigma3; thr=thr3; end
  biS=nan(nSubj,1); bjS=nan(nSubj,1); wSU=nan(nSubj,1); wR=nan(nSubj,1);
  nIn=zeros(nSubj,1); opT=false(nSubj,1);
  for k=1:nSubj
    s_obs=ofcSummary(theta0,Nobs,70000+k,task);
    D=nan(nSU,nR);
    for i=1:nSU, for j=1:nR, D(i,j)=discrepancy(s_obs(ii),Spred(ii,i,j),Sg); end, end
    if k<=3, Dc{v,k}=D; end
    [~,m]=min(D(:)); [bi,bj]=ind2sub(size(D),m); biS(k)=bi; bjS(k)=bj;
    inside=D<=thr; nIn(k)=sum(inside(:)); if nIn(k)==0, continue, end
    [I,J]=find(inside);
    wSU(k)=log10(su_range(max(I))/su_range(min(I)));
    wR(k) =log10(r_range(max(J))/r_range(min(J)));
    opT(k)=any(J==nR);
  end
  Rv(v)=struct('biS',biS,'bjS',bjS,'wSU',wSU,'wR',wR,'nIn',nIn,'opT',opT);
end
save('score30b.mat', 'Rv', 'Dc', 'Sigma5', 'thr5', 'thr3', 'iTrueR', 'iTrueSU', '-v7');
for v=1:2
  R=Rv(v);
  a=sum(R.bjS>iTrueR); e=sum(R.bjS==iTrueR); b=sum(R.bjS<iTrueR);
  as=sum(R.biS>iTrueSU); es=sum(R.biS==iTrueSU); bs=sum(R.biS<iTrueSU);
  n=a+b; p=2*min(binocdf_l(min(a,b),n,0.5),0.5);
  fprintf('=== %s ===\n',lab{v});
  fprintf('  empty %d/30   open-at-top %d/30\n', sum(R.nIn==0), sum(R.opT));
  fprintf('  log-width su %.4f (factor %.3f) | r %.4f (factor %.3f) | ratio %.2f\n', ...
     nanmed(R.wSU),10^nanmed(R.wSU),nanmed(R.wR),10^nanmed(R.wR),nanmed(R.wR)/nanmed(R.wSU));
  fprintf('  best-fit r  index: %d above / %d AT / %d below truth   sign-test p = %.4f (n=%d)\n',a,e,b,p,n);
  fprintf('  best-fit su index: %d above / %d AT / %d below   [CONTROL]\n',as,es,bs);
  fprintf('  median best-fit r = %.3e  (truth %.3e)\n\n', r_range(round(nanmed(R.bjS))), theta0(4));
end
fprintf('=== GATE: subjects 1-8, 5-stat, vs recorded pilot ===\n');
R=Rv(1);
fprintf('  su factor %.3f (rec 1.26) | r factor %.3f (rec 7.20) | ratio %.2f (rec 8.5)\n', ...
   10^nanmed(R.wSU(1:8)),10^nanmed(R.wR(1:8)),nanmed(R.wR(1:8))/nanmed(R.wSU(1:8)));
fprintf('  best-fit su %.3f +/- %.3f (rec 0.398 +/- 0.021)\n', ...
   mean(su_range(R.biS(1:8))), std(su_range(R.biS(1:8))));
fprintf('  median r %.2e (rec 2.78e-03) | r above truth %d of 8 (rec 6) | closed %d of 8 (rec 7)\n', ...
   r_range(round(median(R.bjS(1:8)))), sum(R.bjS(1:8)>iTrueR), sum(~R.opT(1:8)));
fprintf('  same 8 subjects, 3-stat: su factor %.3f | r factor %.3f\n', ...
   10^nanmed(Rv(2).wSU(1:8)), 10^nanmed(Rv(2).wR(1:8)));
fprintf('\n=== D5 - D3 across grid (predicted constant in theta) ===\n');
for k=1:3
  dd=Dc{1,k}(:)-Dc{2,k}(:);
  fprintf('  subj %d: mean %.3f  sd %.4f  range %.4f  (as %%%% of mean: %.2f%%%%)\n', ...
     k, mean(dd), std(dd), max(dd)-min(dd), 100*(max(dd)-min(dd))/mean(dd));
end
