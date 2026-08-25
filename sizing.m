% Local sensitivity at the truth node -> predicted PROFILE widths, tier 1 vs 2.
% J = d E[s] / d theta by central differences on the grid (log coords for the
% multiplicative parameters, so widths come out as log10 widths directly).
% Profile half-width for parameter q = sqrt(chi2inv(0.95,1) * inv(J'*inv(Sig)*J)_qq)
theta0=[0.4;0.4;4;0.002]; Nobs=200; task='sim2';
for tier=1:2
  if tier==1
    G=load('grid4c.mat'); S=statsCovariance(theta0,Nobs,400,60000,task);
    i3=[1 3 4]; S=S(i3,i3); P=G.Spred(i3,:,:,:,:); lab='tier 1 (3 stats)';
  else
    G=load('grid4c_t2.mat'); S=statsCovariance2(theta0,Nobs,400,60000,task,G.tSample);
    P=G.Spred; lab='tier 2 (18 stats)';
  end
  iT=[find(G.su_range==theta0(1)) find(G.ss_range==theta0(2)) ...
      find(G.d_range==theta0(3))  find(G.r_range==theta0(4))];
  ax={log10(G.su_range),log10(G.ss_range),G.d_range,log10(G.r_range)};
  p=size(P,1); J=zeros(p,4);
  for q=1:4
    ip=iT; im=iT; ip(q)=iT(q)+1; im(q)=iT(q)-1;
    sp=P(:,ip(1),ip(2),ip(3),ip(4)); sm=P(:,im(1),im(2),im(3),im(4));
    J(:,q)=(sp(:)-sm(:))/(ax{q}(iT(q)+1)-ax{q}(iT(q)-1));
  end
  F=J'*(S\J); Fi=inv(F); hw=sqrt(3.841458820694124*diag(Fi));
  nm={'sigma_u','sigma_s','delay(steps)','r'};
  fprintf('\n=== %s ===  cond(F) = %.2e\n', lab, cond(F));
  fprintf('%-14s %14s %14s\n','parameter','half-width','full width');
  for q=1:4
    if q==3, fprintf('%-14s %14.3f %12.2f steps\n', nm{q}, hw(q), 2*hw(q));
    else,    fprintf('%-14s %14.4f %12.3fx\n', nm{q}, hw(q), 10^(2*hw(q))); end
  end
end
