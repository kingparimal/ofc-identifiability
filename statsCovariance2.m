function [Sigma, S] = statsCovariance2(theta, N, nRep, seed0, task, tSample)
%STATSCOVARIANCE2  Sampling covariance of the TIER 2 statistics.
%   Same design as statsCovariance: solves the control problem once and rolls
%   out nRep times, so the cost is nRep rollouts and one solve.
%
%   nRep must comfortably exceed the number of statistics. With 18 statistics,
%   nRep = 400 gives a ~6% threshold inflation; nRep = 1000 gives ~2%.
    if nargin<6||isempty(tSample), tSample=[16 23 30 37 44 51]; end
    if nargin<5||isempty(task),  task='sim2'; end
    if nargin<4||isempty(seed0), seed0=10000; end
    if nargin<3||isempty(nRep),  nRep=400;    end

    p = 3*numel(tSample);
    assert(nRep >= 10*p, 'statsCovariance2:nRep', ...
        ['nRep = %d is too few for %d statistics. Use at least %d ' ...
         '(10x), preferably more.'], nRep, p, 10*p);

    opts=struct('sigma_u',theta(1),'sigma_s',theta(2),'r',theta(4));
    M = augmentModel(buildModel2D(task,opts), round(theta(3)));
    [L,K] = solveLQGND(M);
    dp=M.Dpos(1:2)'; nhat=dp/norm(dp); that=[-nhat(2);nhat(1)];
    k=numel(tSample);

    S=zeros(p,nRep);
    for j=1:nRep
        rng(seed0+j);
        [~,Ptraj] = rolloutLQG_traj(M,L,K,theta(1),0,N,tSample);
        mn=zeros(k,1); vn=zeros(k,1); vt=zeros(k,1);
        for q=1:k
            P=Ptraj(:,:,q); pn=nhat'*P; pt=that'*P;
            mn(q)=mean(pn); vn(q)=var(pn); vt(q)=var(pt);
        end
        S(:,j)=[mn;vn;vt];
    end
    Sigma=cov(S');
    c=cond(Sigma);
    if c>1e10
        warning('statsCovariance2:illConditioned', ...
            ['cond(Sigma) = %.2e. Time points are too closely spaced and the ' ...
             'variance profile is nearly redundant across them. Drop points.'], c);
    end
end
