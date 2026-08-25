function [s, info] = ofcSummary2(theta, N, seed, task, tSample)
%OFCSUMMARY2  TIER 2 summary statistics: variance profile along the trajectory.
%
%   [s, info] = OFCSUMMARY2(theta, N, seed, task, tSample)
%
%   theta = [sigma_u; sigma_s; d; r]   as in ofcSummary.
%   tSample : state indices to score (default [16 23 30 37 44 51], i.e.
%             150 200 290 360 430 500 ms at dt = 10 ms).
%
%   Returns, at each sampled time, THREE statistics in task coordinates:
%
%       s = [ mean_n(t1..tk) ; var_n(t1..tk) ; var_t(t1..tk) ]
%
%   so numel(s) = 3*numel(tSample).
%
%   WHY ONLY THREE PER TIME POINT. The axis-swap symmetry of Sim 2 is a
%   property of the PROCESS, not just of the terminal state, so mean_t and
%   cov_nt are exactly zero in expectation at EVERY time step. Including them
%   would add 2k statistics of pure noise and inflate the threshold for
%   nothing -- the same trap tier 1 had with the reported five. gate_tier2.m
%   verifies this along the trajectory rather than assuming it.
%
%   WHY MEAN_N IS INCLUDED. The mean trajectory is where the effort weight r
%   shows itself most directly: a low r gives an aggressive controller that
%   travels fast and decelerates hard, a high r a sluggish one. Endpoints see
%   only the destination. This is the statistic tier 1 could not have.
%
%   SIZE. Six time points gives 18 statistics. At nRep = 400 the Hotelling
%   threshold is ~6% above the chi-square ideal. Do not push much past ~20
%   without raising nRep: at 101 statistics the penalty is 39% and you lose
%   more to estimating Sigma than you gain in information.

    if nargin < 5 || isempty(tSample), tSample = [16 23 30 37 44 51]; end
    if nargin < 4 || isempty(task),    task = 'sim2'; end
    if nargin < 3 || isempty(seed),    seed = 0;      end

    validateattributes(theta, {'numeric'}, {'vector','numel',4,'finite'}, mfilename, 'theta');
    theta = theta(:);
    sigma_u=theta(1); sigma_s=theta(2); d=theta(3); r=theta(4);
    assert(abs(d-round(d))<1e-9, 'ofcSummary2:delayInteger', 'd must be an integer.');
    d = round(d);
    assert(d>=0 && sigma_u>0 && sigma_s>0 && r>0, 'ofcSummary2:positive', ...
           'd >= 0 and sigma_u, sigma_s, r > 0 required.');

    opts = struct('sigma_u',sigma_u,'sigma_s',sigma_s,'r',r);
    M = augmentModel(buildModel2D(task,opts), d);
    assert(max(tSample) <= M.n, 'ofcSummary2:tSample', ...
           'tSample must not exceed n = %d.', M.n);
    [L,K] = solveLQGND(M);

    rng(seed);
    [~, Ptraj] = rolloutLQG_traj(M, L, K, sigma_u, 0, N, tSample);

    dp   = M.Dpos(1:2)';
    nhat = dp/norm(dp);
    that = [-nhat(2); nhat(1)];

    k = numel(tSample);
    mn=zeros(k,1); vn=zeros(k,1); vt=zeros(k,1);
    for q = 1:k
        P  = Ptraj(:,:,q);
        pn = nhat'*P;  pt = that'*P;
        mn(q)=mean(pn); vn(q)=var(pn); vt(q)=var(pt);
    end
    s = [mn; vn; vt];

    if nargout > 1
        info = struct('tSample',tSample,'mean_n',mn,'var_n',vn,'var_t',vt, ...
                      'sd_rel',sqrt(vn),'sd_irrel',sqrt(vt), ...
                      'AR',sqrt(vt./vn),'N',N,'theta',theta,'task',task);
    end
end
