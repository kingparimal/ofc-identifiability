function [Sigma, S] = statsCovariance(theta, N, nRep, seed0, task)
%STATSCOVARIANCE  Sampling covariance of the summary statistics at trial count N.
%
%   [Sigma, S] = STATSCOVARIANCE(theta, N, nRep, seed0, task)
%
%   Repeats the synthetic experiment nRep times with different seeds and
%   returns the 5-by-5 covariance of the resulting summary vectors, plus the
%   raw 5-by-nRep matrix S of replicates.
%
%   This is what turns a raw difference between summary vectors into a
%   calibrated statistical distance. Without it, the discrepancy would weight
%   a 0.001 m difference in mean position the same as a 0.001 m^2 difference
%   in variance -- quantities with different units and wildly different
%   sampling noise.
%
%   REPLICATION, NOT RESAMPLING. Each replicate is an independent run of the
%   forward model, not a bootstrap resample of one dataset. That is more
%   expensive but exact: it measures the true sampling distribution rather
%   than approximating it. Affordable here because the forward model is fast.
%   If it becomes the bottleneck, resampling the N terminal states within a
%   single run is the cheaper approximation.
%
%   COST. nRep solves of the control problem. The gains do not change between
%   replicates (theta is fixed), so this could be sped up by solving once and
%   rolling out nRep times; see statsCovarianceFast below if needed. Kept
%   simple here because correctness at this stage matters more than speed.
%
%   Recommended nRep >= 200. The relative error of an estimated variance from
%   nRep replicates is about sqrt(2/nRep) -- 10% at nRep = 200.

    if nargin < 5 || isempty(task),  task  = 'sim2'; end
    if nargin < 4 || isempty(seed0), seed0 = 10000;  end
    if nargin < 3 || isempty(nRep),  nRep  = 200;    end

    assert(nRep >= 20, 'statsCovariance:nRep', ...
        'nRep = %d is too few to estimate a 5x5 covariance.', nRep);

    % Gains depend only on theta, so solve once and reuse across replicates.
    sigma_u = theta(1);  sigma_s = theta(2);
    d       = round(theta(3));  r = theta(4);

    opts = struct('sigma_u', sigma_u, 'sigma_s', sigma_s, 'r', r);
    M = augmentModel(buildModel2D(task, opts), d);
    [L, K] = solveLQGND(M);

    dp   = M.Dpos(1:2)';
    nhat = dp / norm(dp);
    that = [-nhat(2); nhat(1)];

    S = zeros(5, nRep);
    for k = 1:nRep
        rng(seed0 + k);
        [~, ~, Xf] = rolloutLQG(M, L, K, sigma_u, 0, N);
        P  = Xf(1:2, :);
        pn = nhat' * P;   pt = that' * P;
        C  = cov([pn' pt']);
        S(:,k) = [mean(pn); mean(pt); C(1,1); C(2,2); C(1,2)];
    end

    Sigma = cov(S');

    % ---- conditioning check -------------------------------------------
    % A near-singular Sigma means two summary statistics are (numerically)
    % redundant, and the Mahalanobis distance is not well defined. This is
    % exactly what would happen if aspect ratio or average SD were included
    % alongside the two SDs they are computed from.
    c = cond(Sigma);
    if c > 1e10
        warning('statsCovariance:illConditioned', ...
            ['Sigma has condition number %.2e -- the summary statistics are ' ...
             'close to linearly dependent. Check that no derived quantity ' ...
             'has been added to the vector.'], c);
    end
end
