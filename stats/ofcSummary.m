function [s, info] = ofcSummary(theta, N, seed, task)
%OFCSUMMARY  Endpoint summary statistics for a parameter vector.
%
%   [s, info] = OFCSUMMARY(theta, N, seed, task)
%
%       theta = [sigma_u; sigma_s; d; r]
%           sigma_u  control-dependent (motor) noise scale
%           sigma_s  sensory noise scale
%           d        extra sensory delay in steps (INTEGER, 10 ms each)
%           r        effort weight
%       N     trials
%       seed  RNG seed (reproducibility)
%       task  'sim1' or 'sim2'   (default 'sim2')
%
%   Returns s, a 5-vector of summary statistics in TASK COORDINATES:
%
%       s = [ mean_n ; mean_t ; var_n ; var_t ; cov_nt ]
%
%   where n is the task-relevant direction (gradient of the task error in
%   position space) and t the task-irrelevant direction orthogonal to it.
%
%   WHY FIVE, AND WHY THESE FIVE.
%   The endpoint cloud is 2-D. Treated as Gaussian it is fully described by
%   its mean (2 numbers) and covariance (3 numbers) -- five in total, and no
%   more. Any further statistic computed from the cloud is a function of
%   these and carries no additional information.
%
%   In particular the set usually reported in this literature --
%   bias, SD_rel, SD_irrel, aspect ratio, average SD -- looks like five
%   numbers but is only THREE:
%
%       aspect ratio = SD_irrel / SD_rel        (derived)
%       average SD   = (SD_rel + SD_irrel)/2    (derived)
%
%   This matters twice over. Practically, a covariance matrix built from all
%   five is singular, and the Mahalanobis discrepancy needs its inverse.
%   Substantively, three statistics cannot determine four parameters: there
%   is generically a one-dimensional family of theta producing identical
%   expected statistics. That is STRUCTURAL under-identification, present at
%   infinite N, and independent of any sampling argument.
%
%   The two statistics the reported set discards are mean_t (bias along the
%   task-IRRELEVANT direction) and cov_nt (the correlation, equivalently the
%   tilt of the principal axis away from the task axes). Both are near zero
%   at the nominal parameters -- Gate 2.2 found the principal axis within
%   0.3 deg of the task-irrelevant direction -- so they may carry little
%   information in practice. Whether they do is an empirical question this
%   analysis answers rather than assumes.
%
%   info returns the interpretable quantities (SD_rel, SD_irrel, aspect
%   ratio, bias, average SD, principal axis angle) for reporting. These are
%   NOT used in the discrepancy; they are derived from s.
%
%   Requires buildModel2D, augmentModel, solveLQGND, rolloutLQG on the path.

    if nargin < 4 || isempty(task), task = 'sim2'; end
    if nargin < 3 || isempty(seed), seed = 0;      end

    validateattributes(theta, {'numeric'}, {'vector', 'numel', 4, 'finite'}, ...
        mfilename, 'theta');
    theta = theta(:);

    sigma_u = theta(1);
    sigma_s = theta(2);
    d       = theta(3);
    r       = theta(4);

    assert(abs(d - round(d)) < 1e-9, 'ofcSummary:delayInteger', ...
        'Delay d must be an integer number of steps; got %g.', d);
    d = round(d);
    assert(d >= 0, 'ofcSummary:delayNegative', 'd must be >= 0.');
    assert(sigma_u > 0 && sigma_s > 0 && r > 0, 'ofcSummary:positive', ...
        'sigma_u, sigma_s and r must be positive.');

    % ---- forward model ----
    opts = struct('sigma_u', sigma_u, 'sigma_s', sigma_s, 'r', r);
    M = buildModel2D(task, opts);
    M = augmentModel(M, d);

    [L, K] = solveLQGND(M);

    rng(seed);
    [~, ~, Xf] = rolloutLQG(M, L, K, sigma_u, 0, N);

    % ---- task coordinates, derived from the cost (not hard-coded) ----
    dp   = M.Dpos(1:2)';
    nhat = dp / norm(dp);
    that = [-nhat(2); nhat(1)];

    P = Xf(1:2, :);                       % terminal positions, 2-by-N
    pn = nhat' * P;                       % task-relevant projection
    pt = that' * P;                       % task-irrelevant projection

    C = cov([pn' pt']);                   % 2-by-2

    s = [ mean(pn) ; mean(pt) ; C(1,1) ; C(2,2) ; C(1,2) ];

    if nargout > 1
        sd_n = sqrt(C(1,1));
        sd_t = sqrt(C(2,2));
        info = struct( ...
            'sd_rel',   sd_n, ...
            'sd_irrel', sd_t, ...
            'AR',       sd_t / sd_n, ...
            'avgSD',    (sd_n + sd_t)/2, ...
            'bias',     mean(M.Dpos * Xf), ...
            'angle',    0.5*atan2d(2*C(1,2), C(1,1)-C(2,2)), ...
            'N',        N, ...
            'theta',    theta, ...
            'task',     task, ...
            'nhat',     nhat, ...
            'that',     that);
    end
end
