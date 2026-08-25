function D = discrepancy(s_obs, s_pred, Sigma)
%DISCREPANCY  Calibrated statistical distance between two summary vectors.
%
%   D = DISCREPANCY(s_obs, s_pred, Sigma)
%
%       D = (s_obs - s_pred)' * inv(Sigma) * (s_obs - s_pred)
%
%   If s_obs is a draw from a distribution with mean s_pred and covariance
%   Sigma, and that distribution is approximately Gaussian, then D follows a
%   chi-square distribution with 5 degrees of freedom -- the DIMENSION OF THE
%   SUMMARY VECTOR, not the number of parameters.
%
%   That distinction matters for the threshold. The confidence region
%
%       { theta : D(theta) <= chi2inv(0.95, 5) }
%
%   contains the true theta with probability 0.95. Using the number of
%   parameters (4) instead would give a region that under-covers.
%
%   For a PROFILE of a single parameter -- fixing it and minimising over the
%   rest -- the appropriate threshold is chi2inv(0.95, 1) applied to the
%   increase in D above its global minimum, which is the standard profile
%   likelihood construction. The two thresholds answer different questions
%   and should not be interchanged.
%
%   Solved by mldivide rather than inv() for numerical stability. If Sigma is
%   near-singular this will warn; that indicates redundant summary statistics
%   rather than a numerical accident.

    s_obs  = s_obs(:);
    s_pred = s_pred(:);

    assert(numel(s_obs) == numel(s_pred), 'discrepancy:dim', ...
        'Summary vectors differ in length: %d vs %d.', ...
        numel(s_obs), numel(s_pred));
    assert(isequal(size(Sigma), [numel(s_obs) numel(s_obs)]), ...
        'discrepancy:sigmaDim', 'Sigma must be %d-by-%d.', ...
        numel(s_obs), numel(s_obs));

    e = s_obs - s_pred;
    D = e' * (Sigma \ e);
end
