function [Xf, Ptraj, endpoints] = rolloutLQG_traj(M, L, K, sigma_c, sigma_add, N, tSample, Xi, Om)
%ROLLOUTLQG_TRAJ  rolloutLQG, additionally recording positions along the way.
%
%   [Xf, Ptraj] = rolloutLQG_traj(M, L, K, sigma_c, sigma_add, N, tSample)
%   [Xf, Ptraj] = rolloutLQG_traj(..., tSample, Xi, Om)
%
%   Ptraj : 2-by-N-by-numel(tSample), positions p1,p2 at the requested STATE
%           indices. State t is at time (t-1)*dt; t = 1 is the start, t = n the
%           terminal state.
%
%   With Xi and Om supplied the noise is INJECTED rather than drawn:
%       Xi : nc-by-N-by-(n-1)   control-dependent draws
%       Om : ny-by-N-by-(n-1)   sensory draws
%   That mode exists so the axis-swap symmetry can be tested pathwise and
%   exactly, rather than statistically. See gate_tier2.m G2.
%
%   CRITICAL: in drawing mode this consumes random numbers in exactly the order
%   rolloutLQG does, so with the same seed it reproduces rolloutLQG's terminal
%   state BITWISE. Recording positions adds no draws. gate_tier2.m G1 checks
%   this; if it fails, nothing downstream is trustworthy.

    A=M.A; B=M.B; H=M.H; n=M.n; C=M.C; nc=size(C,3);
    nu=size(B,2); ny=size(H,1);
    sd = sqrt(diag(M.Oomega));
    inject = (nargin >= 9) && ~isempty(Xi) && ~isempty(Om);

    X  = repmat(M.x1,1,N);
    XH = repmat(M.x1,1,N);
    Ptraj = zeros(2, N, numel(tSample));

    for t = 1:n-1
        q = find(tSample==t,1);
        if ~isempty(q), Ptraj(:,:,q) = X(1:2,:); end     % X is state t here

        U = -L(:,:,t) * XH;
        XN = A*X + B*U;
        for j = 1:nc
            if inject, z = Xi(j,:,t); else, z = randn(1,N); end
            XN = XN + (C(:,:,j)*U) .* z;
        end
        if sigma_add > 0
            XN = XN + B * (sigma_add * randn(nu,N));
        end
        if inject, w = Om(:,:,t); else, w = randn(ny,N); end
        Y  = H*X + sd .* w;
        XH = A*XH + B*U + K(:,:,t)*(Y - H*XH);
        X  = XN;
    end
    q = find(tSample==n,1);
    if ~isempty(q), Ptraj(:,:,q) = X(1:2,:); end          % terminal state

    Xf = X;
    endpoints = X(1,:)';
end
