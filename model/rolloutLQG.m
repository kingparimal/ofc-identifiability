function [endpoints, umag, Xf] = rolloutLQG(M, L, K, sigma_c, sigma_add, N)
%ROLLOUTLQG  Monte Carlo of the closed loop, vectorised over trials.
%   Controller acts on the ESTIMATE, never on the true state.
%
%   endpoints : final value of state 1, one per trial   (N-by-1)
%   umag      : |sigma_c * u| at every step             (N*(n-1)-by-1)
%   Xf        : full terminal state, all trials         (nx-by-N)
%
%   All N trials advance together as columns of X. Use Xf for multi-DOF
%   models where "the endpoint" is not a single scalar.
%
%   The sigma_c argument scales the umag calibration output only; the motor
%   noise itself comes entirely from M.C. Keep the two consistent.
%
%   Trials start exactly at M.x1 (M.S1 is not sampled). Fine for Sim 1 and
%   Sim 2, which have Sigma_1 = 0; revisit for Sim 5/6.

    A = M.A;  B = M.B;  H = M.H;  n = M.n;
    C = M.C;  nc = size(C,3);
    nu = size(B,2);
    ny = size(H,1);
    sd = sqrt(diag(M.Oomega));

    X  = repmat(M.x1, 1, N);        % true state,  nx-by-N
    XH = repmat(M.x1, 1, N);        % estimate,    nx-by-N

    umag = zeros(N*(n-1), 1);
    idx  = 0;

    for t = 1:n-1
        U = -L(:,:,t) * XH;                          % nu-by-N, on the estimate

        umag(idx + (1:N)) = sqrt(sum((sigma_c*U).^2, 1))';
        idx = idx + N;

        XN = A*X + B*U;                              % deterministic part
        for j = 1:nc
            XN = XN + (C(:,:,j)*U) .* randn(1, N);   % multiplicative motor noise
        end
        if sigma_add > 0
            XN = XN + B * (sigma_add * randn(nu, N));
        end

        Y  = H*X + sd .* randn(ny, N);               % observe CURRENT state
        XH = A*XH + B*U + K(:,:,t)*(Y - H*XH);
        X  = XN;
    end

    Xf        = X;
    endpoints = X(1,:)';
end
