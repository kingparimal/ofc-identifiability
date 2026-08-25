function [L, K, hist] = solveLQGND(M, maxIter, tol)
%SOLVELQGND  Coordinate descent of Todorov (2005), Eq 4.2 (backward) and 5.2 (forward).
%   [L,K,hist] = solveLQGND(M) with M from buildModel1D / buildModel2D.
%   hist is the expected cost per iteration; it MUST decrease monotonically.
%
%   Dimension-general: state, control and observation sizes are read from
%   A, B and H. Supports c >= 1 control-dependent noise sources: M.C may be
%   nx-by-nu (single source) or nx-by-nu-by-c. Supports a full time-varying
%   state cost M.Q of size nx-by-nx-by-n, so via-point tasks are expressible.

    if nargin < 2, maxIter = 50;    end
    if nargin < 3, tol     = 1e-12; end

    A = M.A;  B = M.B;  H = M.H;  C = M.C;  Q = M.Q;
    Oom = M.Oomega;  Oxi = M.Oxi;  Oeta = M.Oeta;
    x1 = M.x1;  S1 = M.S1;  R = M.R;  n = M.n;

    nx = size(A,1);          % state dimension
    nu = size(B,2);          % control dimension
    ny = size(H,1);          % observation dimension
    nc = size(C,3);          % number of control-dependent noise sources

    assert(size(Q,3) == n, 'M.Q must be nx-by-nx-by-n (one slice per time step)');

    L = zeros(nu, nx, n-1);
    K = zeros(nx, ny, n-1);
    hist = nan(maxIter, 1);

    for iter = 1:maxIter
        % ---- BACKWARD: controller (Eq 4.2) ----
        Sx = Q(:,:,n);  Se = zeros(nx);  s = 0;
        for t = n-1:-1:1
            % sum_i Ci' * (Sx + Se) * Ci
            CSC = zeros(nu);
            for i = 1:nc
                Ci  = C(:,:,i);
                CSC = CSC + Ci'*(Sx + Se)*Ci;
            end

            L(:,:,t) = (R + B'*Sx*B + CSC) \ (B'*Sx*A);
            Kt = K(:,:,t);  Lt = L(:,:,t);
            s = s + trace(Sx*Oxi + Se*(Oxi + Oeta + Kt*Oom*Kt'));
            Sx_new = Q(:,:,t) + A'*Sx*(A - B*Lt);      % <-- Q(:,:,t) was missing
            Se_new = A'*Sx*B*Lt + (A - Kt*H)'*Se*(A - Kt*H);
            Sx = Sx_new;  Se = Se_new;
        end
        hist(iter) = x1'*Sx*x1 + trace((Sx + Se)*S1) + s;

        % ---- FORWARD: estimator (Eq 5.2) ----
        Sef = S1;  Sxh = x1*x1';  Sxe = zeros(nx);
        for t = 1:n-1
            K(:,:,t) = A*Sef*H' / (H*Sef*H' + Oom);
            Kt = K(:,:,t);  Lt = L(:,:,t);  ABL = A - B*Lt;

            % sum_i Ci * Lt * Sxh * Lt' * Ci'
            CLSLC = zeros(nx);
            for i = 1:nc
                Ci    = C(:,:,i);
                CLSLC = CLSLC + Ci*Lt*Sxh*Lt'*Ci';
            end

            % NOTE: trailing factor is A', NOT (A-Kt*H)'.  This is the
            % post-substitution form and only holds at the optimal Kt.
            Sef_n = Oxi + Oeta + (A - Kt*H)*Sef*A' + CLSLC;
            Sxh_n = Oeta + Kt*H*Sef*A' + ABL*Sxh*ABL' ...
                    + ABL*Sxe*H'*Kt' + Kt*H*Sxe'*ABL';
            Sxe_n = ABL*Sxe*(A - Kt*H)' - Oeta;
            Sef = Sef_n;  Sxh = Sxh_n;  Sxe = Sxe_n;
        end

        if iter > 1 && abs(hist(iter) - hist(iter-1)) < tol*abs(hist(iter))
            break
        end
    end
    hist = hist(1:iter);
end