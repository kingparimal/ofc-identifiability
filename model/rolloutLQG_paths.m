function Xf = rolloutLQG_paths(M, L, K, Xi, Om)
%ROLLOUTLQG_PATHS  rolloutLQG with EXTERNALLY SUPPLIED noise.
%   Xi : nc-by-N-by-(n-1) control-dependent draws
%   Om : ny-by-N-by-(n-1) sensory draws
%   Must be bitwise identical to rolloutLQG when fed that function's own draws.
    A=M.A; B=M.B; H=M.H; n=M.n; C=M.C; nc=size(C,3); ny=size(H,1);
    N = size(Xi,2);
    sd = sqrt(diag(M.Oomega));
    X = repmat(M.x1,1,N);  XH = repmat(M.x1,1,N);
    for t = 1:n-1
        U  = -L(:,:,t)*XH;
        XN = A*X + B*U;
        for j = 1:nc
            XN = XN + (C(:,:,j)*U) .* Xi(j,:,t);
        end
        Y  = H*X + sd .* Om(:,:,t);
        XH = A*XH + B*U + K(:,:,t)*(Y - H*XH);
        X  = XN;
    end
    Xf = X;
end
