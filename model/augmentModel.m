function Ma = augmentModel(M, d)
%AUGMENTMODEL  Add d extra steps of sensory delay by state augmentation.
%
%   Ma = augmentModel(M, d)
%
%   The problem formulation already contains ONE implicit step of delay: y_t
%   is received after u_t has been generated. The supplement adds d further
%   steps by carrying delayed sensory copies in the state:
%
%       xtilde_t = [ x_t ; H*x_{t-1} ; H*x_{t-2} ; ... ; H*x_{t-d} ]
%
%   so the augmented dimension is nx + d*ny. The new observation matrix picks
%   out the OLDEST copy, H*x_{t-d}; the new dynamics matrix drops it, shifts
%   the rest along, and inserts H*x_t at the front.
%
%   Total sensorimotor delay is (1 + d)*dt. For Sim 1-6, d = 4 gives 50 ms.
%
%   augmentModel(M, 0) returns M unchanged -- worth testing, it catches
%   every off-by-one in the block indexing.

    if nargin < 2 || d == 0
        Ma = M;
        return
    end

    A = M.A;  B = M.B;  H = M.H;
    nx = size(A,1);  nu = size(B,2);  ny = size(H,1);
    na = nx + d*ny;                            % augmented dimension

    % ---- dynamics ----
    Aa = zeros(na);
    Aa(1:nx, 1:nx) = A;                        % x_{t+1} = A x_t + B u_t
    Aa(nx+(1:ny), 1:nx) = H;                   % first delayed copy = H x_t
    for k = 2:d                                % shift the rest along
        dst = nx + (k-1)*ny + (1:ny);
        src = nx + (k-2)*ny + (1:ny);
        Aa(dst, src) = eye(ny);
    end

    Ba = [B; zeros(d*ny, nu)];

    % ---- observation: the OLDEST copy ----
    Ha = zeros(ny, na);
    Ha(:, nx + (d-1)*ny + (1:ny)) = eye(ny);

    % ---- noise ----
    nc = size(M.C, 3);
    Ca = zeros(na, nu, nc);
    for i = 1:nc
        Ca(1:nx, :, i) = M.C(:,:,i);
    end

    Oxi  = zeros(na);   Oxi(1:nx,1:nx)  = M.Oxi;
    Oeta = zeros(na);   Oeta(1:nx,1:nx) = M.Oeta;
    S1   = zeros(na);   S1(1:nx,1:nx)   = M.S1;

    % ---- cost ----
    n  = M.n;
    Qa = zeros(na, na, n);
    Qa(1:nx, 1:nx, :) = M.Q;                   % delayed copies are NOT penalised

    % ---- initial state: assume the system was at rest before onset ----
    x1a = [M.x1; repmat(H*M.x1, d, 1)];

    % ---- assemble, preserving everything downstream needs ----
    Ma = M;
    Ma.A = Aa;  Ma.B = Ba;  Ma.H = Ha;  Ma.C = Ca;
    Ma.Oxi = Oxi;  Ma.Oeta = Oeta;  Ma.S1 = S1;
    Ma.Q = Qa;  Ma.x1 = x1a;
    Ma.Dpos = [M.Dpos, zeros(1, d*ny)];
    Ma.delay = d;
    Ma.nx_plant = nx;
end
