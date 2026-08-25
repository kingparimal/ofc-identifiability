% Pathwise test of the Sim 2 axis-swap symmetry.  MATLAB and Octave compatible.
% Requires: buildModel2D, augmentModel, solveLQGND, rolloutLQG, rolloutLQG_paths
%
% Two gates:
%   gate  -- rolloutLQG_paths must reproduce rolloutLQG BITWISE (expect 0)
%   sym   -- transformed noise must give exactly -S*x  (expect < 1e-13 for sim2,
%            and visibly nonzero for sim1, which breaks the symmetry via Q)
sw    = [0 1; 1 0];
tasks = {'sim2','sim1'};
fprintf('%-6s %3s %12s %10s %10s %14s\n','task','d','gate(bitwise)','dL','dK','symmetry');
for a = 1:numel(tasks)
    task = tasks{a};
    for d = [0 4]
        M = augmentModel(buildModel2D(task), d);
        [L, K] = solveLQGND(M);
        N  = 64;
        nc = size(M.C,3);
        ny = size(M.H,1);
        Sy = kron(eye(3), sw);
        Su = sw;
        Sa = blkdiag(kron(eye(4),sw), kron(eye(d), Sy));

        % gain equivariance
        eL = 0; eK = 0;
        for t = 1:M.n-1
            eL = max(eL, norm(Su*L(:,:,t) - L(:,:,t)*Sa, 'fro'));
            eK = max(eK, norm(Sa*K(:,:,t) - K(:,:,t)*Sy, 'fro'));
        end

        % reference rollout, then the same draws replayed through the injectable one
        rng(7);
        [~, ~, Xa] = rolloutLQG(M, L, K, M.opts.sigma_u, 0, N);
        rng(7);
        Xi = zeros(nc, N, M.n-1);
        Om = zeros(ny, N, M.n-1);
        for t = 1:M.n-1
            for j = 1:nc
                Xi(j,:,t) = randn(1,N);
            end
            Om(:,:,t) = randn(ny,N);
        end
        Xb   = rolloutLQG_paths(M, L, K, Xi, Om);
        gate = max(abs(Xa(:) - Xb(:)));

        % transformed noise: xi2 -> -xi2,  omega -> -Sy*omega
        Xi2 = Xi;  Xi2(2,:,:) = -Xi2(2,:,:);
        Om2 = zeros(size(Om));
        for t = 1:M.n-1
            Om2(:,:,t) = -Sy * Om(:,:,t);
        end
        Xc  = rolloutLQG_paths(M, L, K, Xi2, Om2);
        sym = max(abs(Xc(:) - reshape(-Sa*Xa, [], 1)));

        fprintf('%-6s %3d %12.1e %10.1e %10.1e %14.2e\n', task, d, gate, eL, eK, sym);
    end
end
