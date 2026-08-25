% Matrix-level equivariance of the model under x -> -S x, u -> -Su u.
% All residuals zero  =>  the terminal distribution is exactly invariant under
% (p1,p2) -> (-p2,-p1), hence mean_t == 0 and cov_nt == 0 for ALL theta.
% MATLAB and Octave compatible.
sw    = [0 1; 1 0];
tasks = {'sim1','sim2'};
for a = 1:numel(tasks)
    task = tasks{a};
    for d = [0 4]
        M  = augmentModel(buildModel2D(task), d);
        Sx = kron(eye(4), sw);
        Sy = kron(eye(3), sw);
        Su = sw;
        Sa = blkdiag(Sx, kron(eye(d), Sy));
        rA  = norm(Sa*M.A - M.A*Sa, 'fro');
        rB  = norm(Sa*M.B - M.B*Su, 'fro');
        rH  = norm(Sy*M.H - M.H*Sa, 'fro');
        rOm = norm(Sy*M.Oomega*Sy - M.Oomega, 'fro');
        rQ  = norm(Sa'*M.Q(:,:,end)*Sa - M.Q(:,:,end), 'fro');
        rR  = norm(Su'*M.R*Su - M.R, 'fro');
        rC1 = norm(Sa*M.C(:,:,1)*Su - M.C(:,:,1), 'fro');   % expect +
        rC2 = norm(Sa*M.C(:,:,2)*Su + M.C(:,:,2), 'fro');   % expect -
        rx1 = norm(Sa*M.x1 + M.x1);                         % need S*x1 = -x1
        fprintf(['%-6s d=%d nx=%3d | A=%.1e B=%.1e H=%.1e Oom=%.1e ' ...
                 'Q=%.1e R=%.1e C1=%.1e C2=%.1e x1=%.1e\n'], ...
                task, d, size(M.A,1), rA, rB, rH, rOm, rQ, rR, rC1, rC2, rx1);
    end
end
