%% gate_tier2.m -- gates that must pass before ANY tier 2 result is trusted.
%
% G1  rolloutLQG_traj reproduces rolloutLQG's terminal state BITWISE.
%     Recording positions adds no random draws, so the streams must be
%     identical. Anything other than 0 means the draw order changed and every
%     tier 2 number would be computed from a different process than tier 1.
%
% G2  The axis-swap symmetry holds ALONG THE TRAJECTORY, not just at the end.
%     If it does, mean_t and cov_nt are exactly zero at every sampled time and
%     are correctly excluded from ofcSummary2. Tested by checking they scale
%     as 1/sqrt(N) (noise) rather than settling on a nonzero value (signal).
%
% G3  Sigma is well conditioned at the chosen time points -- i.e. the variance
%     profile is not numerically redundant across them.
%
% G4  The tier 2 statistics actually MOVE along the escape direction that
%     defeated tier 1 (sigma_u, sigma_s, r all decreasing together). If they
%     do not, tier 2 cannot close the region and there is no point proceeding.

clear; clc;
theta0=[0.4;0.4;4;0.002]; task='sim2'; tS=[16 23 30 37 44 51];
M = augmentModel(buildModel2D(task), 4);
[L,K] = solveLQGND(M);

fprintf('=== G1: bitwise agreement with rolloutLQG ===\n');
for N=[1 64]
    rng(11); [~,~,Xa] = rolloutLQG(M,L,K,0.4,0,N);
    rng(11); Xb = rolloutLQG_traj(M,L,K,0.4,0,N,tS);
    fprintf('  N=%3d  max|Xf_traj - Xf| = %.3e\n', N, max(abs(Xa(:)-Xb(:))));
end

fprintf('\n=== G2: symmetry along the trajectory, PATHWISE and EXACT ===\n');
% The transformation maps the WHOLE trajectory x_t -> -S x_t, so this is a
% machine-precision identity, not a 1/sqrt(N) trend. Expect < 1e-13 for Sim 2.
sw=[0 1;1 0]; Sy=kron(eye(3),sw); nc=size(M.C,3); ny=size(M.H,1); N=64;
rng(21);
Xi=zeros(nc,N,M.n-1); Om=zeros(ny,N,M.n-1);
for t=1:M.n-1
    for j=1:nc, Xi(j,:,t)=randn(1,N); end
    Om(:,:,t)=randn(ny,N);
end
[~,Pa] = rolloutLQG_traj(M,L,K,0.4,0,N,tS,Xi,Om);
Xi2=Xi; Xi2(2,:,:)=-Xi2(2,:,:);
Om2=zeros(size(Om));
for t=1:M.n-1, Om2(:,:,t) = -Sy*Om(:,:,t); end
[~,Pb] = rolloutLQG_traj(M,L,K,0.4,0,N,tS,Xi2,Om2);
err=0;
for q=1:numel(tS)
    err = max(err, max(max(abs(Pb(:,:,q) - (-sw*Pa(:,:,q))))));
end
fprintf('  max |P_transformed - (-S P)| over all %d sampled times = %.3e\n', numel(tS), err);
fprintf('  (Sim 2 must be < 1e-13; this is what justifies scoring 3 stats\n');
fprintf('   per time point rather than 5.)\n');

fprintf('\n=== G3: conditioning of the tier 2 Sigma ===\n');
Sig2 = statsCovariance2(theta0,200,400,60000,task,tS);
p=3*numel(tS);
fprintf('  %d statistics, nRep 400: cond(Sigma) = %.3e\n', p, cond(Sig2));
fprintf('  Hotelling 95%% threshold %.2f  (chi2 ideal %.2f, inflation %.1f%%)\n', ...
        hotelling(p,400), chi2q(0.95,p), 100*(hotelling(p,400)/chi2q(0.95,p)-1));

fprintf('\n=== G4: do the tier 2 statistics move along the escape direction? ===\n');
fprintf('  (walls.m: region escapes with sigma_u, sigma_s and r ALL decreasing)\n');
s0 = ofcSummary2(theta0,20000,12345,task,tS);
fprintf('  %-34s %12s %12s\n','theta','D_tier1','D_tier2');
Sig1 = statsCovariance(theta0,200,400,60000,task);
i3=[1 3 4]; S1=Sig1(i3,i3);
b0 = ofcSummary(theta0,20000,12345,task);
for f=[0.85 0.70 0.55]
    th=[0.4*f; 0.4*f; 4; 0.002*f];
    s1=ofcSummary(th,20000,12345,task);
    s2=ofcSummary2(th,20000,12345,task,tS);
    D1=discrepancy(s1(i3),b0(i3),S1);
    D2=discrepancy(s2,s0,Sig2);
    fprintf('  all x%.2f  (su %.3f ss %.3f r %.1e) %12.2f %12.2f\n', ...
            f, th(1), th(2), th(4), D1, D2);
end
fprintf('\n  thresholds: tier1 %.2f, tier2 %.2f\n', hotelling(3,400), hotelling(p,400));
fprintf('  PASS if D_tier2 rises far above its threshold where D_tier1 does not.\n');
