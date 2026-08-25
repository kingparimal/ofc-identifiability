%% sim1_5stat.m -- ITEM 2, the decisive part. Does Sim 1 gain from its two
%% extra statistics at tier 1?
%
% Sim 1 breaks the axis-swap symmetry (only Q and x1 differ from Sim 2), so
% mean_t and cov_nt are NOT identically zero there. The open question is
% whether they are INFORMATIVE -- derivative nonzero -- or merely nonzero.
% mean_t sat at ~0.119 m and near constant across the theta tried, which would
% make it useless.
%
% rank(J) with all FIVE statistics answers it. If rank is still 3, the extra
% two carry no derivative information and Sim 1 is no better identified than
% Sim 2 at tier 1.

clear; clc;
theta0=[0.4;0.4;4;0.002]; N=200; nRep=400; h=0.02;
for task = {'sim2','sim1'}
    t = task{1};
    % full 5-statistic Jacobian
    f = @(th) ofcSummaryP(th,20000,12345,t,1);
    s0=f(theta0); J=zeros(5,4);
    for q=1:4
        tp=theta0; tm=theta0;
        if q==3, tp(3)=5; tm(3)=3; den=2;
        else, tp(q)=theta0(q)*10^h; tm(q)=theta0(q)*10^(-h); den=2*h; end
        J(:,q)=(f(tp)-f(tm))/den;
    end
    sv = svd(J);
    fprintf('=== %s, tier 1, all 5 statistics ===\n', t);
    fprintf('  singular values of J: %s\n', mat2str(sv',3));
    fprintf('  rank(J) = %d of 4\n', rank(J,1e-8*norm(J)));
    % how much do the two extra rows move, relative to the three informative?
    fprintf('  row norms  mean_n %.3e  mean_t %.3e  var_n %.3e  var_t %.3e  cov_nt %.3e\n', ...
            norm(J(1,:)), norm(J(2,:)), norm(J(3,:)), norm(J(4,:)), norm(J(5,:)));
    fprintf('  |d mean_t/d theta| as %% of |d mean_n/d theta| : %.2f%%\n\n', ...
            100*norm(J(2,:))/norm(J(1,:)));
end
%% ---- seed stability: a real derivative is stable, noise is not ----
fprintf('=== is d(mean_t)/d(theta) SIGNAL or Monte Carlo noise? ===\n');
fprintf('A real derivative reproduces across seeds; noise does not.\n\n');
seeds=[12345 23456 34567];
for task = {'sim2','sim1'}
    t=task{1};
    R=zeros(numel(seeds),4);   % d(mean_t)/d(theta) row, per seed
    R1=zeros(numel(seeds),4);  % d(mean_n)/d(theta) row, per seed (control)
    for a=1:numel(seeds)
        f=@(th) ofcSummaryP(th,20000,seeds(a),t,1);
        for q=1:4
            tp=theta0; tm=theta0;
            if q==3, tp(3)=5; tm(3)=3; den=2;
            else, tp(q)=theta0(q)*10^h; tm(q)=theta0(q)*10^(-h); den=2*h; end
            dd=(f(tp)-f(tm))/den;
            R(a,q)=dd(2); R1(a,q)=dd(1);
        end
    end
    cv  = std(R,0,1)./abs(mean(R,1));
    cv1 = std(R1,0,1)./abs(mean(R1,1));
    fprintf('  %s\n', t);
    fprintf('    mean_t row, coeff of variation across seeds: %s\n', mat2str(cv,2));
    fprintf('    mean_n row (CONTROL)                       : %s\n', mat2str(cv1,2));
    if median(cv) > 0.5
        verdict = 'is NOISE (unstable across seeds)';
    else
        verdict = 'looks like SIGNAL (stable across seeds)';
    end
    fprintf('    -> mean_t %s\n', verdict);
    fprintf('       ratio of medians, mean_t CV / control CV: %.0fx\n\n', ...
            median(cv)/median(cv1));
end
fprintf(['VERDICT. A rank test on a Monte Carlo Jacobian cannot distinguish an\n' ...
         'exactly-zero derivative from noise: Sim 2 returns rank 4 with five\n' ...
         'statistics even though rows 2 and 5 are PROVED zero. Seed stability\n' ...
         'is the test that separates them.\n']);
