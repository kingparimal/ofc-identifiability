%% robustness.m -- ITEMS 2 and 4, via curvature rather than grids.
%
% Both are "does the conclusion survive a different modelling choice". The
% conclusion is about RANK and WIDTHS, and sizing.m gives both from
% F = J' inv(Sigma) J in seconds, box independently. Rebuilding a 34,650-point
% grid per variant would cost ~3 h each and answer the same question.
%
% ITEM 2 -- Sim 1. Sim 1 breaks the axis-swap symmetry (only Q and x1 differ),
% so mean_t and cov_nt are NOT identically zero there and tier 1 could in
% principle have 5 informative statistics instead of 3. The question is whether
% they are INFORMATIVE (derivative nonzero) or merely nonzero: mean_t sat at
% ~0.119 m and near constant across the theta tried, which would make it
% useless. rank(J) answers it directly.
%
% ITEM 4 -- Q scale. buildModel2D takes P as an option and sets the terminal
% cost scale to 1/(P+2). The adopted value is P = 1 -> 1/3; the untested
% alternative in caveat 1 is 1/4, i.e. P = 2. One option change, no code edit.
%
% Reported per variant: cond(F), numerical rank of J, and the four profile
% widths. A variant that changes the RANK changes the study's conclusion; one
% that only shifts widths by a few percent does not.

clear; clc;
theta0 = [0.4;0.4;4;0.002];
N = 200;  nRep = 400;  h = 0.02;        % h = finite-difference step, log10 units
tS = [16 23 30 37 44 51];
nm = {'sigma_u','sigma_s','delay','r'};

variants = { ...
  struct('name','Sim 2, P=1 (ADOPTED)', 'task','sim2', 'P',1), ...
  struct('name','Sim 2, P=2 (Q=1/4)',   'task','sim2', 'P',2), ...
  struct('name','Sim 1, P=1',           'task','sim1', 'P',1)};

for v = 1:numel(variants)
    V = variants{v};
    fprintf('\n======== %s ========\n', V.name);
    for tier = 1:2
        [J, Sig] = jac(theta0, V.task, V.P, N, nRep, tier, tS, h);
        r = rank(J, 1e-8*norm(J));
        F = J'*(Sig\J);
        fprintf('  tier %d: %d stats | rank(J) = %d of 4 | cond(F) = %.2e\n', ...
                tier, size(J,1), r, cond(F));
        if r == 4
            hw = sqrt(3.841458820694124*diag(inv(F)));
            fprintf('          widths:');
            for q=1:4
                if q==3, fprintf('  %s %.2f steps', nm{q}, 2*hw(q));
                else,    fprintf('  %s x%.3f', nm{q}, 10^(2*hw(q))); end
            end
            fprintf('\n');
        else
            fprintf('          RANK DEFICIENT -- widths undefined along the null direction\n');
        end
    end
end
fprintf(['\nRANK is the quantity that matters. rank(J) = 3 at tier 1 is the\n' ...
         'R^4 -> R^3 result; if Sim 1 also gives 3, its two extra statistics\n' ...
         'carry no derivative information and Sim 1 is no better off.\n']);
