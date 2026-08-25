%% compensation_direction.m
% The tier 1 null direction: which change in theta leaves the three
% informative summary statistics unchanged.

theta0 = [0.4; 0.4; 4; 0.002];      % sigma_u, sigma_s, d, r
task   = 'sim2';
Npred  = 200000;                     % high N: we want the expectation
seed   = 12345;                      % common random numbers across points
keep   = [1 3 4];                    % the three informative statistics

% relative step sizes; d is integer so it steps by 1
rel = [0.05; 0.05; NaN; 0.05];

s0 = ofcSummary(theta0, Npred, seed, task);
J  = zeros(3,4);

for k = 1:4
    tp = theta0; tm = theta0;
    if k == 3
        tp(k) = theta0(k) + 1;  tm(k) = theta0(k) - 1;  h = 2;
    else
        step = rel(k)*theta0(k);
        tp(k) = theta0(k) + step;  tm(k) = theta0(k) - step;  h = 2*step;
    end
    sp = ofcSummary(tp, Npred, seed, task);
    sm = ofcSummary(tm, Npred, seed, task);
    J(:,k) = (sp(keep) - sm(keep)) / h;
end

% normalise columns to relative units so the four are comparable
scale = theta0;  %scale(3) = 1;          % d is already in natural units
Jrel  = J .* scale';

fprintf('Jacobian (relative units), rows = [mean_n, var_n, var_t]:\n');
disp(Jrel)
fprintf('rank = %d,  cond = %.3e\n', rank(Jrel), cond(Jrel));

v = null(Jrel);
fprintf('\nnull direction (relative change per unit step):\n');
lbl = {'sigma_u','sigma_s','d','r'};
for k = 1:4
    fprintf('  %-8s %+.4f\n', lbl{k}, v(k));
end
fprintf('\nsigns: %s\n', mat2str(sign(v)'));