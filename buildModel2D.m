function M = buildModel2D(task, opts)
%BUILDMODEL2D  Sim 1 and Sim 2 of Todorov & Jordan (2002), Supplementary Notes.
%
%   M = buildModel2D('sim1')          aiming: end anywhere on a line at -20 deg
%   M = buildModel2D('sim2')          intercept: two 1D masses meet
%   M = buildModel2D('sim1', opts)    override any parameter
%
%   Both simulations share the SAME plant: two independent unit point masses,
%   each driven by one actuator through a second-order filter. Sim 1 calls them
%   the x and y axes of a 2D mass; Sim 2 calls them two separate 1D masses.
%   Mechanically identical -- only Q and x1 differ.
%
%   State (8-D), interleaved so that kron(., eye(2)) builds it directly:
%       x = [p1; p2; v1; v2; f1; f2; g1; g2]
%
%   MUSCLE FILTER -- differs from Todorov (2005).
%   The 2002 supplement writes the force update as
%       f(t+dt) = exp(-dt/tau)*f(t) + u(t) + noise
%   i.e. EXACT exponential decay with a UNIT coefficient on u, so each filter
%   stage has DC gain 1/(1-exp(-dt/tau)) = 4.5208, and the second-order cascade
%   has DC gain 20.438 (newtons of steady force per unit u).
%
%   Todorov (2005) instead uses Euler discretisation, (1-dt/tau) with dt/tau on
%   the input, giving unit DC gain. Both are internally consistent but they are
%   NOT the same plant. Using the 2005 form here makes u ~20x larger for the
%   same force, hence the effort term ~418x too heavy, hence a movement that
%   never completes. Keep the two conventions straight.
%
%   No target is carried in the state: both task errors are homogeneous
%   quadratics, so there is no constant offset to absorb.
%
%   Timing: T = control steps, n = T+1 states, duration = T*dt.

    if nargin < 1 || isempty(task), task = 'sim1'; end
    if nargin < 2, opts = struct(); end

    % ---- defaults (Supplementary Notes sec 2, Sim 1-6) ----
    d.dt      = 0.01;      % s
    d.T       = 50;        % control steps -> 500 ms
    d.m       = 1;         % kg
    d.tau1    = 0.04;      % s, second-order muscle filter
    d.tau2    = 0.04;      % s
    d.wv      = 0.1;       % velocity stopping weight   (1D used 0.2)
    d.wf      = 0.01;      % force stopping weight      (1D used 0.02)
    d.r       = 0.002;     % effort weight
    d.sigma_u = 0.4;       % control-dependent noise scale
    d.sigma_s = 0.4;       % sensory noise scale
    d.P       = 1;         % number of POSITIONAL task constraints -> 1/(P+2)

    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(opts, f{i}), opts.(f{i}) = d.(f{i}); end
    end
    o = opts;

    dt = o.dt;  T = o.T;  n = T + 1;
    M.dt = dt;  M.T = T;  M.n = n;  M.opts = o;  M.task = task;

    % ---- single-axis 4-state block: [p; v; f; g] ----
    %   p(t+dt) = p + dt*v
    %   v(t+dt) = v + dt*f/m
    %   f(t+dt) = exp(-dt/tau2)*f + g          <- unit coefficient on g
    %   g(t+dt) = exp(-dt/tau1)*g + u          <- unit coefficient on u
    a1 = exp(-dt/o.tau1);
    a2 = exp(-dt/o.tau2);
    A1 = [1, dt, 0,  0;
          0, 1,  dt/o.m, 0;
          0, 0,  a2, 1;
          0, 0,  0,  a1];
    B1 = [0; 0; 0; 1];

    M.dcgain = 1 / ((1 - a1) * (1 - a2));    % steady force per unit u

    % ---- two independent copies, interleaved ----
    M.A = kron(A1, eye(2));                 % 8 x 8
    M.B = kron(B1, eye(2));                 % 8 x 2
    M.H = [eye(6), zeros(6,2)];             % p1 p2 v1 v2 f1 f2  (g not sensed)

    % ---- control-dependent noise, c = 2 ----
    % sigma_u * [e1 -e2; e2 e1] * u : circular 2D noise, SD = sigma_u*|u|.
    % NOT two independent per-axis noises -- that would be axis-aligned.
    Jrot = [0 -1; 1 0];
    M.C = zeros(8, 2, 2);
    M.C(:,:,1) = o.sigma_u * M.B;
    M.C(:,:,2) = o.sigma_u * M.B * Jrot;

    % ---- sensory noise ----
    sd = o.sigma_s * [0.01; 0.01; 0.1; 0.1; 1; 1];   % m, m, m/s, m/s, N, N
    M.Oomega = diag(sd.^2);

    M.Oxi  = zeros(8);
    M.Oeta = zeros(8);
    M.S1   = zeros(8);

    % ---- task-specific ----
    switch lower(task)
        case 'sim1'
            M.x1 = [0.2; 0.2; 0; 0; 0; 0; 0; 0];
            dpos = [tand(-20), -1, 0, 0, 0, 0, 0, 0];
        case 'sim2'
            M.x1 = [-0.1; 0.1; 0; 0; 0; 0; 0; 0];
            dpos = [1, -1, 0, 0, 0, 0, 0, 0];
        otherwise
            error('buildModel2D:task', 'task must be ''sim1'' or ''sim2''');
    end
    M.Dpos = dpos;

    % ---- terminal cost ----
    Dvel = o.wv * [0 0 1 0 0 0 0 0;
                   0 0 0 1 0 0 0 0];
    Dfor = o.wf * [0 0 0 0 1 0 0 0;
                   0 0 0 0 0 1 0 0];
    D = [dpos; Dvel; Dfor];

    Qn = (1 / (o.P + 2)) * (D' * D);
    Qn = (Qn + Qn') / 2;

    M.Q = zeros(8, 8, n);
    M.Q(:,:,n) = Qn;
    M.D = D;

    % ---- effort ----
    M.R = (o.r / T) * eye(2);
end