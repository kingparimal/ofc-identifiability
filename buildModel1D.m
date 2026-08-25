function M = buildModel1D(n)
    if nargin < 1, n = 30; end     % default duration: 30 steps (0.29 s)
%BUILDMODEL1D  1D reaching model of Todorov (2005), Neural Computation, sec 6.1.
%   Single source of truth for the plant, cost, and noise. Every script and
%   Simulink block should read its matrices from here, never hard-code them.
%
%   Returns struct M with fields consumed by solveLQG:
%     dt n            time step and horizon
%     A B H           dynamics / control / observation
%     C               control-dependent noise scaling (m x p x c); here c = 1
%     Q               state cost, m x m x n  (zero except final step)
%     R               control cost (already divided by n-1)
%     Oxi Oomega Oeta process / sensory / internal noise covariances
%     x1 S1           initial state mean and covariance

    % ---- timing ----
    dt   = 0.01;          % s
    tau1 = 0.04;          % muscle filter time constants
    tau2 = 0.04;
    m    = 1;             % kg
    %n    = 30;            % steps -> 0.30 s
    M.dt = dt;
    M.n  = n;

    % ---- dynamics.  state x = [p; pdot; f; g; p_target] ----
    M.A = [1, dt, 0,          0,          0;
           0, 1,  dt/m,       0,          0;
           0, 0,  1-dt/tau2,  dt/tau2,    0;
           0, 0,  0,          1-dt/tau1,  0;
           0, 0,  0,          0,          1];
    M.B = [0; 0; 0; dt/tau1; 0];
    M.H = [1 0 0 0 0;
           0 1 0 0 0;
           0 0 1 0 0];

    % ---- noise ----
    sigma_c = 0.5;
    M.C = sigma_c * M.B;                          % control-dependent (c = 1)

    sigma_s = 0.5;
    M.Oomega = diag( (sigma_s*[0.02; 0.2; 1]).^2 );   % square the VECTOR (correct form)

    M.Oxi  = zeros(5);      % no additive process noise
    M.Oeta = zeros(5);      % no internal noise
    M.S1   = zeros(5);      % initial state known exactly

    % ---- initial state ----
    p_target = 0.1;                              % 10 cm
    M.x1 = [0; 0; 0; 0; p_target];

    % ---- cost ----
    w_v = 0.2;  w_f = 0.02;  r = 1e-5;
    p_vec = [1; 0;   0;   0; -1];
    v_vec = [0; w_v; 0;   0;  0];
    f_vec = [0; 0;   w_f; 0;  0];
    Qn = p_vec*p_vec' + v_vec*v_vec' + f_vec*f_vec';

    % solveLQG expects a time-varying stack Q(:,:,t); only the final step is nonzero
    M.Q = zeros(5, 5, n);
    M.Q(:,:,n) = Qn;

    M.R = r / (n-1);                             % effort penalty (paper's 1/(n-1))
end

