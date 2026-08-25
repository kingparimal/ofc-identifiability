%% fig1_endpoint_clouds.m
% Figure 1: the minimum intervention principle, and what removing it costs.
% Rows: Sim 1 (aiming), Sim 2 (interception)
% Cols: optimal controller, desired-state controller
%
% Plotted in TASK COORDINATES, as deviations from each cloud's own mean.
% Task-irrelevant on x, task-relevant on y, so every ellipse is horizontal.
% All four panels share identical axes, so both the elongation and the
% overall size are directly comparable across panels.

clear; clc;
N = 5000;

R = struct();
for k = 1:2
    tk = {'sim1','sim2'}; tk = tk{k};
    M  = buildModel2D(tk);
    o  = M.opts;

    dp   = M.Dpos(1:2)';
    nhat = dp/norm(dp);   that = [-nhat(2); nhat(1)];

    % ---- optimal controller ----
    rng(0);
    [L, Kk] = solveLQGND(M);
    [~,~,Xf] = rolloutLQG(M, L, Kk, o.sigma_u, 0, N);
    Po = Xf(1:2,:);
    R(k).opt.n = nhat'*Po;   R(k).opt.t = that'*Po;
    pstar = mean(Po,2);

    % ---- desired-state controller ----
    Md = struct();
    Md.dt=M.dt; Md.T=M.T; Md.n=M.n; Md.opts=o; Md.task=[tk '_desired'];
    Md.A = blkdiag(M.A,1);
    Md.B = [M.B; 0 0];
    Md.H = [M.H, zeros(6,1)];
    Md.C = zeros(9,2,2);
    for i = 1:2, Md.C(:,:,i) = [M.C(:,:,i); 0 0]; end
    Md.Oomega = M.Oomega;
    Md.Oxi = zeros(9); Md.Oeta = zeros(9); Md.S1 = zeros(9);
    Md.x1  = [M.x1; 1];
    D2 = [1 0 0 0 0 0 0 0 -pstar(1);
          0 1 0 0 0 0 0 0 -pstar(2);
          zeros(1,2), o.wv, 0, zeros(1,4), 0;
          zeros(1,3), o.wv,    zeros(1,4), 0;
          zeros(1,4), o.wf, 0, zeros(1,2), 0;
          zeros(1,5), o.wf,    zeros(1,2), 0];
    Qd = (1/(2+2))*(D2'*D2);
    Md.Q = zeros(9,9,M.n);  Md.Q(:,:,M.n) = (Qd+Qd')/2;
    Md.R = M.R;

    rng(0);
    [Ld, Kd] = solveLQGND(Md);
    [~,~,Xfd] = rolloutLQG(Md, Ld, Kd, o.sigma_u, 0, N);
    Pd = Xfd(1:2,:);
    R(k).des.n = nhat'*Pd;   R(k).des.t = that'*Pd;
    R(k).task = tk;

    fprintf('%s  optimal AR %.3f   desired-state AR %.3f\n', upper(tk), ...
        std(R(k).opt.t)/std(R(k).opt.n), std(R(k).des.t)/std(R(k).des.n));
end

%% ---- figure ------------------------------------------------------------
W = 14; H = 10; FS = 8; LW = 1.4;
cOpt = [0.20 0.40 0.70];   cDes = [0.75 0.35 0.15];

% common scale, set from the largest SPREAD (not the absolute position)
sdT = max([std(R(1).opt.t) std(R(1).des.t) ...
           std(R(2).opt.t) std(R(2).des.t)])*100;
lim = 3.4*sdT;
fprintf('\nlargest task-irrelevant SD %.3f cm  -> axis limit +/-%.2f cm\n', sdT, lim);

close all
fh = figure('Color','w','WindowStyle','normal', ...
            'Units','centimeters','Position',[2 2 W H], ...
            'PaperUnits','centimeters','PaperPositionMode','manual', ...
            'PaperSize',[W H],'PaperPosition',[0 0 W H]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

pan = 0;
for k = 1:2
    for c = 1:2
        pan = pan + 1;
        if c==1, S = R(k).opt; col = cOpt; nm='Optimal';
        else,    S = R(k).des; col = cDes; nm='Desired-state'; end

        % deviations from this cloud's own mean
        t = (S.t - mean(S.t))*100;
        n = (S.n - mean(S.n))*100;

        ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        axis(ax,'equal'); pbaspect(ax,[1.7 1 1]);

        plot(ax, t, n, '.', 'MarkerSize', 2, 'Color', [col 0.30]);
        C = cov([t' n']);
        q = linspace(0,2*pi,200);
        E = 2*chol(C,'lower')*[cos(q); sin(q)];
        plot(ax, E(1,:), E(2,:), '-', 'Color', col, 'LineWidth', LW+0.4);

        xlim(ax,[-lim lim]); ylim(ax,[-lim/1.7 lim/1.7]);
        AR = std(S.t)/std(S.n);
        title(ax, sprintf('%s  %s, %s: AR = %.2f', ...
              char(64+pan), upper(R(k).task), nm, AR), 'FontWeight','normal');
        if k==2, xlabel(ax,'Task-irrelevant (cm)'); end
        if c==1, ylabel(ax,'Task-relevant (cm)');   end
    end
end

set(findall(fh,'-property','FontSize'),'FontSize',FS);
set(findall(fh,'-property','FontName'),'FontName','Helvetica');
set(findall(fh,'Type','axes'),'LineWidth',0.5,'TickDir','out', ...
                              'XTickLabelRotation',0,'YTickLabelRotation',0);

drawnow
arrayfun(@(a) set(a.Toolbar,'Visible','off'), findall(fh,'Type','axes'));

fh.Units='centimeters'; fprintf('size (cm): %s\n', mat2str(round(fh.Position,2)));
print(fh,'fig1_endpoint_clouds.pdf','-dpdf','-vector');
fprintf('written fig1_endpoint_clouds.pdf\n');