%% fig4_compensation.m
% Figure 4: the compensation direction.
% A  endpoint clouds at nominal and displaced -- near identical
% B  the null direction, two independent normalisations

clear; clc;

task   = 'sim2';
theta0 = [0.4; 0.4; 4; 0.002];        % sigma_u, sigma_s, d, r
N      = 5000;

% Null direction, relative units, sign set to the escape direction
% (noise and effort down, delay up).
v = [-0.084; -0.401; +0.906; -0.105];

alpha  = 0.5 / abs(v(2));              % displacement halves sigma_s
theta1 = theta0 .* (1 + alpha*v);
theta1(3) = round(theta1(3));          % delay must be an integer

fprintf('nominal   : su %.4f  ss %.4f  d %d  r %.5f\n', ...
    theta0(1), theta0(2), theta0(3), theta0(4));
fprintf('displaced : su %.4f  ss %.4f  d %d  r %.5f\n\n', ...
    theta1(1), theta1(2), theta1(3), theta1(4));

%% ---- endpoint clouds ---------------------------------------------------
res = struct();
for k = 1:2
    if k==1, th = theta0; nm = 'nom'; else, th = theta1; nm = 'dis'; end

    opts = struct('sigma_u', th(1), 'sigma_s', th(2), 'r', th(4));
    M = augmentModel(buildModel2D(task, opts), round(th(3)));
    [L, K] = solveLQGND(M);

    rng(1);
    [~, ~, Xf] = rolloutLQG(M, L, K, th(1), 0, N);
    dp = M.Dpos(1:2)';  nh = dp/norm(dp);  tt = [-nh(2); nh(1)];
    P = Xf(1:2,:);
    res.(nm).pn = nh'*P;   res.(nm).pt = tt'*P;

    fprintf('%s: SD_rel %.4f cm  SD_irrel %.4f cm  AR %.3f  bias %.4f cm\n', ...
        nm, std(res.(nm).pn)*100, std(res.(nm).pt)*100, ...
        std(res.(nm).pt)/std(res.(nm).pn), mean(res.(nm).pn)*100);
end

fprintf('\nrelative difference: SD_rel %.2f%%, SD_irrel %.2f%%, AR %.2f%%\n', ...
    100*abs(std(res.dis.pn)/std(res.nom.pn)-1), ...
    100*abs(std(res.dis.pt)/std(res.nom.pt)-1), ...
    100*abs((std(res.dis.pt)/std(res.dis.pn))/(std(res.nom.pt)/std(res.nom.pn))-1));

%% ---- figure ------------------------------------------------------------
W = 11.5; H = 5.2; FS = 8; LW = 1.5;
cN = [0.20 0.40 0.70];   cD = [0.85 0.35 0.15];

close all
fh = figure('Color','w','WindowStyle','normal', ...
            'Units','centimeters','Position',[2 2 W H], ...
            'PaperUnits','centimeters','PaperPositionMode','manual', ...
            'PaperSize',[W H],'PaperPosition',[0 0 W H]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% --- A: endpoint clouds ---
ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
axis(ax1,'equal'); pbaspect(ax1,[1.5 1 1]);
plot(ax1, res.nom.pt*100, res.nom.pn*100, '.', 'MarkerSize', 2, 'Color', [cN 0.30]);
plot(ax1, res.dis.pt*100, res.dis.pn*100, '.', 'MarkerSize', 2, 'Color', [cD 0.30]);
hh = gobjects(1,2);
for k = 1:2
    if k==1, nm='nom'; c=cN; else, nm='dis'; c=cD; end
    C  = cov([res.(nm).pt' res.(nm).pn']*100);
    mu = [mean(res.(nm).pt); mean(res.(nm).pn)]*100;
    q  = linspace(0,2*pi,200);
    E  = 2*chol(C,'lower')*[cos(q); sin(q)] + mu;
    hh(k) = plot(ax1, E(1,:), E(2,:), '-', 'Color', c, 'LineWidth', LW);
end
lim = 3.2*max(std(res.nom.pt), std(res.dis.pt))*100;
xlim(ax1,[-lim lim]); ylim(ax1,[-lim/1.5 lim/1.5]);
xlabel(ax1,'Task-irrelevant (cm)'); ylabel(ax1,'Task-relevant (cm)');
title(ax1,'A  Endpoint scatter','FontWeight','normal');
lgA = legend(ax1, hh, {'nominal','displaced'}, 'Location','north', ...
       'Orientation','horizontal','Box','off');

% --- B: null direction, two normalisations ---
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on')
vrel = [-0.084 -0.401 +0.906 -0.105];
vlog = [-0.0969 -0.4846 +0.8618 -0.1142];
b = bar(ax2, [vrel; vlog]', 'grouped', 'LineWidth', 0.4);
b(1).FaceColor = [0.35 0.45 0.60];  b(2).FaceColor = [0.72 0.77 0.84];
yline(ax2, 0, 'k-', 'LineWidth', 0.5);
xticks(ax2,1:4); xticklabels(ax2,{'\sigma_u','\sigma_s','{\itd}','{\itr}'});
ylim(ax2,[-0.7 1.15]); ylabel(ax2,'Null component');
title(ax2,'B  Compensation direction','FontWeight','normal');
lgB = legend(ax2, {'relative','log'}, 'Location','northwest','Box','off');

% --- fonts last, then shrink the legends so they survive ---
set(findall(fh,'-property','FontSize'),'FontSize',FS);
set(findall(fh,'-property','FontName'),'FontName','Helvetica');
set(findall(fh,'Type','axes'),'XTickLabelRotation',0,'YTickLabelRotation',0, ...
                              'LineWidth',0.5,'TickDir','out');
lgA.FontSize = FS - 1.5;
lgB.FontSize = FS - 1.5;

drawnow
fh.Units = 'centimeters';
fprintf('\nfigure size (cm): %s\n', mat2str(round(fh.Position,2)));

print(fh, 'fig4_compensation.pdf', '-dpdf', '-vector');
fprintf('written fig4_compensation.pdf\n');