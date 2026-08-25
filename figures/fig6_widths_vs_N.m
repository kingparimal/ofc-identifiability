%% fig6_widths_vs_N.m
% Figure 6: how the 95% interval width scales with trial count.
% A  the three multiplicative parameters, as log-width
% B  the delay, in steps, with the quantisation floor
%
% NOTE ON UNITS. For sigma_u, sigma_s and r the interval width is a FACTOR
% (hi/lo). It is log10 of that factor -- not the factor itself -- that
% scales as 1/sqrt(N). Plotting the factor on log-log gives a slope near
% -0.1, not -0.5. The delay width is additive in steps and scales directly.
%
% Filled markers: CURVATURE-derived, from the Fisher information.
% Open markers at N = 200: GRID-derived, from the confidence region.
% Grid runs 73-94% of curvature, so curvature is a conservative bound.

clear; clc;

N   = [50 100 200 500 1000 5000];
wsu = [1.248 1.169 1.117 1.072 1.051 1.022];
wss = [1.985 1.624 1.409 1.242 1.166 1.071];
wr  = [1.770 1.497 1.330 1.198 1.136 1.059];
wd  = [3.42  2.42  1.71  1.08  0.76  0.34];      % steps

gsu = 1.11;  gss = 1.33;  gr = 1.23;             % GRID anchors, N = 200

fprintf('fitted slopes of log10(y) vs log10(N), expect -0.5:\n');
sl = @(y) polyfit(log10(N), log10(y), 1);
s1 = sl(log10(wsu));  s2 = sl(log10(wss));  s3 = sl(log10(wr));  s4 = sl(wd);
fprintf('  log-width sigma_u %+.4f\n', s1(1));
fprintf('  log-width sigma_s %+.4f\n', s2(1));
fprintf('  log-width r       %+.4f\n', s3(1));
fprintf('  delay steps       %+.4f\n', s4(1));

fprintf('\ngrid as a fraction of curvature at N = 200 (log-width):\n');
fprintf('  sigma_u %.3f   sigma_s %.3f   r %.3f\n', ...
    log10(gsu)/log10(wsu(3)), log10(gss)/log10(wss(3)), log10(gr)/log10(wr(3)));

%% ---- figure ------------------------------------------------------------
W = 13; H = 5.4; FS = 8; LW = 1.5; MS = 4;
c1 = [0.20 0.40 0.70];   c2 = [0.85 0.35 0.15];   c3 = [0.30 0.55 0.35];
cd = [0.45 0.35 0.60];

close all
fh = figure('Color','w','WindowStyle','normal', ...
            'Units','centimeters','Position',[2 2 W H], ...
            'PaperUnits','centimeters','PaperPositionMode','manual', ...
            'PaperSize',[W H],'PaperPosition',[0 0 W H]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% --- A: multiplicative parameters ---
ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on')
h1 = plot(ax1, N, log10(wsu), 'o-', 'Color', c1, 'MarkerFaceColor', c1, ...
          'LineWidth', LW, 'MarkerSize', MS);
h2 = plot(ax1, N, log10(wss), 's-', 'Color', c2, 'MarkerFaceColor', c2, ...
          'LineWidth', LW, 'MarkerSize', MS);
h3 = plot(ax1, N, log10(wr),  '^-', 'Color', c3, 'MarkerFaceColor', c3, ...
          'LineWidth', LW, 'MarkerSize', MS);

% GRID anchors: open markers, white-filled so they read as distinct
plot(ax1, 200, log10(gsu), 'o', 'Color', c1, 'MarkerFaceColor', 'w', ...
     'LineWidth', LW, 'MarkerSize', MS+2);
plot(ax1, 200, log10(gss), 's', 'Color', c2, 'MarkerFaceColor', 'w', ...
     'LineWidth', LW, 'MarkerSize', MS+2);
plot(ax1, 200, log10(gr),  '^', 'Color', c3, 'MarkerFaceColor', 'w', ...
     'LineWidth', LW, 'MarkerSize', MS+2);

% legend proxy for the grid markers
hg = plot(ax1, NaN, NaN, 'ko', 'MarkerFaceColor','w', ...
          'LineWidth', LW, 'MarkerSize', MS+2);

set(ax1,'XScale','log','YScale','log');
xlim(ax1,[40 6500]); xticks(ax1,[50 100 200 500 1000 5000]);
xlabel(ax1,'Trials {\itN}'); ylabel(ax1,'log_{10} interval width');
title(ax1,'A  Multiplicative parameters','FontWeight','normal');
lgA = legend(ax1, [h1 h2 h3 hg], ...
             {'\sigma_u','\sigma_s','{\itr}','grid, {\itN} = 200'}, ...
             'Location','southwest','Box','off');

% --- B: delay ---
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on')
plot(ax2, N, wd, 'd-', 'Color', cd, 'MarkerFaceColor', cd, ...
     'LineWidth', LW, 'MarkerSize', MS);
yline(ax2, 1, 'k--', 'LineWidth', 1);
text(ax2, 52, 1.06, 'quantisation floor', 'FontSize', FS-1.5, ...
     'VerticalAlignment','bottom', 'HorizontalAlignment','left');
set(ax2,'XScale','log','YScale','log');
xlim(ax2,[40 6500]); xticks(ax2,[50 100 200 500 1000 5000]);
ylim(ax2,[0.25 4.8]); yticks(ax2,[0.5 1 2 4]);
xlabel(ax2,'Trials {\itN}'); ylabel(ax2,'Interval width (steps)');
title(ax2,'B  Sensory delay','FontWeight','normal');

% --- fonts last ---
set(findall(fh,'-property','FontSize'),'FontSize',FS);
set(findall(fh,'-property','FontName'),'FontName','Helvetica');
set(findall(fh,'Type','axes'),'XTickLabelRotation',0,'YTickLabelRotation',0, ...
                              'LineWidth',0.5,'TickDir','out');
lgA.FontSize = FS - 1;

drawnow
fh.Units = 'centimeters';
fprintf('\nfigure size (cm): %s\n', mat2str(round(fh.Position,2)));

print(fh, 'fig6_widths_vs_N.pdf', '-dpdf', '-vector');
fprintf('written fig6_widths_vs_N.pdf\n');