%% fig3_tier1_escape.m
% Tier 1: the accepted region does not close.
% A  (sigma_s, d) -- the dominant pair; the band crosses the full delay range
% B  (sigma_u, r) -- the region exits through the lower corner
%
% Colour is clipped at D = 30. Everything above that is rejected outright,
% so the distinction between "far above threshold" and "very far above"
% carries no information and would only compress the range that matters.

clear; clc;
load thesis_fig_data.mat

CMAX = 30;
lv   = linspace(0, CMAX, 25);

W = 13.5; H = 5.4; FS = 8;
close all
fh = figure('Color','w','WindowStyle','normal', ...
            'Units','centimeters','Position',[2 2 W H], ...
            'PaperUnits','centimeters','PaperPositionMode','manual', ...
            'PaperSize',[W H],'PaperPosition',[0 0 W H]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% --- A: (sigma_s, d) ---
ax1 = nexttile; hold(ax1,'on'); box(ax1,'on')
contourf(ax1, d_axis, log10(ss_axis), min(D_ss_d,CMAX), lv, 'LineColor','none');
contour(ax1, d_axis, log10(ss_axis), D_ss_d, [thr1 thr1], 'w-','LineWidth',2.2);
plot(ax1, 4, log10(0.4), 'w+','MarkerSize',9,'LineWidth',1.8);
clim(ax1,[0 CMAX]);
xlabel(ax1,'delay {\itd} (steps)'); ylabel(ax1,'log_{10} sensory noise \sigma_s');
title(ax1,'A  (\sigma_s, {\itd})','FontWeight','normal');
xlim(ax1,[0 10]); xticks(ax1,0:2:10);
ylim(ax1,[log10(min(ss_axis)) log10(max(ss_axis))]);

% --- B: (sigma_u, r) ---
ax2 = nexttile; hold(ax2,'on'); box(ax2,'on')
contourf(ax2, log10(r_axis), su_axis, min(D_su_r,CMAX), lv, 'LineColor','none');
contour(ax2, log10(r_axis), su_axis, D_su_r, [thr1 thr1], 'w-','LineWidth',2.2);
plot(ax2, log10(0.002), 0.4, 'w+','MarkerSize',9,'LineWidth',1.8);
clim(ax2,[0 CMAX]);
xlabel(ax2,'log_{10} effort weight {\itr}'); ylabel(ax2,'motor noise \sigma_u');
title(ax2,'B  (\sigma_u, {\itr})','FontWeight','normal');
xlim(ax2,[log10(min(r_axis)) log10(max(r_axis))]);
ylim(ax2,[min(su_axis) max(su_axis)]);

cb = colorbar(ax2);
cb.Label.String = 'discrepancy {\itD}';
cb.Ticks = [0 thr1 10 20 30];
cb.TickLabels = {'0','7.9','10','20','\geq30'};
colormap(parula);

set(findall(fh,'-property','FontSize'),'FontSize',FS);
set(findall(fh,'-property','FontName'),'FontName','Helvetica');
set(findall(fh,'Type','axes'),'LineWidth',0.5,'TickDir','out', ...
                              'XTickLabelRotation',0,'YTickLabelRotation',0);
drawnow
fh.Units='centimeters'; fprintf('size (cm): %s\n', mat2str(round(fh.Position,2)));
print(fh,'fig3_tier1_escape.pdf','-dpdf','-vector');
fprintf('written fig3_tier1_escape.pdf\n');