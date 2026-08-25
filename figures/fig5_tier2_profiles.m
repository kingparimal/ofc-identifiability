%% fig5_tier2_profiles.m
% Profile likelihood, tier 1 against tier 2, all four parameters.
%
% NOTE ON THE TIER 2 RANGE. The tier 2 grid was run on a refined range
% around the closed region -- for r it spans 0.6% of the tier 1 range, for
% sigma_u 20%. The blue curves therefore occupy only part of each panel.
% This is not truncation: every tier 2 curve terminates at dD between 17
% and 99, far above the threshold of 3.84, so each is closed with a wide
% margin. The shaded band marks the tier 2 grid extent. The refinement was
% possible only because the region closed; tier 1 could not be refined this
% way, which is itself the result.

clear; clc;
load thesis_fig_data.mat

lbl    = {'\sigma_u','\sigma_s','{\itd} (steps)','{\itr}'};
ax_1   = {su_axis, ss_axis, d_axis, r_axis};
ax_2   = {su_axis_t2, ss_axis_t2, d_axis_t2, r_axis_t2};
tru    = [0.4 0.4 4 0.002];
useLog = [false true false true];

W = 15; H = 9.5; FS = 8; LW = 1.5;
c1 = [0.75 0.35 0.15];   c2 = [0.20 0.40 0.70];
cShade = [0.20 0.40 0.70];
FLOOR = 1e-3;

close all
fh = figure('Color','w','WindowStyle','normal', ...
            'Units','centimeters','Position',[2 2 W H], ...
            'PaperUnits','centimeters','PaperPositionMode','manual', ...
            'PaperSize',[W H],'PaperPosition',[0 0 W H]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

YL = [5e-3 400];

for i = 1:4
    ax = nexttile; hold(ax,'on'); box(ax,'on')

    x1 = ax_1{i}(:);  y1 = max(prof1{i}(:) - minD1, FLOOR);
    x2 = ax_2{i}(:);  y2 = max(prof2{i}(:) - minD2, FLOOR);
    if useLog(i), x1 = log10(x1); x2 = log10(x2); xt = log10(tru(i));
    else,          xt = tru(i);  end

    % --- shaded band marking the tier 2 grid extent ---
    xb = [min(x2) max(x2) max(x2) min(x2)];
    yb = [YL(1) YL(1) YL(2) YL(2)];
    patch(ax, xb, yb, cShade, 'FaceAlpha', 0.10, 'EdgeColor','none');

    grid(ax,'on');

    if i == 3
        h1 = plot(ax, x1, y1, 'o-','Color',c1,'MarkerFaceColor',c1, ...
                  'LineWidth',LW,'MarkerSize',3.5);
        h2 = plot(ax, x2, y2, 's-','Color',c2,'MarkerFaceColor',c2, ...
                  'LineWidth',LW,'MarkerSize',3.5);
    else
        h1 = plot(ax, x1, y1, '-','Color',c1,'LineWidth',LW);
        h2 = plot(ax, x2, y2, '-','Color',c2,'LineWidth',LW);
    end

    yline(ax, chi2_1, 'k--','LineWidth',1);
    xline(ax, xt, 'k:','LineWidth',1);

    set(ax,'YScale','log'); ylim(ax, YL);
    xlim(ax,[min([x1;x2]) max([x1;x2])]);

    if useLog(i), xlabel(ax,['log_{10} ' lbl{i}]); else, xlabel(ax,lbl{i}); end
    ylabel(ax,'\Delta{\itD}');
    title(ax, sprintf('%s  %s', char(64+i), lbl{i}), 'FontWeight','normal');

    if i == 1
        lg = legend(ax,[h1 h2],{'tier 1','tier 2'}, ...
                    'Location','southeast','Box','off');
    end
end

set(findall(fh,'-property','FontSize'),'FontSize',FS);
set(findall(fh,'-property','FontName'),'FontName','Helvetica');
set(findall(fh,'Type','axes'),'LineWidth',0.5,'TickDir','out', ...
                              'XTickLabelRotation',0,'YTickLabelRotation',0);
lg.FontSize = FS - 1;

drawnow
fh.Units='centimeters'; fprintf('size (cm): %s\n', mat2str(round(fh.Position,2)));
print(fh,'fig5_tier2_profiles.pdf','-dpdf','-vector');
fprintf('written fig5_tier2_profiles.pdf\n');