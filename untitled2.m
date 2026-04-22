clc; clear; close all;

%% ===================== FOLDERS =====================
DATAFOLD = 'E:\桌面\乔方昱课题论文\研究内容1\实验11\数据\ratio_samples';
SAVEFOLD = 'E:\桌面\乔方昱课题论文\研究内容1\实验11\图片';
if ~exist(SAVEFOLD,'dir'); mkdir(SAVEFOLD); end

%% ===================== IEEE STYLE (MATCH YOUR HEATMAP SCRIPT) =====================
STYLE.fontName    = 'Times New Roman';
STYLE.fontSize    = 8;
STYLE.figW_in     = 3.5;      % single-column
STYLE.figH_in     = 1.78;
STYLE.dpi         = 600;

%% ===================== SETTINGS =====================
snr_list = [50 45 40 35 30 25];
bus_list = 8:15;
Nb = numel(bus_list);
Ns = numel(snr_list);

% quantiles for spread
qLo = 10;
qMd = 50;
qHi = 90;

outName = sprintf('ratio_bar_Q%d_Q%d_Q%d.emf', qLo, qMd, qHi);
outPath = fullfile(SAVEFOLD, outName);

%% ===================== COLORS (YOUR RGB) =====================
rgb255 = [ ...
   207,216,227;   % SNR 50
   160,180,200;   % SNR 45
   209,221,203;   % SNR 40
   120,146,122;   % SNR 35
   238,226,196;   % SNR 30
   196,155,88];   % SNR 25
C = rgb255/255;

%% ===================== LOAD + COMPUTE METRIC =====================
Qlo = nan(Nb, Ns);
Qmd = nan(Nb, Ns);
Qhi = nan(Nb, Ns);

% 数值稳定：避免 log10(0)
eps_ratio = 1e-12;

for s = 1:Ns
    snr_db = snr_list(s);
    f = fullfile(DATAFOLD, sprintf('ratioSamples_case33bw_h3_bus8to15_SNR%d.mat', snr_db));
    assert(exist(f,'file')==2, 'File not found: %s', f);

    S = load(f);                 % expects ratioSNR: [10000 x 8]
    ratioSNR = S.ratioSNR;
    assert(size(ratioSNR,2)==Nb, 'ratioSNR must be [Nsamp x %d]', Nb);

    for b = 1:Nb
        r = ratioSNR(:,b);

        % ===== 改这里：dB 不平衡（全为正）=====
        a = abs(r);
        a = max(a, eps_ratio);
        % x = abs(20*log10(a));    % |20log10(|ratio|)|  (dB)
        x = log10(a);    % |20log10(|ratio|)|  (dB)

        x = x(isfinite(x) & ~isnan(x));

        Qlo(b,s) = prctile(x, qLo);
        Qmd(b,s) = prctile(x, qMd);
        Qhi(b,s) = prctile(x, qHi);
    end
end

%% ===================== PLOT =====================
fig = figure('Units','inches', ...
             'Position',[1 1 STYLE.figW_in STYLE.figH_in], ...
             'Color','w');

ax = axes(fig);
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
ax.FontName = STYLE.fontName;
ax.FontSize = STYLE.fontSize;
ax.XGrid = 'off';
ax.YGrid = 'on';

% ===== 更紧凑的横坐标：用 1..Nb 代替 8..15，再显示标签为 8..15 =====
xpos = 1:Nb;
hb = bar(ax, xpos, Qmd, 'grouped');   % hb: 1×Ns

% 柱子贴紧 + 颜色
for i = 1:numel(hb)
    hb(i).BarWidth  = 1.0;            % 贴紧
    hb(i).FaceColor = C(i,:);
    hb(i).EdgeColor = [0.2 0.2 0.2];
    hb(i).LineWidth = 0.4;
end

% 收紧 x 轴边界
xlim(ax, [0.5, Nb+0.5]);

% labels
hx = xlabel(ax, 'Node ID', 'FontName', STYLE.fontName, 'FontSize', STYLE.fontSize);
hx.Units = 'normalized';
hx.Position(2) = hx.Position(2);
ylabel(ax, 'Ratio fingerprint feature(log)', 'FontName', STYLE.fontName, 'FontSize', STYLE.fontSize);
ylim([-2,2]);
% ticks 显示成 bus 编号
set(ax, 'XTick', xpos, 'XTickLabel', string(bus_list));

% ===== errorbars (Qlo-Qhi), 不进 legend =====
groupWidth = min(0.95, Ns/(Ns+0.5));  % 让组内更满
for s = 1:Ns
    x = xpos - groupWidth/2 + (2*s-1) * groupWidth / (2*Ns);
    y  = Qmd(:,s);
    lo = y - Qlo(:,s);
    hi = Qhi(:,s) - y;

    eb = errorbar(ax, x, y, lo, hi, 'k', ...
        'LineStyle','none', 'LineWidth', 0.8, 'CapSize', 2);
    eb.HandleVisibility = 'off';
end
%% ===================== INSET: zoom nodes 10-12 =====================
zoomBus = [10 11 12];
idxZoom = ismember(bus_list, zoomBus);   % 对应 bus 10,11,12
xZoomMain = xpos(idxZoom);               % 在主图中的真实位置 = [3 4 5]

% 取局部范围
Qlo_z = Qlo(idxZoom,:);
Qmd_z = Qmd(idxZoom,:);
Qhi_z = Qhi(idxZoom,:);

% -------- 主图上画一个框，标出被放大的区域 --------
yMinBox = min(Qlo_z(:));
yMaxBox = max(Qhi_z(:));
yPad    = 0.05*(yMaxBox - yMinBox + eps);

rectangle(ax, ...
    'Position', [min(xZoomMain)-0.5, yMinBox-yPad-0.30, ...
                 numel(xZoomMain), (yMaxBox-yMinBox)+2*yPad+0.75], ...
    'EdgeColor', [0.80 0.20 0.20] , ...
    'LineStyle', '--', ...
    'LineWidth', 1);

% -------- 右上角插入小图 --------
ax_in = axes('Position',[0.14 0.23 0.30 0.28]);   % [左 下 宽 高]
hold(ax_in,'on'); box(ax_in,'on');

hb_in = bar(ax_in, 1:sum(idxZoom), Qmd_z, 'grouped');

for i = 1:numel(hb_in)
    hb_in(i).BarWidth  = 1.0;
    hb_in(i).FaceColor = C(i,:);
    hb_in(i).EdgeColor = [0.2 0.2 0.2];
    hb_in(i).LineWidth = 0.4;
end

% 误差棒
groupWidth_in = min(0.95, Ns/(Ns+0.5));
for s = 1:Ns
    x_in = (1:sum(idxZoom)) - groupWidth_in/2 + (2*s-1)*groupWidth_in/(2*Ns);
    y_in  = Qmd_z(:,s);
    lo_in = y_in - Qlo_z(:,s);
    hi_in = Qhi_z(:,s) - y_in;

    eb_in = errorbar(ax_in, x_in, y_in, lo_in, hi_in, 'k', ...
        'LineStyle','none', 'LineWidth', 0.6, 'CapSize', 2);
    eb_in.HandleVisibility = 'off';
end

% inset 样式
ax_in.FontName = STYLE.fontName;
ax_in.FontSize = 7;
ax_in.XGrid = 'off';
ax_in.YGrid = 'on';

xlim(ax_in, [0.5, sum(idxZoom)+0.5]);
ylim(ax_in, [yMinBox-yPad, yMaxBox+yPad]);
set(ax_in, 'XTick', 1:sum(idxZoom), 'XTickLabel', string(zoomBus));

% title(ax_in, 'Zoom: Nodes 10-12', ...
%     'FontName', STYLE.fontName, 'FontSize', 7);

% 可选：把 inset 的 y 轴放右边，更紧凑
ax_in.YAxisLocation = 'left';
annotation(fig, 'arrow', ...
    [0.38 0.36], ...
    [0.54 0.46], ...
    'LineWidth', 0.6, ...
    'HeadLength', 6, ...
    'HeadWidth', 6, ...
    'Color', 'k');
%% ===================== LEGEND (2 ROWS, SHORT TOKENS, BOTTOM CENTER) =====================
snrNames = arrayfun(@(x)sprintf('%d dB', x), snr_list, 'UniformOutput', false);

lgd = legend(ax, hb, snrNames, ...
    'Location','south', ...
    'Orientation','horizontal', ...
    'Box','off');

lgd.FontName = STYLE.fontName;
lgd.FontSize = STYLE.fontSize;

% 两排：6项 -> 3列
try
    lgd.NumColumns = 6;
end

% 缩短 legend 色块长度
try
    lgd.ItemTokenSize = [6 6];
catch
end

% 手动控制 legend 长度和位置
lgd.Units = 'normalized';
legendWidth  = 0.78;
legendHeight = 0.12;
lgd.Position = [0.54-legendWidth/2, 0.85, legendWidth, legendHeight];

% 让坐标轴更紧凑
ax.Units = 'normalized';
ax.Position = [0.07, 0.15, 0.9, 0.82];

%% ===================== EXPORT =====================
exportgraphics(fig, outPath, 'ContentType','vector', 'Resolution', STYLE.dpi);
exportgraphics(fig, strrep(outPath,'.emf','.pdf'), 'ContentType','vector', 'Resolution', STYLE.dpi);

fprintf('[SAVED] %s\n', outPath);
