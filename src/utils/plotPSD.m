Uref = U_inf;
Lref = 1;
signals = {SRS_Base_Lift(7000:80000,2)/(q*2), SRS_Control(1:50632,13)/(q*2), Base_lift(:,2)/(q), R75e_7(18420:end,end-5)/(q)};
labels  = { ...
    'uRANS, control off', ...
    'uRANS, control on', ...
    'DDES, control off', ...
    'DDES, control on'};
fs=1/0.4;
plotPSDStrouhalOverlay(signals, fs, Uref, Lref, labels, 'psd_strouhal_overlay.pdf');

function plotPSDStrouhalOverlay(signals, fs, Uref, Lref, labels, exportFile)
    arguments
        signals (1,4) cell
        fs (1,1) double {mustBePositive}
        Uref (1,1) double {mustBePositive}
        Lref (1,1) double {mustBePositive}
        labels (1,4) cell
        exportFile (1,:) char = ''
    end

    % Welch settings
    nSamples = cellfun(@numel, signals);
    minN = min(nSamples);
    winLength = min(4096, floor(minN/4));
    win = hann(winLength, 'periodic');
    noverlap = round(0.5 * winLength);
    nfft = max(2^nextpow2(winLength), winLength);

    fig = figure('Color','w','Units','pixels','Position',[50 50 1100 750]);
    ax = axes(fig);
    hold(ax, 'on');

    % Line styles for publication clarity
    lineStyles = {'-','--','-.',':'};
    lineWidths = [1.8, 1.8, 1.8, 1.8];

    for k = 1:4
        x = signals{k}(:);
        x = x - mean(x, 'omitnan');

        [Pxx, f] = pwelch(x, win, noverlap, nfft, fs, 'psd');

        % Strouhal number and PSD per unit Strouhal
        St = f * Lref / Uref;
        Pxx_St = Pxx * (Uref / Lref);

        % Remove zero frequency for log plotting
        idx = St > 0 & isfinite(St) & isfinite(Pxx_St);

        semilogx(ax, St(idx), Pxx_St(idx), ...
            'LineStyle', lineStyles{k}, ...
            'LineWidth', lineWidths(k), ...
            'Color', 'k');
    end

    grid(ax, 'on');
    box(ax, 'on');
    xlim([1e-6 2])
    set(ax, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 14, ...
        'LineWidth', 1.0, ...
        'TickDir', 'in', ...
        'TickLength', [0.015 0.015], ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on');

    xlabel(ax, 'Strouhal number, $St = fL/U$', ...
        'Interpreter', 'latex', 'FontSize', 16);
    ylabel(ax, 'PSD per Strouhal unit', ...
        'Interpreter', 'latex', 'FontSize', 16);

    legend(ax, labels, ...
        'Location', 'southwest', ...
        'Interpreter', 'none', ...
        'Box', 'off');

    title(ax, 'Power spectral density comparison', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 16, ...
        'FontWeight', 'normal');

    if ~isempty(exportFile)
        exportgraphics(fig, exportFile, 'ContentType', 'vector');
    end
end