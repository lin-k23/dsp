% Analyzes the signal using different Sampling Rates (Fs)
% to demonstrate the effect of aliasing.
%
clear; clc; close all;

% --- 1. Signal Parameters ---
f = [103, 107, 115];
A = [0.8, 1.0, 0.1];
N = 4096; % Fixed N for this analysis

% --- 2. Fs Cases to Test ---
Fs_list = [1000, 500, 250, 220]; % (Hz)
num_cases = length(Fs_list);

figure('Position', [100, 100, 900, 700]);
t = tiledlayout(num_cases, 1, 'TileSpacing', 'compact');
sgtitle('Effect of Sampling Rate ($F_s$) on Spectrum (Aliasing)', 'Interpreter', 'latex');

% --- 3. Loop and Process ---
for i = 1:num_cases
    Fs = Fs_list(i);
    f_nyq = Fs / 2;
    
    % --- Generate Signal ---
    t = (0:N-1)' / Fs;
    s = A(1)*sin(2*pi*f(1)*t) + A(2)*sin(2*pi*f(2)*t) + A(3)*sin(2*pi*f(3)*t);
    
    % --- Process FFT ---
    win = hann(N);
    s_win = s .* win;
    Y = fft(s_win);
    
    P2 = abs(Y(1:N/2+1));
    P_amp = P2 / N;
    P_amp(2:end-1) = 2 * P_amp(2:end-1);
    coherent_gain = sum(win) / N;
    P_amp_corrected = P_amp / coherent_gain;
    
    f_hz = (0:Fs/N:f_nyq)';
    P_db = 20*log10(P_amp_corrected + 1e-10);

    % --- Plot ---
    ax = nexttile;
    plot(f_hz, P_db, 'LineWidth', 1.5);
    grid on;
    title(sprintf('Fs = %.0f Hz  (Nyquist = %.0f Hz)', Fs, f_nyq));
    ylabel('Amplitude (dB)');
    xlim([0,f_nyq])
    
    % Add expected signal lines
    hold on;
    stem(f, 20*log10(A), 'r', 'filled', 'LineStyle', 'none', 'Marker', 'v');
    
    % Add aliased signal lines
    if f(3) > f_nyq
        f_alias = abs(f(3) - Fs); % f_alias = |115 - Fs|
        stem(f_alias, 20*log10(A(3)), 'm:', 'Marker', 'x', 'LineWidth', 2);
    end

    hold off;
end
xlabel('Frequency (Hz)');