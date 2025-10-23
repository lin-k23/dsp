% --- analysis_windows.m ---
%
% Analyzes the signal using four different window functions to
% compare the trade-off between resolution and dynamic range.
%
clear; clc; close all;

% --- 1. Signal & System Parameters ---
Fs = 500;           % Sampling Rate (Hz)
N = 16384;          % FFT Points (2^14), for good resolution
t = (0:N-1)' / Fs;  % Time vector (N x 1)

% Signal components
f = [103, 107, 115];
A = [0.8, 1.0, 0.1];
s = A(1)*sin(2*pi*f(1)*t) + A(2)*sin(2*pi*f(2)*t) + A(3)*sin(2*pi*f(3)*t);

% --- 2. Window Definitions ---
win_rect = rectwin(N);
win_hann = hann(N);
win_hamm = hamming(N);
win_black = blackman(N);

% --- 3. Process FFT for each window ---
[f_hz, P_rect_db] = process_fft(s, win_rect, Fs, N);
[~, P_hann_db]    = process_fft(s, win_hann, Fs, N);
[~, P_hamm_db]    = process_fft(s, win_hamm, Fs, N);
[~, P_black_db]   = process_fft(s, win_black, Fs, N);

% --- 4. Plotting ---
figure('Position', [100, 100, 900, 700]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle(sprintf('Window Function Comparison (N = %d, Fs = %.0f Hz)', N, Fs), 'FontWeight', 'bold');

% --- Plot 1: Linear Scale (Zoomed) ---
% This plot shows the MAIN-LOBE width (Resolution)
ax1 = nexttile;
plot(f_hz, 10.^(P_rect_db/20), 'LineWidth', 1); hold on;
plot(f_hz, 10.^(P_hann_db/20), 'LineWidth', 1.5);
plot(f_hz, 10.^(P_hamm_db/20), 'LineWidth', 1.5);
plot(f_hz, 10.^(P_black_db/20), 'LineWidth', 1.5);
hold off; grid on;
title('Linear Scale (Peak Detail / Resolution)');
ylabel('Amplitude');
legend('Rectangular', 'Hann', 'Hamming', 'Blackman', 'Location', 'northeast');
xlim([100, 120]);
ylim([0, 1.2]);

% --- Plot 2: Log (dB) Scale (Full View) ---
% This plot shows the SIDE-LOBE level (Dynamic Range)
ax2 = nexttile;
plot(f_hz, P_rect_db, 'LineWidth', 1); hold on;
plot(f_hz, P_hann_db, 'LineWidth', 1.5);
plot(f_hz, P_hamm_db, 'LineWidth', 1.5);
plot(f_hz, P_black_db, 'LineWidth', 1.5);
hold off; grid on;
title('Log (dB) Scale (Dynamic Range / Sidelobes)');
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
xlim([0, Fs/2]);
ylim([-120, 10]); % Set dB floor

% --- Local Helper Function ---
function [f, P_db] = process_fft(s, win, Fs, N)
    % 1. Apply window
    s_win = s .* win;
    
    % 2. Calculate FFT
    Y = fft(s_win);
    
    % 3. Get single-sided spectrum and scale by N
    P_amp = abs(Y(1:N/2+1)) / N;
    P_amp(2:end-1) = 2 * P_amp(2:end-1);
    
    % 4. Correct for the window's "coherent gain" to get correct amplitude
    coherent_gain = sum(win) / N;
    P_amp_corrected = P_amp / coherent_gain;
    
    % 5. Convert to dB (with epsilon for log(0))
    P_db = 20*log10(P_amp_corrected + 1e-10);
    f = (0:Fs/N:Fs/2)';
end
