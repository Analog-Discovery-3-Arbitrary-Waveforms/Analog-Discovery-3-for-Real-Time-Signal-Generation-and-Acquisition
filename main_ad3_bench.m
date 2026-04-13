%% Preliminaries & Quick Check - Devices are Connected
clc; clear; close all;
% Before using the bench generator, the following preliminaries must be
% installed:

% 1. Instrument Control Toolbox Support Package for Keysight IO Libraries 
% and VISA Interface (free download)
url = "https://uk.mathworks.com/matlabcentral/fileexchange/52245-instrument" + ...
    "-control-toolbox-support-package-for-keysight-io-libraries-and-visa-interface";
hypertext = "Instrument Control Toolbox Support Package";
fprintf('<a href="matlab: web(''%s'') ">%s</a>\n', url, hypertext);

% 2. NI VISA drivers 
% https://www.keysight.com/us/en/lib/software-detail/computer-software/io-libraries-suite-downloads-2175637.html
url = "https://www.ni.com/en/support/downloads/drivers/download.ni-visa.html#585834";
hypertext = "NI-VISA";
fprintf('<a href="matlab: web(''%s'') ">%s</a>\n', url, hypertext);

% Check if the bench generator is detected:
visadevlist;

% Make sure the path with the digilent toolbox has been added: 
% addpath(['C:\...\digilent'])
addpath(['C:\Users\Kseniia\Dropbox\Apps\Overleaf\(Submitted) ' ...
    'I2MTC 2026 Closed-loop Helmholtz Coil for the Characterisation ' ...
    'of Sensors for Electromagnetic Tracking\Matlab\I2MTC journal\digilent'])
daqlist("digilent")

% Preamble
msg = "Can be helpful to use this page:";
url = "https://uk.mathworks.com/help/daq/transition-your-code-from-session-to-dataacquisition-interface.html";
hypertext = "Matlab implemented functions for Digilent Analog Discovery 3";
fprintf('%s <a href="matlab: web(''%s'') ">%s</a>\n', msg, url, hypertext);

%% AD3 Generator
daqreset; clc; close all;

F_s = 0.5e6;                            % Sampling frequency [Hz]
f_base = 1e3;                           % Base frequency [Hz]
V = 0.1;                                % Signal amplitude [V]
N_awg =  F_s/2;                         % Number of samples (waveform generator)
num_tones = 5;                          % Number of tones

% Arbitrary Waveform Generation
amp = V * ones(num_tones,1);            % Per tone amplitudes
disp(table(amp(:)*1000, 20*log10(amp(:)), ...
    'VariableNames', {'Amplitudes_mV','Amplitudes_dBV'}))

n = (0:N_awg-1).';
waveform = zeros(size(n));

for k = 1:num_tones
    waveform = waveform + amp(k) .* sin(2*pi*k*f_base*n/F_s);
end
waveform = waveform(:);

% Setup DAQ generator
dq_out = daq("digilent");
addoutput(dq_out,"AD3_0","ao0","Voltage");
dq_out.Rate = F_s;

% Preload and run continuously
preload(dq_out, waveform);
start(dq_out, "repeatoutput");

%% Test: Plot Acquired Waveform 
addpath('C:\Matlab projects\exp 1\ad3_bench');
clc; close all;

filename = 'scope_5m_100mV.csv';   
data = readmatrix(filename);
t = data(:,1);                      % time (s)
y = data(:,2);                      % voltage (V)
plot(t,y)

%% Save the Acquired (Bench Scope)
% a file from the benchtop scope was saved using USB stick, and should be
% reconstructed using thsi code:

addpath('C:\Matlab projects\ISSC');
folder = 'C:\Matlab projects\ISSC';

files = "scope_75.csv";
num_tones = 5;

for i = 1:numel(files)

    % Load data
    data = readmatrix(files(i));                    
    t = data(:,1);                              % time axis
    y = data(:,2);                              % voltage axis

    % Remove NaNs from the data
    valid = ~isnan(t) & ~isnan(y) & ~isinf(t) & ~isinf(y);
    t = t(valid);
    y = y(valid);

    dt = diff(t);
    dt = dt(dt > 0);
    F_s = 1/median(diff(t));

    N = 2^floor(log2(length(y)));   
    y = y(1:N);

    % FFT 
    w = flattopwin(N);                          % Flattop window
    cg = sum(w)/N;
    Y = fft(y .* w);
    mag = abs(Y/(N*cg));
    mag = mag(1:N/2+1);
    mag(2:end-1) = 2*mag(2:end-1);
    mag_db = 20*log10(mag);

    % Frequency axis
    f = (0:N/2)' * F_s/N;

    % Skip DC
    f_min = 500;
    valid_idx = f > f_min;

    % Peak detection
    [pks, locs] = findpeaks(mag_db(valid_idx), f(valid_idx), ...
        'MinPeakDistance', 0.8*f_base, ...
        'MinPeakProminence', 6, ...
        'SortStr', 'descend');

    % Keep strongest tones
    n_keep = min(num_tones, numel(pks));
    pks = pks(1:n_keep);
    locs = locs(1:n_keep);

    % Convert to linear
    peak_freq = locs(:);
    peak_amp  = 10.^(pks(:)/20);
    result = [peak_freq, peak_amp];

    % Save
    [~, name, ~] = fileparts(files(i));
    filename = fullfile(folder, name + "_peaks.txt");
    writematrix(result, filename);

    % Plot
    figure('Position',[300 100 700 500]);

    % Time-domain
    subplot(2,1,1)
    idx = 1:floor(length(y)/10);
    plot(t(idx), y(idx)); grid on;
    xlabel("Time (s)"); ylabel("Voltage [V]");

    % FFT
    subplot(2,1,2)
    plot(f, mag_db); hold on;
    plot(peak_freq, pks, 'ro','MarkerFaceColor','r')
    grid on;
    ylim([-120 60]);
    xlim([0 F_s/2]);
    xlabel("Frequency (Hz)"); ylabel("Amplitude (dBV)");

end