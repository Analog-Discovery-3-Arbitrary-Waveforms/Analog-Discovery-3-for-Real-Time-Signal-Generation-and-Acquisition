%% Preliminaries & Quick Check - Devices are Connected
clc; clear; close all;
% Before using the bench generator, the following preliminaries must be
% installed:

% Instrument Control Toolbox Support Package for Keysight IO Libraries 
% and VISA Interface (free download)
url = "https://uk.mathworks.com/matlabcentral/fileexchange/52245-instrument" + ...
    "-control-toolbox-support-package-for-keysight-io-libraries-and-visa-interface";
hypertext = "Instrument Control Toolbox Support Package";
fprintf('<a href="matlab: web(''%s'') ">%s</a>\n', url, hypertext);

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

%% Arbitrary Waveform Generation & Scoping Using AD3 
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

%% AD3 Scope
F_s = 1.5e6;                                % Scope sampling frequncy (~3 multiples of AWG F_s)

% Setup DAQ scope
dq_in = daq("digilent");
addinput(dq_in,'AD3_0','ai0','Voltage'); 

dq_in.Channels.Coupling = 'DC';
dq_in.Rate = F_s;

N_period = F_s / f_base;                    % Number of samples per signal period
T = 100;                                    % Scope T periods
N_fft = round(T * N_period);                % FFT number of samples
num_avg = 10;                               % Number of averages
mag_accum = zeros(N_fft/2 + 1, 1);          % Average signal amplitude
num_spurs = 10;                             % Number of spurios harmonics detected

% FFT setup
window = flattopwin(N_fft);                 % Flattop window to avoid scalloping loss
cg = sum(window)/N_fft;
f = (0:N_fft/2).' * F_s/N_fft;

% Acquisition loop
for i = 1:(num_avg + 1)
    % Acquire data from AD3
    data = read(dq_in, N_fft, "OutputFormat", "Matrix");

    % Skip the first N_fft samples, assuming transient effects may occur in this region
    if i == 1
        continue
    end
    y = data(:,1);                          % y is measured voltage

    % FFT 
    Y = fft(y .* window, N_fft);
    mag = abs(Y/(N_fft*cg));
    mag = mag(1:N_fft/2+1);
    mag(2:end-1) = 2*mag(2:end-1);
    mag_accum = mag_accum + mag;

end

% Average
mag = mag_accum / num_avg;
mag_db = 20*log10(mag);

% Peak & Spur detection
[peak_db, peak_freq] = findpeaks(mag_db(100:end), f(100:end), ...
    'MinPeakHeight', -50, ...
    'NPeaks', num_tones);
[spur_pks, spur_locs] = findpeaks(mag_db(100:end), f(100:end), ...
    'MinPeakHeight', -90, ...
    'NPeaks', num_spurs + num_tones, ...
    'MinPeakDistance', f_base * 0.8);

spur_pks = spur_pks(num_tones + 1:end);
spur_locs = spur_locs(num_tones + 1:end);

n_spur = min(num_spurs, length(spur_pks));
spur_freq = spur_locs(1:n_spur);
spur_db = spur_pks(1:n_spur);

% total harmonic distortion (THD)
Nharm = 10;   
THD_dB = thd(y, F_s, Nharm);
fprintf('THD = %.2f dB\n', THD_dB);

% Plot
t = (0:N_fft-1)'/F_s;
figure('Position',[300 100 700 500]);

% Time-domain plot
subplot(2,1,1)
plot(t(1:end/10), y(1:end/10)); grid on;
xlabel("Time (s)"); ylabel("Amplitude (V)");

% FFT plot 
subplot(2,1,2)
plot(f, mag_db); hold on;
plot(peak_freq, peak_db,'rs','MarkerFaceColor','r')
plot(spur_freq, spur_db,'ks','MarkerFaceColor','b')
grid on; ylim([-120 60]); xlim([0 30*f_base]);
xlabel("Frequency (Hz)"); ylabel("Amplitude (dBV)");

% Labels (peaks)
for k = 1:numel(peak_freq)
    text(peak_freq(k), peak_db(k) + 3, ...
        sprintf('%.2f', peak_db(k)), ...
        'Rotation',60, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontSize',12);
end
% Labels (spurs)
for k = 1:numel(spur_freq)
    text(spur_freq(k), spur_db(k) + 3, ...
        sprintf('%.2f', spur_db(k)), ...
        'Rotation',60, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontSize',12);
end

% AD3 cleanup
clear dq_in;

%% Save Results to File
saveas(gcf, 'exp2.svg');

results = [peak_freq, 10.^(peak_db/20)];
folder = 'C:\Matlab projects\ISSC\exp 3';
filename = fullfile(folder, 'final_exp2_ad3_ad3.txt');
writematrix(results, filename, ...
    'WriteMode','append');

%% Case-study: 2-Scope Operation: AWG + Sensor
F_s = 1.5e6;                                % AD3 scope sampling frequency

% Setup DAQ scopes
dq_in = daq("digilent");
addinput(dq_in,'AD3_0','ai0','Voltage');    % Sensor
addinput(dq_in,'AD3_0','ai1','Voltage');    % AWG

dq_in.Channels(1).Coupling = 'DC';
dq_in.Channels(2).Coupling = 'DC';
dq_in.Rate = F_s;

N_period = F_s / f_base;                    % Number of samples per signal period
T = 100;                                    % Scope T periods
N_fft = round(T * N_period);                % FFT number of samples
num_avg = 10;                               % Number of averages
mag_accum = zeros(N_fft/2 + 1, 1);          % Average signal amplitude
num_spurs = 10;                             % Number of spurios harmonics detected


% FFT setup
window = flattopwin(N_fft);                 % Flattop window
cg = sum(window)/N_fft;
f = (0:N_fft/2).' * F_s/N_fft;

% Accumulated averaged magnitudes
mag_accum_sensor = zeros(N_fft/2 + 1, 1);
mag_accum_awg  = zeros(N_fft/2 + 1, 1);

% Acquisition loop
for i = 1:(num_avg + 1)
    % Acquire data from AD3
    data = read(dq_in, N_fft, "OutputFormat", "Matrix");

    % Skip the first N_fft samples, assuming transient effects may occur in this region
    if i == 1
        continue
    end
    y_sensor = data(:,1);
    y_awg = data(:,2);

    % Sensor FFT 
    Y = fft(y_sensor .* window, N_fft);
    mag = abs(Y/(N_fft*cg));
    mag = mag(1:N_fft/2+1);
    mag(2:end-1) = 2*mag(2:end-1);
    mag_accum_sensor = mag_accum_sensor + mag;

    % AWG FFT
    Y = fft(y_awg .* window, N_fft);
    mag = abs(Y/(N_fft*cg));
    mag = mag(1:N_fft/2+1);
    mag(2:end-1) = 2*mag(2:end-1);
    mag_accum_awg = mag_accum_awg + mag;
end

% Averaged magnitudes
mag_sensor = mag_accum_sensor / num_avg;
mag_awg = mag_accum_awg / num_avg;

mag_db_sensor = 20*log10(mag_sensor);
mag_db_awg = 20*log10(mag_awg);

% Peak & Spur detection (AWG)
[peak_awg_db, peak_awg_freq] = findpeaks(mag_db_awg(50:end), f(50:end), ...
    'MinPeakHeight', -50, ...
    'NPeaks', num_tones);
[spur_awg_pks, spur_awg_locs] = findpeaks(mag_db_awg(50:end), f(50:end), ...
    'MinPeakHeight', -90, ...
    'NPeaks', num_spurs + num_tones, ...
    'MinPeakDistance', f_base * 0.5);

spur_awg_pks = spur_awg_pks(num_tones + 1:end);
spur_awg_locs = spur_awg_locs(num_tones + 1:end);

n_spur = min(num_spurs, length(spur_awg_pks));
spur_awg_freq = spur_awg_locs(1:n_spur);
spur_awg_db = spur_awg_pks(1:n_spur);

% Peak & Spur detection (Sensor)
[peak_sensor_db, peak_sensor_freq] = findpeaks(mag_db_sensor(50:end), f(50:end), ...
    'MinPeakHeight', -50, ...
    'NPeaks', num_tones);
[spur_sensor_pks, spur_sensor_locs] = findpeaks(mag_db_sensor(50:end), f(50:end), ...
    'MinPeakHeight', -80, ...
    'NPeaks', num_spurs + num_tones, ...
    'MinPeakDistance', f_base * 0.5);

spur_sensor_pks = spur_sensor_pks(num_tones + 1:end);
spur_sensor_locs = spur_sensor_locs(num_tones + 1:end);

n_spur = min(num_spurs, length(spur_sensor_pks));
spur_sensor_freq = spur_sensor_locs(1:n_spur);
spur_sensor_db = spur_sensor_pks(1:n_spur);

% Plot
t = (0:N_fft-1)'/F_s;
figure('Position',[300 100 700 500]);

% AWG Time-domain 
subplot(2,1,1)
plot(t(1:end/10), y_awg(1:end/10)); grid on;
xlabel("Time (s)"); ylabel("Amplitude (V)"); ylim([-1.5 1.5]);

% AWG FFT 
subplot(2,1,2)
plot(f, mag_db_awg); hold on;
plot(peak_awg_freq, peak_awg_db,'rs','MarkerFaceColor','r')
plot(spur_awg_freq, spur_awg_db,'ks','MarkerFaceColor','b')
grid on; ylim([-120 40]); xlim([0 30*f_base]);
xlabel("Frequency (Hz)"); ylabel("Amplitude (dBV)");

% Labels (Sensor peaks/spurs)
for k = 1:numel(peak_awg_freq)
    text(peak_awg_freq(k), peak_awg_db(k) + 3, ...
        sprintf('%.2f', peak_awg_db(k)), ...
        'Rotation',60, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontSize',12);
end
for k = 1:numel(spur_awg_freq)
    text(spur_awg_freq(k), spur_awg_db(k) + 3, ...
        sprintf('%.2f', spur_awg_db(k)), ...
        'Rotation',60, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontSize',12);
end
saveas(gcf, 'exp2_awg.svg');

figure('Position',[300 100 700 500]);

% Sensor Time-domain 
subplot(2,1,1)
plot(t(1:end/10), y_sensor(1:end/10)); grid on;
% title("AWG Monitor");
xlabel("Time (s)"); ylabel("Amplitude (V)"); ylim([-1.5 1.5]);

% Sensor FFT 
subplot(2,1,2)
plot(f, mag_db_sensor); hold on;
plot(peak_sensor_freq, peak_sensor_db,'rs','MarkerFaceColor','r')
plot(spur_sensor_freq, spur_sensor_db,'ks','MarkerFaceColor','b')
grid on; ylim([-120 40]); xlim([0 30*f_base]);
xlabel("Frequency (Hz)"); ylabel("Amplitude (dBV)");

% Labels (Sensor peaks/spurs)
for k = 1:numel(peak_sensor_freq)
    text(peak_sensor_freq(k), peak_sensor_db(k) + 3, ...
        sprintf('%.2f', peak_sensor_db(k)), ...
        'Rotation',60, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontSize',12);
end
for k = 1:numel(spur_sensor_freq)
    text(spur_sensor_freq(k), spur_sensor_db(k) + 3, ...
        sprintf('%.2f', spur_sensor_db(k)), ...
        'Rotation',60, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontSize',12);
end
saveas(gcf, 'exp2_sensor.svg');

% AD3 Cleanup
clear dq_in;

%% Save Results to File
results_awg = [peak_awg_freq, 10.^(peak_awg_db/20)];
results_sensor = [peak_sensor_freq, 10.^(peak_sensor_db/20)];

folder = 'C:\Matlab projects\ISSC\exp 3';
filename = fullfile(folder, 'final_bench_ad3_8_peaks_awg.txt');
writematrix(results_awg, filename, ...
    'WriteMode','append');
filename = fullfile(folder, 'final_bench_ad3_8_peaks_sensor.txt');
writematrix(results_sensor, filename, ...
    'WriteMode','append');
