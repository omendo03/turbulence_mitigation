
% This is the most difficult one.
% It seems that we need to re-do the BER calculation.
graphics_toolkit('gnuplot')
clc;
clear variables;
close all;

% Parameters
Pt_dBm = [15, 20, 25]; % Transmitted power in dBm
Pt = db2pow(Pt_dBm) * 1e-3 % Transmitted power in Watts
data_rate_Gbps = [1.5, 3, 4]; % Data rate in Gbps
data_rate = data_rate_Gbps * 1e9; % Data rate in bps
num_channels = 3; % Number of channels assigned with PV code
weight = 2; % PV code weight
code_length = 6; % PV code length
eta_transmitter = 0.99; % Transmitter optical efficiency
eta_receiver = 0.99; % Receiver optical efficiency
area_receiver = 1; % Area of receiver aperture (m^2)
beam_concentrator_gain = 2.1025; % Beam concentrator gain
divergence_angle = 75 * 1e-3; % Beam divergence angle (rad)
responsivity_PD = 1; % Responsivity of PD (A/W)
receiver_noise_temperature = 300; % Receiver noise temperature (K)
receiver_load_resistance = 1030; % Receiver load resistance (ohm)

% Constants
k = 1.38e-23; % Boltzmann constant (J/K)
c = 3e8; % Speed of light (m/s)

% Water properties
waters = {'Pure Sea', 'Clear Ocean', 'Coastal Ocean'};
a = [0.0405, 0.114, 0.179]; % Absorption coefficient
b = [0.0025, 0.037, 0.219]; % Scattering coefficient
c_extinction = [a(1) - b(1), a(2) - b(2), a(3) - b(3)]; % Extinction coefficient

% Other parameters
num_trials = 100; % Number of Monte Carlo trials
L0 = 10; % Outer scale parameter (meters)
lambda = 0.5e-6; % Optical wavelength (meters)
alpha = 0.5; % Alpha parameter of Gamma-Gamma turbulence model
beta = 4; % Beta parameter of Gamma-Gamma turbulence model

% Define distance range
distance_range = linspace(100, 1000, 10); % Distance range from 100 to 1000 meters
%distance_range = linspace(0, 10000, 100); % Distance range from 100 to 1000 meters

% Initialize arrays to store results
snr_avg = zeros(length(data_rate), length(Pt), length(waters), length(distance_range));
ber_avg = zeros(length(data_rate), length(Pt), length(waters), length(distance_range));

for d = 1:length(data_rate)
    for p = 1:length(Pt)
        for w = 1:length(waters)
            for dist = 1:length(distance_range)
                snr_trial = zeros(1, num_trials);
                %ber_trial = zeros(1, num_trials);

                % Kolmogorov turbulence simulation
                %r0 = (0.423 * (lambda / L0)^(6/5))^-1; % Fried parameter
                r0 = (0.423*dist)^(-3/5) * (lambda/(2*pi))^(6/5); % Fried parameter
                phase_distortion = r0 * (randn(num_trials, 1) + 1i * randn(num_trials, 1)); % Phase distortions

                % Generate Rayleigh fading channels for each trial
                h1 = sqrt(0.5) * (randn(num_trials, num_channels) + 1i * randn(num_trials, num_channels));

                for trial = 1:num_trials
                    % Calculate distance-dependent variables
                    distance = distance_range(dist);
                    divergence_angle = asin((lambda / distance) / beam_concentrator_gain); % Review
                    %divergence_angle = atan((lambda / distance) )* beam_concentrator_gain; % Review

                    % Apply Kolmogorov turbulence correction
                    %phase_corrected = exp(1i * angle(h1(trial, :)) - phase_distortion(trial));
                    phase_corrected = exp(1i * angle(h1(trial, :)) - phase_distortion(trial));

                    %%% pow_corrected = abs(Pt(p)*phase_corrected) % IT does nothing of use...


                    % Calculate receiver area (m^2)
                    A_receiver = pi * (divergence_angle * area_receiver)^2;


                    % Calculate receiver noise power (W)
                    electrical_bandwidth = 0.75 * data_rate(d); % Electrical bandwidth (Hz)
                    P_noise = k * receiver_noise_temperature * electrical_bandwidth;

                    % Calculate SNR

                    snr_trial(trial) = (Pt(p) * eta_transmitter * eta_receiver * beam_concentrator_gain^2 * weight * Pt(p) * data_rate(d)) / ...
                        %(c^2 * responsivity_PD^2 * A_receiver * P_noise);
                        (c^2 * divergence_angle^2 * responsivity_PD^2 * A_receiver * P_noise);

                    % Calculate BER % log10?? Seems to help.
                    %ber_trial(trial) = 0.5 * erfc(sqrt(log10(snr_trial(trial))));
                    %ber_trial(trial) = 0.5 * erfc(sqrt(snr_trial(trial)));
                    %erf(sqrt(snr_trial(trial)))
                end

                % Average over trials
                %snr_avg(d, p, w, dist) = mean(snr_trial);
                %ber_avg(d, p, w, dist) = mean(ber_trial);
                snr_avg(d, p, w, dist) = mean(snr_trial);


                ber_avg(d, p, w, dist) = 0.5*erfc(sqrt(snr_avg(d, p, w, dist)));

            end
        end
    end
end

% Plotting SNR
figure;
for w = 1:length(waters)
    for d = 1:length(data_rate)
        for p = 1:length(Pt)
            plot(distance_range, squeeze(snr_avg(d, p, w, :)), '-*', 'linewidth', 1.5); hold on; grid on;
        end
    end
end
legend(waters);
xlabel('Distance (m)');
ylabel('SNR');
title('SNR Comparison for Different Water Types');
legend('Location', 'northwest')
legend('Pure Sea', 'Clear Ocean', 'Coastal Ocean');

size(ber_avg)
length(squeeze(snr_avg(d, p, w, :)))
length(squeeze(ber_avg(d, p, w, :)))
% Plotting BER
figure;
for w = 1:length(waters)
    for d = 1:length(data_rate)
        for p = 1:length(Pt)
            plot(squeeze(pow2db(snr_avg(d, p, w, :))), squeeze(ber_avg(d, p, w, :)), '-*', 'linewidth', 1.5); hold on; grid on;
        end
    end
end
legend(waters);
xlabel('SNR (dm)');
ylabel('BER');
%ylim([0,1e-5])
title('BER Comparison for Different Water Types');
