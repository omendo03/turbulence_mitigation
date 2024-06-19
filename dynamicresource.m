% Figures 3 and 4 (BER modified)
clc; clear variables; close all;

#rng(3420); # Care needed with seeed.

% Distances
d1 = 500; d2 = 200;

% Power allocation coefficients
a1 = 0.75; a2 = 0.25;

N = 5e7
eta = 4; % Path loss exponent

% Rayleigh fading channels
function ray = RayleighFading(eta, d, block)
  ray = sqrt(d^-eta)*(randn(block, 1) + 1i*randn(block, 1))/sqrt(2);
end;

% Channel gains
function g_channel = CalculateG(N, d, eta)
  h_1 = 0.0;
  block_size = 1e6;
  steps = N/block_size
  for hist = 1:steps
    h_1 = h_1 + mean(abs(RayleighFading(eta, d, block_size) + RayleighFading(eta, d, block_size)).^2);
  endfor
  g_channel = h_1/steps
end;

% Transmit power
Pt = -30:5:30; % in dBm
pt = (10^-3)*db2pow(Pt); % linear scale

BW = 10^6;  % bandwidth
% Noise power
No = -174 + 10*log10(BW);   % in dBm
no = (10^-3)*db2pow(No);    % in linear scale

% Turbulence parameters
alpha = 0.5; % Alpha parameter of Gamma-Gamma turbulence model
beta = 4;    % Beta parameter of Gamma-Gamma turbulence model

% Modulation parameters (BPSK)
Eb_No = Pt - 10*log10(BW); % Energy per bit to noise power spectral density ratio (in dB)
Eb_No_linear = 10.^(Eb_No/10); % Convert to linear scale
Pe = 0.5*erfc(sqrt(Eb_No_linear)); % Probability of error for BPSK modulation

ber_theoretical = Pe; % Theoretical BER for BPSK

ber_noma = zeros(2, length(Pt));
ber_oma = zeros(2, length(Pt));
snr_noma = zeros(2, length(Pt));
snr_oma = zeros(2, length(Pt));

g_values = [CalculateG(N, d1, eta) CalculateG(N, d2, eta)]
a_values = [a1 a2];
for g = 1:columns(g_values)
  for u = 1:length(Pt)

     % Dynamic resource allocation for SNR calculations
     % For simplicity, we assume dynamic power allocation here
     dynamic_pt_noma = pt(u) * a_values(g) / (1 + alpha * pt(u) * a_values(g) * g_values(g) / no); % Dynamic power allocation for MIMO-NOMA
     dynamic_pt_oma = pt(u) / (1 + alpha * pt(u) * g_values(g) / no); % Dynamic power allocation for MIMO-OMA

     % SNR calculations with dynamic resource allocation
     snr_noma(g, u) = dynamic_pt_noma * g_values(g) / no; % SNR for MIMO-NOMA with dynamic resource allocation
     snr_oma(g, u) = dynamic_pt_oma * g_values(g) / no; % SNR for MIMO-NOMA with dynamic resource allocation

     % BER calculations (no turbulence model included for simplicity)
     ber_noma(g, u) = 0.5*erfc(sqrt((snr_noma(g, u))));%
     ber_oma(g, u) = 0.5*erfc(sqrt((snr_oma(g, u)))); % BER for MIMO-NOMA
  end
end

% Plotting
figure;
plot(Pt, snr_noma(1,:), '-*b', 'linewidth',1.5); hold on; grid on;
plot(Pt, snr_oma(1, :), '-*r','linewidth',1.5);
plot(Pt, snr_noma(2,:), '-*g', 'linewidth',1.5);
plot(Pt, snr_oma(2, :), '-*o','linewidth',1.5);
legend('Location', 'southeast')
legend('MIMO-NOMA SNR 500 m', 'MIMO-OMA SNR 500 m', 'MIMO-NOMA SNR 200 m', 'MIMO-OMA SNR 200 m');
xlabel('Transmit power (dBm)');
ylabel('Signal-to-Noise Ratio (SNR)');
title('SNR with dynamic resource allocation for optical turbulence mitigation');

figure;

plot(Pt, ber_noma(1, :), '-*b', 'linewidth',1.5); hold on; grid on;
plot(Pt, ber_oma(1, :), '-*r','linewidth',1.5);
plot(Pt, ber_noma(2, :), '-*g', 'linewidth',1.5);
plot(Pt, ber_oma(2, :), '-*o','linewidth',1.5);
legend('MIMO-NOMA BER 500 m', 'MIMO-OMA BER 500 m', 'MIMO-NOMA BER 200 m', 'MIMO-OMA BER 200 m');
xlabel('Transmit power (dBm)');
ylabel('Bit Error Rate (BER)');
title('BER with dynamic resource allocation for optical turbulence mitigation');
