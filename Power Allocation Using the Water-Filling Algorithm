channel_gain = [0.8, 0.5, 0.3, 0.1];  % Channel gains
noise_power = 0.01 * ones(1, 4);  % Noise power for each channel
total_power = 10;  % Total power to be allocated
lambda_val = (total_power + sum(noise_power ./ channel_gain)) / length(channel_gain);
power_allocation = max(0, (1 / lambda_val) - (noise_power ./ channel_gain));

figure;
bar(1:length(channel_gain), power_allocation);
title('Power Allocation using Water-Filling Algorithm');
xlabel('Channel Index');
ylabel('Power Allocated');
grid on;

