function TPM = creattpm(speed_sequence)
acceleration_sequence = [diff(speed_sequence) 0];
speed_bins = floor(min(speed_sequence)):1:ceil(max(speed_sequence));
acceleration_bins = floor(min(acceleration_sequence)):0.1:ceil(max(acceleration_sequence));
% Initialize TPM matrix
TPM = zeros(length(speed_bins) - 1, length(acceleration_bins) - 1);
% Building the TPM Matrix
for i = 1:length(speed_sequence) - 1
    current_speed = speed_sequence(i);
    current_acceleration = acceleration_sequence(i);
    next_speed = speed_sequence(i + 1);
    next_acceleration = acceleration_sequence(i + 1);
    speed_idx = find(current_speed >= speed_bins(1:end-1) & current_speed < speed_bins(2:end));
    acceleration_idx = find(current_acceleration >= acceleration_bins(1:end-1) & current_acceleration < acceleration_bins(2:end));
    next_speed_idx = find(next_speed >= speed_bins(1:end-1) & next_speed < speed_bins(2:end));
    next_acceleration_idx = find(next_acceleration >= acceleration_bins(1:end-1) & next_acceleration < acceleration_bins(2:end));
    TPM(speed_idx, acceleration_idx) = TPM(speed_idx, acceleration_idx) + 1;
end
TPM = TPM ./ sum(TPM, 2);
% Draw the TPM Matrix
% imagesc(TPM);
% colorbar;
% title('Transition Probability Matrix (TPM)');
% xlabel('Acceleration');
% ylabel('Speed');
% xticks(1:length(acceleration_bins) - 1);
% xticklabels(acceleration_bins(1:end-1));
% yticks(1:length(speed_bins) - 1);
% yticklabels(speed_bins(1:end-1));
end