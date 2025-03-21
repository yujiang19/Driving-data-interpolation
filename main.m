% data_trip10s [6:velocity 8:volt 9:current 10:SOC 13: pedal 16:acceleration 17:distance]
% data_trip1s [6:velocity 8:volt 9:current 10:SOC 16:acceleration 17:distance]
clc; clear; close all;
folder_path = fullfile(pwd, 'example_data');
file_info = dir(fullfile(folder_path, '*.mat'));
file_names = {file_info.name};
load (fullfile(pwd, 'example_inpdata'))
for ik = 1:length(D)
    %% Imputation Model
    tStart = tic;
    pro_data = augmentf(D);
    elapsedTime = toc(tStart);
    fprintf('Elapsed time: %.6f seconds\n', elapsedTime);
end

%% validation
aff = valdationpolt(pro_data);
aff.calfea % Calculate characteristic parameter 
aff.vplot % Vehicle speed calibration visualization
aff.socplot % Battery data visualization
aff.vaplot % calculate velocity-acceleration distribution
aff.energyplot % calculate energy consumption




