clear; clc; close all;
addpath('functions')

load('data_treatment_local.mat')
patient_number = 3; %[1:4];

band_power = [4 12]; % [4 8] Hz
folder_path = data_treatment.treatment_neuromod_folder{patient_number};

tremor_freq_right = 5.36; % in [Hz]

%% loading one random .mat file from neuromod assessment in MR with accelerometer:

ls = dir(folder_path);
table_ls = struct2table(ls);
ls = ls(find(table_ls.isdir == 0));
k = ceil(numel(ls)*rand);

filename = ls(k).name;

load(ls(k).name, 'dataBVA', 'self'); % important to understand all the fields in self.
% Closing figure
fig_GUI = findall(groot, 'Type', 'figure', 'Name', 'neuromodulationGUI');
close(fig_GUI)

fprintf('Filename loaded: %s\n', filename)

%% Analyzing data:

[tremor_assessment, left_data, right_data] = analyze_neuromod_data(self, dataBVA, ...
    'imposed_freq_tremor', tremor_freq_right);

%% Plotting

visualize_tremor_assessment_recording(left_data, right_data, tremor_assessment, ...
    'filename', filename)



