addpath('functions')
tremor_freq_right = 5.18; % in Hz

[tremor_assessment, left_data, right_data] = analyze_neuromod_data(self, dataBVA, ...
    'imposed_freq_tremor', tremor_freq_right);

%% Plotting

visualize_tremor_assessment_recording(left_data, right_data, tremor_assessment)