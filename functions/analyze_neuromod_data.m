function [tremor_assessment, left_data, right_data] = analyze_neuromod_data(self, dataBVA, varargin)
    % DESCRIPTION: 
    %   Importing .mat file recorded during the treatment on the PRISMA
    %   stim computer. The .mat file is always structured in the same
    %   fashion: self + neuromod object.
    %
    % WARNING: the fields of neuromod object (self) will change as the
    % script evolves. This function needs to be compatible with all the
    % .mat files saved.
    %
    % SYNTAX:
    %   [tremor_assessment, left_data, right_data] = ...
    %               import_neuromod_file(filename, varargin)
    %
    % INPUTS:
    % - filename : path of the .mat file containing data
    %
    % OPTIONAL INPUTS:
    %
    % OUTPUTS:
    % - tremor assessment: structure containing the following fields:
    %    > tremor_assessment.time_brain_vision
    %    > tremor_assessment.time_tiepie
    %    > tremor_assessment.data_tiepie
    %    > tremor_assessment.insightec_trig_BVA
    %    > tremor_assessment.recording
    %    > tremor_assessment.has_posture_marker
    %    > tremor_assessment.sonication_time_tiepie
    %    > tremor_assessment.sonication_time_brain_vision
    %    > tremor_assessment.posture_interval
    %    > tremor_assessment.filename
    % 
    % - left_data / right data: struc containing the following fields:
    %    > left_data.acceleration
    %    > left_data.spectrum
    %    > left_data.freq_spectrum 
    %    > left_data.time_spectrum 
    %    > left_data.indicators
    %
    %
    % Created March 26, 2021
    % Copyright, Thomas Bancel
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Default parameters %
    %%%%%%%%%%%%%%%%%%%%%%
    
    plot_figure = false;
    threshold_tiepie_insightec = 2;
    threshold_brain_vision_insightec = 200;
    
    normalization = 1450; % 1450mV = 1g = 9,806 m/s2
    g = 9.806; % in [m/s^2]
    band_power = [4 12];
    one_hand = false;
    
    if ~isempty(varargin)
	    for input_index = 1:2:length(varargin)
	        switch varargin{input_index}
                case 'imposed_freq_tremor'
                    imposed_freq_tremor = varargin{input_index + 1};
                case 'plot_figure'
                    plot_figure = varargin{input_index + 1};
                case 'band_power'
                    band_power = varargin{input_index + 1};
            end
        end
    end
    
    if isempty(dataBVA)
        left_data = [];
        right_data = [];
        tremor_assessment = [];
        return
    end 
    
    % Analysing timestamp in filename
    % date_string = filename(1:15);
    % tremor_assessment.start_new_file = datetime(date_string, 'InputFormat', 'yyyyMMdd''T''HHmmss');
    
    % data tiepie
    fs_tiepie = 10e3;
    tiepie = self.scp;
    data_tiepie = self.data_tiepie;
    time_tiepie = (1:size(data_tiepie,1))/fs_tiepie;

    % data brain vision (matlab)
    fs_brain_vision = 1/self.RDA.props.samplingInterval*10^6; % in [Hz]
    time_brain_vision = (1:size(dataBVA,1))/fs_brain_vision;


    channel_names = self.RDA.props.channelNames;
    idx_left = find(contains(channel_names, 'L_ACC_X'));
    idx_right = find(contains(channel_names, 'R_ACC_X')); %strfind(data.label, 'L_ACC_X'); 
    
    if isempty(idx_left)
        one_hand = true;
        idx_right = find(contains(channel_names, 'ACC_X'));
        idx_left = find(contains(channel_names, 'ACC_X'));
    end

    time_recording = time_brain_vision;
    recording = dataBVA;
    raw_acceleration_left = sqrt(sum(dataBVA(:,idx_left:(idx_left+2)).^2, 2)); % raw acceleration
    raw_acceleration_right = sqrt(sum(dataBVA(:,idx_right:(idx_right+2)).^2, 2)); % raw acceleration
    
    left_acceleration = g*raw_acceleration_left/normalization;
    right_acceleration = g*raw_acceleration_right/normalization;

    insightec_trig = dataBVA(:,1);

    if plot_figure
        figure;
        plot(time_brain_vision, dataBVA(:,1)); hold on
        plot(time_tiepie, data_tiepie);
        xlabel('Time (s)'); ylabel('Insightec trig');
        legend('Stream Brain Vision', 'Scope (chan 1)', 'Scope (chan 2)');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Find sonication times %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    sonication_time_tiepie = [];
    if ~isempty(data_tiepie)
        % find tiepie channel containing insightec trig.
        [~, idx_max] = max([max(data_tiepie(:,1)), max(data_tiepie(:,2))]);
        % h = histogram(data_tiepie(:, idx_max));
        
        % find sonication times using tiepie data:
        d = diff(data_tiepie(:, idx_max));
        % limit = 15*std(d);
        limit = threshold_tiepie_insightec; % using absolute value based on data analysis
        a = min(find(d > limit))/fs_tiepie;
        b = max(find(d < -limit))/fs_tiepie;

        if ~isempty(a)
            sonication_time_tiepie = [a b];
        end
    end

    sonication_time_brain_vision = [];
    if ~isempty(insightec_trig)

        d = diff(insightec_trig);
        % limit = 15*std(d);
        limit = threshold_brain_vision_insightec;
        a = min(find(d > limit))/fs_brain_vision;
        b = max(find(d < -limit))/fs_brain_vision;

        if ~isempty(a)
            sonication_time_brain_vision = [a b];
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find posture interval using keyLogger data: %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    kl_data = self.KeyLogger.Data;
    posture_times = cell2mat(kl_data(find(strcmp(kl_data(:,1), 'b')),4));
    relax_times = cell2mat(kl_data(find(strcmp(kl_data(:,1), 'y')),4));

    if isempty(posture_times) || isempty(relax_times)
        has_posture_marker = 0;
    else
        has_posture_marker = 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % setting default posture interval if %
    % relax_times / posture_times are not %
    % working                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isempty(relax_times)
        relax_times = max(time_brain_vision);
    end
    if isempty(posture_times)
        posture_times = min(time_brain_vision);
    end

    posture_interval = [min(posture_times) min(relax_times)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analyzing left and right data %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('==> Left tremor assessment \n')
    [left_acceleration, ~, left_spectrum, left_freq_spectrum, left_t_spectrum, ...
        left_indicators] = compute_tremor_indicators(left_acceleration, fs_brain_vision, ...
        'time_interval', posture_interval, 'plot_figure', plot_figure, ...
        'imposed_freq_tremor', imposed_freq_tremor, 'band_power', band_power);

    fprintf('==> Right tremor assessment \n')
    [right_acceleration, ~, right_spectrum, right_freq_spectrum, right_t_spectrum, ...
        right_indicators] = compute_tremor_indicators(right_acceleration, fs_brain_vision, ...
        'time_interval', posture_interval, 'plot_figure', plot_figure, ...
        'imposed_freq_tremor', imposed_freq_tremor, 'band_power', band_power);

    %%%%%%%%%%%%%%%%%%%
    % Setting outputs %
    %%%%%%%%%%%%%%%%%%%
    
    tremor_assessment.time_brain_vision = time_brain_vision;
    tremor_assessment.time_tiepie = time_tiepie;
    
    tremor_assessment.only_one_hand = one_hand;

    tremor_assessment.data_tiepie = data_tiepie;
    tremor_assessment.insightec_trig_BVA = insightec_trig;
    tremor_assessment.recording = dataBVA; 
    tremor_assessment.has_posture_marker = has_posture_marker;
    tremor_assessment.sonication_time_tiepie = sonication_time_tiepie;
    tremor_assessment.sonication_time_brain_vision = sonication_time_brain_vision;
    tremor_assessment.posture_interval  = posture_interval;
    % tremor_assessment.filename = filename;

    left_data.acceleration = left_acceleration;
    left_data.spectrum = left_spectrum;
    left_data.freq_spectrum = left_freq_spectrum;
    left_data.time_spectrum = left_t_spectrum;
    left_data.indicators = left_indicators;

    right_data.acceleration = right_acceleration;
    right_data.spectrum = right_spectrum;
    right_data.freq_spectrum = right_freq_spectrum;
    right_data.time_spectrum = right_t_spectrum;
    right_data.indicators = right_indicators;

    % For setting the vertical axes in the plot later and having comparable
    % axis
    tremor_assessment.max_power_band = max([max(left_indicators.power_in_band_over_time), ...
        max(right_indicators.power_in_band_over_time)]);

    % for having comparable c_axis for spectrum
    max_left_spectrum = max(abs(left_spectrum(:)));
    min_left_spectrum = min(abs(left_spectrum(:)));
    left_data.spectrum_limit = [min_left_spectrum max_left_spectrum];

    max_right_spectrum = max(abs(right_spectrum(:)));
    min_right_spectrum = min(abs(right_spectrum(:)));
    right_data.spectrum_limit = [min_right_spectrum max_right_spectrum];

end