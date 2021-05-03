function [acceleration, time, spectrum, freq_spectrum, t_spectrum, indicators] = ...
    compute_tremor_indicators(acceleration, fs, varargin)
    % Compute tremor indicators based on acceleration and sampling
    % frequency.
    %
    % SYNTAX:
    %   [acceleration, time, spectrum, freq_spectrum, t_spectrum, indicators] = ...
    %       compute_tremor_indicators(raw_acceleration, fs, varargin)
    % 
    % INPUTS:
    % - raw_acceleration : N x 3 time-vector, where N is the number of data
    % points. Each column is the acceleration in one dimension.
    % - fs : sampling frequency 
    % 
    % OPTIONAL INPUTS:
    % - time_interval : Q x 2 time interval : each row corresponds to a 
    % [min_t max_t] interval in [s] over which the analysis is performed.
    % Default is the whole dataset
    % - freq_tremor : in [Hz] imposed tremor frequency
    % - band_power : [min_f max_f] in [Hz] tremor band where analysis will
    % be performed.
    %
    % OUTPUTS:
    % - acceleration : N x 1 time-vector in [m/s^2], determined with g and
    % normalization value from accelerometer
    % - time : 1 x N time in [s]
    % 
    % - t_spectrum: P x 1 in [s] timesteps at which the spectrum is
    % computed. It is computed with a default 2s window and 1.5s overlap.
    % - freq_spectrum : M x 1 in [Hz] frequencies over which the spectrum
    % is determined in the 2s window. The spectrum is determined over all 
    % fft frequencies and then restricted to [0 50] Hz. The frequency step
    % depends on the number of time values.
    %
    % - spectrum : [M x P] complex double. The spectrum is normalized over 
    % delta_f and fs to be comparable between recordings.
    %
    % - indicators: struct containing different fields:
    %    > indicators.freq_tremor : identified tremor frequency in [Hz]
    %    during the time interval of analysis
    %    > indicators.tremor_width : width of the tremor peak in [Hz]
    %    > indicators.tremor_amplitude : amplitude of tremor in [mm]
    %    > indicators.imposed_freq_tremor : imposed frequency (input)
    %    > indicators.tremor_amplitude_imposed_freq : amplitude at this
    %      specific frequency
    %    > indicators.band : band in [Hz] over which the power is computed
    %    > indicators.power_in_band_over_time : power inside the band
    %    during the time interval
    %    > indicators.mean_power_band : mean power inside the band during
    %    the whole time interval
    
    % WARNING : need to change the organization of the output to be
    % consistent with the fact that there can be several posture intervals.
    
    % Default parameters:
    
    normalization = 1450; % 1450mV = 1g = 9,806 m/s2
    g = 9.806; % in [m/s^2]
    
    time_interval = [];
    max_freq = 30; % in [Hz]
    
    window_time = 2; % in [seconds]
    window_overlap = 1.5; % in [seconds]
    
    min_peak_interval = 3; % in Hz
    
    band_power = [4 12];
    imposed_freq_tremor = [];
    plot_figure = true;
    figure_name = [];
    
    % params for detecting tremor peaks
    tremor_frequency_search_interval = [3 10];
    
    if ~isempty(varargin)
	    for input_index = 1:2:length(varargin)
	        switch varargin{input_index}
                case 'imposed_freq_tremor'
                    imposed_freq_tremor = varargin{input_index + 1};
                case 'time_interval'
                    time_interval = varargin{input_index + 1};
                case 'plot_figure'
                    plot_figure = varargin{input_index + 1};
                case 'figure_name'
                    figure_name = varargin{input_index + 1};
                case 'band_power'
                    band_power = varargin{input_index + 1};
                case 'max_freq'
                    max_freq = varargin{input_index + 1};
                case 'tremor_frequency_search_interval'
                    tremor_frequency_search_interval = varargin{input_index + 1};
            end
        end
    end
   
    % Computing real acceleration (in m/s^2)
    % acceleration = sqrt(sum(raw_acceleration.^2, 2))*g/normalization;
    time = (1:numel(acceleration))/fs;
    
    % Computing spectrogram
    n_window = window_time*fs;
    n_overlap = window_overlap*fs;

    n_fft = 3*n_window;
    freq_spectro = (0:n_fft-1)*fs/(n_fft - 1);
    idx_max_freq = max(find(freq_spectro  < max_freq));
    
    [spectrum, freq_spectrum, t_spectrum] = spectrogram(acceleration, ...
        n_window, n_overlap, freq_spectro(1:idx_max_freq), fs);
    
    % normalizing spectrum:
    delta_f = mean(diff(freq_spectrum));
    spectrum = spectrum*delta_f/fs; % to validate with JF
    
    % Computing mean spectrum in the interesting time interval
    if isempty(time_interval)
        time_interval = [min(time), max(time)];
    end
    
    % Power in specified band
    idx_freq_band = find(freq_spectrum > band_power(1) & freq_spectrum < band_power(2));
    tremor_amplitude_imposed_f_tremor = [];
    
    % Computing data during all posture intervals
    for k = 1:size(time_interval,1)
        
        fprintf('Computing indicators during posture %i (between %2.0fs and  %2.0fs)\n', ...
            [k time_interval(k,1) time_interval(k,2)]);
        
        idx_t_spectro = find(t_spectrum < time_interval(k,2) & t_spectrum > time_interval(k,1));

        % moyenne du spectre pendant le temps d'analyse
        mean_spectrum{k} = mean(abs(spectrum(:, idx_t_spectro)), 2); 

        % Determing tremor frequency (based on peak analysis)
        [peaks{k}, locs{k}, peaks_width] = findpeaks(mean_spectrum{k}, 'MinPeakDistance', floor(min_peak_interval/delta_f));
        
        freq_peaks = freq_spectrum(locs{k});
        [~, idx_sort] = sort(peaks{k}, 'descend');
        
        for jj = 1:numel(idx_sort)
            if freq_peaks(idx_sort(jj)) < tremor_frequency_search_interval(2) && ...
                    freq_peaks(idx_sort(jj)) > tremor_frequency_search_interval(1)
                freq_tremor(k) = freq_peaks(idx_sort(jj));
                idx_freq_tremor = idx_sort(jj);
                amp_tremor_spectrum = peaks{k}(idx_freq_tremor);
                break
            else
                continue
            end       
        end
        
        % [amp_tremor_spectrum, idx_freq_tremor] = max(peaks{k});
        freq_tremor(k) = freq_spectrum(locs{k}(idx_freq_tremor)); % in [Hz]

        % Computing tremor indicators for this determined frequency
        tremor_width(k) = peaks_width(idx_freq_tremor)*delta_f; % in [Hz]
        tremor_amplitude(k) = amp_tremor_spectrum/freq_tremor(k)^2; % in [m]  

        % Computing tremor indicators at imposed frequency
        
        if ~isempty(imposed_freq_tremor)
            idx_imposed_f_tremor = min(find(freq_spectrum > imposed_freq_tremor));
            tremor_amplitude_imposed_f_tremor(k) = mean_spectrum{k}(idx_imposed_f_tremor)/imposed_freq_tremor^2;
        else
            tremor_amplitude_imposed_f_tremor(k) = NaN;
        end

        mean_power_band(k) = sum(abs(mean_spectrum{k}(idx_freq_band)).^2); % in (m/s^2)^2
        
    end

    power_in_band_over_time = sum(abs(spectrum(idx_freq_band,:)).^2, 1);
    
    fprintf("-----------------------------\n");
    fprintf("Mean tremor frequency (peak detection) : %2.2f Hz \n", mean(freq_tremor));
    fprintf("Mean tremor amplitude (peak detection) : %2.2f mm \n", mean(tremor_amplitude)*10^(3));
    fprintf("Mean power in [%i %i] Hz band : %2.4f (m/s^2)^2 \n", [band_power(1) band_power(2) mean(mean_power_band)]);
    % fprintf("Peak width (peak detection) : %2.2f Hz \n", tremor_width);
    if ~isempty(imposed_freq_tremor)
        fprintf("Mean tremor amplitude (imposed freq: %2.2f Hz): %2.2f mm \n", [imposed_freq_tremor, mean(tremor_amplitude_imposed_f_tremor)*10^(3)]);
    end
    fprintf("-----------------------------\n\n");
    
    %%%%%%%%%%%
    % OUTPUTS %
    %%%%%%%%%%%
    
    indicators.time_interval = time_interval;
    
    indicators.freq_tremor = freq_tremor;
    indicators.tremor_width = tremor_width;
    indicators.tremor_amplitude = tremor_amplitude;
    
    indicators.imposed_freq_tremor = imposed_freq_tremor;
    indicators.tremor_amplitude_imposed_freq = tremor_amplitude_imposed_f_tremor;
    
    indicators.band_power = band_power;
    indicators.mean_power_band = mean_power_band;
    indicators.power_in_band_over_time = power_in_band_over_time;
    
    indicators.mean_spectrum = mean_spectrum;
    
    %%%%%%%%%%%%%
    % Plotting: %
    %%%%%%%%%%%%%
    
    if plot_figure
%         f1=figure;
%         f1.Name = sprintf("%s - mean spectrum during 1st posture", figure_name);
%         plot(freq_spectrum, mean_spectrum{1}); hold on;
%         plot(freq_spectrum(locs{1}), peaks{1}, '*');
%         title(sprintf('Average spectrogram between %2.2f s and  %2.2f s', time_interval(1,:)));
%         xlabel('Frequency (Hz)');
%         ylabel('abs(mean(spectrum)) / Hz');
        
%         f2=figure;
%         f2.Name = sprintf("%s - Acceleration during recording", figure_name);
%         plot(time, acceleration);
%         vline(time_interval(1,:), 'r');
%         xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
%         title('Acceleration during recording')
        
        if size(time_interval,1) > 1
            f3=figure;
            f3.Name = sprintf("%s - Evolution of tremor amplitude during different postures", figure_name);
            
            h1 = subplot(311);
            plot(mean_power_band, '-o', 'MarkerFaceColor', 'r')
            title(sprintf('Power in [%i %i] Hz band during different posture', [band_power(1) band_power(2)]));
            ylabel('Power (m/s^2)^2');
            % xlabel('Posture number');
            set(gca,'XTick',(1:size(time_interval,1)))
            
            h2 = subplot(312);
            if ~isempty(imposed_freq_tremor)
                plot(tremor_amplitude_imposed_f_tremor*10^3, '-o', 'MarkerFaceColor', 'r');
                title(sprintf("Tremor amplitude (mm) at imposed freq: %2.2f Hz", imposed_freq_tremor));
            else
                plot(tremor_amplitude*10^3);
                title("Tremor amplitude (mm) at detected peak freq");
            end
            ylabel('Amplitude (mm)');
            % xlabel('Posture number'); 
            xlim([1 size(time_interval,1)])
            set(gca,'XTick',(1:size(time_interval,1)))
            
            h3 = subplot(313);
            plot(indicators.freq_tremor);
            title("Detected tremor peak (Hz)");
            if ~isempty(imposed_freq_tremor)
                hline(imposed_freq_tremor, 'r', 'Imposed freq')
            end
            ylabel('Frequency (Hz)');
            xlabel('Posture number'); 
            xlim([1 size(time_interval,1)])
            set(gca,'XTick',(1:size(time_interval,1)))
            ylim([0 max(indicators.freq_tremor)+5]);
            
            linkaxes([h1 h2 h3], 'x');
            
        end
    end
       
end