function visualize_tremor_assessment_recording(left_data, right_data, tremor_assessment, varargin)
    % DESCRIPTION:
    %
    % SYNTAX:
    %
    %
    % INPUT:

    filename = 'Last sonication';
    
    if ~isempty(varargin)
        for input_index = 1:2:length(varargin)
            switch varargin{input_index}
                case 'filename'
                    filename = varargin{input_index + 1};
            end
        end
    end
    
    n_figures=size(findobj('type','figure'), 1);
    f=figure(n_figures+1); 

    min_freq = 1;
    max_freq = 20;
    idx_min_freq = min(find(right_data.freq_spectrum > min_freq));
    idx_max_freq = max(find(right_data.freq_spectrum < max_freq));

    s_left = abs(left_data.spectrum(idx_min_freq:idx_max_freq, :));
    s_right = abs(right_data.spectrum(idx_min_freq:idx_max_freq, :));

    c_min_axis = min([s_left(:); s_right(:)]);
    c_max_axis = max([s_left(:); s_right(:)]);

    n_figures=size(findobj('type','figure'), 1);
    f=figure(n_figures+1);
    f.Name = filename;

    h1 = subplot(411) 
    plot(tremor_assessment.time_brain_vision, right_data.acceleration); hold on;
    plot(tremor_assessment.time_brain_vision, left_data.acceleration);
    xlabel('Time (s)')
    ylabel('Acceleration (m/s^2)');
    title('Acceleration')
    legend('Right hand', 'Left hand')
    if tremor_assessment.has_posture_marker
        vline(tremor_assessment.posture_interval, 'r', 'P')
    end
    if ~isempty(tremor_assessment.sonication_time_brain_vision)
        vline(tremor_assessment.sonication_time_brain_vision, 'b', 'Son')
    end

    h2 = subplot(412)
    imagesc(left_data.time_spectrum, left_data.freq_spectrum(idx_min_freq:idx_max_freq), ...
         10*log10(abs(left_data.spectrum(idx_min_freq:idx_max_freq, :))));
    caxis(10*log10([c_min_axis c_max_axis]));
    title('Left hand')
    ylabel('Freq [Hz]')
    xlabel('Time (s)')

    if tremor_assessment.has_posture_marker
        vline(tremor_assessment.posture_interval, 'r')
    end
    if ~isempty(tremor_assessment.sonication_time_brain_vision)
        vline(tremor_assessment.sonication_time_brain_vision, 'b')
    end

    h3 = subplot(413)
    imagesc(right_data.time_spectrum, right_data.freq_spectrum(idx_min_freq:idx_max_freq), ...
         10*log10(abs(right_data.spectrum(idx_min_freq:idx_max_freq, :))));
    caxis(10*log10([c_min_axis c_max_axis]));
    title('Right hand')
    ylabel('Freq [Hz]')
    xlabel('Time (s)')

    if tremor_assessment.has_posture_marker
        vline(tremor_assessment.posture_interval, 'r')
    end
    if ~isempty(tremor_assessment.sonication_time_brain_vision)
        vline(tremor_assessment.sonication_time_brain_vision, 'b')
    end

    h4 = subplot(414)
    plot(right_data.time_spectrum, right_data.indicators.power_in_band_over_time); hold on;
    plot(left_data.time_spectrum, left_data.indicators.power_in_band_over_time);
    title('Power in band');
    if tremor_assessment.has_posture_marker
        vline(tremor_assessment.posture_interval, 'r')
    end
    if ~isempty(tremor_assessment.sonication_time_brain_vision)
        vline(tremor_assessment.sonication_time_brain_vision, 'b')
    end
    legend('Right hand', 'Left hand')
    xlabel('Time (s)')
    ylabel(sprintf('Power in band [%i %i] Hz', left_data.indicators.band_power)); 
    
    linkaxes([h1 h2 h3 h4], 'x')
    xlim([min(left_data.time_spectrum) max(left_data.time_spectrum)]);
end