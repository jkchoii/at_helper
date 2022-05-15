classdef plotPackage

    properties(SetAccess = private)
        fft_pack
    end

    methods

        function obj = plotPackage(fft_pack)
            obj.fft_pack = fft_pack;

        end

        function [ handler ] = plotSignalSpectrum(obj, sig_cell, cell_idx, range_src2recv)
            handler = [];

            sig_t = sig_cell{cell_idx};

            handler.h = figure; 
            handler.description = sprintf("r = %.2f, power spectrum", range_src2recv(cell_idx));
            
            plot(obj.fft_pack.fftfreq, 20*log10(abs(fft(sig_t)))); 
            title(['Power spectrum, fs = ' num2str(obj.fft_pack.fs), 'Hz']); 
            ylabel('dB'), xlabel('Hz'); 
            xlim([0, obj.fft_pack.fs/2]); 
            grid on; grid minor;

        end
        
        function [ handler ] = plotPressure(obj, sig_cell, mode_cell, cell_idx, mode_cut, range_src2recv)
            sig_t = sig_cell{cell_idx};
            mode_t = mode_cell{cell_idx};

            handler = repmat(struct("h", 0, "description", ""), [2, 1]);
            
            %%
            handler(1).h = figure;
            handler(1).description = sprintf("r = %.2f, pressure in time series", range_src2recv(cell_idx));
            subplot(211); 
            plot(obj.fft_pack.time, sig_t); 
            title(['P(t,r,z), r = ' num2str(range_src2recv(cell_idx) / 1000), 'km']); 
            xlabel('time(s)'); ylabel('amp'); grid on; grid minor; ylim([-1 1]); %xlim(x_target_lim); 

            subplot(212); 
            plot(obj.fft_pack.time, mode_t); 
            title(['P_m(t,r,z), r = ' num2str(range_src2recv(cell_idx) / 1000), 'km']); 
            xlabel('time(s)'); ylabel('amp'); grid on; grid minor; ylim([-1 1]); %xlim(x_target_lim); 

            %%
            handler(2).h = figure;
            handler(2).description = sprintf("r = %.2f, decomposed pressure in time series", range_src2recv(cell_idx));

            num_modes = size(mode_t, 2);
            if num_modes > mode_cut
                num_modes = mode_cut;
            end
            for idx = 1:num_modes
                subplot(num_modes, 1, idx);
                plot(obj.fft_pack.time, mode_t(:,idx));
                title(sprintf("P_%d(t,r,z), r = %.2fkm", idx, range_src2recv(cell_idx) / 1000));
                xlabel('time(s)'); ylabel('amp'); grid on; grid minor; ylim([-1 1]); % xlim(x_target_lim);
            end
        end

        function [ handler ] = plotSpectrogram(obj, sig_cell, spectro, groups, cell_idx, mode_cut, range_src2recv, stft_pack, is_expand, is_norm)
            
            mode_len = size(spectro(cell_idx).s, 3);
            plot_num_mode = length(groups);
            if plot_num_mode > mode_cut
                plot_num_mode = mode_cut;
            end

            if mode_len > plot_num_mode
                mode_len = plot_num_mode;
            end
            
            interval = 16;
            num_plot = mode_len * interval;


            sig_t = sig_cell{cell_idx};
            s = abs(spectro(cell_idx).s);
            f = spectro(cell_idx).f;
            t = spectro(cell_idx).t;
            if is_norm
                sig_t = sig_t ./ max(sig_t);
                for idx = 1:size(s, 3)
                    s(:,:,idx) = s(:,:,idx) ./ max(max(s(:,:,idx)));
                end
            end



            handler(1).h = figure;
            handler(1).description = sprintf("r = %.2fkm spectrogram", range_src2recv(cell_idx) / 1000);
            c = colororder;
            clen = length(c);
            for idx = 1:mode_len
                c_idx = 1+mod(idx-1, clen);
                subplot(num_plot, 1, 1+(idx-1)*interval:1+(idx-1)*interval+3)
                plot(stft_pack.time(stft_pack.sample_idx(1,cell_idx):stft_pack.sample_idx(2,cell_idx)), sig_t(:,idx), "Color", c(c_idx,:)); 
                %plot(stft_pack.time, sig_t(:,idx), "Color", c(c_idx,:)); 
                ylabel('amp'); ylim([-1 1]); 
                xticks([]);
                grid on; grid minor; 
            
                subplot(num_plot, 1, (1+(idx-1)*interval+4):(1+(idx)*interval-1)); 
                imagesc(t, f, abs(s(:,:,idx))); 
                ylim([min(stft_pack.roi_freq) max(stft_pack.roi_freq)]);  ylabel('frequency(Hz)');
                if is_norm
                    caxis([0 1]);
                else

                    caxis([0 3.5]);
                end
                xt = xticks;
                xticks([]);
                axis xy;
                

                hold on; 
                if mode_len == 1
                    for loop = 1:plot_num_mode
                        %plot(groups(loop).t_group(:,:,cell_idx), groups(loop).f_group, 'LineWidth', 1, "Color", c(loop,:));
                        plot(groups(loop).t_group(:,:,cell_idx), groups(loop).f_group, 'LineWidth', 1.5, "Color", 'red');
                    end
                else
                    plot(groups(idx).t_group(:,:,cell_idx), groups(idx).f_group, 'LineWidth', 1.5, "Color", 'red');
                end
                hold off;
                

            end
            xlabel('time t(s)');
            xticks(xt);
            sgtitle(sprintf("r = %.2fkm, spectrogram", range_src2recv(cell_idx) / 1000));

            if is_expand
                time_expand = 0.5;
                x_arrival = range_src2recv(cell_idx) / 1500;
                x_min = x_arrival - time_expand;
                x_max = x_arrival + time_expand*3;
                if x_min < 0
                    x_min = 0;
                end

                %if x_max > max(stft_pack.time)
                if x_max > max(t)
                    %x_max = max(stft_pack.time)-.1;
                    x_max = max(t)-.1;
                end
                
                for idx = 2:length(handler(1).h.Children)
                    handler(1).h.Children(idx).XLim = [x_min, x_max];
                %handler(1).h.Children(3).XLim = [x_min, x_max];
                end
                xticks('auto');
            end

        end



    end
end

