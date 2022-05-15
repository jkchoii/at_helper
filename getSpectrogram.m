function [spectro, mode_spectro] = getSpectrogram(m_kraken, m_save, sig_t, mode_t, modes, spectro_params)


%%
[ m_helper ] = signalProcessingHelper;
broad_freq = spectro_params.broad_freq;
range_src2recv = spectro_params.range_src2recv;
dir_name = spectro_params.dir_name;
cut_freq_ratio = spectro_params.cut_freq_ratio;

%%

[ groups, mode_t ] = m_kraken.getGroupTravelTimes(modes, mode_t, range_src2recv, cut_freq_ratio);
[ groups, mode_t ] = m_kraken.getRemovedEvanescentWaves(groups, mode_t);



%% spectrograms
stft_pack = m_kraken.fft_pack;
%stft_pack.nfft = 2^(nextpow2(stft_pack.fs*2) - 1);
stft_pack.nfft = 1024;
stft_pack.fftfreq=(0:stft_pack.nfft-1)*stft_pack.fs/stft_pack.nfft;
stft_pack.strides = 1;
%stft_pack.window_length = 31;
%window_ratio = stft_pack.fs/200;
%stft_pack.window_length = 31*window_ratio;
stft_pack.window_length = 31;
stft_pack.pad_length = floor(stft_pack.window_length / 2);
stft_pack.broad_freq = broad_freq;
stft_pack.roi_freq = [0 100];



%% cut samples

if spectro_params.cut_flag == true
    if iscolumn(range_src2recv) == true
        range_src2recv = range_src2recv.';
    end
    fastest_group_speed = max(groups{1}(1).vg);
    slowest_group_speed = min(groups{1}(end).vg);
    time_delay = [range_src2recv - spectro_params.cut_samples.start; range_src2recv + spectro_params.cut_samples.end] / fastest_group_speed;
    for idx = 1:length(range_src2recv)
        time_delay(1,idx) = min(groups{idx}(1).t_group) - 0.2;
    end
    sample_idx = floor(time_delay * m_kraken.fft_pack.fs);
    sample_diff = sample_idx(2,:) - sample_idx(1,:);
    num_sample_desired = max(sample_diff);
    sample_idx(1, :) = sample_idx(1, :) - (num_sample_desired - sample_diff);
    
    for idx = 1:length(sig_t)
        sig_t{idx} = sig_t{idx}(sample_idx(1,idx):sample_idx(2,idx));
        mode_t{idx} = mode_t{idx}(sample_idx(1,idx):sample_idx(2,idx),:);
    
    end

else
    sample_idx = zeros(2, length(range_src2recv));
    sample_idx(1,:) = 1;
    sample_idx(2,:) = length(sig_t{1});

end




[ spectro ] = m_helper.getSpectrogram(sig_t, stft_pack);
[ mode_spectro ] = m_helper.getSpectrogram(mode_t, stft_pack);

%[ groups ] = m_kraken.getGroupRefined(mode_spectro, groups);


stft_pack.sample_idx = sample_idx;
stft_pack.sample_min_time = stft_pack.time(sample_idx(1,:)) - stft_pack.pad_length / stft_pack.fs;

for idx = 1:length(mode_spectro)
    spectro(idx).t = spectro(idx).t + stft_pack.sample_min_time(idx);
    mode_spectro(idx).t = mode_spectro(idx).t + stft_pack.sample_min_time(idx);    
    %stft_pack.time = stft_pack.time(sample_idx(1,:):sample_idx(2,:));
    %stft_pack.time = stft_pack.time(sample_idx(1,:):end);
end


%[rows, cols] = size(spectro(idx).s);
%dd = zeros(rows, cols*2);

%figure('Renderer', 'painters', 'Position', [-900 474 512 512 ]);
%imagesc([]);
%h = gca;

% for idx = 1:length(mode_spectro)
%     t = linspace(spectro(idx).t(1), spectro(idx).t(end), 512);
%     f = spectro(idx).f(1:128);
%     f = linspace(f(1), f(end), 512);
% 
%     [T, F] = meshgrid(spectro(idx).t, spectro(idx).f(1:128));
%     [Tq, Fq] = meshgrid(t, f);
% 
%     s_cut = spectro(idx).s(1:128,:);
%     s_coh = abs(s_cut.^2);
%     mode_s_cut = mode_spectro(idx).s(1:128,:,:);
%     s_incoh = abs(sum(mode_s_cut.^2, 3));
%     s_coh = s_coh ./ max(s_coh(:));
%     s_incoh = s_incoh ./ max(s_incoh(:));
% 
%     s_coh = interp2(T,F, s_coh, Tq, Fq, 'spline');
%     s_incoh = interp2(T,F, s_incoh, Tq, Fq, 'spline');
% 
%     dd = [s_coh, s_incoh];
%     imagesc(dd(257:end,:));
%     immat = getframe(h);
% 
%     imwrite(immat.cdata, ['./incoh/' num2str(range_src2recv(idx)./1000) 'km.png']);
% end



% is_norm = true;
% [ spectro_handler ] = m_plt.plotSpectrogram(sig_t, spectro, groups, cell_idx, mode_cut, range_src2recv, stft_pack, true, is_norm);
% is_norm = true;
% [ mode_spectro_handler] = m_plt.plotSpectrogram(mode_t, mode_spectro, groups, cell_idx, mode_cut, range_src2recv, stft_pack, true, is_norm);


if spectro_params.is_save == true
    [m_save] = m_save.SetStartDataIndex(0);
    target_size = [256, 256];
    
    is_full_modes = false;
    [m_save] = m_save.saveSpectrogramFixedTime(dir_name, mode_t, mode_spectro, groups, stft_pack, [0 100], m_kraken.pos, is_full_modes, target_size );


    is_full_modes = true;
    [m_save] = m_save.saveSpectrogramFixedTime(dir_name, sig_t, spectro, groups, stft_pack, [0 100], m_kraken.pos, is_full_modes, target_size);
    
    
end



end

