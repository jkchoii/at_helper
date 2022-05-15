function [ warp_spectro ] = getWarpedSpectrogram(m_kraken, sig_t, warp_params)

[ m_helper ] = signalProcessingHelper;

range_src2recv = warp_params.range_src2recv;
broad_freq = warp_params.broad_freq;

stft_pack = m_kraken.fft_pack;

sig_warped = cell(length(range_src2recv), 1);
warp_params = cell(length(range_src2recv), 1);
num_samples = [];
for idx = 1:length(range_src2recv) 
    %sig_t{idx} = sig_t{idx}(sample_idx(1,idx):end);
     [ m_warp ] = modeWarping(sig_t{idx}, stft_pack, 1500.0, range_src2recv(idx), num_samples);
    % [ m_warp ] = modeWarping(sig_t{idx}, stft_pack, 1500.0, (stft_pack.sample_min_time) * 1500);
    [ sig_warped{idx}, warp_params{idx} ] = m_warp.process('forward');
end

warp_stft_pack = m_kraken.fft_pack;
warp_stft_pack.fs = warp_params{idx}.warp_fs;
warp_stft_pack.nfft = 2^(nextpow2(warp_stft_pack.nfft) - 1);
warp_stft_pack.fftfreq=(0:warp_stft_pack.nfft-1)*warp_params{idx}.warp_fs/warp_stft_pack.nfft;
warp_stft_pack.strides = 1;
warp_stft_pack.window_length = 301;
warp_stft_pack.broad_freq = broad_freq;
warp_stft_pack.pad_length = floor(warp_stft_pack.window_length/2);

[ warp_spectro ] = m_helper.getSpectrogram(sig_warped, warp_stft_pack);


% [ warp_spectro_handler ] = m_plt.plotSpectrogram(sig_warped, warp_spectro, cell_idx, range_src2recv, warp_stft_pack, true);



end

