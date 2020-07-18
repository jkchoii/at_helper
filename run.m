clear; clc; close all; fclose('all');

addpath(genpath('./at'));


%% signal gen
broad_freq = 25:0.1:75;

dir_path = fullfile(pwd, 'env\ref');
run_dir_path = fullfile(pwd, 'env\run_env');

[ m_kraken ] = krakenProcessor(dir_path, run_dir_path, 'pekeris', 'broadband', broad_freq);

[ m_kraken ] = m_kraken.execute();
[ modes ] = m_kraken.getModeStruct('Propagation');

range_src2recv = 10000; % meter
recv_depth = 100; % meter
%
[ g ] = m_kraken.getGreenSpectrum(modes, range_src2recv, recv_depth);
time_delay = range_src2recv / 1500;
[ sig_t ] = m_kraken.getTimeSeries(g, time_delay);
[ mode_t ] = m_kraken.getModeDecomposition(modes, range_src2recv, recv_depth);


%% mode warping

stft_pack = m_kraken.fft_pack;
stft_pack.nfft = 2^(nextpow2(stft_pack.nfft) - 1);
stft_pack.fftfreq=(0:stft_pack.nfft-1)*stft_pack.fs/stft_pack.nfft;
stft_pack.strides = 1;
stft_pack.window_length = 31;
stft_pack.broad_freq = broad_freq;

m_helper = signalProcessingHelper;

%sig_t = awgn(sig_t, 0, 'measured');


[ STFT ] = m_helper.stft(...
    sig_t, ...
    hamming(stft_pack.window_length), ...
    stft_pack.strides, ...
    stft_pack.nfft, ...
    stft_pack.fs, ...
    'symmetric');

figure; 
subplot(8,1,1); plot(stft_pack.time, sig_t);
subplot(8,1,2:8); imagesc(stft_pack.time, stft_pack.fftfreq, abs(STFT)); 
axis xy; ylim([0 stft_pack.broad_freq(end)])
% figure; 
% subplot(8,1,1); plot(stft_pack.time, sig_t); xlim([6.5 7.5]);
% subplot(8,1,2:8); imagesc(stft_pack.time, stft_pack.fftfreq, abs(STFT)); 
% axis xy;  axis([6.5 7.5 0 stft_pack.broad_freq(end)]);

new_sig = range_src2recv / 1500;


[ m_warp ] = modeWarping(sig_t, stft_pack, 1500.0, range_src2recv);
[ sig_warped, warp_params ] = m_warp.process(sig_t, 'forward');

stft_pack = m_kraken.fft_pack;
stft_pack.nfft = 2^(nextpow2(stft_pack.nfft) - 1);
stft_pack.fftfreq=(0:stft_pack.nfft-1)*stft_pack.fs/stft_pack.nfft;
stft_pack.strides = 1;
stft_pack.window_length = 301;
stft_pack.broad_freq = broad_freq;

[ STFT_warped ] = m_helper.stft(...
    sig_warped, ...
    hamming(stft_pack.window_length), ...
    stft_pack.strides, ...
    stft_pack.nfft, ...
    stft_pack.fs, ...
    'symmetric');

figure; 
subplot(8,1,1); plot(warp_params.t, sig_warped);
subplot(8,1,2:8); imagesc(warp_params.t, stft_pack.fftfreq, abs(STFT_warped)); 
axis xy; ylim([0 stft_pack.broad_freq(end)])

%%
[ sig_unwarped, unwarp_params ] = m_warp.process(sig_warped, 'inverse', warp_params);

figure; plot(unwarp_params.t, sig_t);
hold on; plot(unwarp_params.t, sig_unwarped, 'or')