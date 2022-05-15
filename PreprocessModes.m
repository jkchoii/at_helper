function [m_kraken, sig_t, mode_t] = PreprocessModes(m_kraken, modes, simulation_params)


resample_params = simulation_params.resample_params;
%cut_freq_ratio = simulation_params.cut_freq_ratio;
range_src2recv = simulation_params.range_src2recv;
recv_depth = simulation_params.recv_depth;

[ m_kraken, sig_t, mode_t ] = m_kraken.getSignalDecomposed(modes, range_src2recv, recv_depth, resample_params);

m_kraken.pos = [];
m_kraken.pos = repmat(m_kraken.pos, [length(range_src2recv), 1]);
for idx = 1:length(range_src2recv)
    m_kraken.pos(idx).sd = modes{1}.pos.s.z;
    m_kraken.pos(idx).rd = recv_depth;
    m_kraken.pos(idx).range_src2recv = range_src2recv(idx);
end

%% plot pressure in time series & power spectrum
%[ m_plt ] = plotPackage(m_kraken.fft_pack);
%cell_idx = 1;
%mode_cut = length(mode_t{1});
%[ pressure_handler ] = m_plt.plotPressure(sig_t, mode_t, cell_idx, mode_cut, range_src2recv);
%[ spectrum_handler ] = m_plt.plotSignalSpectrum(sig_t, cell_idx, range_src2recv);

%[ spectrum_handler ] = m_plt.plotSignalSpectrum(mode_t, cell_idx, range_src2recv);




end

