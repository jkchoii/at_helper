clear; clc; close all; fclose('all'); clear read_modes_bin;
cd('D:\Workspace\git\matlab\normal modes\at_helper')

broad_freq = 25:0.1:100;
[m_save] = savePackage('data');

%start_range = [3000; 5000; 7000; 9000; 11000];
%end_range = [5000; 7000; 9000; 11000; 13000];

start_range = [9000];
end_range = [9000];

num_samples = 1;


execute_params = [];
execute_params.bottom_speed_set = linspace(1600, 2000, 100);
execute_params.bottom_depth_set = linspace(100, 200, 100);
execute_params.source_depth_set = linspace(5, 95, 100);
execute_params.ssp_set = linspace(1450, 1550, 100); 
execute_params.rho_set = linspace(1.6, 2.2, 10);
% execute_params.bottom_speed_set = 1600;
% execute_params.bottom_depth_set = 100;
% execute_params.source_depth_set = 25;
% execute_params.ssp_set = 1500;

[range_src2recv, dir_name] = getSampledRange(start_range, end_range, num_samples);

recv_depth = 25; % meter
resample_params = struct('p', 1, 'q', 1);
cut_freq_ratio = 0.95;

simulation_params = [];    
simulation_params.resample_params = resample_params;
simulation_params.recv_depth = recv_depth;
simulation_params.cut_freq_ratio = cut_freq_ratio;
simulation_params.range_src2recv = range_src2recv;

spectro_params = [];
spectro_params.range_src2recv = range_src2recv;
spectro_params.cut_freq_ratio = cut_freq_ratio; 
spectro_params.cut_samples = [];
spectro_params.cut_flag = true;
spectro_params.cut_samples.start = 2000;
spectro_params.cut_samples.end = 2000;
spectro_params.is_save = true;
%spectro_params.is_save = false;
spectro_params.broad_freq = broad_freq;

addpath(genpath('./at'));




%% signal gen
%broad_freq = 50:0.1:100;


dir_path = '.\env\ref';
run_dir_path = '.\env\run_env';

% [ m_kraken ] = krakenProcessor(dir_path, run_dir_path, 'pekeris', 'broadband', broad_freq);
%[ m_kraken ] = krakenProcessor(dir_path, run_dir_path, 'pekeris_test', 'broadband', broad_freq);
[ m_kraken ] = krakenProcessor(dir_path, run_dir_path, 'pekeris_active_modes', 'broadband', broad_freq);
%[ m_kraken ] = m_kraken.execute(execute_params);

mode_info = [];
mode_info.raw_modes = 0;
mode_info.intense_modes = 0;
mode_info = repmat(mode_info, [length(m_kraken.file_list), 1]);

load("num_mode_info.mat");
train_file_list = [];
test_file_list = [];

for idx = 1:length(mode_info)
    if mode_info(idx).raw_modes == 5 && mode_info(idx).intense_modes == 5
        train_file_list = [train_file_list; mode_info(idx).filename];
    end

    if mode_info(idx).raw_modes ~= 5 && mode_info(idx).intense_modes == 5
        test_file_list = [test_file_list; mode_info(idx).filename];
    end
end
generate_file_list = [train_file_list; test_file_list];


tic
for dataset_idx = 1:length(m_kraken.file_list)

    close all;
    mode_info_save_path = split(m_kraken.file_list(dataset_idx), '\');
    mode_info_save_path(4) = [];

    fn = mode_info_save_path(4);
    mode_info_save_path = fullfile(mode_info_save_path(1), mode_info_save_path(2), mode_info_save_path(3), mode_info_save_path(4));
    
    gen = false;
    for file_valid = 1:length(generate_file_list)
        if strcmp(generate_file_list(file_valid), fn) == true
            gen = true;
            break;
        end
    end


    
    
%    fid = fopen([char(mode_info_save_path) '.txt'], 'w');
    if gen == true
    
        [ modes ] = m_kraken.getModeStruct(m_kraken.file_list(dataset_idx), 'PROPAGATION');
        
        [m_kraken, sig_t, mode_t] = PreprocessModes(m_kraken, modes, simulation_params);
        num_raw_mode = size(mode_t{1}, 2); 
    
        
        tokens = split(m_kraken.file_list(dataset_idx), '\');
        
        spectro_params.dir_name = [char(tokens(end)) '_' char(dir_name)];
    
        [sig_spectro, mode_spectro] = getSpectrogram(m_kraken, m_save, sig_t, mode_t, modes, spectro_params);
        
        num_intense_mode = size(mode_spectro(1).s, 3);
    
    %     msg = sprintf("%s %d %d", fn, num_raw_mode, num_intense_mode);
    %     fprintf(fid, "%s\n", msg);
    %     fclose(fid);
        mode_info(dataset_idx).raw_modes = num_raw_mode;
        mode_info(dataset_idx).intense_modes = num_intense_mode;
        mode_info(dataset_idx).filename = fn;

    end


    %figure; imagesc(sig_spectro(1).t, sig_spectro(1).f(1:512), abs(abs(sig_spectro(1).s(1:512,:)).^2)); set(gca, 'YDir', 'normal')
    %xlim([5000/1500-0.2, 5000/1500+0.2])
    %sm = sum(mode_spectro(1).s, 3);
    %figure; imagesc(sig_spectro(1).t, sig_spectro(1).f(1:512), abs(abs(sm(1:512,:)).^2)); set(gca, 'YDir', 'normal')
    %xlim([5000/1500-0.2, 5000/1500+0.2])
end 
toc
[dd, ss] = sort([mode_info.raw_modes]);
mode_info = mode_info(ss);
save("num_mode_info.mat", "mode_info");

warp_params = [];

warp_params.range_src2recv = range_src2recv;
warp_params.broad_freq = broad_freq; 

[ warped_spectro ] = getWarpedSpectrogram(m_kraken, sig_t, warp_params);
figure; imagesc(warped_spectro(1).t, warped_spectro(1).f(1:length(warped_spectro(1).f)/2), abs(warped_spectro(1).s(1:length(warped_spectro(1).f)/2,:)).^2); xlabel('time shifted and time warped'); ylabel('warped frequency'); title([num2str(range_src2recv./1000) 'km coherent']); axis xy;


[ sig_unwarped, unwarp_params ] = m_warp.process(sig_warped, 'inverse', warp_params);


