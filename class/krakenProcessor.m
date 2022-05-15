classdef krakenProcessor
    % KRAKENC.exe/field.exe 실행하고 데이터 뽑아냄
    % 
	%
    % methods
    %   public:
    %       krakenProcessor(*) - 원본 .env, .flp 파일이 위치한 폴더 지정
    %       execute(*) - krakenc.exe/field.exe 실행
    %       getModeStruct(*) - krakenc.exe에서 계산된 모드 정보 및 파라미터
    %       getModeDecomposition(*) - 각각의 Green에서 각각의 mode pressure 뽑아냄
    %       getTimeSeries(*) - Green 함수에서 전체 time series 데이터 뽑아냄
    %       getGreenSpectrum(*) - Green spectrum 구함
    %   private:
    %       modeStructEdit(*) - 허수부 conjugate 및 phi_zs 처리
    %       envParamsReader(*) - .env, .flp 사본 생성
    %       krakenBroadBandProcessor(*) - 광대역 신호 생성을 위한 .env 파일 생성
    %       preprocessModeStructure(*) - Kraken에 의해 계산된 결과 읽고 다루기 쉽게 가공
    %       getFFTPack(*) - sampling theorem에 일치하는 time, freq vector 생성
    %
    
    properties(SetAccess = private)
        file_root
        file_name
        run_file_root
        
        is_broadband
        freq_band
        
        
        is_precision_warning
        file_list
    end

    properties(SetAccess = public)
        pos

    end



    properties(SetAccess = public)
        fft_pack
    end
    
    methods
        function obj = krakenProcessor(file_path, run_path, file_name, varargin)
            %KRAKENPROCESSOR 이 클래스의 인스턴스 생성
            % file_path: 원본 env, flp가 위치한 폴더 경로
            % run_path: 사본 env, flp가 위치할 폴더 경로
            % file_name: env, flp 파일 이름(.env, .flp는 같은 이름이어야 함)
            % varargin: ex)  krakenProcessor(*,*,*,'broadband', 50:100)
            
            env_files = dir(fullfile(file_path, [file_name, '.env']));
            assert(~isempty(env_files), '.env 파일 없음')
            
            flp_files = dir(fullfile(file_path, [file_name, '.flp']));
            assert(~isempty(flp_files), '.flp 파일 없음')
            
            obj.file_root = file_path;
            obj.file_name = file_name;
            
            obj.run_file_root = run_path;
            
            obj.is_broadband = false;
            obj.freq_band = [];
            obj.is_precision_warning = false;
            obj.pos;

            if ~isempty(varargin)
                narginchk(4,5);
                if strcmpi(varargin{1}, 'BROADBAND')
                    obj.is_broadband = true;
                    
                    assert(~isempty(varargin{2}), 'freq band 비었음');
                    assert(isvector(varargin{2}), 'freq band 벡터만 가능');
                    obj.freq_band = varargin{2};
                    [ obj ] = getFFTPack(obj);
                    
                else
                    error('broadband 매개변수 확인');
                end
            end
            obj.file_list = [];
            ext_lists = ["bot_speed"; "depth"; "sd"; "ssp";"rho"];
            for idx = 1:length(ext_lists)
                search_dir = fullfile(run_path, char(ext_lists(idx)));
                lists = dir(fullfile(search_dir, '*.env'));
                for jdx = 1:length(lists)
                    obj.file_list = [obj.file_list; string(fullfile(run_path, char(ext_lists(idx)), lists(jdx).name(1:end-4)))];

                end
            end

        end
        
        
        function [ obj ] = execute(obj, execute_params)
            % 사본 파일이 위치한 곳에서 Kraken 실행
            % env 파일로 .mod 파일 생성
            % flp 파일로 .shd 파일 생성
            
            dir_info = dir([obj.file_root '\*.env']);
            assert(~isempty(dir_info), '잘못된 경로');
            
            
            run_full_path = fullfile(obj.run_file_root, obj.file_name);
            
            if obj.is_broadband
                [obj, exe_file_lists] = krakenBroadBandProcessor(obj, run_full_path, execute_params);
            else % single freq
                [obj, Opt, Pos_flp] = envParamsReader(obj, run_full_path);
            end
 
            tic
            parfor idx = 1:length(exe_file_lists)

                system(['krakenc_jk ' char(exe_file_lists(idx))]);
                system(['field_jk ' char(exe_file_lists(idx))]);

            end
            toc

            obj.file_list = exe_file_lists;
            
            writematrix(exe_file_lists, fullfile(obj.run_file_root, [datestr(now, 'yymmdd-HHMM') '_info.txt']));
%             fid = fopen(fullfile(obj.run_file_root, [datestr(now, 'yymmdd-HHMM') '_info.txt']), 'rt');
%             file_lists = [];
%             while ~feof(fid)
%                 file_lists = [ file_lists; string(fgetl(fid))];
%             end
%             fclose(fid);

            fclose('all'); % adhoc
        end

        function [ groups, mode_t ] = getGroupTravelTimes(obj, modes, mode_t, range_src2recv, cut_ratio)
            
            num_range = length(range_src2recv);
            max_num_modes = length(modes{end}.group_velocity);
            vg = zeros(length(modes), max_num_modes) - 1;
            vp = zeros(length(modes), max_num_modes) - 1;
            f_vg = zeros(length(modes), 1);

            %% merge
            for idx = 1:length(modes)
                
                f_vg(idx) = modes{idx}.freqVec;
                num_modes = length(modes{idx}.group_velocity);
                vg(idx, 1:num_modes) = modes{idx}.group_velocity;
                vp(idx, 1:num_modes) = modes{idx}.phase_velocity;

            end

            t_cell = cell(max_num_modes, num_range);
            pg_ratio_cell = cell(max_num_modes, 1);
            f_cell = cell(max_num_modes, 1);
            
            %% group travel time
            for idx = 1:max_num_modes
                
                valid_modes = vg(:,idx) > 0;
                f_cell{idx} = f_vg(valid_modes);
                pg_ratio_cell{idx} = vg(valid_modes, idx) ./ vp(valid_modes, idx);
                for jdx = 1:num_range
                    t_cell{idx, jdx} = range_src2recv(jdx) ./ vg(valid_modes, idx);
                    %pg_ratio_cell{idx,jdx} = vp(valid_modes, idx) ./ vg(valid_modes, idx);
                end
                
            end


            %% align
            groups = struct('num_mode', 0, 'vg', zeros(size(vg)), 'vp', zeros(size(vg)), ...
                't_group', zeros(length(t_cell{1,1}), 1, num_range),...
                'f_group', zeros(length(f_cell{1}), 1), ...
                'vp_ratio',zeros(length(t_cell{1,1}), 1));
            groups = repmat(groups, [max_num_modes, 1]);

            for idx = 1:max_num_modes
                groups(idx).num_mode = idx;
                groups(idx).vg = vg(vg(:,idx) > 0);
                groups(idx).vp = vp(vp(:,idx) > 0);
                groups(idx).t_group = zeros(length(t_cell{idx,1}), 1, num_range);
                groups(idx).f_group = zeros(length(f_cell{idx}), 1);
                groups(idx).vp_ratio = zeros(length(t_cell{idx,1}), 1);

                groups(idx).vp_ratio = pg_ratio_cell{idx};

                for jdx = 1:num_range
                    groups(idx).t_group(:,:,jdx) = t_cell{idx,jdx};
                    groups(idx).f_group = f_cell{idx};
                end
            end


            max_iter = length(groups);
            
            
            groups = getSurvivedGroups(obj, groups, cut_ratio, max_iter);
%             groups = groups(survived);
% 
%             for idx = 1:length(groups)
%                 [cut_idx] = getCutIndex(obj, groups(idx).t_group(:,:,1));
%                 slice = cut_idx:size(groups(idx).t_group, 1);
%                 groups(idx).t_group = groups(idx).t_group(slice,1,:);
%                 groups(idx).f_group = groups(idx).f_group(slice);
%             end
% 
% 
%             
%             for idx = 1:length(mode_t)
%                 mode_t{idx} = mode_t{idx}(:,survived);
% 
%             end


        end

        function [ groups ] = getGroupRefined(obj, mode_spectro, groups)

            max_val = zeros(1,1,size(mode_spectro(1).s, 3)); 

            max_mode = 0;
            for idx = 1:length(mode_spectro)
                if size(mode_spectro(idx).s, 3) > max_mode
                    max_mode = size(mode_spectro(idx).s, 3);
                end

            end
            s = zeros(size(mode_spectro(1).s, 1), size(mode_spectro(1).s, 2), max_mode);
            for idx = 1:length(mode_spectro)
                current_s = abs((mode_spectro(idx).s)).^2;
                for jdx = 1:size(mode_spectro(idx).s, 3)
                    max_val(1,1,jdx) = max(max(current_s(:,:,jdx)));
                end
                current_s = current_s ./ max_val(1,1, 1:size(current_s, 3));
                s(:,:, 1:size(current_s, 3)) = s(:,:, 1:size(current_s, 3)) + current_s;
                
            end

            s_power = squeeze(sum(s,2));
            s_power = s_power ./ max(s_power);
            nfft = length(mode_spectro(1).f);
            f = mode_spectro(1).f(1:floor(nfft/2));
            s_power = s_power(1:floor(nfft/2), :);
            

            max_pos_set = zeros(size(s_power, 2),1);
            for idx = 1:size(s_power, 2)
                [~, max_pos_set(idx)] = max(s_power(:,idx));
            end


            df = obj.freq_band(2) - obj.freq_band(1);
            for jdx = 1:length(groups)
                for idx = length(groups{jdx}):-1:2
                    if max_pos_set(idx) >= max_pos_set(idx-1)
                        s_check = s_power(:,idx);
                        for kdx = 1:length(s_check)
                            if s_check(kdx) > 0.2
                                cut_pts = kdx;
                                break;
                            end
                        end
                        %[~, cut_pts] = min(abs(s_power(:,idx) - 0.2));
                        [diff_val, cut_freq_pts] = min(abs(f(cut_pts) - groups{jdx}(idx).f_group));
                        if diff_val <= 2*df
                            slice = cut_freq_pts:length(groups{jdx}(idx).f_group);
                                             
                            if length(groups{jdx}(idx).f_group) > length(slice)
                                groups{jdx}(idx).f_group = groups{jdx}(idx).f_group(slice);
                                groups{jdx}(idx).t_group = groups{jdx}(idx).t_group(slice);
                            end
                            
                        end
                    end
                end
            end
            
            for idx = 1:length(groups)
                groups{idx} = rmfield(groups{idx}, ["vp", "vp_ratio"]);

            end
            

        end


        function [cand_idx] = getCutIndex(obj, t_group)
            %% get slope
            check = [];
            cand_idx = 1;
            for idx = (length(t_group)):-1:3
                jdx = idx - 1;
                %slope_sign = (t_group(end) - t_group(idx));
                slope_sign = (t_group(idx) - t_group(idx-1)) * (t_group(jdx) - t_group(jdx-1));
                check = [check; slope_sign];
                
                if slope_sign > 0
                    continue;
                end
                cand_idx = idx;
                break;
            end



        end

        function [obj, sig_cell, mode_cell] = getSignalDecomposed(obj, modes, range_src2recv, recv_depth, resample_params)

            sig_cell = cell(length(range_src2recv), 1);
            mode_cell = cell(length(range_src2recv), 1);
            for idx = 1:length(range_src2recv)
                %
                [ g ] = getGreenSpectrum(obj, modes, range_src2recv(idx), recv_depth);
                time_delay = range_src2recv(idx) / 1500;
                sample_idx = floor(time_delay * obj.fft_pack.fs);
                %[ sig_t ] = m_kraken.getTimeSeries(g, time_delay);
                [ sig_t ] = getTimeSeries(obj, g, 0);
                [ mode_t ] = getModeDecomposition(obj, modes, range_src2recv(idx), recv_depth);
                %sig_t = flipud(sig_t);
                
                sig_cell{idx} = sig_t;
                mode_cell{idx} = mode_t;
             
            end

%             for idx = 1:length(sig_cell)
%                 sig_cell{idx} = filtfilt(resample_params.sos, sig_cell{idx});
%             end
% 
%             for idx = 1:length(mode_cell)
%                 for jdx = 1:size(mode_cell{idx}, 2)
%                     mode_cell{idx}(:,jdx) = filtfilt(resample_params.sos, mode_cell{idx}(:,jdx));
%                 end
%             end


            sig_cell = getResampledSignal(obj, sig_cell, resample_params);
            mode_cell = getResampledSignal(obj, mode_cell, resample_params);
            L = length(sig_cell{1});
            obj.fft_pack.time = (0:(L-1))*resample_params.q/(resample_params.p*obj.fft_pack.fs);
            obj.fft_pack.fs = obj.fft_pack.fs * resample_params.p/resample_params.q;
            df = obj.fft_pack.fs/L;
            obj.fft_pack.fftfreq = 0:df:(df*(L-1));
            

        end 

        function [groups] = getSurvivedGroups(obj, groups, cut_ratio, max_iter, max_nfreq)
%            survived = true(ones(length(groups), 1));
%             sample_cut_threshold = floor((1 - cut_ratio) * max_nfreq);
%             for idx = 1:max_iter
%                 valid_idx = cut_ratio < groups(idx).vp_ratio;
%                 %groups(idx).vg = groups(idx).vg(valid_idx);
%                 groups(idx).t_group = groups(idx).t_group(valid_idx,1,:);
%                 groups(idx).f_group = groups(idx).f_group(valid_idx);
%                 if isempty(groups(idx).t_group) == true
%                     survived(idx) = false;
%                 end
%                 cut_ratio = cut_ratio - 0.05;
%             end
            
            
            for idx = 1:max_iter
                sample_cut_threshold = floor((1.0-cut_ratio) * length(groups(idx).t_group));
%                 num_samples = length(groups(idx).t_group(:,:,1));
%                 if num_samples > sample_cut_threshold+1
                    groups(idx).t_group = groups(idx).t_group(sample_cut_threshold+1:end,:,:);
                    groups(idx).f_group = groups(idx).f_group(sample_cut_threshold+1:end);
%                 else
%                     survived(idx) = false;
%                 end
                
            end

%             for idx = max_iter:-1:1
%                 %if length(groups(idx).vg) > 2.5 * cut_threshold
%                 if survived(idx) == true
%                     if length(groups(idx).vg) < 2 * sample_cut_threshold % 너무 적으면 자름
%                         survived(idx) = false;
%                     %else
%                         %start_pts = floor(sample_cut_threshold/10) + 1;
%                         
%                         %groups(idx).t_group = groups(idx).t_group(start_pts:end,:,:);
%                         %groups(idx).f_group = groups(idx).f_group(start_pts:end,:,:);
%                     end
%                 end
%                 
%             end
        end

        
        function [ modes ] = getModeStruct(obj, filename, varargin)
            % 계산된 .mod 파일에서 모드 반환
            % varargin: ex) getModeStruct('propagation')
            %           - propagation modes vs leaky modes
            %           default: leaky modes
            
            
            narginchk(2,3);
            prop_mode_flag = false;
            if ~isempty(varargin) % exclude leaky mode
                narginchk(3,3);    

                if strcmpi(varargin{1}, 'PROPAGATION') 
                    prop_mode_flag = true;
                else
                    prop_mode_flag = false;
                    %error('Propagation 매개변수 확인');
                end
            end
            
            %filename = fullfile(obj.run_file_root, obj.file_name);
            filename = char(filename);
            modes = cell(length(obj.freq_band),1);
            
            parfor idx = 1:length(obj.freq_band)
            %for idx = 1:length(obj.freq_band)
                
                [modes{idx}] = preprocessModeStructure(obj, filename, obj.freq_band(idx), prop_mode_flag);
                
                modes{idx}.freqVec = obj.freq_band(idx);
                modes{idx}.Nfreq = idx;
            end

            [modes] = getModeGroupVelocities(obj, [filename '.grv'], modes);

        end
        
        function [ sig_t, norm_val ] = getTimeSeries(obj, g, time_delay)
            % Green 함수에서 time series 생성
            % g: Green spectrum
            % time_delay: 시간 지연 ex) range(m)/1500(m/s)
            % sig_t: time series
            
            narginchk(3,3);

            fftpack = obj.fft_pack;
            min_freq = min(obj.freq_band);
            [~, insert_idx] = min(abs(min_freq - fftpack.fftfreq));

            %g_pad = zeros(fftpack.nfft, size(g, 2));
            g_pad = zeros(length(fftpack.fftfreq), size(g, 2));

            %gg = fft(hanning(31), length(g));
            %g = conv(gg,g);

            g_pad(insert_idx:insert_idx+size(g,1)-1,:) = g;


            %g_pad(1:size(g,1),:) = g;
            %g_pad(1:size(g,1),:) = g;
            %sig_f = g_pad.*exp(1j.*2*pi*obj.fft_pack.fftfreq'.*time_delay);
            sig_f = g_pad;
            %sig_t = ifft(sig_f*exp(-1j.*2*pi*[50:0.1:100]), 'symmetric');
            sig_t = ifft(sig_f, 'symmetric');
            
            %sig_t = flipud(sig_t);
            norm_val = max(abs(sig_t));
            sig_t = sig_t / norm_val;
            
        end
        
        function [ mode_t ] = getModeDecomposition(obj, modes, r, varargin)
            % 각각의 Green에서 각각의 mode pressure 뽑아냄
            % modes: getModeStruct(*)에 의해 얻어진 모드 구조체
            % r: src-recv간 range(m)
            % varargin: recv 수심 <- 해당 수심에서의 mode-pressure 얻어냄
            %           default: 모든 수심(구현 안한듯. 기억이 안남)
            % mode_t: mode1 ~ modeN 각각의 pressure
            
            
            ref_mode = modes{end};
            if ~isempty(varargin)
                recv_depth = varargin{1};
                [~, recv_idx] = min(abs(recv_depth - ref_mode.z));
            else
                recv_idx = 1:length(ref_mode.z); % 미구현됨
            end
            
            %% processing
            freq_phi_zs = zeros(length(modes), length(modes{end}.k));
            freq_phi_zr = zeros(length(modes), length(modes{end}.k));
            freq_eigenvals = zeros(length(modes), length(modes{end}.k));
            for idx = 1:length(modes)
                for jdx = 1:length(modes{idx}.k)
                    freq_phi_zs(idx, jdx) = modes{idx}.phi_zs(jdx);
                    freq_phi_zr(idx, jdx) = modes{idx}.phi(recv_idx, jdx);
                    freq_eigenvals(idx, jdx) = modes{idx}.k(jdx);
                end
            end
            
            mode_g = zeros(length(modes), length(modes{end}.k));
            for idx = 1:length(modes{end}.k)
                valid_idx = find(freq_eigenvals(:,idx));
%                 mode_g(valid_idx,idx) = ...
%                     freq_phi_zr(valid_idx,idx) .* freq_phi_zs(valid_idx,idx) ...
%                     .*besselh(0, freq_eigenvals(valid_idx,idx)*r);

                mode_g(valid_idx, idx) = freq_phi_zr(valid_idx,idx) .* freq_phi_zs(valid_idx,idx) ...
                                .* exp(-1j*freq_eigenvals(valid_idx,idx)*r) .* sqrt(2./(pi*freq_eigenvals(valid_idx,idx)*r));

            end
            rho_w = ref_mode.Top.rho;
            mode_g = 1j ./ rho_w*(sqrt(8*pi*r)) * exp(1j*pi/4) * mode_g;
           


            fftpack = obj.fft_pack;
            min_freq = min(obj.freq_band);
            [~, insert_idx] = min(abs(min_freq - fftpack.fftfreq));

            mode_f = zeros(fftpack.nfft,length(modes{end}.k));
            mode_f(insert_idx:insert_idx+length(obj.freq_band)-1,:) = mode_g;
            
            mode_t = ifft(mode_f, 'symmetric');
            %norm_val = norm(sum(mode_t, 2));
            norm_val = max(abs(sum(mode_t, 2)));
            mode_t = mode_t ./ norm_val;
            %mode_t = mode_t ./ (norm_val*ones(size(mode_t)));
            %mode_t = flipud(mode_t);
        end

        function [groups, mode_t] = getRemovedEvanescentWaves(obj, groups, mode_t)

            rem_threshold = 0.01; % 1% 미만 제거
            
            num_range = length(mode_t);
            new_groups = cell(1, 1);
            new_groups{1} = groups;
            for idx = 1:length(new_groups{1})
                new_groups{1}(idx).valid_flag = true;
            end
            new_groups = repmat(new_groups, [num_range, 1]);

            check = [];
            for idx = 1:num_range
                
                sig_m_norm = mode_t{idx} ./ max(mode_t{idx}(:));
                energy_ratio = max(sig_m_norm);
                %check = [check; energy_ratio];
                survived_idx = energy_ratio > rem_threshold;
                

                for jdx = 1:length(survived_idx)
                    new_groups{idx}(jdx).t_group = groups(jdx).t_group(:,:,idx);
                    
                    if survived_idx(jdx) == false
                        new_groups{idx}(jdx).valid_flag = false;
                        %disp(num2str(idx));
                    end
                    
                end
                mode_t{idx} = mode_t{idx}(:,survived_idx);

                %survived_idx = true(length(energy_ratio), 1);

                %mode_t{idx} = mode_t{idx}(:,survived_idx);
                
                %new_groups{idx} = groups(survived_idx);
                %for jdx = 1:length(new_groups{idx})
                %    new_groups{idx}(jdx).t_group = groups(jdx).t_group(:,:,idx);
                %end

            end

            groups = new_groups;


        end
        
        
        function [g] = getGreenSpectrum(obj, modes, r, varargin)
            % Green spectrum 구함
            % modes: getModeStruct(*)에서 얻어진 모드 구조체
            % r: src-recv간 거리(m)
            % varargin: recv 수심 <- 해당 수심에서의 mode-pressure 얻어냄
            %           default: 모든 수심(구현 안한듯. 기억이 안남)
            % g: Green spectrum
            % ref: Computational Ocean Acoustics, p.351, Eq(5.63)
            
            ref_mode = modes{end};
            
            if ~isempty(varargin)
               recv_depth = varargin{1};
               [~, recv_idx] = min(abs(recv_depth - ref_mode.z));
            else
                recv_idx = 1:length(ref_mode.z);
            end
            
            g = zeros(length(obj.freq_band), length(recv_idx));
            rho_w = modes{1}.Top.rho;
            
            for idx = 1:length(obj.freq_band) %length(modes)
                
                phi_z = modes{idx}.phi(recv_idx,:);
                phi_zs = modes{idx}.phi_zs;
                k = modes{idx}.k;
                %g(idx,:) = sum(phi_z.*phi_zs .* besselh(0, k*r).', 2);
                

                %asymptotic form
                %asympt_form = phi_z .* phi_zs ./ sqrt(k*r).' .* (sqrt(2./(pi*k*r)) .* exp(1j.*(k*r - pi/4))).';
                %hankel = (sqrt(2./(pi*k*r)) .* exp(1j.*(k*r - pi/4))).';
                asympt_form = phi_z .* phi_zs .* exp(-1j*k*r).' .* sqrt(2./(pi*k*r)).';
                %asympt_form = phi_z .* phi_zs .* exp(-1j*k*r).' ./ sqrt(2./(pi*k*r)).' .* exp(1j*pi/4);
                %grn = asympt_form;
                %grn = asympt_form .* hankel;
                %grn = grn ./ vecnorm(grn);
                g(idx,:) = sum(asympt_form);
              
            end
            
            %g = 1j ./ (rho_w * sqrt(8*pi*r)) * exp(-1j*pi/4) * g;
            %g = 1j ./ (4*rho_w) * exp(1j*pi/4) * g;
            %g = 1j ./ (4*rho_w) * exp(-1j*pi/4) * g;
            %g = 1j ./ rho_w*(sqrt(8*pi*r)) * exp(-1j*pi/4) * g;
            g = 1j ./ (4*rho_w) * exp(1j*pi/4) * g;
            
            g_pack = [];
            %g_pack.grn = g ./ vecnorm(g);
            g_pack.freqs = obj.freq_band;
            %g_pack.
        end
        
    end
    
    %% private region
    methods(Access = private)
        
        function [ obj, Opt, Pos_flp, TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax ] ...
                = envParamsReader(obj, run_full_path)
            
            root_path = fullfile(obj.file_root, obj.file_name);
            [ TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, ~] = read_env( root_path, 'kraken' );
            [ ~, Opt, ~, ~, ~, ~, R, Pos_flp ] = read_flp( root_path ) ; 
            Pos_flp.r.range = R/1000.0;
            %obj.pos = Pos_flp;
            obj.freq_band = freq;
            write_env( run_full_path, 'kraken', TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax);
            
        end
        
        function [new_phi, phi_zs, z, zs] = modeStructEdit(obj, phi, z, pos)
        % 소스 위치와 리시버 위치가 불일치할 경우 고유 함수가 src 위치에서도 나옴
        % zs, z의 고유 함수를 따로 관리
        % phi: [z, zs]에서의 고유 함수 (num recv+불일치 소스 수 x 모드 수)
        % new_phi: z에서의 고유 함수 (num recv x 모드 수)
        % phi_zs: zs에서의 고유함수

            recv_depth = pos.r.z;
            src_depth = pos.s.z;
            assert((length(src_depth) == 1), 'source가 여러개면 phi 구성 다시하는거 새로 짜야함')

            for idx = 1:length(src_depth)
                [min_val, min_idx] = min(abs(src_depth(idx) - recv_depth));
                if min_val < 1e-3 % recv depth와 src_depth가 일치
                    phi_zs = phi(min_idx,:);
                    new_phi = phi;
                else
                    % recv depth와 src_depth가 불일치
                    [~, src_idx] = min(abs(z - src_depth(idx)));
                    logic_idx = 1:length(z);
                    phi_zs = phi(src_idx, :);
                    new_phi = phi(setdiff(logic_idx, src_idx), :);
                end

            end
            z = recv_depth;
            zs = src_depth;

            if size(phi_zs, 2) == 1
                phi_zs = phi_zs.';
            end
        end

        function [file_list] = WriteVariousEnvironmentFiles(obj, run_full_path, BB_freq, execute_params)
            total_set_count = length(execute_params.bottom_depth_set) ...
                + length(execute_params.bottom_speed_set) ...
                + length(execute_params.source_depth_set) ...
                + length(execute_params.ssp_set) ...
                + length(execute_params.rho_set);

            list_count = 1;
            file_list = string(zeros(total_set_count, 1));

            path_list = WriteEnvironmentFile(obj, 'rho', BB_freq, execute_params);
            file_list(list_count:(list_count + path_list.count - 1)) = path_list.file_name;
            list_count = list_count + path_list.count;


            path_list = WriteEnvironmentFile(obj, 'sd', BB_freq, execute_params);
            file_list(list_count:(list_count + path_list.count - 1)) = path_list.file_name;
            list_count = list_count + path_list.count;

            path_list = WriteEnvironmentFile(obj, 'ssp', BB_freq, execute_params);
            file_list(list_count:(list_count + path_list.count - 1)) = path_list.file_name;
            list_count = list_count + path_list.count;

            path_list = WriteEnvironmentFile(obj, 'depth', BB_freq, execute_params);
            file_list(list_count:(list_count + path_list.count - 1)) = path_list.file_name;
            list_count = list_count + path_list.count;

            path_list = WriteEnvironmentFile(obj, 'bot_speed', BB_freq, execute_params);
            file_list(list_count:(list_count + path_list.count - 1)) = path_list.file_name;

            


        end

        function [path_list] = WriteEnvironmentFile(obj, tag, BB_freq, var_env)
            dir_path = fullfile(obj.run_file_root, tag);
            is_exist = ~isempty(dir(dir_path));
            if is_exist == false
                mkdir(dir_path);
            end
            
            [ path_list ] = WriteCurrentEnvironment(obj, tag, BB_freq, var_env);

        end

        function [path_list] = WriteCurrentEnvironment(obj, tag, BB_freq, var_env)
            path_list = [];
            path_list.count = 0;
            path_list.file_name = [];


            [ obj, Opt, Pos_flp, TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax ] ...
                = envParamsReader(obj, fullfile(obj.run_file_root, obj.file_name ));
            path = '';

            if numel(Opt) <= 4
                opt_tmp = char(zeros(1,4));
                opt_tmp(1:numel(Opt)) = Opt;
                opt_tmp(numel(Opt)+1:end) = ' ';
            end
            Opt = opt_tmp;

            if strcmp(tag, 'sd') == true
                source_depth_set = var_env.source_depth_set;
                path = 'sd';
                path_list.file_name = string(zeros(length(source_depth_set), 1));
                
                for idx = 1:length(source_depth_set)

                    Pos.s.z = source_depth_set(idx);
                    Pos_flp.s.z = source_depth_set(idx);


                    this_file_name = [obj.file_name '_' path '_' num2str(source_depth_set(idx))];
                    write_path = fullfile(obj.run_file_root, path, this_file_name);

                    write_env_broadband( write_path, 'kraken', TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, BB_freq);
                    write_fieldflp(write_path, Opt, Pos_flp);

                    path_list.count = path_list.count + 1;
                    path_list.file_name(idx) = string(write_path);
                end
                

            elseif strcmp(tag, 'ssp') == true
                ssp_set = var_env.ssp_set;
                path = 'ssp';
                path_list.file_name = string(zeros(length(ssp_set), 1));

                for idx = 1:length(ssp_set)

                    SSP.c(:) = ssp_set(idx);
                    SSP.raw.alphaR(:) = ssp_set(idx);

                    this_file_name = [obj.file_name '_' path '_' num2str(ssp_set(idx))];
                    write_path = fullfile(obj.run_file_root, path, this_file_name);

                    write_env_broadband( write_path, 'kraken', TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, BB_freq);
                    write_fieldflp(write_path, Opt, Pos_flp);
                    
                    path_list.count = path_list.count + 1;
                    path_list.file_name(idx) = string(write_path);
                end
                

            elseif strcmp(tag, 'depth') == true
                bottom_depth_set = var_env.bottom_depth_set;
                path = 'depth';
                
                path_list.file_name = string(zeros(length(bottom_depth_set), 1));

                for idx = 1:length(bottom_depth_set)
                    Bdry.Bot.depth = bottom_depth_set(idx);
                    SSP.depth(end) = bottom_depth_set(idx);
                    SSP.raw.z(end) = bottom_depth_set(idx);
                    SSP.z(end) = bottom_depth_set(idx);
                    Pos.r.z = linspace(0, bottom_depth_set(idx), length(Pos.r.z));
                    
                    this_file_name = [obj.file_name '_' path '_' num2str(bottom_depth_set(idx))];
                    write_path = fullfile(obj.run_file_root, path, this_file_name);

                    write_env_broadband( write_path, 'kraken', TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, BB_freq);
                    write_fieldflp(write_path, Opt, Pos_flp);

                    path_list.count = path_list.count + 1;
                    path_list.file_name(idx) = string(write_path);

                end
                

            elseif strcmp(tag, 'bot_speed') == true
                bottom_speed_set = var_env.bottom_speed_set;
                path = 'bot_speed';

                path_list.file_name = string(zeros(length(bottom_speed_set), 1));

                for idx = 1:length(bottom_speed_set)
                    Bdry.Bot.cp = bottom_speed_set(idx);
                    Bdry.Bot.HS.alphaR = bottom_speed_set(idx);

                    this_file_name = [obj.file_name '_' path '_' num2str(bottom_speed_set(idx))];
                    write_path = fullfile(obj.run_file_root, path, this_file_name);

                    write_env_broadband( write_path, 'kraken', TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, BB_freq);
                    write_fieldflp(write_path, Opt, Pos_flp);

                    path_list.count = path_list.count + 1;
                    path_list.file_name(idx) = string(write_path);
                end

            elseif strcmp(tag, 'rho') == true
                rho_set = var_env.rho_set;
                path = 'rho';

                path_list.file_name = string(zeros(length(rho_set), 1));
                for idx = 1:length(rho_set)
                    Bdry.Bot.HS.rho = rho_set(idx);
                    Bdry.Bot.rho = rho_set(idx);
                    

                    this_file_name = [obj.file_name '_' path '_' num2str(rho_set(idx))];
                    write_path = fullfile(obj.run_file_root, path, this_file_name);

                    write_env_broadband( write_path, 'kraken', TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, BB_freq);
                    write_fieldflp(write_path, Opt, Pos_flp);

                    path_list.count = path_list.count + 1;
                    path_list.file_name(idx) = string(write_path);
                end


                
            end

        end
        
        
        function [obj, file_list] = krakenBroadBandProcessor(obj, run_full_path, execute_params)
            % broadband 환경 파일 작성
            
            BB_freq = obj.freq_band;
            if size(BB_freq, 1) > 1
                BB_freq = BB_freq';
            end
            assert(~(min(BB_freq) < 0), '양의 freq만 가능');

            obj.is_broadband = true;
            obj.freq_band = BB_freq;

            [file_list] = WriteVariousEnvironmentFiles(obj, run_full_path, BB_freq, execute_params);
            
        end
        
        function [modes] = preprocessModeStructure(obj, filename, freq, prop_mode_flag)
            % 김경섭 박사 논문과 같은 허수부 부호
            % 플래그에 따라 propagation mode만 구할 수 있음
            
            [ ~, ~, ~, ~, ~, pos, ~ ] = read_shd_bin( [filename '.shd'], freq );
            [modes] = read_modes( [filename '.mod'], freq );
                
            
            %modes.k = conj(modes.k);
            %modes.raw_phi = conj(modes.phi);
            modes.raw_phi = modes.phi;

            [modes.phi, modes.phi_zs, modes.z, modes.zs] = modeStructEdit(obj, modes.raw_phi, modes.z, pos);
            modes.pos = pos;
            

            if prop_mode_flag % exclude leaky modes
                %prop_modes = imag(modes.k) < eps;
                k2 = 2 * pi * freq ./ real(modes.Bot.cp);
                prop_modes = real(modes.k) > k2;
                modes.prop_modes = prop_modes;
                modes.k = modes.k(prop_modes);
                modes.phi = modes.phi(:,prop_modes);
                modes.raw_phi = modes.raw_phi(:,prop_modes);
                modes.phi_zs = modes.phi_zs(:,prop_modes);

            else
                prop_modes = true(length(modes.k), 1);
                modes.prop_modes = prop_modes;

            end
        end

        function [ modes ] = getModeGroupVelocities(obj, filename, modes)

            fid = fopen(filename, "rb");
            
            fseek(fid, 4, "cof");
            n_packets = fread(fid, 1, "int32");
            fseek(fid, 4, "cof");
            
            grpv = cell(n_packets, 1);
            
            for idx = 1:n_packets
            
                fseek(fid, 4, "cof");
                freq = fread(fid, 1, "double");
                fseek(fid, 4, "cof");
                
                fseek(fid, 4, "cof");
                mode_len = fread(fid, 1, "int32");
                fseek(fid, 4, "cof");
            
                grpv{idx} = zeros(mode_len, 1);
                
                packet_size = 16;
                packet = fread(fid, packet_size*mode_len, "uint8=>uint8");
                
                slice_map = reshape(1:(mode_len*packet_size), [packet_size, mode_len]);
                slice_idx = 5:12;
                slicer = slice_map(slice_idx,:);
                
                grpv{idx} = typecast(packet(slicer(:)), "double");
                modes{idx}.group_velocity = grpv{idx}(modes{idx}.prop_modes);
                modes{idx}.phase_velocity = 2*pi*modes{idx}.freqVec ./ real(modes{idx}.k);
            end
            fclose(fid);

        end

        
        %function [ obj ] = getFFTPack(obj, fs)
        function [ obj ] = getFFTPack(obj)
            % sampling theory에 맞게 time, freq 벡터 구함
            % green function 관련 계산에서 써먹기 위함
            % Note: 다른 fs로 time series를 구하고 싶으면 이 함수를 수정해야 함
            % ref: https://community.sw.siemens.com/s/article/digital-signal-processing-sampling-rates-bandwidth-spectral-lines-and-more
            obj.fft_pack = [];
            
            df = abs(obj.freq_band(2) - obj.freq_band(1));
            obj.fft_pack.fs = 2 * max(obj.freq_band);
            %obj.fft_pack.fs = fs;
            %T = 1;
            %obj.fft_pack.nfft = obj.fft_pack.fs * T;
            %obj.fft_pack.nfft = 2^nextpow2(round(obj.fft_pack.fs / df));
            obj.fft_pack.nfft = obj.fft_pack.fs / df;

            if ~isinteger(obj.fft_pack.nfft)
                obj.fft_pack.nfft = round(obj.fft_pack.fs / df);
                if ~obj.is_precision_warning
                    warning('nfft 반올림 됨 / 정밀도 오차인지 확인')
                    obj.is_precision_warning = true;
                end
            end
            
            %obj.fft_pack.fftfreq = linspace(0, obj.fft_pack.fs, obj.fft_pack.nfft);
            obj.fft_pack.fftfreq = (0:obj.fft_pack.nfft-1)*df; 
            obj.fft_pack.time = (0:obj.fft_pack.nfft-1)/obj.fft_pack.fs;

            
        end
        
        function [ s ] = getResampledSignal(obj, s, resampled_params)
            
            for idx = 1:size(s,1)
                for jdx = 1:size(s,2)
                    s{idx,jdx} = resample(s{idx,jdx}, resampled_params.p, resampled_params.q);
                end
            end
        
        end
    end
    
end


