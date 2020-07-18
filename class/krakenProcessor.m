classdef krakenProcessor
    % KRAKENC.exe/field.exe 실행하고 데이터 뽑아냄
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
        
        fft_pack
        is_precision_warning
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
        end
        
        
        function [ obj ] = execute(obj)
            % 사본 파일이 위치한 곳에서 Kraken 실행
            % env 파일로 .mod 파일 생성
            % flp 파일로 .shd 파일 생성
            
            dir_info = dir([obj.file_root '\*.env']);
            assert(~isempty(dir_info), '잘못된 경로');
            
            
            run_full_path = fullfile(obj.run_file_root, obj.file_name);
            
            if obj.is_broadband
                [obj, Opt, Pos_flp] = krakenBroadBandProcessor(obj, run_full_path);
            else % single freq
                [obj, Opt, Pos_flp] = envParamsReader(obj, run_full_path);
            end
 
            write_fieldflp(run_full_path, Opt, Pos_flp);
            
            system(['krakenc ' run_full_path]);
            system(['field ' run_full_path]);
            %krakenc(run_full_path);
            %field(run_full_path);
            
            fclose('all'); % adhoc
        end
        
        
        function [ modes ] = getModeStruct(obj, varargin)
            % 계산된 .mod 파일에서 모드 반환
            % varargin: ex) getModeStruct('propagation')
            %           - propagation modes vs leaky modes
            %           default: leaky modes
            
            
            narginchk(1,2);
            prop_mode_flag = false;
            if ~isempty(varargin) % exclude leaky mode
                narginchk(2,2);    

                if strcmpi(varargin{1}, 'PROPAGATION') 
                    prop_mode_flag = true;
                else
                    error('Propagation 매개변수 확인');
                end
            end
            
            filename = fullfile(obj.run_file_root, obj.file_name);
            modes = cell(length(obj.freq_band),1);
            
            parfor idx = 1:length(obj.freq_band)
                
                [modes{idx}] = preprocessModeStructure(obj, filename, obj.freq_band(idx), prop_mode_flag);
                
                modes{idx}.freqVec = obj.freq_band(idx);
                modes{idx}.Nfreq = idx;
            end

        end
        
        function [ sig_t ] = getTimeSeries(obj, g, time_delay)
            % Green 함수에서 time series 생성
            % g: Green spectrum
            % time_delay: 시간 지연 ex) range(m)/1500(m/s)
            % sig_t: time series
            
            narginchk(3,3);

            fftpack = obj.fft_pack;
            %sig_f = zeros(fftpack.nfft, size(g, 2));
            g_pad = zeros(fftpack.nfft, size(g, 2));
            
            g = flipud(g);
            g_pad(1:size(g,1),:) = g;
            sig_f = g_pad.*exp(1j.*2*pi*obj.fft_pack.fftfreq'.*time_delay);
            sig_t = ifft(sig_f, 'symmetric');
            
            sig_t = sig_t / vecnorm(sig_t);
            
            %sig_t = flipud(sig_t);
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
                mode_g(valid_idx,idx) = ...
                    freq_phi_zr(valid_idx,idx) .* freq_phi_zs(valid_idx,idx) ...
                    ./ sqrt(freq_eigenvals(valid_idx, idx)*r) .*besselh(0, freq_eigenvals(valid_idx,idx)*r);
            end
            
            fftpack = obj.fft_pack;
            
            mode_f = zeros(fftpack.nfft,length(modes{end}.k));
            mode_f(1:length(obj.freq_band),:) = mode_g;
            
            mode_t = ifft(mode_f, 'symmetric');
            mode_t = flipud(mode_t);
        end
        
        
        function [g] = getGreenSpectrum(obj, modes, r, varargin)
            % Green spectrum 구함
            % modes: getModeStruct(*)에서 얻어진 모드 구조체
            % r: src-recv간 거리(m)
            % varargin: recv 수심 <- 해당 수심에서의 mode-pressure 얻어냄
            %           default: 모든 수심(구현 안한듯. 기억이 안남)
            % g: Green spectrum
            %
            
            ref_mode = modes{end};
            
            if ~isempty(varargin)
               recv_depth = varargin{1};
               [~, recv_idx] = min(abs(recv_depth - ref_mode.z));
            else
                recv_idx = 1:length(ref_mode.z);
            end
            
            g = zeros(length(obj.freq_band), length(recv_idx));
            
            for idx = 1:length(modes)
                
                phi_z = modes{idx}.phi(recv_idx,:);
                phi_zs = modes{idx}.phi_zs;
                k = modes{idx}.k;

                g(idx,:) = sum(phi_z.*phi_zs ./ sqrt(k*r).' .* besselh(0, k*r).', 2);

            %     %asymptotic form
            %     
            %     asympt_form = phi_z .* phi_zs ./ sqrt(k*r).' .* (sqrt(2./(pi*k*r)) .* exp(1j.*(k*r - pi/4))).';
            %     hankel = (sqrt(2./(pi*k*r)) .* exp(1j.*(k*r - pi/4))).';
            %
            %     grn = asympt .* hankel;
            %     %grn = grn ./ vecnorm(grn);
            %     all_g(idx,:) = sum(grn, 2);
            end
            
            g = g ./ vecnorm(g);
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
        
        function [new_phi, phi_zs] = modeStructEdit(obj, phi, z, pos)
        % 소스 위치와 리시버 위치가 불일치할 경우 고유 함수가 src 위치에서도 나옴
        % zs, z의 고유 함수를 따로 관리
        % phi: [z, zs]에서의 고유 함수 (num recv+불일치 소스 수 x 모드 수)
        % new_phi: z에서의 고유 함수 (num recv x 모드 수)
        % phi_zs: zs에서의 고유함수

            recv_depth = pos.r.depth;
            src_depth = pos.s.depth;

            z_redundant = setdiff(z, recv_depth);
            if(any(z_redundant)) % src pos != recv pos
                recv_idx = z_redundant == z;
                phi_zs = phi(recv_idx, :);
                new_phi = phi(~recv_idx, :);
            else % src pos == recv pos
                src_idx = src_depth == recv_depth;
                phi_zs = phi(src_idx,:);
                new_phi = phi;
            end

            if size(phi_zs, 2) == 1
                phi_zs = phi_zs.';
            end
        end
        
        
        function [obj, Opt, Pos_flp] = krakenBroadBandProcessor(obj, run_full_path)
            % broadband 환경 파일 작성
            
            BB_freq = obj.freq_band;
            if size(BB_freq, 1) > 1
                BB_freq = BB_freq';
            end
            assert(~(min(BB_freq) < 0), '양의 freq만 가능');

            [ obj, Opt, Pos_flp, TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax ] ...
                = envParamsReader(obj, run_full_path);
            
            obj.is_broadband = true;
            obj.freq_band = BB_freq;
            
            write_env_broadband( run_full_path, 'kraken', ...
                TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, BB_freq)

        end
        
        function [modes] = preprocessModeStructure(obj, filename, freq, prop_mode_flag)
            % 김경섭 박사 논문과 같은 허수부 부호
            % 플래그에 따라 propagation mode만 구할 수 있음
            
            [ ~, ~, ~, ~, pos, ~ ] = read_shd_bin( [filename '.shd'], freq );
            [modes] = read_modes( [filename '.mod'], freq );
                
            modes.raw_phi = double(conj(modes.phi));
            modes.raw_phi = modes.raw_phi ./ vecnorm(modes.raw_phi);
            
            modes.k = double(conj(modes.k));

            [modes.phi, modes.phi_zs] = modeStructEdit(obj, modes.raw_phi, modes.z, pos);

            if prop_mode_flag % exclude leaky modes
                prop_modes = imag(modes.k) < eps;
                modes.k = modes.k(prop_modes);
                modes.phi = modes.phi(:,prop_modes);
                modes.raw_phi = modes.raw_phi(:,prop_modes);
                modes.phi_zs = modes.phi_zs(:,prop_modes);

            end
        
        end
        
        function [ obj ] = getFFTPack(obj)
            % sampling theory에 맞게 time, freq 벡터 구함
            % green function 관련 계산에서 써먹기 위함
            % Note: 다른 fs로 time series를 구하고 싶으면 이 함수를 수정해야 함
            
            obj.fft_pack = [];
            
            df = abs(obj.freq_band(2) - obj.freq_band(1));
            obj.fft_pack.fs = 2 * max(obj.freq_band);
            obj.fft_pack.nfft = obj.fft_pack.fs / df;
            if ~isinteger(obj.fft_pack.nfft)
                obj.fft_pack.nfft = round(obj.fft_pack.fs / df);
                if ~obj.is_precision_warning
                    warning('nfft 반올림 됨 / 정밀도 오차인지 확인')
                    obj.is_precision_warning = true;
                end
            end
            obj.fft_pack.fftfreq = (0:obj.fft_pack.nfft-1)*df; 
            obj.fft_pack.time = (0:obj.fft_pack.nfft-1)/obj.fft_pack.fs;
        end
        
        
        
    end
    
end


