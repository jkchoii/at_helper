classdef modeWarping
    %MODEWARPING 이 클래스의 요약 설명 위치
    %   자세한 설명 위치
    
    properties
        tmin
        tmax
        c
        r
        fft_pack
        dt

    end
    
    methods
        function obj = modeWarping(sig, fft_pack, c, r)
            if iscolumn(sig)
                sig = sig.';
            end
            
            sig_length = length(sig);
            obj.dt = 1/fft_pack.fs;
            obj.tmin = r/c + obj.dt;
            obj.tmax = sig_length/fft_pack.fs + r/c;
            
            obj.r = r;
            obj.c = c;
            obj.fft_pack = fft_pack;
            
        end
        
        function [ sig, params ] = process(obj, signal, flag, varargin)
            narginchk(3,4);
            
            if iscolumn(signal)
                signal = signal.';
            end
            
            if strcmpi(flag, 'FORWARD')
                [ sig, params ] = processForward(obj, signal);

            elseif strcmpi(flag, 'INVERSE')
                narginchk(4,4);
                if ~isempty(varargin)
                    params = varargin{1};
                    [ sig, params ] = processInverse(obj, signal, params);
                end
            else
                error('flag 확인')
            end
            
        end
    end
    
    methods(Access = private)
        function [t_warped] = warp_t(obj, t)
            t_warped = sqrt(t.^2 + (obj.r/obj.c).^2);
        end
        
        function [t_iwarped] = iwarp_t(obj, t)
            t_iwarped = sqrt(t.^2 - (obj.r/obj.c).^2);
        end
        
        function [ sig_warp, warp_params ] = processForward(obj, signal)
            
            %dt_warped = iwarp_t(obj.tmax) - iwarp_t(obj.tmax - obj.dt);
            dt_warp = 1/obj.fft_pack.fs * obj.tmax / iwarp_t(obj, obj.tmax);
            freq_warp = 1/dt_warp;
            fs_warp = 2 * freq_warp;

            t_warp_max = iwarp_t(obj, obj.tmax);
            M = ceil(t_warp_max * fs_warp);
            
            % uniform sampling
            t_warp_homo = (0:M-1)/fs_warp;
            
            % non-uniform sampling
            t_warp_inhomo = warp_t(obj, t_warp_homo);
            
            coeff = sqrt(t_warp_homo ./ t_warp_inhomo);
            
            %% sinc interp
            [Ts, T] = ndgrid(t_warp_inhomo ,obj.fft_pack.time + obj.tmin);
            sig_warp = coeff.' .* sinc(obj.fft_pack.fs*(Ts-T))*signal.';
            
            
            %%
            warp_params = [];
            warp_params.t = t_warp_homo;
            warp_params.warp_fs = fs_warp;
            warp_params.origin_fs = obj.fft_pack.fs;
            warp_params.origin_length = length(signal);
            warp_params.warp_length = length(sig_warp);
        end
        
        function [ sig_iwarp, params ] = processInverse(obj, signal, params)
            
            t = (1:params.origin_length)/params.origin_fs + obj.r/obj.c;
            t_iwarp = iwarp_t(obj, t);
            coeff = sqrt(t ./ t_iwarp);
            
            %% sinc interp
            
            [Ts, T] = ndgrid(t_iwarp, (0:params.warp_length-1)/params.warp_fs);
            sig_iwarp = real(coeff.' .* sinc(params.warp_fs*(Ts-T))*signal.');
            
            %%
            params.t = t;
        end
        
    end
end

