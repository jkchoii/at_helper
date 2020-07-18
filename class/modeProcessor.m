classdef modeProcessor
    % mode coefficient 계산 및 추정을 위한 modeFilter 내부 함수 쿼리
    % 
    % methods
    %   public:
    %       modeProcessor(*)
    %       getExactModeCoefficients(*) - 정확히 계산된 모드 계수 반환
    %       getEstimatedModeCoefficients(*) - 모드 필터로 추정된 모드 계수 반환
    %   private:
    %   
    
    properties(SetAccess = private)
        pressure
        eig_func
        eig_val
        eig_func_zs
        mode_filter
        
    end
    
    methods
        function obj = modeProcessor(pressure, phi, phi_zs, k)
            %MODEPROCESSOR 이 클래스의 인스턴스 생성
            % pressure: ( #recv x 1 )
            % phi: 모드 함수 ( #recv x #modes )
            % phi_zs: src depth에서의 모드 함수 ( 1 x #modes )
            % k: eigenvalues ( #modes x 1 )
            
            if(size(pressure, 2) ~= 1)
                obj.pressure = pressure.';
            else
                obj.pressure = pressure;
            end
            
            num_recv = [size(pressure, 1), size(phi, 1)];
            assert(1 == length(unique(num_recv)), '리시버 수 안맞음');
            
            num_modes = [size(phi, 2), length(k), size(phi_zs,2)];
            assert(1 == length(unique(num_modes)), '모드 수 안맞음');
            
            obj.eig_func_zs = phi_zs;
            obj.eig_func = phi;
            obj.eig_val = k;
            
        end
        
        function [m_coeff] = getExactModeCoefficients(obj, src_rho, range)
            % 계산을 통한 정확한 모드 반환
            % src_rho: rho(Zs)
            % range: src - recv 간 거리
            % m_coeff: 정규화된 모드 계수
            
            assert(src_rho > eps || range > eps, '잘못된 매개변수')
            
            const_coeff = (sqrt(2*pi) * exp(1j*(pi/4))) / src_rho;
            mode_term = exp(1j*obj.eig_val*range) ./ sqrt(obj.eig_val*range) .* obj.eig_func_zs.';
            m_coeff = const_coeff .* mode_term;
            m_coeff = m_coeff ./ max(m_coeff);
            
        end
        
        function [m_coeff] = getEstimatedModeCoefficients(obj, method)
            % 모드 필터로 추정한 모드 반환
            % method: mode filter 선택
            % m_coeff: 정규화된 모드 계수
            
            obj.mode_filter = modeFilter(obj.eig_func);

            m_coeff = obj.mode_filter.process(obj.pressure, method);
            m_coeff = m_coeff ./ max(m_coeff);
            
        end
        
    end
end

