classdef modeProcessor
    % mode coefficient ��� �� ������ ���� modeFilter ���� �Լ� ����
    % 
    % methods
    %   public:
    %       modeProcessor(*)
    %       getExactModeCoefficients(*) - ��Ȯ�� ���� ��� ��� ��ȯ
    %       getEstimatedModeCoefficients(*) - ��� ���ͷ� ������ ��� ��� ��ȯ
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
            %MODEPROCESSOR �� Ŭ������ �ν��Ͻ� ����
            % pressure: ( #recv x 1 )
            % phi: ��� �Լ� ( #recv x #modes )
            % phi_zs: src depth������ ��� �Լ� ( 1 x #modes )
            % k: eigenvalues ( #modes x 1 )
            
            if(size(pressure, 2) ~= 1)
                obj.pressure = pressure.';
            else
                obj.pressure = pressure;
            end
            
            num_recv = [size(pressure, 1), size(phi, 1)];
            assert(1 == length(unique(num_recv)), '���ù� �� �ȸ���');
            
            num_modes = [size(phi, 2), length(k), size(phi_zs,2)];
            assert(1 == length(unique(num_modes)), '��� �� �ȸ���');
            
            obj.eig_func_zs = phi_zs;
            obj.eig_func = phi;
            obj.eig_val = k;
            
        end
        
        function [m_coeff] = getExactModeCoefficients(obj, src_rho, range)
            % ����� ���� ��Ȯ�� ��� ��ȯ
            % src_rho: rho(Zs)
            % range: src - recv �� �Ÿ�
            % m_coeff: ����ȭ�� ��� ���
            
            assert(src_rho > eps || range > eps, '�߸��� �Ű�����')
            
            const_coeff = (sqrt(2*pi) * exp(1j*(pi/4))) / src_rho;
            mode_term = exp(1j*obj.eig_val*range) ./ sqrt(obj.eig_val*range) .* obj.eig_func_zs.';
            m_coeff = const_coeff .* mode_term;
            m_coeff = m_coeff ./ max(m_coeff);
            
        end
        
        function [m_coeff] = getEstimatedModeCoefficients(obj, method)
            % ��� ���ͷ� ������ ��� ��ȯ
            % method: mode filter ����
            % m_coeff: ����ȭ�� ��� ���
            
            obj.mode_filter = modeFilter(obj.eig_func);

            m_coeff = obj.mode_filter.process(obj.pressure, method);
            m_coeff = m_coeff ./ max(m_coeff);
            
        end
        
    end
end

