classdef modeFilter
    % mode filter ���� �� mode coefficients ���
    % ����� ��漷�ڻ� ���� �� ����
    %
	%
    % methods
    %   public:
    %       modeFilter(*)
    %       process(*)
    %   private:
    %       sampledModeShape(*): H = E^H
    %       pseudoInverseNaive(*): H = (E^H * E)^-1 * E^H
    %       pseudoInverseEVD(*): H = hat( (E^H * E)^-1 ) * E^H
    %       diagonalWeighting(*): H = (E^H*E + beta*I)^-1 * E^H
    %
    % TODO:
    %       recursive calc for pseudoInverseEVD
    %       maximum a posterior method
    
    
    properties(SetAccess = private)
        H
        phi

    end
    
    methods(Access = public)
        
        function obj = modeFilter(phi)
            
            obj.H = [];
            obj.phi = phi;
            
        end
        
        function [m_coeff] = process(obj, p, method, beta)
            
            assert(isstring(method) || ischar(method), '������ �Ǵ� ��Ʈ���̾�� ��')
            
            if(strcmp(method, 'sampledModeShape'))
                obj = sampledModeShape(obj.phi);
                
            elseif(strcmp(method, 'pseudoInverseNaive'))
                obj = pseudoInverseNaive(obj, obj.phi);
                
            elseif(strcmp(method, 'pseudoInverseEVD'))
                obj = pseudoInverseEVD(obj, obj.phi);
                
            elseif(strcmp(method, 'diagonalWeighting'))
                assert(beta>eps, 'beta �ʹ� ����')
                obj = diagonalWeighting(obj, obj.phi, beta);
                
            else
                assert(false, '���� ����');
            end
            
            m_coeff = obj.H*p;
            
        end
        
    end
    
    methods(Access = private)
        function [obj] = sampledModeShape(obj, phi)
            obj.H = phi';
        end
        
        function [obj] = pseudoInverseNaive(obj, phi)
            %obj.H = pinv(phi);
            obj.H = (phi'*phi)\phi';
        end
        
        function [obj] = pseudoInverseEVD(obj, phi)
            % truncated EVD( ���Ǽ� ~1% �̸� ����)
            
            [e, k] = eig(phi'*phi);
            k = diag(k);

            valid_eig = (k./max(k))/100 < k./max(k); %  more 1%

            e = e(:,valid_eig);
            k = k(valid_eig);

            reconstruct = zeros(length(e), length(e), length(k));
            for idx = 1:length(k)
                reconstruct(:,:,idx) = 1./sqrt(k(idx)) * (e(:,idx) * e(:,idx)');
            end
            
            obj.H = sum(reconstruct, 3) * phi'; % bar inv(E^H*E) * E^H

        end
        
        function [obj] = diagonalWeighting(obj, phi, beta)
            obj.H = phi'*phi+beta*ones(size(phi, 2), 1)\phi';
        end
        
    end
end

